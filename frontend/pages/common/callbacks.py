import json
import math
from logging import getLogger
import time
from functools import partial

import requests
import dash
from dash import dcc, html, dash_table
import pandas as pd
import plotly.graph_objects as go

from config import backend_root
from pages.common.pairplot import (
    get_individual_transforms_df,
    get_group_by_fragment_df,
    get_group_by_fragment_table,
    get_prop_labels,
    rename_column_headers,
    map_colors_to_ids,
    get_range_filter_args,
    get_styles,
    generate_pair_tables,
    df_round
)

logger = getLogger(__name__)


def get_individual_transforms_table(create_id, df, default_x_label, default_y_label):

    initial_data = df.to_dict('records')
    transform_table = dash_table.DataTable(
        id=create_id('table_transforms'),
        # Dash seems to only wrap text at spaces or - marks
        # But property column names can be long with multiple _, and Dash is setting widths of these columns to be narrow sometimes
        columns=[dict(name=i.replace('_', ' '), id=i, type='text', presentation='markdown') for i in df.columns if i not in ['id']],
        data=initial_data,
        fixed_rows={'headers': True},
        style_header={
            'whiteSpace': 'normal',
            'height': 'auto',
            'backgroundColor': '#86818a',
            'color': '#FFFFFF',
            'fontWeight': 'bold',
            'textAlign': 'left',
            'fontSize': 15,
        },
        editable=False,
        filter_action="native",
        sort_action="native",
        sort_mode='multi',
        row_selectable='multi',
        row_deletable=False,
        selected_rows=[],
        page_action='native',
        page_current=0,
        page_size=15,
        markdown_options={'link_target': '_blank', 'html': True},
        style_data_conditional=(get_styles(df, statistic='median')),
        style_table={'height': 500, 'maxHeight': 500, 'overflowY': 'auto'},
        active_cell={'row': 0, 'column': 0, 'row_id': initial_data[0]['id'], 'column_id': 'FROM'} if len(df.index) > 0 else None
    )
    return transform_table


def get_group_by_fragment_table(create_id, df, default_x_label, default_y_label):
    initial_data_dict = df.to_dict('records')
    transform_table = dash_table.DataTable(
        id=create_id('table_transforms'),
        # Dash seems to only wrap text at spaces or - marks
        # But property column names can be long with multiple _, and Dash is setting widths of these columns to be narrow sometimes
        columns=[dict(name=i.replace('_', ' '), id=i, type='text', presentation='markdown') for i in df.columns if i not in ['id', 'rule_id_array']],
        data=initial_data_dict,
        # fixed_columns={'headers': True, 'data': 4},
        fixed_rows={'headers': True},
        editable=False,
        filter_action="native",
        sort_action="native",
        sort_mode='multi',
        row_selectable='multi',
        row_deletable=False,
        selected_rows=[],
        page_action='native',
        page_current=0,
        page_size=20,
        markdown_options={'link_target': '_blank', 'html': True},
        style_data_conditional=(get_styles(df, statistic='median median')),
        style_table={'height': 500, 'maxHeight': 500, 'overflowY': 'auto', 'minWidth': '100%'},
        style_header={
            'whiteSpace': 'normal',
            'height': 'auto',
            'backgroundColor': '#86818a',
            'color': '#FFFFFF',
            'fontWeight': 'bold',
            'textAlign': 'left',
            'fontSize': 12,
        },
        style_cell={
            'textAlign': 'left',
            'minWidth': '100px'  # , 'width': '180px', 'maxWidth': '180px',
        },
        active_cell={'row': 0, 'column': 0, 'row_id': initial_data_dict[0]['id'], 'column_id': 'TO'}
    )
    return transform_table


def run_persistent_query(input_data):
    input_data_dict = json.loads(input_data)
    # To load a new example query snapshot to the database upon DB build, uncomment below line, create a shareable link after running a query, then load that link
    # Then copy and paste the below dictionary that is printed out, into matcher/backend/initialize_db/example_queries.json,
    #  taking care to match the format of the other JSON lines that are in example_queries.json (e.g. delete extra fields such as schema and query_id)
    # logging.info(input_data_dict)
    if input_data_dict.get('query_id') != -1:
        # This means we loaded a snapshot that already has saved results in the DB
        # Therefore, instead of running that query again and duplicating those results, we just skip the query and reference those previous results
        return [input_data, json.dumps({'query_id': input_data_dict['query_id']})]

    # Start the query and get an ID which we will use to poll the DB for existence of results over intervals
    schema = input_data_dict['schema']
    started_query = requests.post(backend_root + f'/start_query?schema={schema}', data=input_data)
    started_query = json.loads(started_query.json())

    # If there were problems with input from UI, terminate here
    if 'observations' in started_query:
        return [input_data, started_query]

    query_id = started_query['query_id']

    # Set the time cutoff, in seconds, where we should give up on waiting for the query to finish
    timeout = 3600
    elapsed = 0
    start_time = time.perf_counter()
    # Set the time, in seconds, for how frequently we poll the DB for existence of results that indicate query completion
    # We start with 2 seconds between each check, then increase interval as query becomes longer
    intervals = [
        {'interval_time': 2, 'start_next_interval_at': 14},
        {'interval_time': 5, 'start_next_interval_at': 90},
        # Purposefully use a time longer than the timeout so that this last interval is used until the end
        {'interval_time': 10, 'start_next_interval_at': timeout * 2},
    ]
    current_interval = 0

    while elapsed < timeout:

        time.sleep(intervals[current_interval]['interval_time'])
        query_status = requests.get(backend_root + f'/check_query/{query_id}?schema={schema}')
        query_status = json.loads(query_status.json())

        if query_status['finished'] == "True":
            if query_status.get('finished_with_no_results') == "True":
                started_query['observations'] = {'no_results': "True"}
                return [input_data, started_query]
            elif query_status.get('finished_with_exception') == "True":
                started_query['observations'] = {'exception': "True"}
                return [input_data, started_query]
            else:
                #  The only case leading to display of result data, if we started a query
                # Connect output data (via query_id, which is used in the DB to represent the unique query to which a results set belongs) to the input state,
                # for use in snapshot generation
                input_data_dict['query_id'] = query_id
                input_data = json.dumps(input_data_dict)

                # If we loaded a snapshot that is not referencing a query_id in the DB, we need to associate snapquery_id and query_id here
                if input_data_dict['snapquery_id'] != '':
                    response = requests.get(backend_root + f"/bind_snapquery_to_results/{input_data_dict['snapquery_id']}?query_id={query_id}&schema={schema}")
                return [input_data, json.dumps({'query_id': query_id})]

        elapsed = time.perf_counter() - start_time
        if elapsed > intervals[current_interval]['start_next_interval_at']:
            current_interval += 1

    started_query['observations']['query_timeout'] = "True"
    return [input_data, json.dumps(started_query)]


def instantiate_output_elements(create_id, pair_data, query_data, property_metadata, display_name_to_property_name):
    pair_data = json.loads(pair_data)
    if "observations" in pair_data:
        # This means that something went wrong with the query, so we are returning an empty layout
        # We need to return a pair_tables div, because the iframe resizing callback depends on this div as an input
        return [html.Div(
            children=[
                html.Div(
                    ['no_output'],
                    id=create_id('pair_tables'),
                    style={'display': 'none'}
                )]
        )]

    query_id = pair_data["query_id"]

    snapquery_dict = json.loads(query_data)
    aggregation_type = snapquery_dict['advanced_options']['aggregation_type']

    """
    # identifier associated with rows in the transform table
    if aggregation_type == 'individual_transforms':
        identifier = 'rule_id'
    elif aggregation_type == 'group_by_fragment':
        identifier = 'TO'
    else:
        raise ValueError("invalid query_type, query_type must be 'individual_transforms' or 'group_by_fragment'")
    """

    if snapquery_dict['snapfilter_string'] != '':
        snapfilter_exists = True
    else:
        snapfilter_exists = False
        # Initialize range filters key because we reference it during range_filters initialization, even if there's no snapfilter
        snapfilter = {
            'range_filters': {},
            'transform_row_ids': [],
        }

    rf_args = []
    if snapfilter_exists:
        snapfilter_list = snapquery_dict['snapfilter_string'].split(';')

        transform_row_ids = snapfilter_list[4].split(',')
        # We migrated from Oracle to Postgres in Nov 2021, causing snapfilters with id < 322 to have transform_row_ids that are outdated
        # Therefore we need to not use these ids, if snapfilter id is less than 322
        if transform_row_ids != [''] and int(snapquery_dict['snapfilter_id']) > 321:
            if aggregation_type == 'individual_transforms' and '_' not in transform_row_ids[0]:
                transform_row_ids = [int(x) for x in transform_row_ids]
        else:
            transform_row_ids = []

        # Example range_filters_list is 'A_Papp,None,20@B_Papp,20,None'
        range_filters_list = snapfilter_list[5].split('@')
        range_filters_dict = {}
        if range_filters_list != ['']:
            rf_names, min_rfs, max_rfs = [0] * len(range_filters_list), [0] * len(range_filters_list), [0] * len(range_filters_list)
            for i, rf_name_min_max in enumerate(range_filters_list):
                rf_names[i], min_rfs[i], max_rfs[i] = rf_name_min_max.split(',')
                # for applying filter values to UI elements
                range_filters_dict[rf_names[i]] = {
                    'min': float(min_rfs[i]) if min_rfs[i] != 'None' else None,
                    'max': float(max_rfs[i]) if max_rfs[i] != 'None' else None,
                }
            rf_args = get_range_filter_args(rf_names, min_rfs, max_rfs, display_name_to_property_name)

        snapfilter = {
            'default_x_label': snapfilter_list[0],
            'default_y_label': snapfilter_list[1],
            'default_x_type': snapfilter_list[2],
            'default_y_type': snapfilter_list[3],
            'transform_row_ids': transform_row_ids,
            'range_filters': range_filters_dict
        }

    separator = '-'

    # Use a set in case the user selects the same properties between required and optional properties
    all_props = set(snapquery_dict["REQUIRED_properties"].split(",") + snapquery_dict["OPTIONAL_properties"].split(","))
    prop_labels = get_prop_labels(all_props, property_metadata=property_metadata)

    # Set the default axis labels on the plot, which are populated in the two dropdown menus
    if snapfilter_exists:
        default_x_label, default_y_label, default_x_type, default_y_type = (
            snapfilter['default_x_label'], snapfilter['default_y_label'], snapfilter['default_x_type'], snapfilter['default_y_type']
        )
        initial_agg_prop_labels = []
        for prop_label in prop_labels:
            if default_x_label in prop_label.values() or default_y_label in prop_label.values():
                initial_agg_prop_labels.append(prop_label)
    else:
        initial_agg_prop_labels = prop_labels[0:2]
        if len(prop_labels) == 1:
            default_x_label = prop_labels[0]['A']
            default_x_type = prop_labels[0]['default_A_axis_type']
            default_y_label = prop_labels[0]['change_label']
            default_y_type = prop_labels[0]['default_change_axis_type']
        elif len(prop_labels) > 1:
            default_x_label = prop_labels[0]['change_label']
            default_x_type = prop_labels[0]['default_change_axis_type']
            default_y_label = prop_labels[1]['change_label']
            default_y_type = prop_labels[1]['default_change_axis_type']

    # On load, we aggregate statistics based on the first 2 user-selected properties (or 1, if only 1 was selected)
    agg_args = {
        'query_id': query_id,
        'aggregation_type': aggregation_type,
        'range_filters': rf_args,
        'statistics': [
            {
                'statistic': 'median',
                'property_name': prop['prop'],
                'change_type': prop['change_type'],
                'base': prop['base'],
                'units': prop['units'],
            } for prop in initial_agg_prop_labels
        ]
    }
    schema = snapquery_dict['schema']
    table_data = requests.post(backend_root + f'/aggregate_transforms?schema={schema}', data=json.dumps(agg_args))
    table_data = table_data.json()
    minmax = table_data['minmax']
    # The minmax columns were named within postgres, which restricts certain characters
    # We need the minmax column names to match the names in the x/y dropdowns, therefore we need to rename some columns
    renamed_minmax_headers = rename_column_headers(list(minmax.keys()), minmax=True, property_metadata=property_metadata)
    minmax = dict(zip(renamed_minmax_headers, minmax.values()))
    # Detect whether the environment was included in the aggregation, which is only needed for multicut queries
    environment_included = False
    for column in table_data['column_headers']:
        if 'smiles_env' in column:
            environment_included = True
            break
    columns = rename_column_headers(table_data['column_headers'], property_metadata=property_metadata)
    table_data = table_data['rows']
    df = pd.DataFrame(table_data, columns=columns)

    if environment_included:
        if aggregation_type == 'individual_transforms':
            df['id'] = df['FROM'].apply(lambda x: str(x) + '_') + df['TO']
        elif aggregation_type == 'group_by_fragment':
            df['id'] = df['TO'].copy()
    else:
        # By default for individual_transforms with environment_included == False, the API call will return a column with 'rule_id' header
        # which gets changed to 'id' header by rename_column_headers
        if aggregation_type == 'group_by_fragment':
            df['id'] = df['TO'].copy()

    colors_mapped_to_ids = map_colors_to_ids(aggregation_type, df, environment_included)

    # Approach with dropdowns, radiobuttons, and callback scatter plotting adapted from https://dash.plotly.com/interactive-graphing
    # UI Component definitions

    x_y_dropdown_options = []
    for label in ['change_label', 'A', 'B', 'average_label']:
        for prop in prop_labels:
            x_y_dropdown_options.append({'label': prop[label], 'value': prop[label]})

    x_Dropdown = dcc.Dropdown(
        id=create_id('crossfilter-xaxis-column'),
        options=x_y_dropdown_options,
        value=default_x_label
    )

    x_RadioItems = dcc.RadioItems(
        id=create_id('crossfilter-xaxis-type'),
        options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
        value=default_x_type,
        labelStyle={'display': 'inline-block', 'marginTop': '5px'}
    )

    y_Dropdown = dcc.Dropdown(
        id=create_id('crossfilter-yaxis-column'),
        options=x_y_dropdown_options,
        value=default_y_label
    )

    y_RadioItems = dcc.RadioItems(
        id=create_id('crossfilter-yaxis-type'),
        options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
        value=default_y_type,
        labelStyle={'display': 'inline-block', 'marginTop': '5px'}
    )

    main_scatterplot = dcc.Graph(
        id=create_id('pairPlot')
    )

    range_filters = []
    range_filter_options = []
    for label in ['A', 'B']:
        for prop in prop_labels:
            range_filter_options.append(prop[label])

    for column in range_filter_options:
        range_filters.append(
            (
                column,
                dcc.Input(
                    id={'type': create_id('min_range_filter'), 'index': column},
                    type='number',
                    debounce=True,
                    style={'width': '45%', 'display': 'inline-block'},
                    value=snapfilter['range_filters'][column]['min'] if column in snapfilter['range_filters'] else None,
                ),
                dcc.Input(
                    id={'type': create_id('max_range_filter'), 'index': column},
                    type='number',
                    debounce=True,
                    style={'width': '45%', 'display': 'inline-block'},
                    value=snapfilter['range_filters'][column]['max'] if column in snapfilter['range_filters'] else None,
                )
            )
        )
    range_filter_names = json.dumps(range_filter_options)

    range_filters_dropdown = dcc.Dropdown(
        id=create_id('range_filters_dropdown'),
        options=[{'label': name, 'value': name} for name in range_filter_options],
        value=range_filters[0][0]
    )

    if aggregation_type == 'individual_transforms':
        initialize_transform_table = get_individual_transforms_table
        get_updated_transform_table_df = get_individual_transforms_df
    elif aggregation_type == 'group_by_fragment':
        initialize_transform_table = get_group_by_fragment_table
        get_updated_transform_table_df = get_group_by_fragment_df
    else:
        raise ValueError("aggregation_type has invalid value, must be 'individual_transforms' or 'group_by_fragment'")

    df = get_updated_transform_table_df(df)
    transform_table = partial(initialize_transform_table, create_id)(df, default_x_label, default_y_label)

    initial_agg_data = df.to_dict('records')
    total_row_count = len(initial_agg_data)

    # Above we calculated num_colors based on having a unique color for each ID in the transform_table
    # This number will shrink if we are looking at any combination of properties other than the single most common property in the query
    initial_total_row_count = df['id'].nunique()
    num_rows_displayed_Plaintext = html.Plaintext(
        id=create_id('num_rows_displayed_plaintext'),
        children='Displaying ' + str(initial_total_row_count) + ' out of ' + str(initial_total_row_count) + ' rows',
        style={'font-size': '15px', 'font-weight': 'bold'}
    )

    children = [
        dcc.Store(id=create_id('snapfilter'), data=json.dumps(snapfilter)),
        dcc.Store(id=create_id('colors_mapped_to_ids'), data=[json.dumps(colors_mapped_to_ids)]),
        # Use in the callback that updates transform_table, to enforce initialization of row filtration, using snapfilter['transform_row_ids']
        # dcc.Store(id='snapfilter_applied', data= [{'snapfilter_applied': False if snapfilter_exists else True}]),
        dcc.Store(id=create_id('snapfilter_applied'), data=False if snapfilter_exists else True),
        html.Button("not_intended_to_display", id=create_id("start_highlight_first_row"), n_clicks=0, style={'display': 'none'}),
        html.Button("not_intended_to_display", id=create_id("finish_highlight_first_row"), n_clicks=0, style={'display': 'none'}),
        dcc.Store(id=create_id('identifier'), data='id'),
        dcc.Store(id=create_id('agg_data'), data=json.dumps(initial_agg_data)),
        dcc.Store(id=create_id('minmax'), data=json.dumps(minmax)),

        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                num_rows_displayed_Plaintext,
                                dcc.Store(id=create_id('displayed_row_count'), data=json.dumps(total_row_count)),
                                dcc.Store(id=create_id('total_row_count'), data=json.dumps(total_row_count)),
                            ],
                            style={'display': 'inline-block', 'margin-right': '10px'}
                        ),
                        html.Div(
                            [
                                html.Button("Reset All Filters", id="reset_all_filters_button"),
                            ],
                            style={'display': 'inline-block'}
                        )
                    ],
                    style={'display': 'inline-block'},
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Plaintext(
                                    'Copy shareable link:',
                                    style={'font-size': '15px', 'color': '#8334eb', 'font-weight': 'bold', 'display': 'inline-block', 'margin-right': '10px'}
                                ),
                                dcc.Loading(
                                    [
                                        dcc.Clipboard(
                                            id=create_id("copy_link_clipboard"), content='Error occurred during copying, try once more', className='button',
                                            style={'color': '#ffffff', 'background-color': '#8334eb', 'border': '1px solid #8334eb', 'margin-right': '10px'}
                                        ),
                                        html.Button(
                                            "Enumerate Checked Rows",
                                            id=create_id("enumerate_button"),
                                            style={'display': 'inline-block', 'margin-right': '10px', 'backgroundColor': '#eb8334', 'color': '#ffffff', 'border-style': 'none'}),
                                        dcc.Download(id=create_id("download_enumerations")), dcc.Store(id=create_id("selected_row_data"), data=['{}']),
                                        html.Button(
                                            "Download Raw Data",
                                            id=create_id("download_button"),
                                            style={'display': 'inline-block'}
                                        ),
                                        dcc.Download(id=create_id("download-dataframe-csv")),
                                    ],
                                    type='circle',
                                    style={'display': 'inline-block'},
                                    parent_style={'display': 'inline-block'},
                                ),
                            ],
                            style={'display': 'inline-block'}
                        )
                    ],
                    style={'display': 'inline-block'},
                ),
            ],
            style={'display': 'flex', 'justify-content': 'space-between', 'margin': '0', 'padding': '0'}
        ),

        html.Div(
            [
                html.Plaintext('Select rows:', style={'font-size': '15px', 'font-weight': 'bold', 'display': 'inline-block', 'margin-right': '10px'}),
                html.Button("Select All", id=create_id("select_all_button"), style={'display': 'inline-block'}),
                html.Button("Deselect All", id=create_id("deselect_all_button"), style={'display': 'inline-block'})
            ],
        ),

        html.Div(
            id=create_id('range_filters_div'),
            children=[
                html.Div([
                    html.Div(
                        [html.Plaintext('Filter rows:', style={'font-size': '15px', 'font-weight': 'bold'})],
                        style={'display': 'inline-block', 'margin-right': '10px'}
                    ),
                    html.Div(
                        [html.Button("Filter selected", id=create_id("filter_rows_button"))],
                        style={'display': 'inline-block'}
                    ),
                    html.Div(
                        [html.Button("Reset", id=create_id("reset_rows_filter_button"))],
                        style={'display': 'inline-block'}
                    ),
                ], style={'display': 'inline-block'}),

                html.Div([
                    html.Div(
                        [html.Plaintext(
                            'Range Filters:',
                            style={'display': 'inline-block', 'vertical-align': 'middle', 'font-size': '15px', 'font-weight': 'bold', 'margin-right': '10px'})],
                        style={'width': '15%', 'display': 'inline-block'}
                    ),
                    html.Div(
                        [range_filters_dropdown],
                        style={'width': '30%', 'display': 'inline-block', 'vertical-align': 'middle', 'margin-right': '10px'}
                    ),
                ] + [
                    html.Div(
                        id={'type': create_id('range_filter_div'), 'index': range_filters_name},
                        style={'width': '20%', 'display': 'inline-block', 'margin-right': '10px'} if range_filters_name == range_filters[0][0] else {'display': 'none'},
                        children=[range_filters_min, html.Plaintext(' - ', style={'width': '10%', 'font-weight': 'bold', 'display': 'inline-block'}), range_filters_max]
                    ) for (range_filters_name, range_filters_min, range_filters_max) in range_filters
                ] + [
                    html.Div(
                        [html.Button("Apply", id=create_id("apply_filters_button"))],
                        style={'display': 'inline-block'}
                    ),
                    html.Div(
                        [html.Button("Reset All", id=create_id("reset_filters_button"))],
                        style={'display': 'inline-block'}
                    ),
                ], style={'display': 'inline-block'})
            ],
            style={'display': 'flex', 'justify-content': 'space-between', 'margin': '0', 'padding': '0'}
        ),

        dcc.Store(id=create_id('row_selection_filter'), data=[]),
        dcc.Store(id=create_id('original_row_mappings'), data={}),
        dcc.Store(id=create_id('range_filter_names'), data=range_filter_names),

        html.Div(
            [transform_table],
        ),

        html.Br(), html.Br(),

        html.Div(
            [
                html.Div(
                    [

                        html.Div(
                            [
                                html.Div(
                                    [html.Plaintext(
                                        'X-Axis',
                                        style={'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)', 'font-size': '15px', 'font-weight': 'bold'}
                                    )],
                                    style={'width': '8%', 'display': 'inline-block'}
                                ),
                                html.Div(
                                    [x_Dropdown],
                                    style={'width': '24%', 'display': 'inline-block'}
                                ),
                                html.Div(
                                    [x_RadioItems],
                                    style={'width': '16%', 'display': 'inline-block', 'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)'}
                                ),
                                html.Div(
                                    [],
                                    style={'width': '2%', 'display': 'inline-block'}
                                ),
                                html.Div(
                                    [html.Plaintext(
                                        'Y-Axis',
                                        style={'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)', 'font-size': '15px', 'font-weight': 'bold'}
                                    )],
                                    style={'width': '8%', 'display': 'inline-block'}
                                ),
                                html.Div(
                                    [y_Dropdown],
                                    style={'width': '24%', 'display': 'inline-block'}
                                ),
                                html.Div(
                                    [y_RadioItems],
                                    style={'width': '16%', 'display': 'inline-block', 'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)'}
                                )
                            ],
                        ),

                        html.Div(
                            [main_scatterplot],
                            style={'width': '98%', 'display': 'inline-block'}
                        )
                    ],
                    style={'width': '49%', 'display': 'inline-block', 'vertical-align': 'top'},
                ),

                html.Div(
                    [
                        html.Div(
                            [html.Plaintext(
                                'To see matched pairs, left click on plot and drag-select points',
                                style={'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)', 'font-size': '15px', 'font-weight': 'bold'}
                            )],
                            style={'width': '8%', 'display': 'inline-block'}
                        ),
                        html.Div(
                            [html.Button('Clear', id=create_id('clearButton'), n_clicks=0)]
                        ),

                        html.Div(
                            [],
                            id=create_id('pair_tables')
                        )

                    ],
                    style={'width': '49%', 'display': 'inline-block', 'vertical-align': 'top'}
                )
            ],
            style={'display': 'flex', 'justify-content': 'start', 'margin': '0', 'padding': '0'}
        ),
    ]
    return [html.Div(children=children)]


def aggregate_statistics_by_rule(x_name, y_name, min_rfs, max_rfs, rf_names, pair_data, query_data, property_metadata, display_name_to_property_name):
    pair_data = json.loads(pair_data)
    query_id = pair_data['query_id']
    rf_names = json.loads(rf_names)

    # The xaxis and yaxis dropdown menu properties were generated using the get_prop_labels function. We can use this function again,
    # to generate a dictionary which lets us lookup which properties are selected
    snapquery_dict = json.loads(query_data)
    all_props = snapquery_dict["REQUIRED_properties"].split(",") + snapquery_dict["OPTIONAL_properties"].split(",")
    aggregation_type = snapquery_dict['advanced_options']['aggregation_type']
    prop_dicts = get_prop_labels(all_props, property_metadata=property_metadata)
    prop_dicts = {prop_dict['prop']: prop_dict for prop_dict in prop_dicts}
    selected_props = set()
    for prop in prop_dicts:
        if x_name in prop_dicts[prop].values() or y_name in prop_dicts[prop].values():
            selected_props.add(prop)

    rf_args = get_range_filter_args(rf_names, min_rfs, max_rfs, display_name_to_property_name)

    agg_args = {
        'query_id': query_id,
        'aggregation_type': aggregation_type,
        'range_filters': rf_args,
        'statistics': [
            {
                'statistic': 'median',
                'property_name': prop,
                'change_type': prop_dicts[prop]['change_type'],
                'base': prop_dicts[prop]['base'],
                'units': prop_dicts[prop]['units'],
            } for prop in selected_props
        ]
    }

    schema = snapquery_dict['schema']
    table_data = requests.post(backend_root + f'/aggregate_transforms?schema={schema}', data=json.dumps(agg_args))
    table_data = table_data.json()
    minmax = table_data['minmax']
    # The minmax columns were named within postgres, which restricts certain characters
    # We need the minmax column names to match the names in the x/y dropdowns, therefore we need to rename some columns
    renamed_minmax_headers = rename_column_headers(list(minmax.keys()), minmax=True, property_metadata=property_metadata)
    minmax = dict(zip(renamed_minmax_headers, minmax.values()))
    # Detect whether the environment was included in the aggregation, which is only needed for multicut queries
    environment_included = False
    for column in table_data['column_headers']:
        if 'smiles_env' in column:
            environment_included = True
            break
    columns = rename_column_headers(table_data['column_headers'], property_metadata=property_metadata)
    table_data = table_data['rows']
    df = pd.DataFrame(table_data, columns=columns)
    new_table_styles = (get_styles(df))

    if environment_included:
        if aggregation_type == 'individual_transforms':
            df['id'] = df['FROM'].apply(lambda x: str(x) + '_') + df['TO']
        elif aggregation_type == 'group_by_fragment':
            df['id'] = df['TO'].copy()
    else:
        # By default for individual_transforms with environment_included == False, the API call will return a column with 'rule_id' header which gets changed to 'id'
        # header by rename_column_headers
        if aggregation_type == 'group_by_fragment':
            df['id'] = df['TO'].copy()

    colors_mapped_to_ids = map_colors_to_ids(aggregation_type, df, environment_included)
    total_row_count = df['id'].nunique()

    if aggregation_type == 'individual_transforms':
        df = get_individual_transforms_df(df)
    elif aggregation_type == 'group_by_fragment':
        df = get_group_by_fragment_df(df)

    # Dash seems to only wrap text at spaces or - marks
    # But property column names can be long with multiple _, and Dash is setting widths of these columns to be narrow sometimes
    new_table_columns = [dict(name=i.replace('_', ' '), id=i, type='text', presentation='markdown') for i in df.columns if i not in ['id', 'rule_id_array']]
    data = df.to_dict('records')

    return json.dumps(data), new_table_columns, new_table_styles, json.dumps(minmax), json.dumps(total_row_count), [json.dumps(colors_mapped_to_ids)]


def selected_point_to_pair_tables(clickData, selectedData, children, n_clicks, pair_data, x_name, y_name, query_data, property_metadata):
    query_id = json.loads(pair_data)['query_id']

    snapquery_dict = json.loads(query_data)
    schema = snapquery_dict['schema']
    all_props = snapquery_dict["REQUIRED_properties"].split(",") + snapquery_dict["OPTIONAL_properties"].split(",")
    prop_dicts = get_prop_labels(all_props, property_metadata=property_metadata)
    prop_dicts = {prop_dict['prop']: prop_dict for prop_dict in prop_dicts}
    selected_props = set()
    for prop in prop_dicts:
        if x_name in prop_dicts[prop].values() or y_name in prop_dicts[prop].values():
            selected_props.add(prop)

    data = {
        'prop_changes': [
            {
                'property_name': prop,
                'change_type': prop_dicts[prop]['change_type'],
                'base': prop_dicts[prop]['base'],
                'units': prop_dicts[prop]['units']
            } for prop in selected_props
        ],
    }

    if n_clicks > 0:
        return [], 0, None
    elif selectedData:
        pairs = [point['customdata'] for point in selectedData['points']]
        # Limit number of displayed structures to maximum 50 for now, for performance
        if len(pairs) > 50:
            pairs = pairs[0:50]

        data['construct_id_pairs'] = pairs
        pair_data = requests.post(backend_root + f'/get_pair_data?schema={schema}', data=json.dumps(data))
        pair_data = pair_data.json()
        columns = rename_column_headers(pair_data['column_headers'], property_metadata=property_metadata)
        df = pd.DataFrame(pair_data['rows'], columns=columns)
        for column in df.columns:
            if column not in ['a', 'a_id', 'b', 'b_id', 'constant']:
                df[column] = df[column].apply(lambda x: df_round(x))

        multi_pair_table = generate_pair_tables(df.to_dict('records'))
        return multi_pair_table, 0, None
    elif clickData:
        single_pair = clickData['points'][0]['customdata']

        data['construct_id_pairs'] = [single_pair]
        pair_data = requests.post(backend_root + f'/get_pair_data?schema={schema}', data=json.dumps(data))
        pair_data = pair_data.json()
        columns = rename_column_headers(pair_data['column_headers'], property_metadata=property_metadata)
        df = pd.DataFrame(pair_data['rows'], columns=columns)
        for column in df.columns:
            if column not in ['a', 'a_id', 'b', 'b_id', 'constant']:
                df[column] = df[column].apply(lambda x: df_round(x))

        single_pair_table = generate_pair_tables(df.to_dict('records'))
        return single_pair_table + children, 0, None
    else:
        return children, 0, None


def map_axis_labels_to_API_args(label, property_metadata={}):
    if label[0:2] in ('A_', 'B_'):
        change_type = label[0]
        display_name = label[2:]
    elif label[-12:] == '_fold-change':
        change_type = 'fold_change'
        display_name = label[0:-12]
    elif label[-6:] == '_delta':
        change_type = 'delta'
        display_name = label[0:-6]
    # Could be 'Average ' or 'Average(log) '
    elif label[0:7] == 'Average':
        change_type = 'average'
        display_name = label.split(' ')[-1]

    for prop in property_metadata:
        if property_metadata[prop]['display_name'] == display_name:
            property_name = prop
            break

    axis_args = {'change_type': change_type, 'property_name': property_name}

    if property_metadata[property_name].get('display_unit') is not None:
        axis_args['units'] = property_metadata[property_name]['display_unit']

    return axis_args


def update_graph(
    colors_mapped_to_ids, active_cell, derived_viewport_row_ids, selected_row_ids,
    xaxis_radio_type, yaxis_radio_type, rf_names, min_rfs, max_rfs, minmax, pair_data, start_highlight_first_row, finish_highlight_first_row,
    columns, input_selected_cells, query_data, xaxis_column_name, yaxis_column_name, property_metadata, display_name_to_property_name
):
    pair_data = json.loads(pair_data)
    query_id = pair_data['query_id']
    rf_names = json.loads(rf_names)
    colors_mapped_to_ids = json.loads(colors_mapped_to_ids[0])

    triggered = dash.callback_context.triggered
    if triggered:
        for trig in triggered:
            # In a separate callback where we change the number of rows in the table, we want to activate crossfiltering based on the first row,
            # out of the newly displayed set of rows
            # Because a given id/value pair can only be the output of one callback, and we don't want to merge this callback with that other one,
            # we need two separate flags to accomplish this behavior, which are below:
            # We increment start_highlight_first_row in the first callback to trigger this behavior,
            # then we increment finish_highlight_first_row below, to turn off this behavior
            #
            if trig['prop_id'] == "table_transforms.derived_viewport_row_ids" and finish_highlight_first_row < start_highlight_first_row:
                if len(derived_viewport_row_ids) > 0:
                    active_cell = {'row': 0, 'column': 0, 'column_id': columns[0]['id'], 'row_id': derived_viewport_row_ids[0]}
                # Deactivate this trigger, until we activate it again in the separate callback by incrementing start_highlight_first_row
                finish_highlight_first_row = start_highlight_first_row

    output_selected_cells = []
    # When user clicks on a row, highlight all cells in that row, and get the active_row_id which we will use to crossfilter to scatterplot
    if active_cell is not None:
        active_row_id = active_cell['row_id']
        for i in range(len(columns)):
            cell = {
                'row': active_cell['row'],
                'column': i,
                'column_id': columns[i]['name'],
                'row_id': active_cell['row_id']
            }
            output_selected_cells.append(cell)
    # If we didn't click on any row to trigger this callback, then preserve selected cells, and associated crossfiltered points in scatterplot
    elif input_selected_cells != []:
        output_selected_cells = input_selected_cells
        active_row_id = input_selected_cells[0]['row_id']

    else:
        active_row_id = None

    # The table's row_ids will either be a single integer per row (if aggregation_type == 'individual_transforms'),
    # or a json string of list of integers per row (if aggregation_type == 'group_by_fragment')
    # Either way, we need to give the API call the same argument: a single list of all rule_ids being looked up
    snapquery_dict = json.loads(query_data)
    aggregation_type = snapquery_dict['advanced_options']['aggregation_type']

    grouped_by_environment = False
    query_row_ids = []
    collected_row_ids = []
    if selected_row_ids:
        collected_row_ids += selected_row_ids
    if active_row_id:
        collected_row_ids += [active_row_id]
    if aggregation_type == 'individual_transforms':
        # row_id is either an integer (rule_id) if we did not include env in results,
        # or (for multicut queries with environment where we save the environments) a string with (from_smiles_env)_(to_smiles_env)
        query_row_ids = collected_row_ids
        if query_row_ids:
            grouped_by_environment = True if '_' in str(query_row_ids[0]) else False
    elif aggregation_type == 'group_by_fragment':
        if collected_row_ids:
            if 'rule_id_array' in colors_mapped_to_ids[collected_row_ids[0]]:
                # Only evaluates true for 1-cut many-to-many queries
                for row_id in collected_row_ids:
                    query_row_ids += colors_mapped_to_ids[row_id]['rule_id_array']
            else:
                # This means we did include environment in results (should only happen for multicut queries)
                grouped_by_environment = True
                query_row_ids = collected_row_ids

    if len(query_row_ids) == 0:
        # The callback was triggered without intent to crossfilter, therefore abort the callback
        # This can happen because various actions, such as sorting the table, will trigger derived_viewport_row_ids or active_cell
        # In this case, we clear the scatterplot, so the user is not confused by a scatterplot remaining that is not tied to any crossfiltered row(s) in the graph
        fig = go.Figure(
            data=[],
            layout=go.Layout(
                margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                hovermode='closest',
                dragmode='lasso',
                plot_bgcolor='#e8effa',
                showlegend=False,
            )
        )
    else:
        rf_args = get_range_filter_args(rf_names, min_rfs, max_rfs, display_name_to_property_name)

        plot_args = {
            'query_id': query_id,
            'aggregation_type': aggregation_type,
            # 'grouped_by_environment': True if 'environment_smarts' in (column['name'] for column in columns) else False,
            'grouped_by_environment': grouped_by_environment,
            'ids': query_row_ids,
            'range_filters': rf_args,
            'x_data': map_axis_labels_to_API_args(xaxis_column_name, property_metadata=property_metadata),
            'y_data': map_axis_labels_to_API_args(yaxis_column_name, property_metadata=property_metadata)
        }

        schema = snapquery_dict['schema']
        plot_data = requests.post(backend_root + f'/get_plot_data?schema={schema}', data=json.dumps(plot_args))
        plot_data = plot_data.json()

        columns = rename_column_headers(plot_data['column_headers'], property_metadata=property_metadata)
        columns = columns[0:-2] + ['x', 'y']
        plot_data = plot_data['rows']
        df = pd.DataFrame(plot_data, columns=columns)

        color_default = '#000000'
        colors, opacities, sizes, border_widths, border_colors = [], [], [], [], []
        selected_id_set = set(selected_row_ids or [])

        if aggregation_type == 'individual_transforms':
            id_column = 'id'
            if grouped_by_environment is True:
                df['id'] = df['FROM'].apply(lambda x: x + '_') + df['TO']
        elif aggregation_type == 'group_by_fragment':
            id_column = 'TO'

        for row_id in df[id_column]:
            if row_id == active_row_id:
                colors.append(colors_mapped_to_ids[str(row_id)]['color'])
                opacities.append(1)
                sizes.append(8)
                border_widths.append(2)
            elif row_id in selected_id_set:
                colors.append(colors_mapped_to_ids[str(row_id)]['color'])
                opacities.append(0.7)
                sizes.append(6)
                border_widths.append(0)
            else:
                colors.append(color_default)
                opacities.append(0.8)
                sizes.append(5)
                border_widths.append(0)
        df['colors'] = colors
        df['opacities'] = opacities
        df['sizes'] = sizes
        df['border_widths'] = border_widths
        df['border_colors'] = ['#000000'] * df['border_widths'].size

        rows_selected = df['colors'] != color_default
        df_selected = df

        #  For enumerating building blocks, and/or targets
        # Get list of unique fragments in 'TO' column, which we can take into a jupyter notebook for now

        fig = go.Figure(
            data=[
                go.Scattergl(
                    x=df_selected['x'],
                    y=df_selected['y'],
                    customdata=list(df_selected[['from_construct_id', 'to_construct_id']].to_records(index=False)),
                    mode='markers',
                    marker={
                        'size': df_selected['sizes'],
                        'color': df_selected['colors'],
                        'opacity': df_selected['opacities'],
                        'line': {'width': df_selected['border_widths'], 'color': df_selected['border_colors']},
                    }
                )
            ],

            layout=go.Layout(
                margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                hovermode='closest',
                dragmode='lasso',
                plot_bgcolor='#e8effa',
                showlegend=False,
            )
        )

    # HANDLE PLOT AXIS SCALING
    plot_x_args = {
        'title': xaxis_column_name,
    }

    plot_y_args = {
        'title': yaxis_column_name,
    }

    xaxis_log_disabled = False
    yaxis_log_disabled = False

    # We store min and max field values in this dictionary, minmax, which was created during statistics aggregation
    minmax = json.loads(minmax)
    # Make minmax names consistent with UI names
    minmax = {key.replace('average_log_', 'average(log) ').replace('average_', 'average '): value for key, value in minmax.items()}

    # In these cases, we force linear scale (can't take log of negative numbers).
    if minmax['max_' + xaxis_column_name.lower()] <= 0 or minmax['min_' + xaxis_column_name.lower()] <= 0:
        xaxis_log_disabled = True
        xaxis_radio_type = 'Linear'
    if minmax['max_' + yaxis_column_name.lower()] <= 0 or minmax['min_' + yaxis_column_name.lower()] <= 0:
        yaxis_log_disabled = True
        yaxis_radio_type = 'Linear'

    xaxis_linearLog_radioRuttons_options = [
        {'label': 'Linear', 'value': 'Linear', 'disabled': False},
        {'label': 'Log', 'value': 'Log', 'disabled': xaxis_log_disabled}
    ]
    yaxis_linearLog_radioRuttons_options = [
        {'label': 'Linear', 'value': 'Linear', 'disabled': False},
        {'label': 'Log', 'value': 'Log', 'disabled': yaxis_log_disabled}
    ]

    if xaxis_radio_type == 'Linear':
        plot_x_args['type'] = 'linear'
        if 'delta' in xaxis_column_name:
            x_max = 1.05 * max(abs(minmax['max_' + xaxis_column_name.lower()]), abs(minmax['min_' + xaxis_column_name.lower()]))
            x_min = None
            vline_x_value = 0
        else:
            x_max = minmax['max_' + xaxis_column_name.lower()]
            x_max += 0.05 * (abs(x_max))
            x_min = minmax['min_' + xaxis_column_name.lower()]
            x_min -= 0.05 * (abs(x_max))
    elif xaxis_radio_type == 'Log':
        plot_x_args['type'] = 'log'
        plot_x_args['dtick'] = 'D2'
        if 'fold-change' in xaxis_column_name:
            x_max = 1.05 * max(abs(math.log(minmax['max_' + xaxis_column_name.lower()], 10)), abs(math.log(minmax['min_' + xaxis_column_name.lower()], 10)))
            x_min = None
            vline_x_value = 1
        else:
            x_max = math.log(minmax['max_' + xaxis_column_name.lower()], 10)
            x_max += 0.05 * (abs(x_max))
            x_min = math.log(minmax['min_' + xaxis_column_name.lower()], 10)
            x_min -= 0.05 * (abs(x_max))

    if yaxis_radio_type == 'Linear':
        plot_y_args['type'] = 'linear'
        if 'delta' in yaxis_column_name:
            y_max = 1.05 * max(abs(minmax['max_' + yaxis_column_name.lower()]), abs(minmax['min_' + yaxis_column_name.lower()]))
            y_min = None
            hline_y_value = 0
        else:
            y_max = minmax['max_' + yaxis_column_name.lower()]
            y_max += 0.05 * (abs(y_max))
            y_min = minmax['min_' + yaxis_column_name.lower()]
            y_min -= 0.05 * (abs(y_max))
    elif yaxis_radio_type == 'Log':
        plot_y_args['type'] = 'log'
        plot_y_args['dtick'] = 'D2'
        if 'fold-change' in yaxis_column_name:
            y_max = 1.05 * max(abs(math.log(minmax['max_' + yaxis_column_name.lower()], 10)), abs(math.log(minmax['min_' + yaxis_column_name.lower()], 10)))
            y_min = None
            hline_y_value = 1
        else:
            y_max = math.log(minmax['max_' + yaxis_column_name.lower()], 10)
            y_max += 0.05 * (abs(y_max))
            y_min = math.log(minmax['min_' + yaxis_column_name.lower()], 10)
            y_min -= 0.05 * (abs(y_max))

    if x_min is None:
        plot_x_args['range'] = [-1 * x_max, x_max]
    else:
        plot_x_args['range'] = [x_min, x_max]
    if y_min is None:
        plot_y_args['range'] = [-1 * y_max, y_max]
    else:
        plot_y_args['range'] = [y_min, y_max]

    fig.update_xaxes(**plot_x_args)
    fig.update_yaxes(**plot_y_args)
    if x_min is None:
        fig.add_vline(x=vline_x_value)
    if y_min is None:
        fig.add_hline(y=hline_y_value)

    # We only wanted the active cell to know what the user clicked on; once we select the entire row, we no longer need information about the active cell.
    # If we did leave the active cell activated, then it causes annoying behavior, where clicking a row checkbox away from the active cell causes the focus
    # to jump to the active cell
    active_cell = None

    return [
        active_cell, finish_highlight_first_row, output_selected_cells, fig, None, xaxis_radio_type, yaxis_radio_type, xaxis_linearLog_radioRuttons_options,
        yaxis_linearLog_radioRuttons_options
    ]
