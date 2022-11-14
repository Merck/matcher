import dash
from dash import Dash, dash_table, dcc, html
from dash.dependencies import Input, Output, State, MATCH, ALL
from frontend_api import app as Flask_app
import plotly.graph_objects as go
import pandas as pd
import urllib
import math
import json
import time
import os
import requests
import logging
logging.basicConfig(level=logging.INFO)
from datetime import datetime
import sys

backend_root = os.getenv("MATCHER_BACKEND_URL")
external_frontend_root = os.getenv("EXTERNAL_FRONTEND_URL")
external_backend_root = os.getenv("EXTERNAL_BACKEND_URL")

def df_round(value):
    if type(value) != type(5.5): return value

    if value > 100:
        value = int(value)
    elif value > 10:
        value = round(value, 1)
    elif value > 0.1:
        value = round(value, 2)
    elif value > 0.01:
        value = round(value, 3)
    elif value > 0.001:
        value = round(value, 4)
    elif value < 0.0001:
        value = round(value, 7)
    return value


def smiles_to_image_link(smiles, x, y, scaleImage='False'):
    return '![{}]({}/smilesToImage?smiles={}&x={}&y={}&scaleImage={})'.format(smiles, external_backend_root, urllib.parse.quote(smiles, safe=''), x, y, scaleImage)


def pair_to_aligned_image_links(A_smiles, B_smiles, pattern, x, y):
    return (
        '![{}]({}/smilesToImage_aligned?smiles={}&pattern={}&x={}&y={})'.format(A_smiles, external_backend_root, urllib.parse.quote(A_smiles, safe=''), urllib.parse.quote(pattern, safe=''), x, y),
        '![{}]({}/smilesToImage_aligned?smiles={}&pattern={}&x={}&y={})'.format(B_smiles, external_backend_root, urllib.parse.quote(B_smiles, safe=''), urllib.parse.quote(pattern, safe=''), x, y)
    )


def get_individual_transforms_df(df):

    # Image will get bigger if > 10 atoms in molecule, up to a limit of 2x the size defined below
    df['FROM'] = df['FROM'].apply(lambda x: smiles_to_image_link(x, 100, 50, scaleImage='True'))
    df['TO'] = df['TO'].apply(lambda x: smiles_to_image_link(x, 100, 50, scaleImage='True'))
    for column in df.columns:
        if '_fold-change' in column or '_delta' in column:
            df[column] = df[column].apply(lambda x: df_round(x))
    df.sort_values(by='pairs', inplace=True, ascending=False)

    return df


def get_individual_transforms_table(df, default_x_label, default_y_label):

    initial_data = df.to_dict('records')
    transform_table = dash_table.DataTable(
                id='table_transforms',
                # Dash seems to only wrap text at spaces or - marks
                # But property column names can be long with multiple _, and Dash is setting widths of these columns to be narrow sometimes
                columns=[dict(name=i.replace('_', ' '), id=i, type='text', presentation='markdown') for i in df.columns if i not in ['id']],
                data = initial_data,
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

def get_group_by_fragment_df(df):

    for column in df.columns:
        if '_fold-change' in column or '_delta' in column or r'% transforms' in column:
            df[column] = df[column].apply(lambda x: df_round(x))
    
    columns = df.columns.tolist()

    # Image will get bigger if > 10 atoms in molecule, up to a limit of 2x the size defined below
    df['TO'] = df['TO'].apply(lambda x: smiles_to_image_link(x, 100, 50, scaleImage='True'))
    df['FROM'] = df['FROM'].apply(lambda x: ' '.join(list(map(lambda y: smiles_to_image_link(y,100,50,scaleImage='True'), x))))

    stats_columns = []
    for column in columns:
        if 'median' in column or r'%' in column:
            stats_columns.append(column)

    # Control the order, from left to right, in which the columns will be displayed in the transform table
    # 'id' column will not be displayed, because it is specifically omitted during dash datatable generation
    new_df = df[['id', 'TO', 'transforms', 'pairs'] + stats_columns + ['FROM']]
    new_df.sort_values(by='transforms', inplace=True, ascending=False)
    return new_df

def get_group_by_fragment_table(df, default_x_label, default_y_label):

    initial_data_dict = df.to_dict('records')
    transform_table = dash_table.DataTable(
                id='table_transforms',
                # Dash seems to only wrap text at spaces or - marks
                # But property column names can be long with multiple _, and Dash is setting widths of these columns to be narrow sometimes
                columns=[dict(name=i.replace('_', ' '), id=i, type='text', presentation='markdown') for i in df.columns if i not in ['id', 'rule_id_array']],
                data = initial_data_dict,
                #fixed_columns={'headers': True, 'data': 4},
                fixed_rows={'headers': True},
                editable=False,
                filter_action="native",
                sort_action="native",
                sort_mode='multi',
                row_selectable='multi',
                row_deletable=False,
                selected_rows=[],
                page_action='native',
                page_current= 0,
                page_size= 20,
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
                    'minWidth': '100px'#, 'width': '180px', 'maxWidth': '180px',
                },
                active_cell={'row' : 0, 'column' : 0, 'row_id' : initial_data_dict[0]['id'], 'column_id' : 'TO'}
    )
    return transform_table

def get_styles(df, statistic = ''):
    styles = []
    for column in df.columns:
        prop_name = None
        if 'fold-change' in column:
            styles += data_bars_diverging(df, column)
            prop_name = column[:-12]
        elif 'delta' in column:
            prop_name = column[:-6]
            if 'EPSA' in column:
                styles += data_bars_diverging(df, column, log_or_linear='linear', col_min=-30, col_max=30, midpoint=0)
            elif 'PXR' in column:
                styles += data_bars_diverging(df, column, log_or_linear='linear', col_min=-50, col_max=50, midpoint=0)
            elif 'logD' in column:
                styles += data_bars_diverging(df, column, log_or_linear='linear', col_min=-3, col_max=3, midpoint=0)
            else:
                styles += data_bars_diverging(df, column, log_or_linear='linear', col_min=-10, col_max=10, midpoint=0)
        elif r'% transforms' in column:
            styles += data_bars_diverging(df, column, log_or_linear='linear', col_min = 0, col_max = 100, midpoint = 50)
        if prop_name is not None:
            styles += data_bars_diverging(df, r'% transforms (+' + prop_name + ')', log_or_linear='linear', col_min=0, col_max=100, midpoint=50)
            #styles += data_bars_diverging(df, r'% transforms (+' + prop_name + ')', log_or_linear='linear', col_min=0, col_max=100, midpoint=50, color_above='#0343ab', color_below='#ff9900')
    return styles

# Makes the statistics appear with red / green bars in the transform table
# Adapted from https://dash.plotly.com/datatable/conditional-formatting
# Made the ranges increase exponentially so that fold-change data is appropriately colored / scaled
def data_bars_diverging(df, column, log_or_linear='log', col_min=0.1, col_max=10, midpoint=1, color_above='#3D9970', color_below='#FF4136'):
    n_bins = 100
    bounds = [i * (1.0 / n_bins) for i in range(n_bins + 1)]
    col_max = col_max
    col_min = col_min
    # For linear scale, which we intend to apply to deltas, midpoint of 0 means nothing changes: B - A = 0
    # For log scale, which we intend to apply to fold-changes, midpoint of 1 means nothing changes: B / A = 1
    midpoint = midpoint
    ranges = [col_min]
    num = col_min
    if log_or_linear == 'log':
        # Exponential growth factor, such that compounding col_min by the factor 100 times will equal the col_max
        # If our midpoint is halfway (on the log scale) between col_min and col_max, then the above factor should give us 50 bins below the midpoint, and 50 bins above, spaced logarithmically
        factor = (col_max / col_min) ** (1/100)
        # For example, if we consider a 0.1x fold-change to be a full red bar, and a 10x change to be a full green bar, then we use the below settings:
        # For col_min = 0.1, col_max = 10, midpoint = 1: then the factor = 1.047128
        for i in range(1, 101):
            num = num * factor
            ranges.append(num)
    elif log_or_linear == 'linear':
        # For example: if we're looking at EPSA, and we consider +30 units (or -30 units) to be an off the chart change,
        # resulting in a full green bar, or full red bar, respectively, then the total_span is 60 units
        increment = (col_max - col_min) / n_bins
        for i in range(1, 101):
            num += increment
            ranges.append(num)
    else:
        raise ValueError("log_or_linear must be 'log' or 'linear'")

    len_bounds = len(bounds)
    styles = []
    for i in range(1, len_bounds + 1):
        if i < len_bounds:
            min_bound = ranges[i - 1]
            max_bound = ranges[i]
            min_bound_percentage = bounds[i - 1] * 100
            max_bound_percentage = bounds[i] * 100
            style = {
                'if': {
                    'filter_query': (
                        '{{{column}}} >= {min_bound}' +
                        (' && {{{column}}} < {max_bound}' if (i < len_bounds - 1) else '')
                    ).format(column=column, min_bound=min_bound, max_bound=max_bound),
                    'column_id': column
                },
                'paddingBottom': 2,
                'paddingTop': 2
            }
        # Create a single style for the case where the row value is less than the lower threshold, coloring the bar as -100% red
        elif i == len_bounds:
            min_bound, max_bound = ranges[0], ranges[0]
            min_bound_percentage, max_bound_percentage = 0, 0
            style = {
                'if': {
                    'filter_query': '{{{column}}} < {min_bound}'.format(column=column, min_bound=min_bound),
                    'column_id': column
                },
                'paddingBottom': 2,
                'paddingTop': 2
            }

        if max_bound > midpoint:
            background = (
                """
                    linear-gradient(90deg,
                    white 0%,
                    white 50%,
                    {color_above} 50%,
                    {color_above} {max_bound_percentage}%,
                    white {max_bound_percentage}%,
                    white 100%)
                """.format(
                    max_bound_percentage=max_bound_percentage,
                    color_above=color_above
                )
            )
        else:
            background = (
                """
                    linear-gradient(90deg,
                    white 0%,
                    white {min_bound_percentage}%,
                    {color_below} {min_bound_percentage}%,
                    {color_below} 50%,
                    white 50%,
                    white 100%)
                """.format(
                    min_bound_percentage=min_bound_percentage,
                    color_below=color_below
                )
            )
        style['background'] = background
        styles.append(style)

    return styles

# Takes in a row from the dataframe containing MMPs (1 pair per row)
# Turns that row into a Dash DataTable containing the structures, IDs, properties, and fold-changes for those properties
def generate_pair_tables(rows: list):
# keys:
# A, B, constant, A_ID, B_ID, A_prop, B_prop, prop_fold-change, prop_delta
    pair_tables = []
    for row in rows:

        try:
            A_image, B_image = pair_to_aligned_image_links(row['a'], row['b'], row['constant'], 300, 150)
        except Exception as exc:
            print(exc)
            A_image = smiles_to_image_link(row['a'], 300, 150)
            B_image = smiles_to_image_link(row['b'], 300, 150)

        _column = ['Structure', 'ID']
        A_column = [A_image, row['a_id']]
        B_column = [B_image, row['b_id']]
        __column = ['', 'Fold-changes']
        # For some properties we calculated deltas and not fold-changes
        # First look for fold-changes, and add them, while keeping track of seen delta props; then add deltas underneath
        delta_props = set()
        for key in row.keys():
            # Infer property names, e.g. if key = A_Papp, then key[2:] = Papp is a property
            if key[0:2] == 'a_' and key != 'a_id':
                if key[2:] + '_delta' in row.keys():
                    delta_props.add(key[2:])
                    continue
                _column.append(key[2:])
                A_column.append(row[key])
                B_column.append(row['b_' + key[2:]])
                __column.append(row[key[2:] + '_fold-change'])
        
        # If no fold-change props were found, change the column header, otherwise append a blank row with a new header for Deltas
        if len(__column) == 2:
            __column[1] = 'Deltas'
        elif len(delta_props) > 0:
            _column.append('')
            A_column.append('')
            B_column.append('')
            __column.append('Deltas')
            
        for prop in delta_props:
            _column.append(prop)
            A_column.append(row['a_' + prop])
            B_column.append(row['b_' + prop])
            __column.append(row[prop + '_delta'])

        pair_df = pd.DataFrame({
            '_' : _column,
            'A' : A_column,
            'B' : B_column,
            '__' : __column
        })

        pair_table = dash_table.DataTable(
                    id=str(len(pair_tables)),
                    columns=[dict(name=i, id=i, type='text', presentation='markdown') for i in ['_', 'A', 'B', '__']],
                    data = pair_df.dropna().to_dict('records'),
                    style_header = {'display': 'none'},
                    style_cell={
                        'whiteSpace': 'normal',
                        'height': 'auto',
                    },
                    row_selectable=False,
                    markdown_options={'link_target': '_blank', 'html': True},
                    style_table={'overflowX' : 'auto'},
                    css=[{
                            'selector': 'tr:first-child',
                            'rule': 'display: none',
                        }]
        )
        pair_tables.append(pair_table)

    return pair_tables

def base10_to_color_hex(number):
    hex_number = format(number, 'X')
    hex_number_for_color = '#' + '0'*(6 - len(hex_number)) + hex_number
    return hex_number_for_color

def get_prop_labels(props, property_metadata={}):
    prop_labels = []
    for prop in props:
        if not prop:
            continue

        # display_name is the name for the property that the user will see
        # In contrast, prop represents the name of the property in the database, property_name.name column
        display_name = property_metadata[prop]['display_name']
        default_change_axis_type = default_A_axis_type = 'Log'
        # database_base indicates whether the prop is stored as log, negative_log, or raw in the DB compound_property.value column
        database_base = property_metadata[prop]['base']
        average_label = 'Average(log) ' + display_name if database_base in ('log', 'negative_log') else 'Average ' + display_name

        frontend_base = property_metadata[prop]['display_base']
        units = property_metadata[prop].get('display_unit')
        change_type = property_metadata[prop]['change_displayed']
        change_label = display_name + '_' + change_type
        # Below is necessary to prevent SQL from attempting subtraction
        change_type = change_type.replace('-', '_')

        """
        # EPSA, logD_HPLC_pH7, and PXR should not be transformed from their stored values
        elif prop in ['EPSA', 'logD_HPLC_pH7', 'PXR']:
            average_label = 'Average ' + prop
            change_type = 'delta'
            change_label = prop + '_delta'
            default_change_axis_type = 'Linear'
            if prop in ['logD_HPLC_pH7', 'PXR']:
                default_A_axis_type = 'Linear'
        """

        labels = {
            'prop': prop,
            'A': 'A_' + display_name,
            'B': 'B_' + display_name,
            'average_label': average_label,
            'change_type': change_type,
            'change_label': change_label,
            'default_change_axis_type': default_change_axis_type,
            'default_A_axis_type': default_A_axis_type,
            'base': frontend_base,
            'units': units,
        }
        prop_labels.append(labels)
    return prop_labels

def rename_column_headers(headers, minmax=False, property_metadata={}):
    # We can't completely control case in the column headers retrieved from asyncpg via the API, nor use dashes (-) in normal postgres names, therefore we rename headers as desired here
    renamed = []
    rename_map = {'rule_id': 'id', 'from_smiles':'FROM', 'from_smiles_env':'FROM', 'from_smiles_array':'FROM', 'from_smiles_env_array':'FROM',
        'to_smiles':'TO', 'to_smiles_env':'TO', 'pair_count':'pairs', 'transform_count': 'transforms'}
    original_case_names = sorted(property_metadata.keys(), key=len, reverse=True)

    for header in headers:
        for name in original_case_names:
            if name.lower() in header:
                # Use the user-defined display_name rather than the name from the DB property_name.name column
                if minmax == False:
                    header = header.replace(name.lower(), property_metadata[name]['display_name'])
                elif minmax == True:
                    header = header.replace(name.lower(), property_metadata[name]['display_name'].lower())

                # Only use the longest name replacement, prevent further replacement by substrings of the longest name
                #  above, we sorted original_case_names to have the longest names first
                break
        if header in rename_map:
            renamed.append(rename_map[header])
        elif 'fold_change' in header:
            renamed.append(header.replace('fold_change', 'fold-change'))
        elif 'percent_increased_' in header:
            prop = header.split('_', 2)[2]
            renamed.append(f"% transforms (+{prop})")
        else:
            renamed.append(header)
    return renamed

def map_colors_to_ids(aggregation_type, df, environment_included):

    # Generate randomly distributed colors for plot
    # Convert hex number to integer: this is the largest hex color value we want to use, because any larger hex color values are often too light-colored
    six_hex_base10_limit = int('9aaaaa', 16)
    num_colors = df['id'].nunique()
    if num_colors == 0:
        return {}
    step = int(six_hex_base10_limit / num_colors)
    all_colors = set(range(0, six_hex_base10_limit, step))
    # Tuple of color hex codes, ('#000000', ... , '#XXXXXX')
    all_colors = list(map(base10_to_color_hex, all_colors))

    colors_mapped_to_ids = {}
    current_color = 0
    if aggregation_type == 'individual_transforms' and not environment_included:
        # ID will either be a Numpy int64 integer (rule_id), or a from_smiles_env + '_' + to_smiles_env string
        type_convert = lambda x: int(x)
    else:
        type_convert = lambda x: x

    for row_idx in range(len(df.index)):
        row_id = type_convert(df.iloc[row_idx]['id'])
        if row_id not in colors_mapped_to_ids:
            colors_mapped_to_ids[row_id] = {'color': all_colors[current_color]}
            current_color += 1
            if aggregation_type == 'group_by_fragment':
                # It also benefits us to store the mappings from TO_smiles to rule_id_array in a concise format (as opposed to having to load all data from the transform table), which will be used in scatterplot generation
                if not environment_included:
                    colors_mapped_to_ids[row_id]['rule_id_array'] = df.iloc[row_idx]['rule_id_array']
                else:
                    colors_mapped_to_ids[row_id]['from_env_array'] = df.iloc[row_idx]['FROM']
    return colors_mapped_to_ids

def get_range_filter_args(rf_names, min_rfs, max_rfs, display_name_to_property_name):
    rf_args = []
    for name, rf_min, rf_max in zip(rf_names, min_rfs, max_rfs):
        assert name[0:2] in ('A_', 'B_')
        compound, display_name_prop = name.split('_', 1)
        # e.g. converting display_name_prop='CYP3A4_IC50_uM' to value in DB column property_name.name='CYP3A4'
        prop = display_name_to_property_name[display_name_prop]
        if rf_min not in [None, 'None']:
            rf_args.append({'compound': compound, 'property_name': prop, 'operator': '>=', 'value': rf_min})
        if rf_max not in [None, 'None']:
            rf_args.append({'compound': compound, 'property_name': prop, 'operator': '<=', 'value': rf_max})
    return rf_args


def create_dash_app(dataframe = None, requests_pathname_prefix='/dash/', aggregation_type=None, snapquery_dict=None, snapfilter_dict=None, flask_server=None):

    if flask_server is not None:
        # USING routes_pathname_prefix IS KEY FOR THE DASH INTERFACE TO LOAD FROM FLASK WITH THIS SETUP
        app_dash = Dash(__name__, assets_folder="css/dash/" , routes_pathname_prefix=requests_pathname_prefix, prevent_initial_callbacks=True, server=flask_server, serve_locally=False)
    else:
        raise Exception("A Flask server must be provided as an argument")

    # The initially loaded layout will only show submit button.
    # Clicking on one of these two buttons will obtain data, then instantiate and initialize the rest of the elements
    # Subsequent clicking on one of these two buttons will overwrite all the previously instantiated elements
    def generate_layout():

        layout = html.Div([

            dcc.Store(id='resize_iframe_dummy_output', data=['{}']),
            dcc.Store(id='error_message_dummy_output', data=['{}']),
            dcc.Store(id='property_metadata', data={}),
            dcc.Store(id='display_name_to_property_name', data={}),

            # When the query is running, or data is loading, show a loading icon in place of the submit button, so user doesn't spam queries by accident
            dcc.Loading(
                [
                    dcc.Store(id='input_data', data=['{}']),
                    dcc.Store(id='query_data', data=['{}']),
                    dcc.Store(id='pair_data', data=['{}']),

                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Button("Submit query", id="submit_button", style={'backgroundColor': '#7a0dd9', 'color':'#ffffff', 'border-style': 'none'}),
                                ],
                                style={'display': 'inline-block'}
                            ),

                        ],
                        style={'display': 'flex', 'justify-content': 'center', 'margin': '0', 'padding': '0'}
                    ),
                ],
                type="circle"
            ),
            html.Br(), html.Br(), html.Br(),

            # Here is where we will place the rest of the Dash app elements, after the pair_data is loaded
            html.Div(id='output_div',
                children=[]
            ),
        ])
        return layout

    # pair_data flows in from submit button (via DB query)
    # Once we have the data, we can use the data to initialize all of the output elements
    @app_dash.callback(
        [Output('output_div', 'children')],
        [Input('pair_data', 'data')],
        [State('query_data', 'data'),
        State('property_metadata', 'data'),
        State('display_name_to_property_name', 'data'),]
    )
    def instantiate_output_elements(pair_data, query_data, property_metadata, display_name_to_property_name):

        pair_data = json.loads(pair_data)
        if "observations" in pair_data:
            # This means that something went wrong with the query, so we are returning an empty layout
            # We need to return a pair_tables div, because the iframe resizing callback depends on this div as an input
            return [html.Div(children=[
                            html.Div(
                                ['no_output'],
                                id = 'pair_tables',
                                style={'display': 'none'}
                            )])
            ]

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
                rf_names, min_rfs, max_rfs = [0]*len(range_filters_list), [0]*len(range_filters_list), [0]*len(range_filters_list)
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
            default_x_label, default_y_label, default_x_type, default_y_type = (snapfilter['default_x_label'], snapfilter['default_y_label'], snapfilter['default_x_type'], snapfilter['default_y_type'])
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
            # By default for individual_transforms with environment_included == False, the API call will return a column with 'rule_id' header which gets changed to 'id' header by rename_column_headers
            if aggregation_type == 'group_by_fragment':
                df['id'] = df['TO'].copy()

        colors_mapped_to_ids = map_colors_to_ids(aggregation_type, df, environment_included)

        # Approach with dropdowns, radiobuttons, and callback scatter plotting adapted from https://dash.plotly.com/interactive-graphing
        #### UI Component definitions

        x_y_dropdown_options = []
        for label in ['change_label', 'A', 'B', 'average_label']:
            for prop in prop_labels:
                x_y_dropdown_options.append({'label': prop[label], 'value': prop[label]})

        x_Dropdown = dcc.Dropdown(
            id='crossfilter-xaxis-column',
            options = x_y_dropdown_options,
            value=default_x_label
        )

        x_RadioItems = dcc.RadioItems(
            id='crossfilter-xaxis-type',
            options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
            value=default_x_type,
            labelStyle={'display': 'inline-block', 'marginTop': '5px'}
        )

        y_Dropdown = dcc.Dropdown(
            id='crossfilter-yaxis-column',
            options = x_y_dropdown_options,
            value=default_y_label
        )

        y_RadioItems = dcc.RadioItems(
            id='crossfilter-yaxis-type',
            options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
            value=default_y_type,
            labelStyle={'display': 'inline-block', 'marginTop': '5px'}
        )

        main_scatterplot = dcc.Graph(
            id='pairPlot'
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
                        id={'type': 'min_range_filter', 'index': column},
                        type='number',
                        debounce=True,
                        style={'width': '45%', 'display': 'inline-block'},
                        value=snapfilter['range_filters'][column]['min'] if column in snapfilter['range_filters'] else None,
                    ),
                    dcc.Input(
                        id={'type': 'max_range_filter', 'index': column},
                        type='number',
                        debounce=True,
                        style={'width':'45%', 'display': 'inline-block'},
                        value=snapfilter['range_filters'][column]['max'] if column in snapfilter['range_filters'] else None,
                    )
                )
            )
        range_filter_names = json.dumps(range_filter_options)

        range_filters_dropdown = dcc.Dropdown(
            id='range_filters_dropdown',
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
        transform_table = initialize_transform_table(df, default_x_label, default_y_label)

        initial_agg_data = df.to_dict('records')
        total_row_count = len(initial_agg_data)
    
        # Above we calculated num_colors based on having a unique color for each ID in the transform_table
        # This number will shrink if we are looking at any combination of properties other than the single most common property in the query
        initial_total_row_count = df['id'].nunique()
        num_rows_displayed_Plaintext = html.Plaintext(id='num_rows_displayed_plaintext', children='Displaying ' + str(initial_total_row_count) + ' out of ' + str(initial_total_row_count) + ' rows', style={'font-size': '15px', 'font-weight': 'bold'})

        children = [
            dcc.Store(id='snapfilter', data=json.dumps(snapfilter)),
            dcc.Store(id='colors_mapped_to_ids', data=[json.dumps(colors_mapped_to_ids)]),
            # Use in the callback that updates transform_table, to enforce initialization of row filtration, using snapfilter['transform_row_ids']
            #dcc.Store(id='snapfilter_applied', data= [{'snapfilter_applied': False if snapfilter_exists else True}]),
            dcc.Store(id='snapfilter_applied', data= False if snapfilter_exists else True),
            html.Button("not_intended_to_display", id="start_highlight_first_row", n_clicks=0, style={'display' : 'none'}),
            html.Button("not_intended_to_display", id="finish_highlight_first_row", n_clicks=0, style={'display' : 'none'}),
            dcc.Store(id='identifier', data='id'),
            dcc.Store(id='agg_data', data=json.dumps(initial_agg_data)),
            dcc.Store(id='minmax', data=json.dumps(minmax)),

            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    num_rows_displayed_Plaintext,
                                    dcc.Store(id='displayed_row_count', data=json.dumps(total_row_count)),
                                    dcc.Store(id='total_row_count', data=json.dumps(total_row_count)),
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
                                    html.Plaintext('Copy shareable link:', style={'font-size': '15px', 'color': '#8334eb', 'font-weight': 'bold', 'display': 'inline-block', 'margin-right': '10px'}),
                                    dcc.Loading(
                                        [
                                            dcc.Clipboard(
                                                id="copy_link_clipboard", content='Error occurred during copying, try once more', className='button',
                                                style={'color': '#ffffff', 'background-color': '#8334eb', 'border': '1px solid #8334eb', 'margin-right': '10px'}
                                            ),
                                            html.Button("Enumerate Checked Rows", id="enumerate_button", style={'display': 'inline-block', 'margin-right': '10px', 'backgroundColor': '#eb8334', 'color':'#ffffff', 'border-style':'none'}),
                                            dcc.Download(id="download_enumerations"), dcc.Store(id="selected_row_data", data=['{}']),
                                            html.Button("Download Raw Data", id="download_button", style={'display': 'inline-block'}), dcc.Download(id="download-dataframe-csv"),
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
                    html.Button("Select All", id="select_all_button", style={'display': 'inline-block'}),
                    html.Button("Deselect All", id="deselect_all_button", style={'display': 'inline-block'})
                ],
            ),

            html.Div(id='range_filters_div',
                children=[
                    html.Div([
                        html.Div(
                            [html.Plaintext('Filter rows:', style={'font-size': '15px', 'font-weight': 'bold'})],
                            style={'display': 'inline-block', 'margin-right': '10px'}
                        ),
                        html.Div(
                            [html.Button("Filter selected", id="filter_rows_button")],
                            style={'display': 'inline-block'}
                        ),
                        html.Div(
                            [html.Button("Reset", id="reset_rows_filter_button")],
                            style={'display': 'inline-block'}
                        ),
                    ], style={'display': 'inline-block'}),

                    html.Div([
                        html.Div(
                            [html.Plaintext('Range Filters:', style = {'display': 'inline-block', 'vertical-align': 'middle', 'font-size': '15px', 'font-weight': 'bold', 'margin-right': '10px'})],
                            style={'width': '15%', 'display': 'inline-block'}
                        ),
                        html.Div(
                            [range_filters_dropdown],
                            style={'width': '30%', 'display': 'inline-block', 'vertical-align': 'middle', 'margin-right' : '10px'}
                        ),
                    ] + [
                        html.Div(
                            id={'type': 'range_filter_div', 'index': range_filters_name},
                            style={'width': '20%', 'display': 'inline-block', 'margin-right': '10px'} if range_filters_name == range_filters[0][0] else {'display': 'none'},
                            children=[range_filters_min, html.Plaintext(' - ', style={'width': '10%', 'font-weight': 'bold', 'display': 'inline-block'}), range_filters_max]
                        ) for (range_filters_name, range_filters_min, range_filters_max) in range_filters
                    ] + [
                        html.Div(
                            [html.Button("Apply", id="apply_filters_button")],
                            style={'display': 'inline-block'}
                        ),
                        html.Div(
                            [html.Button("Reset All", id="reset_filters_button")],
                            style={'display': 'inline-block'}
                        ),
                    ], style={'display': 'inline-block'})
                ],
                style={'display': 'flex', 'justify-content': 'space-between', 'margin': '0', 'padding': '0'}
            ),

            dcc.Store(id='row_selection_filter', data=[]),
            dcc.Store(id='original_row_mappings', data={}),
            dcc.Store(id='range_filter_names', data=range_filter_names),

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
                                        [html.Plaintext('X-Axis', style = {'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)', 'font-size': '15px', 'font-weight': 'bold'})],
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
                                        [html.Plaintext('Y-Axis', style = {'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)', 'font-size': '15px', 'font-weight': 'bold'})],
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
                        style={'width': '49%', 'display': 'inline-block', 'vertical-align' : 'top'},
                    ),

                    html.Div(
                        [
                            html.Div(
                                    [html.Plaintext('To see matched pairs, left click on plot and drag-select points', style = {'position': 'relative', 'top': '50%', 'transform': 'translateY(-50%)', 'font-size': '15px', 'font-weight': 'bold'})],
                                    style={'width': '8%', 'display': 'inline-block'}
                            ), 
                            html.Div(
                                [html.Button('Clear', id='clearButton', n_clicks=0)]
                            ),

                            html.Div(
                                [],
                                id = 'pair_tables'
                            )

                        ], 
                        style={'width': '49%', 'display': 'inline-block', 'vertical-align' : 'top'}
                    )
                ],
                style={'display': 'flex', 'justify-content': 'start', 'margin': '0', 'padding': '0'}
            ),
        ]
        return [html.Div(children=children)]

    app_dash.layout = generate_layout

    # Collect user-entered data from matcher.html form
    app_dash.clientside_callback(
    """
    function (n_clicks) {

        // This callback will always fire on page load, in order to support automatic query firing of snapshot-based queries
        // However, if we loaded the page without any snapshot-based query, then the fields will be empty; in that case, don't confuse the user with an error message
        // This block disables itself after the initial page load, so subsequently the user will get error messages if they try to submit an incomplete form
        if (parent.suppress_initial_errors === true) {
            parent.suppress_initial_errors = false;

            // Only exit the callback if there was no snapshot query loaded in via template, in which case the variable will be an empty string
            if (parent.snapquery_id === "") {
                return {};
            }
        }

        //console.log('advanced options data');
        //parent.getAdvancedOptions();

        var highlights = parent.highlights;
        // const backend_root = parent.backend_root;

        // Get the sketched structure
        //const sketched_content = parent.document.getElementById('sketched_content').getAttribute('value');
        const sketched_content = parent.getSketchedContent();

        // Get radio button selections
        var query_type_radios = parent.document.getElementsByName('query_type');
        var query_type;
        for (var i=0; i < query_type_radios.length; i++) {
            if (query_type_radios[i].checked) {
                query_type = query_type_radios[i].value;
            }
        }

        var transform_order_radios = parent.document.getElementsByName('transform_order');
        var transform_order;
        for (var i=0; i < transform_order_radios.length; i++) {
            if (transform_order_radios[i].checked) {
                transform_order = transform_order_radios[i].value;
            }
        }

        // Handle cases where user does not provide all required information
        var missing = "<br/>Cannot submit:<br/><br/>";
        var missing_fields = missing;
        if ((parent.getProps('req') === "") && (parent.getProps('opt') === "")) {
            missing_fields += "<strong/>Please select at least 1 property</strong/><br/><br/>";
        }
        if (highlights.mol1.variable.atoms.length === 0 && highlights.mol2.variable.atoms.length === 0) {
            missing_fields += "<strong/>Please select atoms in the sketcher, and click the green 'Set Variable Atoms' button</strong/><br/><br/>";
        }
        if (missing_fields !== missing) {
            parent.document.getElementById("query_progress_div").innerHTML = missing_fields;
            parent.document.documentElement.style.setProperty('--query_progress_div_visibility', 'block');
            return {};
        }

        parent.document.documentElement.style.setProperty('--query_progress_div_visibility', 'none');

        data = {
            snapquery_id: parent.snapquery_id,
            query_type: query_type,
            query_id: parent.query_id,
            transform_order: transform_order,
            sketched_content: sketched_content,
            // variable_atoms: (JSON.stringify(highlights[0]) !== '{}') ? highlights[0].indexes.atoms.join(",") : '',
            // variable_bonds: (JSON.stringify(highlights[0]) !== '{}') ? highlights[0].indexes.bonds.join(",") : '',
            // environment_atoms: (JSON.stringify(highlights[1]) !== '{}') ? highlights[1].indexes.atoms.join(",") : '',
            // environment_bonds: (JSON.stringify(highlights[1]) !== '{}') ? highlights[1].indexes.bonds.join(",") : '',
            REQUIRED_properties: parent.getProps('req'),
            OPTIONAL_properties: parent.getProps('opt'),
            advanced_options: parent.getAdvancedOptions(),
            snapfilter_id: parent.snapfilter_id,
            snapfilter_string: parent.snapfilter_string,
            schema: parent.schema
        }

        output = JSON.stringify(data);

        // We need to clear the below variables from the input form,
        // so that when the user does a subsequent query, the UI knows to create new IDs for these form data
        parent.snapquery_id = "";
        parent.query_id = -1;
        parent.snapfilter_id = "";
        parent.snapfilter_string = "";

        return [output, parent.property_metadata, parent.display_name_to_property_name];
    }
    """,
    [Output('input_data', 'data'),
    Output('property_metadata', 'data'),
    Output('display_name_to_property_name', 'data'),],
    [Input('submit_button', 'n_clicks')],
    prevent_initial_call=False
    )

    # Use this function to run the query, if your production environment does not kill long requests in a way outside of your control
    # TODO instead, use Celery or similar service
    """
    @app_dash.callback(
        [Output('query_data', 'data'),
        Output('pair_data', 'data')],
        [Input('input_data', 'data')]
    )
    def run_query(input_data):

        pair_data = requests.post(backend_root + '/async_query_v3_backend', data = input_data)
        pair_data = pair_data.json()

        return [input_data, pair_data]
    """

    # Use this function to run the query, if your production environment kills long requests in a way that you cannot prevent
    # TODO instead, use Celery or similar service
    @app_dash.callback(
        [Output('query_data', 'data'),
        Output('pair_data', 'data')],
        [Input('input_data', 'data')]
    )
    def run_persistent_query(input_data):

        input_data_dict = json.loads(input_data)
        if input_data_dict.get('query_id') != -1:
            # This means we loaded a snapshot that already has saved results in the DB
            # Therefore, instead of running that query again and duplicating those results, we just skip the query and reference those previous results
            return [input_data, json.dumps({'query_id': input_data_dict['query_id']})]

        # Start the query and get an ID which we will use to poll the DB for existence of results over intervals
        schema = input_data_dict['schema']
        started_query = requests.post(backend_root + f'/start_query?schema={schema}', data = input_data)
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
                    ### The only case leading to display of result data, if we started a query
                    # Connect output data (via query_id, which is used in the DB to represent the unique query to which a results set belongs) to the input state, for use in snapshot generation
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

    # Display various error messages to the user
    app_dash.clientside_callback(
        """
        function(data) {
            //data = JSON.parse(data);
            let error_div = parent.document.getElementById("query_progress_div");

            // If something went wrong during the query, the data["observations"] field will have values
            if (data["observations"] === undefined) {

                // This means the query went well. Adjust the position of the window to look at the results
                parent.document.getElementById("jump_to_results_button").click();
            
                return {};
            }

            parent.document.documentElement.style.setProperty('--query_progress_div_visibility', 'block');

            if (data["observations"]["all_atoms_selected_as_variable"] === "True") {
                error_div.innerHTML = "<br/><strong/>Cannot submit:</strong/><br/><br/><strong/>The variable fragment can't be the entire molecule.</strong/><br/><br/>At least 1 rotatable bond must be attached to the variable fragment.<br/><br/>";

            } else if (data["observations"]["variable_selection_over_limit"] === "True") {
                let max_variable_atoms = data["observations"]["max_variable_atoms"].toString();
                //console.log(max_variable_atoms);
                error_div.innerHTML = "<br/><strong/>Cannot submit:</strong/><br/><br/>The variable fragment can't be larger than " + max_variable_atoms + " heavy atoms.<br/><br/>Select a variable fragment that is " + max_variable_atoms + " atoms or less.<br/>";

            } else if (data["observations"]["missing_constant"] === "True") {
                error_div.innerHTML = "<br/><strong/>Cannot submit:</strong/><br/><br/>If you select constant (environment) atoms on one side of the reaction arrow, you must also select constant atoms on the other side of the arrow.<br/>";
            } else if (data["observations"]["query_too_general"] === "no_variable_no_constant") {
                error_div.innerHTML = "<br/><strong/>Cannot submit:</strong/><br/><br/>Your input structure is too general, you're basically looking for everything.<br/><br/>This would return too much data for our system to support at this stage.<br/>";
            } else if (data["observations"]["query_too_general"] === "many_to_many_too_general") {
                error_div.innerHTML = "<br/><strong/>Cannot submit:</strong/><br/><br/>Your input structure is too general.<br/><br/>This would return too much data for our system to support at this stage.<br/>Please make the red or green selections larger.<br/>";
            } else if (data["observations"]["no_results"] === "True") {
                error_div.innerHTML = "<br/><strong/>Search complete</strong/><br/><br/>No results matching the search.<br/><br/>Try making the search more general.<br/>";
            } else if (data["observations"]["exception"] === "True") {
                error_div.innerHTML = "<br/><strong/>Error occurred</strong/><br/><br/>Please consider sending a screenshot of your input to the application owner so they can fix this bug.<br/>";
            } else if (data["observations"]["query_timeout"] === "True") {
                error_div.innerHTML = "<br/><strong/>The query timed out after 1 hour</strong/><br/><br/>Let the application owner know what you were looking for; it may be possible to increase the timeout limit.<br/>";
            } else if (data["error_message"] !== undefined) {
                error_div.innerHTML = "<br/><strong/>Error: Something went wrong, sorry about that!</strong/><br/><br/>Please ask the application owner to fix the problem, and include the below details along with what your search was:<br/><br/>" + String(data["observations"]) + "<br/>";
            }
            return {};
        }
        """,
        [Output('error_message_dummy_output', 'data')],
        [Input('pair_data', 'data')]
    )

    snapshot_states = [
        State('crossfilter-xaxis-column', 'value'),
        State('crossfilter-yaxis-column', 'value'),
        State('crossfilter-xaxis-type', 'value'),
        State('crossfilter-yaxis-type', 'value'),
        State('row_selection_filter', 'data'),
        State('range_filter_names', 'data'),
        State({'type': 'min_range_filter', 'index': ALL}, 'value'),
        State({'type': 'max_range_filter', 'index': ALL}, 'value'),
        State('query_data', 'data')
    ]

    def get_new_snapshot_id(inputs):
        snapquery_dict = json.loads(inputs[9])

        # Take snapshot of output UI state and store it in snapfilter_string
        row_selection_filter = inputs[5]
        if row_selection_filter:
            transform_row_ids = ','.join([str(x) for x in row_selection_filter])
        else:
            transform_row_ids = ''

        range_filter_names = json.loads(inputs[6])
        nonempty_range_filters = []

        for name, range_min, range_max in list(zip(range_filter_names, inputs[7], inputs[8])):
            if (range_min is not None) or (range_max is not None):
                range_min = str(range_min) if range_min is not None else 'None'
                range_max = str(range_max) if range_max is not None else 'None'
                nonempty_range_filters.append(','.join([name, range_min, range_max]))

        nonempty_range_filters = '@'.join(nonempty_range_filters)

        # Add two blank strings at end, as placeholders for table-level sorting and filtering strings, which can be implemented later
        snapfilter_string = ';'.join(list(inputs[1:5]) + [
            transform_row_ids, nonempty_range_filters, '', ''
        ])

        # Include the original query in the snapshot, but modify the snapfilter_string with the newly harvested filters
        new_snapshot = snapquery_dict.copy()
        new_snapshot['snapfilter_string'] = snapfilter_string
        # In theory we could track if the snapfilter_id changed after snap page load, and reuse the id if not
        new_snapshot['snapfilter_id'] = ''

        schema = snapquery_dict['schema']
        new_snapshot_id = requests.post(backend_root + f'/snap_write?schema={schema}', data=json.dumps(new_snapshot))
        
        # new_snapshot contains all the information needed to reproduce a query
        #logging.info(json.dumps(new_snapshot))
        # To write this query to DB and associate with the /snap endpoint, if needed outside of the normal application (e.g. during DB initialization),
        #   remove the 'query_id' key/value from new_snapshot and `await backend_api.snap_write(models.ExampleQuery.parse_obj(new_snapshot))`

        new_snapshot_id = new_snapshot_id.json()['new_snapshot_id']

        return new_snapshot_id

    @app_dash.callback(
        [Output('copy_link_clipboard', 'content')],
        [Input('copy_link_clipboard', 'n_clicks')],
        snapshot_states
    )
    def spawn_and_copy_snapshot_link(*inputs):

        new_snapshot_id = get_new_snapshot_id(inputs)
        schema = json.loads(inputs[-1])['schema']
        query_schema = f"?schema={schema}" if schema != "None" else ""
        content = external_frontend_root + '/snap/' + new_snapshot_id + query_schema

        # dcc.Clipboard expects a list for some reason, but will copy the single element in this list to the clipboard
        return [content]

    
    @app_dash.callback(
        [Output("download-dataframe-csv", "data")],
        [Input("download_button", "n_clicks")],
        [State('pair_data', 'data'),
        State('query_data', 'data')],
        prevent_initial_call=True,
    )
    def download_csv(n_clicks, pair_data, query_data):
        query_id = json.loads(pair_data)["query_id"]
        data = {
            'query_id': query_id,
            'query': json.loads(query_data)
        }
        schema = data['query']['schema']
        all_raw_data = requests.post(backend_root + f'/get_all_raw_data?schema={schema}', data=json.dumps(data))

        all_raw_data = all_raw_data.json()
        columns = all_raw_data['column_headers']
        all_raw_data = all_raw_data['rows']
        df = pd.DataFrame(all_raw_data, columns=columns)
        return [dcc.send_data_frame(df.to_csv, "matcher_dataset.csv")]
    
    @app_dash.callback(
        [Output({'type': 'range_filter_div', 'index': ALL}, 'style')],
        [Input('range_filters_dropdown', 'value')],
        [State('range_filter_names', 'data')]
    )
    def display_filters_from_dropdown(dropdown_value, range_filter_names):
        div_display_settings = []
        range_filter_names = json.loads(range_filter_names)
        # The below assignment relies on the order being preserved (since initialization) for both range_filter_names, and the list of range_filter_div outputs
        # Those orders must match
        for name in range_filter_names:
            if name == dropdown_value:
                div_display_settings.append(
                    {'width': '20%', 'display': 'inline-block', 'margin-right': '10px'}
                )
            else:
                div_display_settings.append(
                    {'display': 'none'}
                )
        return [div_display_settings]

    # Clear all range filters
    @app_dash.callback(
        [Output({'type': 'min_range_filter', 'index': ALL}, 'value'),
        Output({'type': 'max_range_filter', 'index': ALL}, 'value')],
        [Input('reset_filters_button', 'n_clicks')],
        [State('range_filter_names', 'data')]
    )
    def reset_filters(n_clicks, range_filter_names):
        range_filter_names = json.loads(range_filter_names)
        return [[None]*len(range_filter_names)]*2

    def select_all_rows(derived_virtual_data):
        if derived_virtual_data is None:
            selected_rows = []
            selected_row_ids = []
        else:
            selected_rows = [row_num for row_num in range(len(derived_virtual_data))]
            selected_row_ids = [row['id'] for row in derived_virtual_data]
        return selected_rows, selected_row_ids

    def deselect_all_rows():
        selected_rows = []
        selected_row_ids = []
        return selected_rows, selected_row_ids

    @app_dash.callback(
        [Output('num_rows_displayed_plaintext', 'children'),
        Output('num_rows_displayed_plaintext', 'style'),],
        [Input('displayed_row_count', 'data'),
        Input('total_row_count', 'data')],
    )
    def show_number_of_displayed_rows(displayed_row_count, total_row_count):

        displayed_row_count = json.loads(displayed_row_count)
        total_row_count = json.loads(total_row_count)
        if displayed_row_count != total_row_count:
            text_color = "#fc0303"
        else:
            text_color = "#000000"
        output_text = 'Displaying ' + str(displayed_row_count) + ' out of ' + str(total_row_count) + ' rows'
        output_style = {'font-size': '15px', 'font-weight': 'bold', 'color' : text_color}

        return output_text, output_style

    @app_dash.callback(
        [Output('reset_filters_button', 'n_clicks'),
        Output('reset_rows_filter_button', 'n_clicks')],
        Input('reset_all_filters_button', 'n_clicks'),
        [State('reset_filters_button', 'n_clicks'),
        State('reset_rows_filter_button', 'n_clicks')]
    )
    def reset_all_filters(n_clicks_A, n_clicks_B, n_clicks_C):
        # Simulate pressing of all of the reset filters buttons
        if n_clicks_B is None:
            n_clicks_B = 0
        if n_clicks_C is None:
            n_clicks_C = 0
        return n_clicks_B + 1, n_clicks_C + 1

    # Update aggregated data when user selects different properties or range filters
    # We output agg_data, which is the parent data after aggregation, but which can undergo subsequent filtration by other clientside callbacks, before being loaded to the table
    @app_dash.callback(
        [Output('agg_data', 'data'),
        Output('table_transforms', 'columns'),
        Output('table_transforms', 'style_data_conditional'),
        Output('minmax', 'data'),
        Output('total_row_count', 'data'),
        Output('colors_mapped_to_ids', 'data')],
        [Input('crossfilter-xaxis-column', 'value'),
        Input('crossfilter-yaxis-column', 'value'),
        Input({'type': 'min_range_filter', 'index': ALL}, 'value'),
        Input({'type': 'max_range_filter', 'index': ALL}, 'value')],
        [State('range_filter_names', 'data'),
        State('pair_data', 'data'),
        State('query_data', 'data'),
        State('property_metadata', 'data'),
        State('display_name_to_property_name', 'data'),]
    )
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
            'statistics' : [
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
            # By default for individual_transforms with environment_included == False, the API call will return a column with 'rule_id' header which gets changed to 'id' header by rename_column_headers
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

    """
    This callback controls which rows in the transform table are displayed vs. hidden
    This callback also controls which rows are checked, when we change which rows are displayed
    We are writing this callback to run clientside, because we don't want to pass the (potentially large) data from the table back and forth with server,
    and these are fairly trivial functions to write in javascript, requiring no special python libraries etc.

    agg_data (list of dictionaries) is the originally aggregated statistics data, before any clientside filtering of the rows.
    We never modify agg_data aside from completely replacing agg_data from scratch using an API endpoint.

    table_data (list of dictionaries) was initialized to be identical to agg_data. table_data can be further filtered in this callback.
    If we want to undo all filtering, we again set table_data = agg_data

    selected_rows (list of int) are the integer indices (starting at 0) of rows that have a checked box next to them, and are a subset of all displayed rows.
    selected_rows lets us keep track of the *positions in the UI* of the selected rows

    selected_row_ids (list of int, or list of str): each id is unique to data in a row, e.g. the rule_id for a MMP tranformation
    selected_row_ids lets us keep track of the *content* of the selected rows
    """
    app_dash.clientside_callback(
        """
        function filter_rows (
            a,b,c,d, agg_data, snapfilter_applied, snapfilter, table_data, derived_virtual_data, 
            selected_row_ids, selected_rows, original_row_mappings, row_selection_filter, 
            page_current, start_highlight_first_row
        ) {

            agg_data = JSON.parse(agg_data);

            // Preserve the 1:1 relationship between the unique index of the checked box, and the unique ID of the data in the row, in the unfiltered data
            // Then we can return to this state if the user wants to undo the hiding of rows
            function set_original_row_mappings (original_row_mappings, selected_row_ids, selected_rows) {
                if (Object.keys(original_row_mappings).length === 0) {
                    if (selected_row_ids.length > 0) {
                        for (let a = 0; a < selected_row_ids.length; a++) {
                            original_row_mappings[selected_row_ids[a]] = selected_rows[a]
                        }
                    }
                }
                return original_row_mappings;
            }

            // Set the checked box indexes to [0 ... k] where k is the number of rows after filtering
            function select_all_filtered_rows (row_selection_filter) {
                if (row_selection_filter.length > 0) {
                    let selected_rows = [];
                    for (let k = 0; k < row_selection_filter.length; k++) {
                        selected_rows.push(k);
                    }
                    return selected_rows;
                } else {
                    return undefined;
                }
            }

            let triggered = dash_clientside.callback_context.triggered;
            snapfilter = JSON.parse(snapfilter);

            if (snapfilter_applied == false) {
                // This is only meant to fire on page load, when loaded via a saved state (snapshot / snapquery / snapfilter)
                selected_row_ids = snapfilter['transform_row_ids'];
                if (selected_row_ids.length == 0) {
                    snapfilter_applied = true; 
                }
                else {
                    row_selection_filter = selected_row_ids;
                    var filter_set = new Set(row_selection_filter);
                    selected_rows = [];
                    for (let m = 0; m < agg_data.length; m++) {
                        if (filter_set.has(agg_data[m]['id'])) {
                            selected_rows.push(m);
                        }
                    }
                    original_row_mappings = set_original_row_mappings(original_row_mappings, selected_row_ids, selected_rows);
                    page_current = 0;
                    selected_rows = select_all_filtered_rows(row_selection_filter);
                    snapfilter_applied = true;
                    start_highlight_first_row += 1;
                }
            }
            else {
                for (let i = 0; i < triggered.length; i++) {
                    if (triggered[i]['prop_id'] == 'select_all_button.n_clicks') {
                        selected_rows = [];
                        selected_row_ids = [];
                        if (Object.keys(derived_virtual_data).length > 0) {
                            // derived_virtual_data is all data after clientside filtering
                            for (let j = 0; j < derived_virtual_data.length; j++) {
                                selected_rows.push(j);
                                selected_row_ids.push(derived_virtual_data[j]['id']);
                            }
                        }
                    }
                    else if (triggered[i]['prop_id'] == 'deselect_all_button.n_clicks') {
                        selected_rows = [];
                        selected_row_ids = [];
                    }
                    else if (triggered[i]['prop_id'] == 'filter_rows_button.n_clicks') {
                        row_selection_filter = selected_row_ids;
                        var filter_set = new Set(row_selection_filter);
                        original_row_mappings = set_original_row_mappings(original_row_mappings, selected_row_ids, selected_rows);
                        page_current = 0;
                        selected_rows = select_all_filtered_rows(row_selection_filter);
                        start_highlight_first_row += 1;
                    }
                    else if (triggered[i]['prop_id'] == 'reset_rows_filter_button.n_clicks') {
                        if (Object.keys(original_row_mappings).length > 0) {
                            selected_rows = [];
                            for (let x = 0; x < selected_row_ids.length; x++) {
                                let row_id = String(selected_row_ids[x]);
                                selected_rows.push(original_row_mappings[row_id]);
                            }
                            original_row_mappings = {};
                        }
                        row_selection_filter = [];
                        page_current = 0;
                        start_highlight_first_row += 1;
                        // We reset the displayable rows to include all of the originally obtained agg_data
                        table_data = agg_data;
                    }
                    else if (triggered[i]['prop_id'] == 'agg_data.data') {
                        // This means that the user triggered a reaggregation of all the data that goes into the table.
                        // In this case, we need to wipe the slate clean, start from scratch with loading the new agg_data into the table.

                        table_data = agg_data
                        selected_rows = [];
                        selected_row_ids = [];
                        row_selection_filter = [];
                        page_current = 0;
                        original_row_mappings = {};
                    }
                }
            }

            // We defined the variable filter_set above, if and only if we need to filter out rows from the currently
            // displayed table_data
            if (filter_set != undefined) {
                var filtered_table_data = [];
                for (let a = 0; a < table_data.length; a++) {
                    if (filter_set.has(table_data[a]['id'])) {
                        filtered_table_data.push(table_data[a]);
                    }
                }
                table_data = filtered_table_data;
            }

            let displayed_row_count = JSON.stringify(table_data.length);

            let return_values = [
                table_data, snapfilter_applied, selected_rows, selected_row_ids, row_selection_filter, page_current,
                displayed_row_count, start_highlight_first_row, original_row_mappings
            ];

            return return_values;
        }
        """,
        [Output('table_transforms', 'data'),
        Output('snapfilter_applied', 'data'),
        Output('table_transforms', 'selected_rows'),
        Output('table_transforms', 'selected_row_ids'),
        Output('row_selection_filter', 'data'),
        Output('table_transforms', 'page_current'),
        Output('displayed_row_count', 'data'),
        Output('start_highlight_first_row', 'n_clicks'),
        Output('original_row_mappings', 'data')],
        [Input('select_all_button', 'n_clicks'),
        Input('deselect_all_button', 'n_clicks'),
        Input('filter_rows_button', 'n_clicks'),
        Input('reset_rows_filter_button', 'n_clicks'),
        Input('agg_data', 'data')],
        [State('snapfilter_applied', 'data'),
        State('snapfilter', 'data'),
        State('table_transforms', 'data'),
        State('table_transforms', 'derived_virtual_data'),
        State('table_transforms', 'selected_row_ids'),
        State('table_transforms', 'selected_rows'),
        State('original_row_mappings', 'data'),
        State('row_selection_filter', 'data'),
        State('table_transforms', 'page_current'),
        State('start_highlight_first_row', 'n_clicks')],
        prevent_initial_call = False
    )

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

    ### GENERATE FIGURE IN GRAPH
    ### CROSS-FILTER FIGURE BASED ON SELECTED ROWS IN TABLE_TRANSFORMS
    # Approach to connect selected row to elements in plot adapted from https://dash.plotly.com/datatable/interactivity
    @app_dash.callback(
        [Output('table_transforms', 'active_cell'),
        Output('finish_highlight_first_row', 'n_clicks'),
        Output('table_transforms', 'selected_cells'),
        Output('pairPlot', 'figure'),
        Output('pairPlot', 'selectedData'),
        Output('crossfilter-xaxis-type', 'value'),
        Output('crossfilter-yaxis-type', 'value'),
        Output('crossfilter-xaxis-type', 'options'),
        Output('crossfilter-yaxis-type', 'options')],
        [Input('colors_mapped_to_ids', 'data'),
        Input('table_transforms', 'active_cell'),
        Input('table_transforms', 'derived_viewport_row_ids'),
        Input('table_transforms', 'selected_row_ids'),
        Input('crossfilter-xaxis-type', 'value'),
        Input('crossfilter-yaxis-type', 'value'),
        Input('range_filter_names', 'data'),
        Input({'type': 'min_range_filter', 'index': ALL}, 'value'),
        Input({'type': 'max_range_filter', 'index': ALL}, 'value'),
        # Changing the x-axis or y-axis values will trigger reaggregation of properties, and calculation of a new minmax field, which then triggers this callback
        Input('minmax', 'data')],
        [State('pair_data', 'data'),
        State('start_highlight_first_row', 'n_clicks'),
        State('finish_highlight_first_row', 'n_clicks'),
        State('table_transforms', 'columns'),
        State('table_transforms', 'selected_cells'),
        State('query_data', 'data'),
        State('crossfilter-xaxis-column', 'value'),
        State('crossfilter-yaxis-column', 'value'),
        State('property_metadata', 'data'),
        State('display_name_to_property_name', 'data'),],
        prevent_initial_call = True
    )
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
                # In a separate callback where we change the number of rows in the table, we want to activate crossfiltering based on the first row, out of the newly displayed set of rows
                # Because a given id/value pair can only be the output of one callback, and we don't want to merge this callback with that other one, we need two separate flags to accomplish this behavior, which are below:
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
                    'row' : active_cell['row'],
                    'column' : i,
                    'column_id' : columns[i]['name'],
                    'row_id' : active_cell['row_id']
                }
                output_selected_cells.append(cell)
        # If we didn't click on any row to trigger this callback, then preserve selected cells, and associated crossfiltered points in scatterplot
        elif input_selected_cells != []:
            output_selected_cells = input_selected_cells
            active_row_id = input_selected_cells[0]['row_id']

        else:
            active_row_id = None

        # The table's row_ids will either be a single integer per row (if aggregation_type == 'individual_transforms'), or a json string of list of integers per row (if aggregation_type == 'group_by_fragment')
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
            # row_id is either an integer (rule_id) if we did not include env in results, or (for multicut queries with environment where we save the environments) a string with (from_smiles_env)_(to_smiles_env)
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
            fig = go.Figure(data = [], layout = go.Layout(
                    margin = {'l': 40, 'b': 40, 't': 10, 'r': 0},
                    hovermode='closest',
                    dragmode='lasso',
                    plot_bgcolor='#e8effa',
                    showlegend=False,
                ))
        else:
            rf_args = get_range_filter_args(rf_names, min_rfs, max_rfs, display_name_to_property_name)

            plot_args = {
                'query_id': query_id,
                'aggregation_type': aggregation_type,
                #'grouped_by_environment': True if 'environment_smarts' in (column['name'] for column in columns) else False,
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
                if grouped_by_environment == True:
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
            df['border_colors'] = ['#000000']*df['border_widths'].size

            rows_selected = df['colors'] != color_default
            df_selected = df

            ############# For enumerating building blocks, and/or targets
            # Get list of unique fragments in 'TO' column, which we can take into a jupyter notebook for now

            fig = go.Figure(
                data = [
                    go.Scattergl(
                        x = df_selected['x'],
                        y = df_selected['y'],
                        customdata = list(df_selected[['from_construct_id','to_construct_id']].to_records(index=False)),
                        mode = 'markers',
                        marker = {
                            'size' : df_selected['sizes'],
                            'color' : df_selected['colors'],
                            'opacity' : df_selected['opacities'],
                            'line' : {'width': df_selected['border_widths'], 'color': df_selected['border_colors']},
                        }
                    )      
                ],

                layout = go.Layout(
                    margin = {'l': 40, 'b': 40, 't': 10, 'r': 0},
                    hovermode='closest',
                    dragmode='lasso',
                    plot_bgcolor='#e8effa',
                    showlegend=False,
                )
            )

        ### HANDLE PLOT AXIS SCALING
        plot_x_args = {
            'title' : xaxis_column_name,
        }

        plot_y_args = {
            'title' : yaxis_column_name,
        }

        xaxis_log_disabled = False
        yaxis_log_disabled = False

        # We store min and max field values in this dictionary, minmax, which was created during statistics aggregation
        minmax = json.loads(minmax)
        # Make minmax names consistent with UI names
        minmax = {key.replace('average_log_', 'average(log) ').replace('average_', 'average ') : value for key, value in minmax.items()}

        # In these cases, we force linear scale (can't take log of negative numbers).
        if minmax['max_' + xaxis_column_name.lower()] <= 0 or minmax['min_' + xaxis_column_name.lower()] <= 0:
            xaxis_log_disabled = True
            xaxis_radio_type = 'Linear'
        if minmax['max_' + yaxis_column_name.lower()] <= 0 or minmax['min_' + yaxis_column_name.lower()] <= 0:
            yaxis_log_disabled = True
            yaxis_radio_type = 'Linear'

        xaxis_linearLog_radioRuttons_options = [
            {'label': 'Linear', 'value': 'Linear', 'disabled' : False},
            {'label': 'Log', 'value': 'Log', 'disabled': xaxis_log_disabled}
        ]
        yaxis_linearLog_radioRuttons_options = [
            {'label': 'Linear', 'value': 'Linear', 'disabled' : False},
            {'label': 'Log', 'value': 'Log', 'disabled': yaxis_log_disabled}
        ]

        if xaxis_radio_type == 'Linear':
            plot_x_args['type'] = 'linear'
            if 'delta' in xaxis_column_name:
                x_max = 1.05*max(abs(minmax['max_' + xaxis_column_name.lower()]), abs(minmax['min_' + xaxis_column_name.lower()]))
                x_min = None
                vline_x_value = 0
            else:
                x_max = minmax['max_' + xaxis_column_name.lower()]
                x_max += 0.05*(abs(x_max))
                x_min = minmax['min_' + xaxis_column_name.lower()]
                x_min -= 0.05*(abs(x_max))
        elif xaxis_radio_type == 'Log':
            plot_x_args['type'] = 'log'
            plot_x_args['dtick'] = 'D2'
            if 'fold-change' in xaxis_column_name:
                x_max = 1.05*max(abs(math.log(minmax['max_' + xaxis_column_name.lower()], 10)), abs(math.log(minmax['min_' + xaxis_column_name.lower()], 10)))       
                x_min = None
                vline_x_value = 1
            else:
                x_max = math.log(minmax['max_' + xaxis_column_name.lower()], 10)
                x_max += 0.05*(abs(x_max))
                x_min = math.log(minmax['min_' + xaxis_column_name.lower()], 10)
                x_min -= 0.05*(abs(x_max))

        if yaxis_radio_type == 'Linear':
            plot_y_args['type'] = 'linear'
            if 'delta' in yaxis_column_name:
                y_max = 1.05*max(abs(minmax['max_' + yaxis_column_name.lower()]), abs(minmax['min_' + yaxis_column_name.lower()]))
                y_min = None
                hline_y_value = 0
            else:
                y_max = minmax['max_' + yaxis_column_name.lower()]
                y_max += 0.05*(abs(y_max))
                y_min = minmax['min_' + yaxis_column_name.lower()]
                y_min -= 0.05*(abs(y_max))
        elif yaxis_radio_type == 'Log':
            plot_y_args['type'] = 'log'
            plot_y_args['dtick'] = 'D2'
            if 'fold-change' in yaxis_column_name:
                y_max = 1.05*max(abs(math.log(minmax['max_' + yaxis_column_name.lower()], 10)), abs(math.log(minmax['min_' + yaxis_column_name.lower()], 10)))
                y_min = None
                hline_y_value = 1
            else:
                y_max = math.log(minmax['max_' + yaxis_column_name.lower()], 10)
                y_max += 0.05*(abs(y_max))
                y_min = math.log(minmax['min_' + yaxis_column_name.lower()], 10)
                y_min -= 0.05*(abs(y_max))

        if x_min is None:
            plot_x_args['range'] = [-1*x_max, x_max]
        else:
            plot_x_args['range'] = [x_min, x_max]
        if y_min is None:
            plot_y_args['range'] = [-1*y_max, y_max]
        else:
            plot_y_args['range'] = [y_min, y_max]
            
        fig.update_xaxes(**plot_x_args)
        fig.update_yaxes(**plot_y_args)
        if x_min is None:
            fig.add_vline(x=vline_x_value)      
        if y_min is None:
            fig.add_hline(y=hline_y_value)

        # We only wanted the active cell to know what the user clicked on; once we select the entire row, we no longer need information about the active cell.
        # If we did leave the active cell activated, then it causes annoying behavior, where clicking a row checkbox away from the active cell causes the focus to jump to the active cell
        active_cell = None

        return [
            active_cell, finish_highlight_first_row, output_selected_cells, fig, None, xaxis_radio_type, yaxis_radio_type, xaxis_linearLog_radioRuttons_options,
            yaxis_linearLog_radioRuttons_options
        ]

    # When user clicks on a point in pair plot, push that pair to the top of the pair_tables
    @app_dash.callback(
        [Output('pair_tables', 'children'),
        Output('clearButton', 'n_clicks'),
        Output('pairPlot', 'clickData')],
        [Input('pairPlot', 'clickData'),
        Input('pairPlot', 'selectedData'),
        Input('pair_tables', 'children'),
        Input('clearButton', 'n_clicks')],
        [State('pair_data', 'data'),
        State('crossfilter-xaxis-column', 'value'),
        State('crossfilter-yaxis-column', 'value'),
        State('query_data', 'data'),
        State('property_metadata', 'data')]
    )
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

    # Writing this clientside so we don't have to send entire table data to the server
    # Instead we send only the specific rows we want
    app_dash.clientside_callback(
        """
        function get_enumeration_transforms (n_clicks, selected_row_ids, table_data) {

            let extracted = [];
            let selected_set = new Set(selected_row_ids);

            function extract_smiles (markdown_text) {

                // Beginning of example input string has this format:
                // let markdown_text = '![[*:1][H]](http://)';
                // Objective of this function is to extract the smiles within the ![] before the (http://)

                let match = /\]\(http/.exec(markdown_text);
                let index = match.index
                return markdown_text.slice(2, index);
            }

            // Extract rows based on selected row ids
            for (let i=0; i < table_data.length; i++) {
                if (selected_set.has(table_data[i]['id'])) {
                    // modified_row = JSON.parse(JSON.stringify(table_data[i]));
                    // Note for group_by_fragment table data, only the first 'FROM' smiles will be extracted
                    // table_data[i]['FROM'] = extract_smiles(table_data[i]['FROM']);

                    // Just extract everything
                    extracted.push(table_data[i]);
                }
            }

            return [extracted];
        }
        """,
        [Output('selected_row_data', 'data')],
        [Input('enumerate_button', 'n_clicks')],
        [State('table_transforms', 'selected_row_ids'),
        # the parent table_transforms.data represents the table's unsorted data
        # By using derived_virtual_data, we are accessing the rows after they have been sorted
        State('table_transforms', 'derived_virtual_data')]
    )

    @app_dash.callback(
        [Output("download_enumerations", "data")],
        [Input('selected_row_data', 'data')],
        snapshot_states + [
            State('query_data', 'data'),
            State('pair_data', 'data')
        ]
    )
    def enumerate_designs(*inputs):

        selected_row_data, query_data, pair_data = inputs[0], inputs[-2], inputs[-1]

        pair_data = json.loads(pair_data)
        query_id = pair_data['query_id']
        if len(selected_row_data) == 0:
            return [{'content':"No rows were selected. Click the checkboxes next to desired rows.", 'filename': "no_rows_selected_error.txt"}]

        new_snapshot_id = get_new_snapshot_id(inputs)

        schema = json.loads(inputs[-2])['schema']
        query_schema = f"?schema={schema}" if schema != "None" else ""
        link_address = external_frontend_root + '/snap/' + new_snapshot_id + query_schema

        data = {
            'query': json.loads(query_data),
            'query_id': query_id,
            'row_data': selected_row_data,
            'link_address': link_address,
        }
        sdf_string = requests.post(backend_root + f'/enumerate?schema={schema}', data=json.dumps(data)).json()

        if sdf_string == 'wildcards_in_core_error':
            return [{'content':"Wildcard atoms are not currently supported in the enumerated core (nonvariable part), because they cause problems during enumeration", 'filename': "wildcards_in_core_error.txt"}]

        current_time = datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
        return [{'content': sdf_string, 'filename': f"matcher_enumeration_{current_time}.sdf"}]

    # Callbacks must be defined before the server starts,
    # but our callbacks reference some elements that are not loaded on server start, which gives errors in the browser console
    # to suppress these errors, uncomment the below line:
    app_dash.config.suppress_callback_exceptions = True

    # Resize the Dash app iframe when its vertical size changes
    app_dash.clientside_callback(
        """
        function (output_div_children, pair_tables_children) {
            if (pair_tables_children.length === 1) {
                if (pair_tables_children[0] === "no_output") {
                    return {};
                }
            }
            parent.resize_iframe(parent.document.getElementById("outputPlot"));
            return {};
        }
        """,
        [Output('resize_iframe_dummy_output', 'data')],
        [Input('output_div', 'children'),
        Input('pair_tables', 'children')]
    )

    return app_dash

app = create_dash_app(flask_server=Flask_app)
server = app.server


if __name__ == '__main__':
    app.run_server(debug=True)
