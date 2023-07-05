from dash import dash_table
import pandas as pd
import urllib
import logging
from config import external_backend_root

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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


def get_styles(df, statistic=''):
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
    hex_number_for_color = '#' + '0' * (6 - len(hex_number)) + hex_number
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
