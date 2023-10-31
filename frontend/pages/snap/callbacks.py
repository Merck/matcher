from datetime import datetime
import json
import logging
from functools import partial

from dash import dcc, callback, clientside_callback
from dash.dependencies import Input, Output, State, ALL
import pandas as pd
import requests

from frontend.pages.common.callbacks import (run_persistent_query, instantiate_output_elements, aggregate_statistics_by_rule, selected_point_to_pair_tables, update_graph)
from frontend.config import backend_root, external_frontend_root
from frontend.pages.snap.constants import create_id

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# pair_data flows in from submit button (via DB query)
# Once we have the data, we can use the data to initialize all of the output elements
snap_instantiate_output_elements = callback(
    [Output('output_div', 'children')],
    [Input('pair_data', 'data')],
    [
        State('query_data', 'data'),
        State('property_metadata', 'data'),
        State('display_name_to_property_name', 'data'),
    ]
)((partial(instantiate_output_elements, create_id)))


# Collect user-entered data from matcher.html form
clientside_callback(
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
    [
        Output('input_data', 'data'),
        Output('property_metadata', 'data'),
        Output('display_name_to_property_name', 'data'),
    ],
    [Input('submit_button', 'n_clicks')],
    prevent_initial_call=False
)

# Use this function to run the query, if your production environment does not kill long requests in a way outside of your control
# TODO instead, use Celery or similar service
"""
@callback(
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
snap_run_persistent_query = callback(
    [
        Output('query_data', 'data'),
        Output('pair_data', 'data')
    ],
    [Input('input_data', 'data')]
)(run_persistent_query)


# Display various error messages to the user
clientside_callback(
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
    # logging.info(json.dumps(new_snapshot))
    # To write this query to DB and associate with the /snap endpoint, if needed outside of the normal application (e.g. during DB initialization),
    #   remove the 'query_id' key/value from new_snapshot and `await backend_api.snap_write(models.ExampleQuery.parse_obj(new_snapshot))`

    new_snapshot_id = new_snapshot_id.json()['new_snapshot_id']

    return new_snapshot_id


@callback(
    [Output('copy_link_clipboard', 'content')],
    [Input('copy_link_clipboard', 'n_clicks')],
    snapshot_states,
    prevent_initial_call=True
)
def spawn_and_copy_snapshot_link(*inputs):

    new_snapshot_id = get_new_snapshot_id(inputs)
    schema = json.loads(inputs[-1])['schema']
    query_schema = f"?schema={schema}" if schema != "None" else ""
    content = external_frontend_root + '/snap/' + new_snapshot_id + query_schema

    # dcc.Clipboard expects a list for some reason, but will copy the single element in this list to the clipboard
    return [content]


@callback(
    [Output("download-dataframe-csv", "data")],
    [Input("download_button", "n_clicks")],
    [
        State('pair_data', 'data'),
        State('query_data', 'data')
    ],
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


@callback(
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
@callback(
    [
        Output({'type': 'min_range_filter', 'index': ALL}, 'value'),
        Output({'type': 'max_range_filter', 'index': ALL}, 'value')
    ],
    [Input('reset_filters_button', 'n_clicks')],
    [State('range_filter_names', 'data')]
)
def reset_filters(n_clicks, range_filter_names):
    range_filter_names = json.loads(range_filter_names)
    return [[None] * len(range_filter_names)] * 2


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


@callback(
    [
        Output('num_rows_displayed_plaintext', 'children'),
        Output('num_rows_displayed_plaintext', 'style'),
    ],
    [
        Input('displayed_row_count', 'data'),
        Input('total_row_count', 'data')
    ],
)
def show_number_of_displayed_rows(displayed_row_count, total_row_count):

    displayed_row_count = json.loads(displayed_row_count)
    total_row_count = json.loads(total_row_count)
    if displayed_row_count != total_row_count:
        text_color = "#fc0303"
    else:
        text_color = "#000000"
    output_text = 'Displaying ' + str(displayed_row_count) + ' out of ' + str(total_row_count) + ' rows'
    output_style = {'font-size': '15px', 'font-weight': 'bold', 'color': text_color}

    return output_text, output_style


@callback(
    [
        Output('reset_filters_button', 'n_clicks'),
        Output('reset_rows_filter_button', 'n_clicks')
    ],
    Input('reset_all_filters_button', 'n_clicks'),
    [
        State('reset_filters_button', 'n_clicks'),
        State('reset_rows_filter_button', 'n_clicks')
    ]
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
snap_aggregate_statistics_by_rule = callback(
    [
        Output('agg_data', 'data'),
        Output('table_transforms', 'columns'),
        Output('table_transforms', 'style_data_conditional'),
        Output('minmax', 'data'),
        Output('total_row_count', 'data'),
        Output('colors_mapped_to_ids', 'data')
    ],
    [
        Input('crossfilter-xaxis-column', 'value'),
        Input('crossfilter-yaxis-column', 'value'),
        Input({'type': 'min_range_filter', 'index': ALL}, 'value'),
        Input({'type': 'max_range_filter', 'index': ALL}, 'value')
    ],
    [
        State('range_filter_names', 'data'),
        State('pair_data', 'data'),
        State('query_data', 'data'),
        State('property_metadata', 'data'),
        State('display_name_to_property_name', 'data'),
    ]
)(aggregate_statistics_by_rule)

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
clientside_callback(
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
    [
        Output('table_transforms', 'data'),
        Output('snapfilter_applied', 'data'),
        Output('table_transforms', 'selected_rows'),
        Output('table_transforms', 'selected_row_ids'),
        Output('row_selection_filter', 'data'),
        Output('table_transforms', 'page_current'),
        Output('displayed_row_count', 'data'),
        Output('start_highlight_first_row', 'n_clicks'),
        Output('original_row_mappings', 'data')
    ],
    [
        Input('select_all_button', 'n_clicks'),
        Input('deselect_all_button', 'n_clicks'),
        Input('filter_rows_button', 'n_clicks'),
        Input('reset_rows_filter_button', 'n_clicks'),
        Input('agg_data', 'data')
    ],
    [
        State('snapfilter_applied', 'data'),
        State('snapfilter', 'data'),
        State('table_transforms', 'data'),
        State('table_transforms', 'derived_virtual_data'),
        State('table_transforms', 'selected_row_ids'),
        State('table_transforms', 'selected_rows'),
        State('original_row_mappings', 'data'),
        State('row_selection_filter', 'data'),
        State('table_transforms', 'page_current'),
        State('start_highlight_first_row', 'n_clicks')
    ],
    prevent_initial_call=False
)


# GENERATE FIGURE IN GRAPH
# CROSS-FILTER FIGURE BASED ON SELECTED ROWS IN TABLE_TRANSFORMS
# Approach to connect selected row to elements in plot adapted from https://dash.plotly.com/datatable/interactivity
snap_update_graph = callback(
    [
        Output('table_transforms', 'active_cell'),
        Output('finish_highlight_first_row', 'n_clicks'),
        Output('table_transforms', 'selected_cells'),
        Output('pairPlot', 'figure'),
        Output('pairPlot', 'selectedData'),
        Output('crossfilter-xaxis-type', 'value'),
        Output('crossfilter-yaxis-type', 'value'),
        Output('crossfilter-xaxis-type', 'options'),
        Output('crossfilter-yaxis-type', 'options')
    ],
    [
        Input('colors_mapped_to_ids', 'data'),
        Input('table_transforms', 'active_cell'),
        Input('table_transforms', 'derived_viewport_row_ids'),
        Input('table_transforms', 'selected_row_ids'),
        Input('crossfilter-xaxis-type', 'value'),
        Input('crossfilter-yaxis-type', 'value'),
        Input('range_filter_names', 'data'),
        Input({'type': 'min_range_filter', 'index': ALL}, 'value'),
        Input({'type': 'max_range_filter', 'index': ALL}, 'value'),
        # Changing the x-axis or y-axis values will trigger reaggregation of properties, and calculation of a new minmax field, which then triggers this callback
        Input('minmax', 'data')
    ],
    [
        State('pair_data', 'data'),
        State('start_highlight_first_row', 'n_clicks'),
        State('finish_highlight_first_row', 'n_clicks'),
        State('table_transforms', 'columns'),
        State('table_transforms', 'selected_cells'),
        State('query_data', 'data'),
        State('crossfilter-xaxis-column', 'value'),
        State('crossfilter-yaxis-column', 'value'),
        State('property_metadata', 'data'),
        State('display_name_to_property_name', 'data'),
    ],
    prevent_initial_call=True
)(update_graph)

# When user clicks on a point in pair plot, push that pair to the top of the pair_tables
snap_selected_point_to_pair_tables = callback(
    [
        Output('pair_tables', 'children'),
        Output('clearButton', 'n_clicks'),
        Output('pairPlot', 'clickData')
    ],
    [
        Input('pairPlot', 'clickData'),
        Input('pairPlot', 'selectedData'),
        Input('pair_tables', 'children'),
        Input('clearButton', 'n_clicks')
    ],
    [
        State('pair_data', 'data'),
        State('crossfilter-xaxis-column', 'value'),
        State('crossfilter-yaxis-column', 'value'),
        State('query_data', 'data'),
        State('property_metadata', 'data')
    ]
)(selected_point_to_pair_tables)


# Writing this clientside so we don't have to send entire table data to the server
# Instead we send only the specific rows we want
clientside_callback(
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
    [
        State('table_transforms', 'selected_row_ids'),
        # the parent table_transforms.data represents the table's unsorted data
        # By using derived_virtual_data, we are accessing the rows after they have been sorted
        State('table_transforms', 'derived_virtual_data')
    ],
    prevent_initial_call=True
)


@callback(
    [Output("download_enumerations", "data")],
    [Input('selected_row_data', 'data')],
    snapshot_states + [
        State('query_data', 'data'),
        State('pair_data', 'data')
    ],
    prevent_initial_call=True
)
def enumerate_designs(*inputs):
    selected_row_data, query_data, pair_data = inputs[0], inputs[-2], inputs[-1]

    pair_data = json.loads(pair_data)
    query_id = pair_data['query_id']
    if len(selected_row_data) == 0:
        return [{'content': "No rows were selected. Click the checkboxes next to desired rows.", 'filename': "no_rows_selected_error.txt"}]

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
        return [{
            'content': "Wildcard atoms are not currently supported in the enumerated core (nonvariable part), because they cause problems during enumeration",
            'filename': "wildcards_in_core_error.txt"
        }]

    current_time = datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    return [{'content': sdf_string, 'filename': f"matcher_enumeration_{current_time}.sdf"}]

# Callbacks must be defined before the server starts,
# but our callbacks reference some elements that are not loaded on server start, which gives errors in the browser console
# to suppress these errors, uncomment the below line:
# TODO app_dash.config.suppress_callback_exceptions = True


# Resize the Dash app iframe when its vertical size changes
clientside_callback(
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
    [
        Input('output_div', 'children'),
        Input('pair_tables', 'children')
    ]
)
