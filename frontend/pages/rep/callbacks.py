import logging
import time
import json
import requests
from functools import partial

from dash import callback, clientside_callback
from dash.dependencies import Input, Output, State, ALL

from config import backend_root
from pages.rep.constants import create_id
from pages.common.callbacks import (instantiate_output_elements, aggregate_statistics_by_rule, selected_point_to_pair_tables, update_graph)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


clientside_callback(
    """
    function (n_clicks) {
        // This callback will always fire on page load, in order to support automatic query firing of snapshot-based queries
        // However, if we loaded the page without any snapshot-based query, then the fields will be empty; in that case, don't confuse the user with an error message
        // This block disables itself after the initial page load, so subsequently the user will get error messages if they try to submit an incomplete form
        input_data = {
            schema: parent.schema,
            rule_environment_statistics_id: parent.rule_environment_statistics_id,
            rep_query_id: parent.rep_query_id,
            OPTIONAL_properties: parent.rep_property,

            // Fixed values for rule/environment/property (rep) view
            REQUIRED_properties: "",
            advanced_options: {'aggregation_type': "individual_transforms"},
            snapfilter_string: ""
        }
        input_data = JSON.stringify(input_data);
        return [input_data, parent.property_metadata, parent.display_name_to_property_name];
    }
    """,
    [
        Output(create_id('input_data'), 'data'),
        Output(create_id('property_metadata'), 'data'),
        Output(create_id('display_name_to_property_name'), 'data'),
    ],
    [Input(create_id('submit_button'), 'n_clicks')],
    prevent_initial_call=False
)

# Initiate a query and poll the DB for results over intervals
def run_persistent_query(input_data):
    input_data = json.loads(input_data)

    if input_data['rep_query_id'] != -1:
        # This means there are already stored results in the DB
        # Therefore, instead of running that query again and duplicating those results, we just skip the query and reference those previous results
        return [json.dumps(input_data), json.dumps({'query_id': input_data['rep_query_id']})]

    # Start the query and get an ID which we will use to poll the DB for existence of results over intervals
    schema = input_data['schema']
    started_query = requests.post(backend_root + f'/start_rep_query?schema={schema}', json={'rule_environment_statistics_id': input_data['rule_environment_statistics_id']})
    started_query = json.loads(started_query.json())
    query_id = started_query['query_id']

    # Set the time cutoff, in seconds, where we should give up on waiting for the query to finish
    timeout = 120
    elapsed = 0
    start_time = time.perf_counter()
    # Set the time, in seconds, for how frequently we poll the DB for existence of results that indicate query completion
    # We start with 2 seconds between each check, then increase interval as query becomes longer
    intervals = [
        {'interval_time': 2, 'start_next_interval_at': 14},
        # Purposefully use a time longer than the timeout so that this last interval is used until the end
        {'interval_time': 5, 'start_next_interval_at': timeout * 2},
    ]
    current_interval = 0

    while elapsed < timeout:

        time.sleep(intervals[current_interval]['interval_time'])
        query_status = requests.get(backend_root + f'/check_query/{query_id}?schema={schema}')
        query_status = json.loads(query_status.json())

        if query_status['finished'] == "True":
            #  The only case leading to display of result data, if we started a query
            # Connect output data (via query_id, which is used in the DB to represent the unique query to which a results set belongs) to the input state,
            # for use in snapshot generation
            input_data['query_id'] = query_id
            break

        elapsed = time.perf_counter() - start_time
        if elapsed > intervals[current_interval]['start_next_interval_at']:
            current_interval += 1

    return [json.dumps(input_data), json.dumps(started_query)]

rep_run_persistent_query = callback(
    [
        Output(create_id('query_data'), 'data'),
        Output(create_id('pair_data'), 'data')
    ],
    [Input(create_id('input_data'), 'data')]
)(run_persistent_query)


# pair_data flows in from submit button (via DB query)
# Once we have the data, we can use the data to initialize all of the output elements
rep_instantiate_output_elements = callback(
    [Output(create_id('output_div'), 'children')],
    [Input(create_id('pair_data'), 'data')],
    [
        State(create_id('query_data'), 'data'),
        State(create_id('property_metadata'), 'data'),
        State(create_id('display_name_to_property_name'), 'data'),
    ]
)(partial(instantiate_output_elements, create_id))


# Update aggregated data when user selects different properties or range filters
# We output agg_data, which is the parent data after aggregation, but which can undergo subsequent filtration by other clientside callbacks, before being loaded to the table
rep_aggregate_statistics_by_rule = callback(
    [
        Output(create_id('agg_data'), 'data'),
        Output(create_id('table_transforms'), 'columns'),
        Output(create_id('table_transforms'), 'style_data_conditional'),
        Output(create_id('minmax'), 'data'),
        Output(create_id('total_row_count'), 'data'),
        Output(create_id('colors_mapped_to_ids'), 'data')
    ],
    [
        Input(create_id('crossfilter-xaxis-column'), 'value'),
        Input(create_id('crossfilter-yaxis-column'), 'value'),
        Input({'type': create_id('min_range_filter'), 'index': ALL}, 'value'),
        Input({'type': create_id('max_range_filter'), 'index': ALL}, 'value')
    ],
    [
        State(create_id('range_filter_names'), 'data'),
        State(create_id('pair_data'), 'data'),
        State(create_id('query_data'), 'data'),
        State(create_id('property_metadata'), 'data'),
        State(create_id('display_name_to_property_name'), 'data'),
    ]
)(aggregate_statistics_by_rule)


# When user clicks on a point in pair plot, push that pair to the top of the pair_tables
rep_selected_point_to_pair_tables = callback(
    [
        Output(create_id('pair_tables'), 'children'),
        Output(create_id('clearButton'), 'n_clicks'),
        Output(create_id('pairPlot'), 'clickData')
    ],
    [
        Input(create_id('pairPlot'), 'clickData'),
        Input(create_id('pairPlot'), 'selectedData'),
        Input(create_id('pair_tables'), 'children'),
        Input(create_id('clearButton'), 'n_clicks')
    ],
    [
        State(create_id('pair_data'), 'data'),
        State(create_id('crossfilter-xaxis-column'), 'value'),
        State(create_id('crossfilter-yaxis-column'), 'value'),
        State(create_id('query_data'), 'data'),
        State(create_id('property_metadata'), 'data')
    ]
)(selected_point_to_pair_tables)


rep_update_graph = callback(
    [
        Output(create_id('table_transforms'), 'active_cell'),
        Output(create_id('finish_highlight_first_row'), 'n_clicks'),
        Output(create_id('table_transforms'), 'selected_cells'),
        Output(create_id('pairPlot'), 'figure'),
        Output(create_id('pairPlot'), 'selectedData'),
        Output(create_id('crossfilter-xaxis-type'), 'value'),
        Output(create_id('crossfilter-yaxis-type'), 'value'),
        Output(create_id('crossfilter-xaxis-type'), 'options'),
        Output(create_id('crossfilter-yaxis-type'), 'options')
    ],
    [
        Input(create_id('colors_mapped_to_ids'), 'data'),
        Input(create_id('table_transforms'), 'active_cell'),
        Input(create_id('table_transforms'), 'derived_viewport_row_ids'),
        Input(create_id('table_transforms'), 'selected_row_ids'),
        Input(create_id('crossfilter-xaxis-type'), 'value'),
        Input(create_id('crossfilter-yaxis-type'), 'value'),
        Input(create_id('range_filter_names'), 'data'),
        Input({'type': create_id('min_range_filter'), 'index': ALL}, 'value'),
        Input({'type': create_id('max_range_filter'), 'index': ALL}, 'value'),
        # Changing the x-axis or y-axis values will trigger reaggregation of properties, and calculation of a new minmax field, which then triggers this callback
        Input(create_id('minmax'), 'data')
    ],
    [
        State(create_id('pair_data'), 'data'),
        State(create_id('start_highlight_first_row'), 'n_clicks'),
        State(create_id('finish_highlight_first_row'), 'n_clicks'),
        State(create_id('table_transforms'), 'columns'),
        State(create_id('table_transforms'), 'selected_cells'),
        State(create_id('query_data'), 'data'),
        State(create_id('crossfilter-xaxis-column'), 'value'),
        State(create_id('crossfilter-yaxis-column'), 'value'),
        State(create_id('property_metadata'), 'data'),
        State(create_id('display_name_to_property_name'), 'data'),
    ],
    prevent_initial_call=True
)(update_graph)
