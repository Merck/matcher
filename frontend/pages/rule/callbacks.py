import logging
from functools import partial

from dash import callback, clientside_callback
from dash.dependencies import Input, Output, State, ALL

from pages.rule.constants import create_id
from pages.common.callbacks import (run_persistent_query, instantiate_output_elements, aggregate_statistics_by_rule, selected_point_to_pair_tables, update_graph)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


clientside_callback(
    """
    function (n_clicks) {
        // This callback will always fire on page load, in order to support automatic query firing of snapshot-based queries
        // However, if we loaded the page without any snapshot-based query, then the fields will be empty; in that case, don't confuse the user with an error message
        // This block disables itself after the initial page load, so subsequently the user will get error messages if they try to submit an incomplete form
        data = {
            snapquery_id: parent.document.getElementById('rule_snapquery_id').value,
            query_type: parent.document.getElementById('rule_query_type').value,
            query_id: parent.document.getElementById('rule_query_id').value,
            transform_order: parent.document.getElementById('rule_transform_order').value,
            REQUIRED_properties: parent.document.getElementById('rule_REQUIRED_properties').value,
            OPTIONAL_properties: parent.document.getElementById('rule_OPTIONAL_properties').value,
            advanced_options: {
                'aggregation_type': parent.document.getElementById('rule_aggregation_type').value,
            },
            snapfilter_id: parent.document.getElementById('rule_snapfilter_id').value,
            snapfilter_string: parent.document.getElementById('rule_snapfilter_string').value,
            schema: parent.schema
        }
        output = JSON.stringify(data);
        return [output, parent.property_metadata, parent.display_name_to_property_name];
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


# Use this function to run the query, if your production environment kills long requests in a way that you cannot prevent
# TODO instead, use Celery or similar service
rule_run_persistent_query = callback(
    [
        Output(create_id('query_data'), 'data'),
        Output(create_id('pair_data'), 'data')
    ],
    [Input(create_id('input_data'), 'data')]
)(run_persistent_query)


# pair_data flows in from submit button (via DB query)
# Once we have the data, we can use the data to initialize all of the output elements
rule_instantiate_output_elements = callback(
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
rule_aggregate_statistics_by_rule = callback(
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
rule_selected_point_to_pair_tables = callback(
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


rule_update_graph = callback(
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
