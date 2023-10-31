import dash
from dash import dcc, html

import logging
from frontend.pages.rep.constants import DASH_PAGE_PREFIX, create_id

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


dash.register_page(__name__, path=f'/{DASH_PAGE_PREFIX}/')


# The initially loaded layout will only show submit button.
# Clicking on one of these two buttons will obtain data, then instantiate and initialize the rest of the elements
# Subsequent clicking on one of these two buttons will overwrite all the previously instantiated elements
def generate_layout():

    layout = html.Div([

        dcc.Store(id=create_id('resize_iframe_dummy_output'), data=['{}']),
        dcc.Store(id=create_id('error_message_dummy_output'), data=['{}']),
        dcc.Store(id=create_id('property_metadata'), data={}),
        dcc.Store(id=create_id('display_name_to_property_name'), data={}),

        # When the query is running, or data is loading, show a loading icon in place of the submit button, so user doesn't spam queries by accident
        dcc.Loading(
            [
                dcc.Store(id=create_id('input_data'), data=['{}']),
                dcc.Store(id=create_id('query_data'), data=['{}']),
                dcc.Store(id=create_id('pair_data'), data=['{}']),
                dcc.Store(id='rule_environment_statistics_id', data=0),
                dcc.Store(id='rep_property', data=''),

                html.Div(
                    [
                        html.Div(
                            [
                                html.Button(
                                    "Submit query",
                                    id=create_id('submit_button'),
                                    style={'backgroundColor': '#7a0dd9', 'color': '#ffffff', 'border-style': 'none'}),
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
        html.Div(
            id=create_id('output_div'),
            children=[]
        ),
    ])
    return layout


layout = generate_layout
