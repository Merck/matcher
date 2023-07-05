import dash
from dash import dcc, html

from pages.snap.constants import DASH_PAGE_PREFIX

dash.register_page(__name__, path=f'/{DASH_PAGE_PREFIX}/')


# The initially loaded layout will only show submit button.
# Clicking on one of these two buttons will obtain data, then instantiate and initialize the rest of the elements
# Subsequent clicking on one of these two buttons will overwrite all the previously instantiated elements
def generate_layout():

    return html.Div([

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
                                html.Button("Submit query", id="submit_button", style={'backgroundColor': '#7a0dd9', 'color': '#ffffff', 'border-style': 'none'}),
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
            id='output_div',
            children=[]
        ),
    ])


layout = generate_layout
