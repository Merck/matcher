from dash import Dash

from frontend.frontend_api import app as Flask_app

# imports callbacks required to register callbacks with the app
import frontend.pages.rep.callbacks as rep_callbacks  # noqa
import frontend.pages.snap.callbacks as snap_callbacks  # noqa


app_dash = Dash(
    __name__,
    assets_folder="css/dash/",
    routes_pathname_prefix='/dash/',
    prevent_initial_callbacks=True,
    server=Flask_app,
    serve_locally=True,
    use_pages=True,
    suppress_callback_exceptions=True
)

server = app_dash.server
