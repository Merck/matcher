from gevent import monkey
monkey.patch_all()

from dash_app import server  # noqa
