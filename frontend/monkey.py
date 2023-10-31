from gevent import monkey
monkey.patch_all()

from frontend.dash_app import server  # noqa
