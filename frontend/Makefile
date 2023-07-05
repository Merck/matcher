include .env

install:
	pip install -r requirements.txt
	pip install -r requirements.dev.txt

run:
	dotenv run gunicorn --reload --worker-class gevent --workers 1 --bind localhost:8000 monkey:server