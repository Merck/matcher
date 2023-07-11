# This Makefile is intended to run frontend server in isolation
include .env

venv: ## Create virtual environment and install packages, for frontend server
	python3 -m venv venv
	. venv/bin/activate && pip install -r requirements.txt && pip install -r requirements.dev.txt

run: ## Run frontend server
	. venv/bin/activate && dotenv run gunicorn --reload --worker-class gevent --workers 1 --bind localhost:8000 monkey:server
