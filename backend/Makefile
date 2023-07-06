# This Makefile is intended to be use to run backend server in isolation with a hosted database server
# This .env file has credentials for database
include .env

venv: ## Create virtual environment and install packages, for backend server
	python3 -m venv venv
	. venv/bin/activate && pip install -r requirements.txt && pip install -r requirements.dev.txt

run: ## Run the backend server
	. venv/bin/activate && dotenv run uvicorn backend_api:app --reload --host 0.0.0.0 --port 8001
