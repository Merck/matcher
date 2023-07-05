include .env

install:
	pip install -r requirements.txt
	pip install -r requirements.dev.txt

run:
	dotenv run uvicorn backend_api:app --reload --host 0.0.0.0 --port 8001