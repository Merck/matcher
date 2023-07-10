from flask import Flask, render_template, request
import os
import requests

app = Flask(__name__, static_url_path='', static_folder='', template_folder='')

backend_root = os.getenv("MATCHER_BACKEND_URL")
data = {
    'external_frontend_root': os.getenv("EXTERNAL_FRONTEND_URL"),
    'external_backend_root': os.getenv("EXTERNAL_BACKEND_URL"),
}


@app.route("/")
def matcher_homepage():
    data['schema'] = request.args.get('schema')
    return render_template("templates/matcher.html", **data)


@app.route("/examples")
def examples():
    data['schema'] = request.args.get('schema')
    return render_template("./examples/examples.html", **data)


@app.route("/snap/<snapshot_id>")
def snapshot_query(snapshot_id):
    schema = request.args.get('schema')
    query_schema = f"?schema={schema}" if schema is not None else ""
    snapshot_data = requests.get(backend_root + f"/snap_read/{snapshot_id}{query_schema}")
    snapshot_data = snapshot_data.json()
    snapshot_data = {**snapshot_data, **data}
    snapshot_data['schema'] = schema
    return render_template("templates/matcher.html", **snapshot_data)


@app.route("/rep/<rule_environment_id>")
def rep(rule_environment_id):
    # Load all data connected to a specific rule, environment, property (rep) triplet
    schema = request.args.get('schema')
    query_schema = f"?schema={schema}" if schema is not None else ""
    rep_property = request.args.get('prop')
    rep_data = requests.get(backend_root + f"/rep_read/{rule_environment_id}?prop={rep_property}").json()
    rep_data['schema'] = schema
    rep_data['rep_property'] = rep_property
    rep_data['rule_environment_id'] = rule_environment_id
    rep_data = {**rep_data, **data}
    return render_template("templates/rule.html", **rep_data)
