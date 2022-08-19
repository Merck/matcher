from flask import Flask, render_template
import os
import requests

app = Flask(__name__, static_url_path='', static_folder='', template_folder='')

backend_root = os.getenv("MATCHER_BACKEND_URL")

@app.route("/")
def matcher_homepage():
    return render_template("matcher.html")

@app.route("/examples")
def examples():
    return render_template("./examples/examples.html")

@app.route("/snap/<snapshot_id>")
def snapshot_query(snapshot_id):
    snapshot_data = requests.get(backend_root + f"/snap_read/{snapshot_id}")
    snapshot_data = snapshot_data.json()
    return render_template("matcher.html", **snapshot_data)
