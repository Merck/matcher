FROM condaforge/miniforge3

WORKDIR /opt

COPY environment.yml requirements.txt ./

RUN conda env create -f /opt/environment.yml

COPY frontend_api.py monkey.py config.py dash_app.py ./
COPY examples/ ./examples/
COPY css/ ./css/
COPY js/* ./js/
COPY ketcher/ ./ketcher/
COPY pages/ ./pages/
COPY templates/ ./templates/

EXPOSE 8000

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "matcher-frontend", "gunicorn", "--worker-class", "gevent", "--workers", "2", "--bind", "0.0.0.0:8000", "monkey:server"]
