FROM condaforge/miniforge3

WORKDIR /opt

COPY ./frontend/environment.yml ./frontend/requirements.txt ./frontend/

RUN conda env create -f ./frontend/environment.yml

COPY ./frontend ./frontend
COPY ./scripts ./scripts

EXPOSE 8000

ENV PYTHONPATH=/opt

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "matcher-frontend", "gunicorn", "--worker-class", "gevent", "--workers", "2", "--bind", "0.0.0.0:8000", "frontend.monkey:server"]
