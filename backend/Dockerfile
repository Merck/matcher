FROM condaforge/miniforge3:23.3.1-1

USER root

RUN apt-get update -y && apt-get install libxrender1 -y

WORKDIR /opt

COPY ./backend/environment.yml ./backend/requirements.txt ./backend/

RUN conda env create -f ./backend/environment.yml

COPY ./backend ./backend
COPY ./scripts ./scripts

EXPOSE 8001

ENV PYTHONPATH=/opt

ENTRYPOINT ["bash", "-c", "echo $RUN_TESTS && chmod +x ./backend/initialize_db/test_start.sh && ./backend/initialize_db/test_start.sh"]
