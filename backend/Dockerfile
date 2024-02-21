FROM condaforge/miniforge3

USER root

RUN apt-get update -y && apt-get install libxrender1 -y

WORKDIR /opt

COPY ./backend/environment.yml ./backend/requirements.txt ./backend/

RUN conda env create -f ./backend/environment.yml

COPY ./backend ./backend
COPY ./scripts ./scripts

EXPOSE 8001

ENV PYTHONPATH=/opt

RUN chmod +x ./backend/entrypoint.sh
ENTRYPOINT ["./backend/entrypoint.sh"]
