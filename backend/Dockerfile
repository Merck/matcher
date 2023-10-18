FROM condaforge/miniforge3

USER root

RUN apt-get update -y && apt-get install libxrender1 -y

WORKDIR /opt

COPY environment.yml requirements.txt ./

RUN conda env create -f /opt/environment.yml

ADD mmpdb ./mmpdb
ADD initialize_db ./initialize_db

COPY backend_api.py models.py search_algorithm.py smiles_syntax.py ss_select.py ./

EXPOSE 8001

ENTRYPOINT ["bash", "-c", "chmod +x ./initialize_db/test_start.sh && ./initialize_db/test_start.sh"]
