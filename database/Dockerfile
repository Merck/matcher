FROM postgres:13

# Install "rdkit" PostgreSQL extension
RUN apt-get update -y && apt-get install postgresql-13-rdkit -y

# Copy in the load-extensions script
COPY load-extensions.sh /docker-entrypoint-initdb.d/
