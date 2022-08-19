#!/bin/bash

structures=./initialize_db/ChEMBL_CYP3A4_hERG_structures.smi
fragments=./initialize_db/ChEMBL_CYP3A4_hERG_structures.fragments
properties=./initialize_db/ChEMBL_CYP3A4_hERG_props.txt
metadata=./initialize_db/ChEMBL_CYP3A4_hERG_metadata.csv
example_queries=./initialize_db/example_queries.json

COMPLETION_FILE=./mmpdb_build_complete
FAILURE_FILE=./mmpdb_build_failed
if [ -f $COMPLETION_FILE ]; then
    echo "mmpdb build was done in a previous container startup, skipping..."
elif [ -f $FAILURE_FILE ]; then
    echo "There was an error with mmpdb database initialization during the first startup of backend_api container. Please recreate the containers with docker-compose up --force-recreate"
    exit 1
else
    # Populate mmpdb database from scratch
    {
        # Standard mmpdb command for generating fragments, except we defined new 'matcher_alpha' fragmentation criteria
        conda run --no-capture-output -n matcher-api python ./mmpdb/mmpdb.py fragment "${structures}" -o "${fragments}" --cut-smarts 'matcher_alpha' && \
        # Standard mmpdb command for identifying MMPs and loading data to DB, except we introduced postgres support, and extended the data model
        conda run --no-capture-output -n matcher-api python ./mmpdb/mmpdb.py index "${fragments}" -o 'database$5432$postgres' && \
        # Standard mmpdb command for loading property data to DB, except we introduced postgres support and ability to add property metadata
        conda run --no-capture-output -n matcher-api python ./mmpdb/mmpdb.py loadprops -p "${properties}" --metadata "${metadata}" 'database$5432$postgres' && \

        # Here we load JSON representations of query input states, to a table in the database
        # Each input query will have a unique integer snapshot_id assigned, starting from 1, going up to the n total number of example query inputs that we write to DB
        # The query can be run by calling a frontend API endpoint using the this snapshot_id as an argument, for example to run the first query, http://localhost:8000/snap/1
        # We embed these example queries in the "Run Example Query" page at http://localhost:8000/examples
        # To generate the required JSON for example_queries.json, print out the JSON that is passed to backend_api.snap_write upon clicking the "Copy Shareable Link" button in Matcher frontend,
        #   then either delete the 'query_id' key/value, or change the 'query_id' value to 'NULL'
        conda run --no-capture-output -n matcher-api python ./initialize_db/load_example_queries.py "${example_queries}" && \
        touch "${COMPLETION_FILE}"
    } || {
        touch "${FAILURE_FILE}"
        exit 1
    }
fi

conda run --no-capture-output -n matcher-api uvicorn backend_api:app --host 0.0.0.0 --port 8001
