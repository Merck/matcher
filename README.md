Matcher is a tool for understanding how chemical structure optimization problems have been solved.

Matcher enables deep control over searching structure/activity relationships (SAR) derived from large datasets, and takes the form of an accessible web application with simple deployment.

Matcher is built around the [mmpdb](https://github.com/rdkit/mmpdb) platform for matched molecular pair (MMP) analysis. Matcher extends the mmpdb data model, introduces new search algorithms, provides a backend API for executing queries and fetching results, and provides a frontend user interface.

# Table of Contents

1. [Quick Start](#quick_start)
2. [Run Example Query](#run_example_query)
3. [Data Included](#data_included)
4. [Use Different Data](#use_different_data)
5. [Metadata Information](#metadata_info)
6. [OpenAPI for Backend](#backend_OpenAPI)
7. [Using mmpdb Commands](#mmpdb_commands)

# Quick Start <a id="quick_start"></a>

Clone this repository, and the mmpdb submodule (with `git submodule init` and `git submodule update`).

Navigate to the parent matcher directory, then execute:

```
docker-compose up
```

Three containers will be launched:

1. database container, which contains the PostgreSQL database
2. backend container, which is controlled through a FastAPI accessible at localhost:8001
3. frontend container, which hosts a Dash app through a Flask API, accessible at localhost:8000

The rate-determining step is for the backend container to use mmpdb commands to process input data, and write resulting MMP data to the database. This process is complete when the below output is observed in the docker-compose console:

```
INFO:uvicorn.error:Uvicorn running on http://0.0.0.0:8001 (Press CTRL+C to quit)
```

After the above is observed, in a web browser (Google Chrome or Microsoft Edge), navigate to localhost:8000. The Matcher query interface should load after a few seconds.

# Run Example Query <a id="run_example_query"></a>

Matcher query logic can be learned by example: click the "Run Example Search" link at the top of the matcher query frontend page (localhost:8000), which will direct to localhost:8000/examples.

Several example inputs/outputs will be displayed. Upon clicking on an example, a new tab will open with the live Matcher interface, containing appropriately-populated input. Output results will load below the input.

Important: The example queries are only guaranteed to work with the example data provided in this package, because these queries specify properties that are present in the example data. If new data is used which has property names differing from the example data's property names, an exception will be thrown in the client js layer, and the example queries will not execute when loaded.

# Data Included <a id="data_included"></a>

Data is present in the backend/initialize_db directory.

<strong>Quick Start data (default)</strong>: Data filenames begin with "quick". Contains 1078 ChEMBL compounds, the minimum to fully reproduce queries described in our publication (TODO: Add hyperlink here).

<strong>Rapidly test/debug the deployment</strong>: Data filenames begin with "test". Contains 14 ChEMBL compounds, a subset of the Quick Start 1078 compounds, for the purpose of rapid testing during development or troubleshooting. All example queries work, but return only a few results.

<strong>Full ChEMBL dataset</strong>: Data filenames begin with "ChEMBL_CYP3A4_hERG". Contains 20267 ChEMBL compounds having CYP3A4 inhibition and/or hERG inhibition data, which were included with the [mmpdb publication](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00173). A superset of the Quick Start data.

# Use Different Data <a id="use_different_data"></a>

The default input compound/property dataset is intentionally very small, so that the containers will initialize quickly for demo purposes.

To use arbitrary data, follow the below steps.

As an example, we illustrate how to use a "medium-size" dataset containing 20,267 compounds taken from the [mmpdb publication](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00173), and referenced in our publication (TODO: Add hyperlink here).

1. Add raw data to the matcher/backend/initialize_db directory. Two files are required, a third file is optional:
    * **Required**: File containing compound SMILES and compound IDs.
        * For this example, ChEMBL_CYP3A4_hERG_structures.smi is already included.<br></br>
    * **Required**: File containing compound IDs and property values.
        * For this example, ChEMBL_CYP3A4_hERG_props.txt is already included.<br></br>
    * **Optional**: File containing metadata about the compound property data (whether the data is log transformed, the units, and how the data should be displayed to users).
        * For this example, ChEMBL_CYP3A4_hERG_metadata.csv is already included. If no metadata file is provided, then by default, property labels and data will be displayed to users exactly as provided in the above property value file, and changes between two properties will be treated as differences (B - A).
<br></br>

2. Edit matcher/backend/initialize_db/quick_start.sh to reference the new data.
    * For this example, we have already performed the editing, and the resulting file is matcher/backend/initialize_db/ChEMBL_CYP3A4_hERG_start.sh.
        * To use this example file, either rename matcher/backend/initialize_db/ChEMBL_CYP3A4_hERG_start.sh to matcher/backend/initialize_db/quick_start.sh (thus replacing the original file with same name), or edit matcher/backend/Dockerfile by changing the two occurrences of quick_start.sh to ChEMBL_CYP3A4_hERG_start.sh.
<br></br>

3. Recreate the containers:
    * Navigate to the parent directory (matcher/), and execute the following command:

```
docker-compose up --force-recreate
```

This time, around 20 minutes will be required to build the database (depending on computer), due to the larger number of compounds in the input data as compared to the original Quick Start data.

# Metadata Information <a id="metadata_info"></a>

Optionally, property metadata can be passed to the mmpdb loadprops command, for 
the purpose of customizing how data is displayed to end users in matcher's web UI.

If metadata is not provided, then all data will be displayed as it exists in
the database, and all changes will be calculated and displayed as deltas (B-A).

Here is an example property metadata table (values should be tab separated):

```
  property_name base  unit  display_name  display_base  display_unit  change_displayed
  hERG_pIC50  negative_log  M hERG_IC50_uM  raw uM  fold-change
```

Example:

  % mmpdb loadprops --properties hERG.csv --metadata hERG_medatadata.csv 'database$5432$postgres'

In the above case, hERG data is provided in the --properties file as a 
negative log of molarity, but we provide additional metadata
which causes the data to appear to users in the matcher web UI as micromolar IC50 values.
This is the most common base/unit conversion use case; not all conversions are supported.

Metadata is stored in columns within the property_name table, and therefore
is easy to modify (albeit manually with SQL) even after the database is built.

Supported metadata options are below. All other values, or * value, will be
converted to default. None means that the column value will be NULL in property_name table.

```
  base [default=raw]: raw, log, negative_log 
  unit [default=None]: M, uM
  display_name [default=property_name]: characters other than *
  display_base [default=raw]: raw
  display_unit [default=None]: uM, M
  change_displayed [default=delta]: delta, fold-change
```

# OpenAPI for Backend <a id="backend_OpenAPI"></a>

Matcher's backend API can be used for querying and gathering results, independently of the frontend, if desired.

The backend API endpoints are documented in backend/openapi.json. This documentation can be viewed when the matcher application is running, at localhost:8001/docs

# Using mmpdb Commands <a id="mmpdb_commands"></a>

The matcher database is an extended mmpdb database, and is reverse-compatible with mmpdb commands.

For example, to run `mmpdb transform` with matcher's database, outputting results to `results.csv` within your local directory:

First launch matcher as described in [Quick Start](#quick_start), then run this command:

```
docker exec -it \
"$(docker ps | grep 'matcher_backend' | awk '{ print $1 }')" \
conda run -n matcher-api \
python /opt/mmpdb/mmpdb.py transform --smiles 'O=C1NC2=C(C=NC(OC)=N2)N=C1' 'public$postgres' > results.csv
```
