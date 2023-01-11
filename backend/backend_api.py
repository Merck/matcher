import uvicorn
import asyncpg
from fastapi import FastAPI, HTTPException, Body, Request, Response, File, UploadFile, BackgroundTasks
from fastapi.responses import ORJSONResponse, HTMLResponse, FileResponse
from fastapi.routing import APIRoute
from fastapi.middleware.cors import CORSMiddleware
from starlette.responses import StreamingResponse
from typing import Callable, List
from models import QueryInput, ValidateSelectionInput, AggregationParameters, PlotParameters, PairsCondensed, \
    GetAllRawData, EnumerationData

# Modules for chemistry operations
# Processes user-defined structure information in the Input UI
import ss_select
# Processes Matcher queries, and does post-processing on query results
import search_algorithm

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor, TemplateAlign
import io
import json
import time
import os
import logging
import copy
from collections import OrderedDict
logging.basicConfig(level=logging.INFO)


username = os.getenv("POSTGRES_USER")
password = os.getenv("POSTGRES_PASSWORD")
database_name = os.getenv("POSTGRES_DB")
database_hostname = os.getenv("POSTGRES_HOST")
database_port = os.getenv("POSTGRES_PORT")


async def get_matcher_conn(schema='public', timeout=3600):

    conn = await asyncpg.connect(
        host=database_hostname,
        port=database_port,
        database=database_name,
        user=username,
        password=password,
        timeout=timeout
    )
    if schema not in ["public", "None", ""]:
        # Enforce lowercase schema names. This prevents headaches involving always needing to use quotes in subsequent SQL
        schema = schema.lower()

        # asyncpg doesn't let us use parameters in SET statements, so we need another mechanism to prevent SQL injection
        # Therefore, here we validate the schema (which comes from an API query parameter),
        #  against schema in the DB
        valid_schemas = await conn.fetch("""
SELECT schema_name
FROM information_schema.schemata
WHERE schema_name NOT IN ('information_schema', 'pg_catalog')
""")
        valid_schemas = set(row['schema_name'] for row in valid_schemas)
        assert schema in valid_schemas

        # Need to include public schema in search_path for rdkit cartridge to work
        await conn.execute(f"SET search_path={schema}, public")
    return conn


class ValidationErrorLoggingRoute(APIRoute):
    def get_route_handler(self) -> Callable:
        original_route_handler = super().get_route_handler()

        async def custom_route_handler(request: Request) -> Response:
            try:
                return await original_route_handler(request)
            except Exception as exc:
                """
                body = await request.body()
                detail = {"errors": exc.errors(), "body": body.decode()}
                raise HTTPException(status_code=422, detail=detail)
                """

                # Insert desired exceptions to raise here
                # For now, just raise the exception as would otherwise occur
                raise exc

                print(exc)

                detail = {"error_message": str(exc.args)}
                raise HTTPException(status_code=500, detail=detail)

        return custom_route_handler

description = """
This API handles core cheminformatics logic of the matcher application.

This backend API can be used to run queries and interrogate results, independently of the matcher frontend.

## Running a query asynchronously

Follow the below steps to run a matcher-style query against the database containing MMP data.

Queries are run asynchronously to prevent blocking by longer queries.

0. (Optional) Validate the variable/environment atom/bond selections within input structure(s) using the /validateSelection endpoint

1. Pass query input to the /start_query endpoint, which quickly returns a query_id. The API then triggers the query asynchronously. <strong>Use the returned query_id for all subsequent endpoints below</strong>

2. Monitor query progress using the /check_query endpoint. The query is complete when the returned JSON "finished" key has a value of "True"

3. Results (if any are found by the query) are stored in the database query_result table, where query_result.query_id is the query_id returned by the /start_query endpoint above

## Get query results

<strong>To dump all results at once:</strong> when the query is finished, as indicated by /check_query, pass the query_id obtained from /start_query to /get_all_raw_data to obtain a list of MMPs (matched molecular pairs) found by the query, along with all associated structure/activity/transform data for each MMP.

<strong>What the matcher application does to get results:</strong>

The below endpoints were developed for more selective, iterated, scalable results gathering. For example, if a query finds 1 million MMPs, it's impractical to pass all that data over the network and/or store it in memory, especially because the user often won't care about most of the results. Instead, we can look selectively at interesting subsets of the results.

Run the below endpoints in sequence to go from high-level transform data to low-level MMP data:

<strong>/aggregate_transforms --> /get_plot_data --> /get_pair_data</strong>

* Use /aggregate_transforms to obtain MMP transform-level data, including statistics, e.g. median change in a specified property for each MMP transform found by the query. The matcher frontend uses this data to populate the table of transforms in the frontend results.

* Use /get_plot_data to get data for all MMPs belonging to specific transform(s). As input, the necessary ids for the transforms are provided by /aggregate_transforms. The matcher frontend uses this data to generate a scatterplot, e.g. change in property A vs. change in property B, where each point is an MMP belonging to specified transform(s).

* Use /get_pair_data to get data for specific MMPs. As input, the necessary ids for the pairs are provided by /get_plot_data. The matcher frontend uses this data to display data for individual MMPs. 
"""

app = FastAPI(
    title="matcher backend",
    description=description,
)
origins = [
    "http://127.0.0.1:8000",
    "http://localhost:8000",
    "http://frontend:8000"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.router.route_class = ValidationErrorLoggingRoute


@app.post("/validateSelection/", response_class=ORJSONResponse)
async def validateSelection(data: ValidateSelectionInput):
    """
    Checks whether specific atom/bond selections will lead to an illegal query, or will conflict with already selected atoms/bonds

    Typical use case: set selected atoms/bonds as either variable or environment
    
    "Autocorrects" illegal selections to prevent most such selections that would never lead to results.

    For example, an exact search for a variable fragment that is a subsection of a ring: this search would not find any results, for typical use cases where MMPs were defined with 'fragment and index' algorithms that do not cut cyclic bonds. Therefore, validateSelection would automatically grow the variable selection to be the entire ring, instead of just the subsection of the ring.

    Another example: When trying to set selected atoms as variable atoms, but if some of these atoms were previously set as environment atoms, all environment labels will be cleared due to the overlap.

    In addition, for "variable" selections, this endpoint will automatically select single neighboring atom(s) (relative to "variable" selection) as an "environment" selection (radius 1).

    The matcher frontend buttons labeled "Set variable atoms" (i.e. variable_atoms) and "Set constant atoms" (i.e. environment_atoms) call this endpoint

    Parameters:

    (note that all atom/bond indices are provided as strings of comma-separated integers, e.g. "1,2,7")

    - **selection_type**: "variable" or "environment", the category to which the selected atom indices will be added in sketched_content, after successful validation
    - **sketched_content**: chemical structure constraints of left-hand-side and right-hand-side portions of the MMP transform and environment. The variable/environment atoms/bonds refer to PREVIOUSLY validated selections. **sketched_content** is an object with below items:
        - **mol1_molfile**: MDL molfile string of the 'starting', 'left-hand-side' structure in the query (if any)
        - **mol1_variable_atoms**: indices of atoms in the variable fragment (which change as part of the transform). The order of atoms in the molfile defines the atom indices, starting from 0
        - **mol1_variable_bonds**: indices of bonds in the variable fragment. The order of bonds in the molfile defines the bond indices, starting from 0
        - **mol1_environment_atoms**: indices of atoms in the environment fragment (which do not change as part of the transform)
        - **mol1_environment_bonds**: indices of bonds in the environment fragment
        - **mol2_molfile**: MDL molfile string of the 'ending', 'right-hand-side' structure in the query (if any)
        - **mol2_variable_atoms**: indices of atoms in the variable fragment (which change as part of the transform). The order of atoms in the molfile defines the atom indices, starting from 0
        - **mol2_variable_bonds**: indices of bonds in the variable fragment. The order of bonds in the molfile defines the bond indices, starting from 0
        - **mol2_environment_atoms**: indices of atoms in the environment fragment (which do not change as part of the transform)
        - **mol2_environment_bonds**: indices of bonds in the environment fragment
    - **mol1_selected_atoms**: indices of atoms in mol1 selected by the user, which are about to undergo validation
    - **mol1_selected_bonds**: indices of bonds in mol1 selected by the user, which are about to undergo validation
    - **mol2_selected_atoms**: indices of atoms in mol2 selected by the user, which are about to undergo validation
    - **mol2_selected_bonds**: indices of bonds in mol2 selected by the user, which are about to undergo validation

    Returns JSON with validated variable and environment atom/bond indices:

    - **mol1_variable_atoms**: array of integers
    - **mol1_variable_bonds**: array of integers
    - **mol1_environment_atoms**: array of integers
    - **mol1_environment_bonds**: array of integers
    - **mol1_entire_molecule_selected**: "True" or "False": Used to signal whether the entire structure was selected as a variable fragment. If "True", the query will return no results, because matcher relies on detection of borders between variable and non-variable atoms to determine point of attachment of the variable atoms.
    - **mol2_variable_atoms**: array of integers
    - **mol2_variable_bonds**: array of integers
    - **mol2_environment_atoms**: array of integers
    - **mol2_environment_bonds**: array of integers
    - **mol2_entire_molecule_selected**: "True" or "False": Used to signal whether the entire structure was selected as a variable fragment. If "True", the query will return no results, because matcher relies on detection of borders between variable and non-variable atoms to determine point of attachment of the variable atoms.
    """
    keys = ('mol1_variable_atoms', 'mol1_variable_bonds', 'mol1_environment_atoms', 'mol1_environment_bonds', 'mol1_entire_molecule_selected',
            'mol2_variable_atoms', 'mol2_variable_bonds', 'mol2_environment_atoms', 'mol2_environment_bonds', 'mol2_entire_molecule_selected',)
    values = ss_select.validate_selection(data.selection_type, data.sketched_content, data.mol1_selected_atoms, data.mol1_selected_bonds, data.mol2_selected_atoms, data.mol2_selected_bonds)
    return_value = dict(zip(keys, values))
    return return_value


def parse_props(required, optional):
    if required != '':
        req_props = required.split(',')
    else:
        req_props = []

    if optional != '':
        opt_props = optional.split(',')
        req_props_set = set(req_props)
        unique_opt_props = []
        for prop in opt_props:
            if prop in req_props_set:
                continue
            unique_opt_props.append(prop)
        opt_props = unique_opt_props
    else:
        opt_props = []

    return [req_props, opt_props]


def parse_query_input(query_input: QueryInput):

    req_props, opt_props = parse_props(query_input.REQUIRED_properties, query_input.OPTIONAL_properties)
    # attach_maps should work with sketched_content as either molfile or rxnfile, assuming variable/environment atom indices correctly correspond to the ordering in the atom block(s)
    # mapped_molfile = ss_select.attach_maps(query_input.sketched_content, query_input.variable_atoms, query_input.environment_atoms)

    mapped_mol1_molfile, mapped_mol2_molfile = None, None
    assert query_input.sketched_content.mol1_variable_atoms != '' or query_input.sketched_content.mol2_variable_atoms != ''

    if query_input.sketched_content.mol1_variable_atoms != '':
        mapped_mol1_molfile = ss_select.attach_maps(query_input.sketched_content.mol1_molfile, query_input.sketched_content.mol1_variable_atoms, query_input.sketched_content.mol1_environment_atoms)
    if query_input.sketched_content.mol2_variable_atoms != '':
        mapped_mol2_molfile = ss_select.attach_maps(query_input.sketched_content.mol2_molfile, query_input.sketched_content.mol2_variable_atoms, query_input.sketched_content.mol2_environment_atoms)

    variable1, variable2, constant, observations = ss_select.extract_variable_constant(mapped_mol1_molfile, mapped_mol2_molfile)

    parsed_query_input = {
        'query_type': query_input.query_type,
        'transform_order': query_input.transform_order,
        'REQUIRED_properties': req_props,
        'OPTIONAL_properties': opt_props,
        'constant': constant,
        'variable1': variable1,
        'variable2': variable2,
        'compound1': None,
        'compound2': None,
        'observations': observations,
        'variable_min_heavies': query_input.advanced_options.variable_min_heavies,
        'variable_max_heavies': query_input.advanced_options.variable_max_heavies,
        'compound_min_heavies': query_input.advanced_options.compound_min_heavies,
        'compound_max_heavies': query_input.advanced_options.compound_max_heavies,
        'aggregation_type': query_input.advanced_options.aggregation_type,
    }
    return parsed_query_input


@app.post("/async_query_v3_backend", deprecated=True)
async def async_query_v3_backend(query: QueryInput, schema: str = 'public'):
    """
    Deprecated, replaced with /start_query and /check_query
    """

    parsed_query = parse_query_input(query)
    observations = parsed_query["observations"]
    
    # Prevent query if it will be WAY too general, returning many many GB: only happens in extreme cases
    if parsed_query['constant'] is None and parsed_query['variable1'] is None:
        if "variable_selection_over_limit" not in observations:
            observations["query_too_general"] = 'no_variable_no_constant'
    if parsed_query['query_type'] == 'substructure':
        if parsed_query['constant'] is None or 'has_one_atom_constant' in observations:
            if parsed_query['variable1'] is None or 'has_one_atom_variable' in observations:
                observations["query_too_general"] = 'many_to_many_too_general'

    if "all_atoms_selected_as_variable" in observations or "variable_selection_over_limit" in observations \
            or "query_too_general" in observations:
        return json.dumps({
            "url": "fatal",
            "observations": observations,
        })
    observations["no_results"] = "False"

    req, opt = parsed_query["REQUIRED_properties"], parsed_query["OPTIONAL_properties"]
    num_props = len(req) + len(opt)

    conn = await get_matcher_conn(schema=schema)

    # Make sure we close the connection even if an exception is raised
    try:
        # Commit when we finish the async with conn.transaction() block, rollback if exception arises before then
        async with conn.transaction():
            # Control the join order carefully, because the postgres optimizer seems to plan queries poorly when a lot of compound_property joins are included (e.g. 100x slower with 3+ property_names)
            # The critical concept is that the from_construct and to_construct should be filtered, by constant/rule_smiles structures, before they are joined together
            setparam = await conn.fetch("SET join_collapse_limit=1")
            setparam = await conn.fetch("SET from_collapse_limit=1")

            new_query_id = await conn.fetch("INSERT INTO query (inserted_at) VALUES (clock_timestamp()) RETURNING id")
            new_query_id = new_query_id[0][0]
            parsed_query["query_id"] = new_query_id

            statement, mol_pkls = search_algorithm.query_v3_statement(parsed_query)
            logging.info('Beginning query execution')
            logging.info(statement)
            logging.info(mol_pkls)
            t1 = time.perf_counter()
            await conn.execute(statement, *mol_pkls)
            t2 = time.perf_counter()

        logging.info("Time in seconds for search and saving results: " + str(round(t2-t1, 6)))

        num_pairs = await conn.fetch(f"SELECT count(*) FROM query_result WHERE query_result.query_id = {new_query_id}")
        logging.info(f"Number of pairs: {num_pairs[0][0]}")

    finally:
        await conn.close()

    if num_pairs[0][0] == 0:
        observations["no_results"] = "True"
        return json.dumps(
            {
                "url": "fatal",
                "observations": observations,
            }
        )

    return json.dumps({'query_id': new_query_id})


@app.post("/start_query")
async def start_query(query: QueryInput, background_tasks: BackgroundTasks, schema: str = 'public'):
    """
    Input structure / property constraints, and receive a query_id that will be used for getting results. The query is then initiated asynchronously.

    Then, proceed to /check_query to find out when the query is done.

    Body Parameters:

    When using this API in isolation from the matcher application, "**use default value**" is recommended for some parameters.

    - **snapquery_id**: **use default value**, matcher uses this to avoid duplicating query input state information for snapshotted queries
    - **query_id**: **use default value**, matcher uses this to directly load results and avoid rerunning queries, for snapshotted queries
    - **query_type**: "exact" or "substructure": specify what kind of search is run on the variable atoms (the atoms that change as part of the MMP transform)
    - **transform_order**: **use default value**, currently only "first_order" transform queries are implemented, "second_order" could be implemented in the future
    - **sketched_content**: These parameters refer to structures, and highlighted atoms, that would appear in sketcher(s) in the matcher frontend.
        - See documentation for **sketched_content** under the /validateSelection endpoint. **sketched_content** provides all structural constraints on returned transforms and MMPs. At least one of either mol1 or mol2 below must be filled out. For each mol, a molfile is required, and variable atom(s) indices are required. The environment atom indices are optional.
        - Not all atom/bond selections, as defined in **sketched_content**, will result in legal queries. The /validateSelection endpoint is designed to detect and correct atom/bond selections that will result in fruitless queries.
    - **OPTIONAL_properties**: comma-separated string of property_name values. Require that both compounds in each MMP result have at least one of these OPTIONAL_properties
    - **REQUIRED_properties**: comma-separated string of property_name values. Require that all compounds belonging to all MMPs in the results have every one of the REQUIRED_properties
        - At least one property must be provided, as either a REQUIRED or OPTIONAL property
    - **advanced_options**: (OPTIONAL):
        - **variable_min_heavies**: integer, minimum number of heavy (non-H) atoms in BOTH variable fragments that comprise a transform, for every transform found by the query
        - **variable_max_heavies**: integer, maximum number of heavy (non-H) atoms in BOTH variable fragments that comprise a transform, for every transform found by the query
        - **compound_min_heavies**: integer, minimum number of heavy (non-H) atoms in BOTH compounds, in every MMP found by the query
        - **compound_max_heavies**: integer, maximum number of heavy (non-H) atoms in BOTH compounds, in every MMP found by the query
        - **aggregation_type**: **use default value**: matcher telegraphs this parameter to subsequent endpoints that control how results are displayed in the frontend (starting with /aggregate_transforms), but this parameter should not affect the initial query itself. Accepted values are "individual_transforms" or "group_by_fragment"
        - **snapfilter_id**: **use default value**, can be used to avoid duplication of output filters for snapshotted queries
        - **snapfilter_string**: **use default value**, used to filter output for snapshotted queries

    Returns **query_id**, a positive integer
    """

    parsed_query = parse_query_input(query)
    observations = parsed_query["observations"]
    
    # Prevent query if it will be WAY too general, returning many many GB: only happens in extreme cases
    if parsed_query['constant'] is None and parsed_query['variable1'] is None:
        if "variable_selection_over_limit" not in observations:
            observations["query_too_general"] = 'no_variable_no_constant'
    if parsed_query['query_type'] == 'substructure':
        if parsed_query['constant'] is None or 'has_one_atom_constant' in observations:
            if parsed_query['variable1'] is None or 'has_one_atom_variable' in observations:
                observations["query_too_general"] = 'many_to_many_too_general'

    if ("all_atoms_selected_as_variable" in observations or "variable_selection_over_limit" in observations or "query_too_general" in observations \
            or "missing_constant" in observations):
        return json.dumps({
            "url": "fatal",
            "observations": observations,
        })
    observations["no_results"] = "False"

    conn = await get_matcher_conn(schema=schema)
    try:
        async with conn.transaction():
            new_query_id = await conn.fetch("INSERT INTO query (inserted_at) VALUES (clock_timestamp()) RETURNING id")
            new_query_id = new_query_id[0][0]
    finally:
        await conn.close()

    parsed_query["query_id"] = new_query_id
    background_tasks.add_task(run_query, parsed_query, schema)
    return json.dumps({'query_id': new_query_id})


async def run_query(parsed_query, schema):

    req, opt = parsed_query["REQUIRED_properties"], parsed_query["OPTIONAL_properties"]
    num_props = len(req) + len(opt)
    new_query_id = parsed_query["query_id"]
    # If num_pairs==1 when we enter the finally block, then we know an exception occurred, and can communicate that back to the client
    num_pairs = -1

    conn = await get_matcher_conn(schema=schema)
    # Make sure we close the connection even if an exception is raised
    try:
        # Commit when we finish the async with conn.transaction() block, rollback if exception arises before then
        # We need to do the initial insert, and all follow-up editing, within the same transaction so that client does not try to aggregate an intermediate result
        async with conn.transaction():
            # Control the join order carefully, because the postgres optimizer seems to plan queries poorly when a lot of compound_property joins are included (e.g. 100x slower with 3+ property_names)
            # The critical concept is that the from_construct and to_construct should be filtered, by constant/rule_smiles structures, before they are joined together
            setparam = await conn.fetch("SET join_collapse_limit=1")
            setparam = await conn.fetch("SET from_collapse_limit=1")

            find_and_save_pairs_statement, mol_pkls, parsed_query = search_algorithm.query_v3_statement(parsed_query)
            logging.info('Beginning query execution')
            logging.info(find_and_save_pairs_statement)
            logging.info(mol_pkls)
            t1 = time.perf_counter()
            await conn.execute(find_and_save_pairs_statement, *mol_pkls)
            num_pairs = await conn.fetch(f"SELECT count(*) FROM query_result WHERE query_result.query_id = {new_query_id}")
            num_pairs = num_pairs[0][0]
            logging.info(f"Number of pairs: {num_pairs}")

            if parsed_query['insert_environment_smarts'] == True and num_pairs > 0:
                # Only relevant for multicut queries where the user specified a constant environment
                # Consolidate chemically identical transformations that (due to attachment point semantics) have nonidentical from/to smiles and environment_smarts

                unglued_ingredients = await conn.fetch(search_algorithm.select_variables_and_env(new_query_id))
                update_with_glued_fragments_statement = search_algorithm.glue_variables_and_env(unglued_ingredients)
                await conn.execute(update_with_glued_fragments_statement)
            t2 = time.perf_counter()

        logging.info("Time in seconds for search and saving results: " + str(round(t2-t1, 6)))

        # If no pairs were found, insert a dummy row to signal to the check_query endpoint that the query completed with no pairs found
        # The row being inserted should not correspond to any real matched pair that would be returned from a legit completed search result
        if num_pairs == 0:
            async with conn.transaction():
                await conn.execute(f"INSERT INTO query_result (query_id, rule_id, use_original_direction, from_construct_id, to_construct_id) VALUES ({new_query_id}, 0, TRUE, 1, 1)")
    finally:
        if num_pairs == -1:
            # This means an exception occurred in the try block; therefore we insert a signal that can be interpreted by the check_query endpoint that is typically called with an interval polling function from frontend
            async with conn.transaction():
                await conn.execute(f"INSERT INTO query_result (query_id, rule_id, use_original_direction, from_construct_id, to_construct_id) VALUES ({new_query_id}, 0, FALSE, 1, 1)")

        await conn.close()
    return


@app.get("/check_query/{query_id}")
async def check_query(query_id: int, schema: str = 'public'):
    """
    Monitor query progress.

    Parameters:

    - **query_id**: The query_id returned by /start_query

    Returns JSON with following items:
    
    - **finished**: will evaluate to "True" when the query is finished, otherwise "False"
    - **finished_with_no_results**: will evaluate to "True" if no results were found and saved in DB by the query, otherwise "False"
    """

    conn = await get_matcher_conn(schema=schema)
    try:
        sample_rows = await conn.fetch(f"SELECT rule_id, use_original_direction, from_construct_id, to_construct_id FROM query_result WHERE query_id={query_id} LIMIT 1")
    finally:
        await conn.close()

    if len(sample_rows) == 0:
        # We use str and not bool to make comparisons consistent in both python and JS layers
        query_status = {'finished': "False"}
    else:
        query_status = {'finished': "True"}
        if list(sample_rows[0]) == [0, True, 1, 1]:
            query_status['finished_with_no_results'] = "True"
        elif list(sample_rows[0]) == [0, False, 1, 1]:
            query_status['finished_with_exception'] = "True"
        else:
            query_status['finished_with_no_results'] = "False"

    return json.dumps(query_status)


@app.post("/get_all_raw_data")
async def get_all_raw_data(inputs: GetAllRawData, schema: str = 'public'):
    """
    Retrieve all results for a specific query, after query has finished as indicated by /check_query.

    Parameters:

    - **query_id**: Number returned by /start_query endpoint to run the query of interest
    - **query**: The same data that was passed to /start_query to generate above query_id

    Returns JSON with items:
    - **column_headers**: Array of strings (headers) associated with data in **rows**, can be used as columns for a table: smiles, compound_ids, property data, etc.
    - **rows**: Array of arrays (i.e., rows), whose elements match the above column_headers, where each row represents a unique MMP.
    """

    req, opt = parse_props(inputs.query.REQUIRED_properties, inputs.query.OPTIONAL_properties)
    aggregation_type = inputs.query.advanced_options.aggregation_type

    statement = search_algorithm.get_all_raw_data_statement(inputs, req, opt)
    logging.info(statement)

    conn = await get_matcher_conn(schema=schema)

    try:
        num_pairs = await conn.fetch(f"SELECT count(*) FROM query_result WHERE query_result.query_id = {inputs.query_id}")
        num_pairs = num_pairs[0][0]
        if num_pairs < 200000:
            async with conn.transaction():
                setparam = await conn.fetch("SET join_collapse_limit=100")
                setparam = await conn.fetch("SET from_collapse_limit=100")

                logging.info('Fetching all raw data')
                t1 = time.perf_counter()
                all_raw_data = await conn.fetch(statement)
                t2 = time.perf_counter()
            logging.info("Time in seconds for query and download: " + str(round(t2-t1, 6)))
        else:
            logging.info("Number of pairs exceeds 200000, not downloading due to risk of crashing API server and frontend server")
    finally:
        await conn.close()

    """
    # Optionally prevent network/server overload
    if num_pairs >= 200000:
        return {
            "column_headers": ['Number of pairs exceeded 200000', 'Did not download due to risk of crashing servers'],
            "rows": ["Please contact application owner for alternative means of download", '']
        }
    """
    
    rows = [list(row.values()) for row in all_raw_data]
    if aggregation_type == 'many_to_many':
        rows += search_algorithm.flip_pairs(rows)

    output = {
        "column_headers" : list(all_raw_data[0].keys()),
        "rows" : rows
    }

    return output


@app.post("/aggregate_transforms")
async def aggregate_transforms(filters: AggregationParameters, schema: str = 'public'):
    """
    Obtain MMP transform-level data, including statistics, e.g. median change in a specified property for each MMP transform found by the query. The matcher frontend uses this data to populate the table of transforms in the frontend results.

    Body Parameters:

    - **query_id**: Number returned by /start_query endpoint to run the query of interest
    - **aggregation_type**: "individual_transforms" or "group_by_fragment": "individual_transforms" is the default and will return statistics on a per-transform basis. Alternatively, "group_by_fragment" is useful for queries where mol1_variable and mol2_variable fragments are similar or identical, for example searching for "privileged" fragments with a specific substructure. Each row returned by "group_by_fragment" will show how a specific fragment compares to all other fragments found in the query to which that fragment is related.
    - **statistics**: list of statistics of interest, with below parameters:
        - **statistic**: currently only "median" is implemented, and is implemented at the database level for best performance. For example, obtain the median change in a specific property for all transforms found by the query.
        - **property_name**: must match one of the property_names in the database property_name table, which were specified with data that was input to the mmpdb loadprops command during database initialization
        - **change_type**: "fold_change" or "delta", "fold_change" returns a quotient (B / A), "delta" returns a difference (B - A), where B is the property value of the right-hand-side compound in the MMP, and A is the property value of the left-hand-side compound in the MMP.
        - **base**: OPTIONAL: "raw" or "log": "raw" is the default and should be used, unless metadata was provided to the mmpdb loadprops command during database initialization. In the latter case, this parameter enables customization of what data statistics are performed upon. Example: a property is stored in the DB as "log" units (as defined in the metadata), but we wish to perform statistics on the non-log, "raw" units, then specify "raw". If no metadata was provided, "raw" will use the raw data that was loaded with mmpdb loadprops.
        - **units**: OPTIONAL: omit this field to use data as-is. Other options are "M" or "uM". Similar to **base** above, if metadata was provided, this parameter can be used to transform values e.g. from molar values stored in the DB, to micromolar values that is more digestible by certain users.
    - **range_filters**: OPTIONAL list of range_filter objects containing the below items, used to filter MMPs prior to statistical aggregation, based on the below criteria, e.g. every 'A' (i.e. left-hand-side) compound in an MMP must have a property_value <= x
        - **compound**: either 'A' or 'B', referring to LHS or RHS of the MMP
        - **property_name**: same description as property_name above
        - **operator**: '<', '>', '<=', '>=', '=', '!='
        - **base**: OPTIONAL: same description as base above
        - **units**: OPTIONAL: same description as units above
        - **value**: value to filter by

    Returns JSON with following items:

    - **minmax**: minimum and maximum values for properties and their statistics.
    - **column_headers**: Array of strings (i.e., headers) matching elements in below **rows**. The specific headers depend on **aggregation_type**. For example, when **aggregation_type** = "individual_transforms", below are the **column_headers** (details about what these headers refer to is provided):
        - **"rule_id"**: id for the transform that defines the difference between the two compounds in MMPs
        - **"from_smiles"**: SMILES for the left-hand-side fragment in the transform
        - **"to_smiles"**: SMILES for the right-hand-side fragment in the transform
        - **"pair_count"**: number of MMPs belonging to this transform
        - **label of the form statistic + property_name + change_type**
    - **rows**: array of arrays (i.e., rows), with elements matching column_headers, where each row correspond to either a transform (when **aggregation_type** = "individual_transforms"), or a fragment (when **aggregation_type** = "group_by_fragment") 
    - **grouped_by_environment**: false or true: whether or not the transforms are grouped by environment. Use this value as input to /get_plot_data endpoint. Only relevant for some queries with multiple points of attachment between variable and environment atoms, depending on symmetry
    """

    conn = await get_matcher_conn(schema=schema)

    try:
        use_environment = await conn.fetch(f"SELECT count(query_result.from_smiles_env) AS count FROM query_result WHERE query_result.query_id = {filters.query_id} AND query_result.from_smiles_env IS NOT NULL")
        use_environment = True if use_environment[0]['count'] > 0 else False
        property_metadata = await get_property_metadata(schema=schema, conn=conn, close_conn=False)
        statement = search_algorithm.aggregate_transforms_v2(filters, use_environment, property_metadata=property_metadata)
        logging.info('aggregation statement:')
        logging.info(statement)

        # async with conn.transaction():
        setparam = await conn.fetch("SET join_collapse_limit=100")
        setparam = await conn.fetch("SET from_collapse_limit=100")

        logging.info('Beginning aggregation')
        t1 = time.perf_counter()
        table_data = await conn.fetch(statement)
        t2 = time.perf_counter()
        logging.info("Time in seconds for query and download: " + str(round(t2-t1, 6)))

    finally:
        await conn.close()

    if filters.aggregation_type == 'individual_transforms':
        # When we include environment in from/to smiles, we don't return rule_id as it isn't needed
        non_stat_fields = 3 if use_environment else 4
        num_table_fields = non_stat_fields + len(filters.statistics)
    elif filters.aggregation_type == 'group_by_fragment':
        num_table_fields = 4 + len(filters.statistics)
        # This is the case for group_by_fragment aggregation without environment constraints
        if 'rule_id_array' in table_data[0].keys():
            num_table_fields += 1

    minmax_headers = list(table_data[0].keys())[num_table_fields : ]
    minmax_data = list(table_data[0].values())[num_table_fields : ]
    column_headers = list(table_data[0].keys())[0: num_table_fields]
    rows = [list(row.values())[0: num_table_fields] for row in table_data[1:]]

    output = {
        'minmax': dict(zip(minmax_headers, minmax_data)),
        'column_headers': column_headers,
        'rows': rows,
        'grouped_by_environment': use_environment,
    }

    return output


@app.post("/get_plot_data")
async def get_plot_data(filters: PlotParameters, schema: str = 'public'):
    """
    Get data for all MMPs associated with specific transform(s) or variable fragment(s).
    
    Before calling this endpoint, first obtain output from /aggregate_transforms.
    
    As input for this endpoint, the necessary ids for the transforms are provided by /aggregate_transforms. The matcher frontend uses this data to generate a scatterplot, e.g. change in property A vs. change in property B, where each point is an MMP belonging to specified transform(s).

    Body Parameters:

    - **query_id**: Number returned by /start_query endpoint to run the query of interest
    - **aggregation_type**: "individual_transforms" or "group_by_fragment": use the same value that was input to the /aggregate_transforms endpoint. see the description provided with /aggregate_transforms endpoint. This parameter affects the ID paradigm used for organizing groups of MMPs ("ids" parameter below)
    - **grouped_by_environment**: false or true: use the value returned by the /aggregate_transforms endpoint.
    - **ids**: list of ids for transforms or fragments: the format of the id depends on **aggregation_type** and **grouped_by_environment**:
        - if **aggregation_type** = "individual_transforms" and **grouped_by_environment** = false, then each id is an integer rule_id returned by /aggregate_transforms
        - if **aggregation_type** = "individual_transforms" and **grouped_by_environment** = true, then each id is a string with the form: from_smiles_env + "_" + to_smiles_env, where from_smiles_env and to_smiles_env are values returned by /aggregate_transforms
        - if **aggregation_type** = "group_by_fragment" and **grouped_by_environment** = false, then each id is an integer rule_id, from the rule_id_array (all rule_id associated with a fragment), returned by /aggregate_transforms
        - if **aggregation_type** = "group_by_fragment" and **grouped_by_environment** = true, then each id is a string of a variable fragment plus its environment: to_smiles_env, returned by /aggregate_transforms
    - **range_filters**: OPTIONAL list of range_filter objects, please refer to the documentation for range_filters under /aggregate_transforms endpoint

    **x_data** and **y_data**: This is just plain data about the MMPs. The reason we are calling the data "x_data" and "y_data" is because this endpoint was originally designed to get data for 2d graphs. This just means that you can get two types of data for MMPs with each call to this endpoint. x_data and y_data arguments can be identical.

    **x_data** and **y_data** are JSON objects with the following parameters:

    - **property_name**: This property_name should match one of the property_names that was passed in the previous call to /aggregate_transforms.
    - **change_type**: "fold_change", "delta" are two of the options, which respectively are B / A and B - A, where B is the property value of the RHS compound of the MMP, and A is the property value of the LHS compound of the MMP. Additional options (which do not actually represent changes in properties from one MMP to another) are "A" and "B", for the LHS and RHS compound property values, respectively, and "average", the average of A and B.
    - **base**: OPTIONAL: see description for base in /aggregate_transforms
    - **units**: OPTIONAL: see description for units in /aggregate_transforms

    Returns JSON with items:

    - **column_headers**: Array of following strings (details about what these labels refer to is provided):
        - **"from_construct_id"**: id for the left-hand-side compound in the MMP (more technically, id for the data structure that ties the compound to its variable/constant fragments)
        - **"to_construct_id"**: id for the right-hand-side compound in the MMP
        - **"rule_id"**: id for the transform that defines the difference between the two compounds in the MMP
        - Depending on the above input parameters "aggregation_type" and "grouped_by_environment", additional ids associated with such descriptors could be included
        - **name for x_axis generated based on choice of property_name and change_type**
        - **name for y_axis generated based on choice of property_name and change_type**
    - **rows**: array of arrays (i.e., rows) with elements matching the above column headers, one row per MMP
    """

    conn = await get_matcher_conn(schema=schema)

    property_metadata = await get_property_metadata(schema=schema, conn=conn, close_conn=False)
    statement = search_algorithm.get_plot_data_statement(filters, property_metadata=property_metadata)
    logging.info(statement)

    try:
        async with conn.transaction():
            setparam = await conn.fetch("SET join_collapse_limit=100")
            setparam = await conn.fetch("SET from_collapse_limit=100")

            logging.info('Fetching plot data')
            t1 = time.perf_counter()
            plot_data = await conn.fetch(statement)
            t2 = time.perf_counter()
        logging.info("Time in seconds for query and download: " + str(round(t2-t1, 6)))

    finally:
        await conn.close()

    plot_data = {
        "column_headers": list(plot_data[0].keys()),
        "rows": [list(row.values()) for row in plot_data]
    }

    return plot_data


@app.post("/get_pair_data")
async def get_pair_data(pairs: PairsCondensed, schema: str = 'public'):
    """
    Get data for specific MMPs

    Body Parameters:

    - **prop_changes**: Defines what property data is retrieved for each MMP. List of JSON objects with below parameters:
        - **property_name**: This property_name should match one of the property_names that was passed in the previous calls to /aggregate_transforms and /get_plot_data.
        - **change_type**: Either "fold_change" or "delta", see description for change type under /get_plot_data
        - **base**: OPTIONAL: see description for base in /aggregate_transforms
        - **units**: OPTIONAL: see description for units in /aggregate_transforms
    - **construct_id_pairs**: Array of arrays that define the MMPs: [from_construct_id, to_construct_id], which were returned by a previous call to /get_plot_data.

    Returns JSON with items:

    - **column_headers**: Array of the following strings labels (details about what these labels refer to is provided):
        - **"a_id"**: original identifier for compound "A" (left-hand-side compound of MMP), e.g. a ChEMBL number
        - **"a"**: SMILES for compound "A"
        - **"b_id"**: original identifier for compound "B" (right-hand-side compound of MMP), e.g. a ChEMBL number
        - **"b"**: SMILES for compound "B"
        - **"constant"**: SMILES for the unchanging common substructure of both compounds in the MMP, i.e. the portion of the MMP compounds' structures that is not affected by the transform
        - **label with the form a_property_name, depending on property name**: property value for compound A
        - **label with the form b_property_name, depending on property name**: property value for compound B
        - **label with the form property_name + "_" + change_type**: either a fold-change, or delta, of B's property value relative to A's property value
    - **rows**: array of arrays (i.e., rows) with elements matching the above column headers, one row per MMP
    """

    conn = await get_matcher_conn(schema=schema)

    property_metadata = await get_property_metadata(schema=schema, conn=conn, close_conn=False)
    statement = search_algorithm.get_pair_data_statement(pairs, property_metadata=property_metadata)
    logging.info(statement)

    try:
        async with conn.transaction():
            setparam = await conn.fetch("SET join_collapse_limit=100")
            setparam = await conn.fetch("SET from_collapse_limit=100")

            logging.info('Fetching pair data')
            t1 = time.perf_counter()
            pair_data = await conn.fetch(statement)
            t2 = time.perf_counter()
        logging.info("Time in seconds for query and download: " + str(round(t2-t1, 6)))

    finally:
        await conn.close()

    pair_data = {
        "column_headers": list(pair_data[0].keys()),
        "rows": [list(row.values()) for row in pair_data]
    }

    return pair_data


@app.post("/enumerate")
async def enumerate_designs(data: EnumerationData, schema: str = 'public'):
    """
    Apply specific transforms from query results to the query input structures, to enumerate new designs based on transforms of interest.

    Currently depends on data passed from the table of transform results displayed in the frontend.
    
    In brief, the necessary transform data for this endpoint can be obtained from the /aggregate_transforms endpoint.
    """

    mol1_molfile, mol2_molfile = data.query.sketched_content.mol1_molfile, data.query.sketched_content.mol2_molfile
    assert mol1_molfile != '' or mol2_molfile != ''
    molfile = mol1_molfile if mol1_molfile else mol2_molfile
    variable_atoms = data.query.sketched_content.mol1_variable_atoms if mol1_molfile else data.query.sketched_content.mol2_variable_atoms
    environment_atoms = data.query.sketched_content.mol1_environment_atoms if mol1_molfile else data.query.sketched_content.mol2_environment_atoms

    aggregation_type = data.query.advanced_options.aggregation_type
    mapped_molfile = ss_select.attach_maps(molfile, variable_atoms, environment_atoms)
    if ss_select.wildcards_in_molecule(mapped_molfile):
        return 'wildcards_in_core_error'
    variable, core = ss_select.extract_variable_and_enumeration_core(mapped_molfile)
    num_cuts = ss_select.get_num_cuts(variable)
    assert num_cuts > 0 and num_cuts < 4 and len(data.row_data) > 0

    # Order and uniquify the ids so that we can output them in the same order that the user viewed them in the table
    transform_ids = OrderedDict((row['id'], None) for row in data.row_data)
    use_environment = True
    if '[*:1]' in str(data.row_data[0]['id']) or type(data.row_data[0]['id']) == type(123):
        # For single-cut many_to_many queries, and multicut many_to_many queries with no user-specified environment, ID will be to_smiles containing at least [*:1]
        # For single-cut one_to_all queries, and multicut one_to_all queries with no user-specified environment, ID will be integer rule_id encoding from_smiles >> to_smiles
        #   (and the rule_id will be negative if we are displaying the transformation in the opposite direction to how it is encoded in the DB)
        use_environment = False
        # Otherwise the id will be a single glued-together variable/environment smiles (many_to_many), or two such smiles (from and to) with a _ inbetween

    if aggregation_type == 'group_by_fragment' and use_environment == False:
        if num_cuts == 1:
            # The most trivial enumeration case, we already have the 'to' fragment as the id, ready to be glued onto the core
            transform_data = [{'transform_id': transform_id, 'from_smiles': None, 'to_smiles': transform_id, 'environment_smarts': None} for transform_id in transform_ids]
        else:
            transform_data = [{'transform_id': transform_id, 'from_smiles': None, 'to_smiles': transform_id, 'environment_smarts' : None} for transform_id in transform_ids]
    else:
        # We need to lookup from/to smiles and potentially environments
        conn = await get_matcher_conn(schema=schema)
        try:
            setparam = await conn.fetch("SET join_collapse_limit=100")
            setparam = await conn.fetch("SET from_collapse_limit=100")
            t1 = time.perf_counter()
            transform_data = await conn.fetch(search_algorithm.ids_to_transforms(transform_ids, data.query_id, aggregation_type, use_environment))
            t2 = time.perf_counter()
        finally:
            await conn.close()

    enumerations = ss_select.apply_transforms_to_input(variable, core, transform_data, aggregation_type, num_cuts)

    # Prevent awkward drawings
    rdDepictor.SetPreferCoordGen(True)

    # Cleanup the core for substructure matching to product mols
    alignment_template = copy.deepcopy(core)
    for atom in alignment_template.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetIsotope(0)
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1)
    alignment_template = Chem.RemoveHs(alignment_template)

    for transform_id in enumerations:
        for product in enumerations[transform_id]['products']:
            # Modifies product in place, to align the atom coordinates to the template
            # The user doesn't want the products to be drawn in unfamiliar orientations.
            # Therefore align the products as best as possible to the input they drew themselves in the sketcher
            if product != 'error':
                align_mols(product, alignment_template, pattern_format='rdkit_mol', match_whole_pattern=False)

    sdf_string = io.StringIO()
    writer = Chem.rdmolfiles.SDWriter(sdf_string)

    # Attach data from the data table to the enumerated molecules, so the data will show up next to each enumerated molecule wherever we load it
    keys = list(key for key in data.row_data[0].keys() if 'transforms' in key or 'median' in key or key in ('pairs')) + ['Data Link', 'Notes']

    # Start the output with the input query molecule, for the user's reference
    whole_molecule = Chem.MolFromMolBlock(molfile)
    none_list = [whole_molecule.SetProp(str(key), '') for key in keys]
    whole_molecule.SetProp('Data Link', data.link_address)
    whole_molecule.SetProp('Notes', 'Starting Molecule')
    whole_molecule.SetProp('ID', 'Matcher-0')
    writer.write(whole_molecule)

    for row in data.row_data:
        transform_id = row['id']
        row['Data Link'] = data.link_address
        row['Notes'] = ''
        if transform_id in enumerations:
            for key in keys:
                for product in enumerations[transform_id]['products']:
                    if product != 'error':
                        product.SetProp(str(key), str(row[key]))

    matcher_id = 1
    # Write the molecules in the same order as the transformations in the user's table
    for transform_id in transform_ids:
        row = enumerations.get(transform_id)
        if row is not None:
            for product in row['products']:
                if product != 'error':
                    product.SetProp('ID', f"Matcher-{matcher_id}")
                    matcher_id += 1
                    writer.write(product)
                else:
                    product = Chem.MolFromSmiles('OOPS')
                    none_list = [product.SetProp(str(key), '') for key in keys]
                    product.SetProp('Data Link', data.link_address)
                    product.SetProp('Notes', 'Error in molecule construction')
                    product.SetProp('ID', f"Matcher-{matcher_id}")
                    matcher_id += 1
                    writer.write(product)

    writer.close()
    return sdf_string.getvalue()


@app.get("/snap_read/{snapshot_id}", response_class=ORJSONResponse)
async def snap_read(snapshot_id: str, schema: str = 'public'):
    """
    Read snapshot from database. The snapshot refers to input state and output filter state for a query. With this data, we can recreate input + output for a query of interest.

    This way, users can save and share results simply by passing around a unique link containing the {snapshot_id} (in combination with functionality provided by the frontend API)

    The snapshot would have previously been written using the snap_write endpoint.
    """

    try:
        conn = await get_matcher_conn(schema=schema)

        statement = f"""
SELECT 
    snapquery.id, snapquery.version_id, snapquery.query_id, snapquery.query_type, snapquery.transform_order,
    snapquery.REQUIRED_properties, snapquery.OPTIONAL_properties, snapfilter.id, snapfilter.snapfilter_string,

    mol1_sketched_content.molfile, mol1_sketched_content.variable_atoms, mol1_sketched_content.variable_bonds,
    mol1_sketched_content.environment_atoms, mol1_sketched_content.environment_bonds,
    mol2_sketched_content.molfile, mol2_sketched_content.variable_atoms, mol2_sketched_content.variable_bonds,
    mol2_sketched_content.environment_atoms, mol2_sketched_content.environment_bonds,

    snapquery.variable_min_heavies, snapquery.variable_max_heavies, snapquery.compound_min_heavies, snapquery.compound_max_heavies, snapquery.aggregation_type
FROM
    snapshot
    INNER JOIN snapquery ON snapquery.id = snapshot.snapquery_id
    LEFT OUTER JOIN sketched_content mol1_sketched_content ON mol1_sketched_content.id = snapquery.mol1_sketched_content_id 
    LEFT OUTER JOIN sketched_content mol2_sketched_content ON mol2_sketched_content.id = snapquery.mol2_sketched_content_id 
    INNER JOIN snapfilter ON snapfilter.id = snapshot.snapfilter_id
WHERE
    snapshot.id = {snapshot_id}
"""

        snapshot_data = await conn.fetch(statement)
        snapshot_data = list(snapshot_data[0])
    finally:
        await conn.close()

    for i, value in enumerate(snapshot_data):
        if value is None:
            snapshot_data[i] = ''

    snapshot_data_keys = (
        "snapshot_id", "snapquery_id", "version_id", "query_id", "query_type", "transform_order", "REQUIRED_properties", "OPTIONAL_properties", "snapfilter_id", "snapfilter_string",
        "mol1_molfile", "mol1_variable_atoms", "mol1_variable_bonds", "mol1_environment_atoms", "mol1_environment_bonds",
        "mol2_molfile", "mol2_variable_atoms", "mol2_variable_bonds", "mol2_environment_atoms", "mol2_environment_bonds",
        "variable_min_heavies", "variable_max_heavies", "compound_min_heavies", "compound_max_heavies", "aggregation_type"
    )

    snapshot = dict(zip(snapshot_data_keys, [snapshot_id] + snapshot_data))

    return snapshot


@app.post("/snap_write", response_class=ORJSONResponse)
async def snap_write(s: QueryInput, schema: str = 'public'):
    """
    Write snapshot to database. The snapshot refers to input state and output filter state for a query. With this data, we can recreate input + output for a query of interest.

    This way, users can save and share results simply by passing around a unique link containing the {snapshot_id} (in combination with functionality provided by the frontend API)
    """

    conn = await get_matcher_conn(schema)

    try:
        async with conn.transaction():

            setparam = await conn.fetch("SET standard_conforming_strings = OFF")

            # arbitrarily setting version_id = 1 for the initial implementation of this feature
            version_id = 1

            # If we didn't get to the output results from a previous snapshot, then we need to store the input query state, otherwise we reuse the input query state
            if s.snapquery_id == '':

                ids = {}
                for mol in ('mol1', 'mol2'):
                    skc = dict(s.sketched_content)
                    if skc[mol+'_molfile'] == '':
                        ids[mol] = 'NULL'
                    else:
                        data = str((skc[mol+'_molfile'], skc[mol+'_variable_atoms'], skc[mol+'_variable_bonds'], skc[mol+'_environment_atoms'], skc[mol+'_environment_bonds']))
                        statement = f"""
INSERT INTO sketched_content (molfile, variable_atoms, variable_bonds, environment_atoms, environment_bonds)
VALUES {data}
RETURNING id
"""
                        result = await conn.fetch(statement)
                        ids[mol] = result[0][0]

                snapquery_write_variables = str((version_id, s.query_type, s.query_id, s.transform_order, ids['mol1'], ids['mol2'],
                s.REQUIRED_properties, s.OPTIONAL_properties, s.advanced_options.variable_min_heavies, s.advanced_options.variable_max_heavies,
                s.advanced_options.compound_min_heavies, s.advanced_options.compound_max_heavies, s.advanced_options.aggregation_type)).replace("'NULL'", "NULL")

                statement = f"""
INSERT INTO snapquery (
    version_id, query_type, query_id, transform_order, mol1_sketched_content_id, mol2_sketched_content_id,
    REQUIRED_properties, OPTIONAL_properties, variable_min_heavies, variable_max_heavies,
    compound_min_heavies, compound_max_heavies, aggregation_type
)
VALUES {snapquery_write_variables}
RETURNING id
"""
                result = await conn.fetch(statement)
                snapquery_id = result[0][0]
            else:
                snapquery_id = int(s.snapquery_id)
    
            # Store the output results filters
            snapfilter_write_variables = str((version_id, s.snapfilter_string))

            statement = f"""
INSERT INTO snapfilter (version_id, snapfilter_string)
VALUES {snapfilter_write_variables}
RETURNING id
"""
            result = await conn.fetch(statement)
            snapfilter_id = result[0][0]

            # Combine the unique ids obtained by the above two queries, into one snapshot
            snapshot_write_variables = str((snapquery_id, snapfilter_id))

            statement = f"""
INSERT INTO snapshot (snapquery_id, snapfilter_id)
VALUES {snapshot_write_variables}
RETURNING id
"""
            result = await conn.fetch(statement)
            snapshot_id = result[0][0]

    finally:
        await conn.close()

    return {"new_snapshot_id": str(snapshot_id)}




@app.get("/bind_snapquery_to_results/{snapquery_id}")
async def bind_snapquery_to_results(snapquery_id: int, query_id: int, schema: str = 'public'):
    """
    When the same query is run multiple times, it is wasteful to repeat the actual query, both in terms of time and disk space.

    For queries saved with the snapshot feature, this endpoint allows us to avoid rerunning the query.

    Here we are just associating the snapshot data with the query_id of results in the database. When the frontend next sees a snapshot with this snapquery_id, the frontend will skip the query and instead load cached results.
    """

    conn = await get_matcher_conn(schema=schema)

    try:
        async with conn.transaction():
            await conn.execute(f"UPDATE snapquery SET query_id = {query_id} WHERE id = {snapquery_id}")
    finally:
        await conn.close()

    return {'success': True}

@app.get("/propertyNames/", response_class=ORJSONResponse)
async def propertyNames(schema: str = 'public'):
    """
    This endpoint is used to populate the properties menu in the frontend input form.

    Returns JSON with items:

    - **properties**: Array of objects with items:
        - **property_name**: string, name of property provided in input file to mmpdb loadprops command
        - **display_name**: string, name of property to be displayed to users in frontend (if different from property_name, as defined in property metadata)
    """

    conn = await get_matcher_conn(schema=schema)
    try:
        names = await conn.fetch("SELECT property_name.name, property_name.display_name FROM property_name")
    finally:
        await conn.close()

    props = {'properties': [dict(zip(('property_name', 'display_name'), row.values())) for row in names]}

    return props

@app.get("/propertyMetadata/", response_class=ORJSONResponse)
async def get_property_metadata(schema: str = 'public', conn=None, close_conn=True):
    """
    This function is needed for obtaining property metadata that was written to the database using the mmpdb loadprops command.

    This metadata is depended upon by some of this API's endpoints that handle compound property data.

    Returns JSON with the following (for each property saved in the property_name table in the database):

    - **property_name**: object with below items:
        - **base**: see /aggregate_transforms documentation for more details
        - **unit**: see /aggregate_transforms documentation for more details
        - **display_name**: alternate name of the property to display to users in frontend
        - **display_unit**: preferred unit to display to users in frontend
        - **change_displayed**: preferred change type, for compound property values of compounds in an MMP, to display to users in frontend
    """

    # When the frontend is calling this endpoint, conn should be None
    # When the backend is calling this endpoint, for current use cases, we already have a conn with a defined schema name, so avoid making another conn
    if conn is None:
        conn = await get_matcher_conn(schema=schema)
    try:
        results = await conn.fetch("SELECT name, base, unit, display_name, display_base, display_unit, change_displayed FROM property_name")
    finally:
        if close_conn == True:
            await conn.close()

    metadata = {row['name'] : dict(row) for row in results}
    return metadata

@app.get("/")
async def root():
    return {"message": "Greetings from the matcher backend API"}


@app.get("/smilesToImage")
async def smilesToImage(smiles: str, x: int, y: int, scaleImage: bool = False):
    """
    Utility endpoint for generating images in the frontend.

    Parameters:

    - **smiles**: SMILES string of structure from which to generate an image
    - **x**: integer width of image
    - **y**: integer height of image
    - **scaleImage**: boolean controlling whether we scale the size of the image based on how large the molecule is
    """
    return _smilesToImage(smiles, x, y, scaleImage)


def _smilesToImage(smiles: str, x: int, y: int, scaleImage=False):

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        # set sanitize=False so that funny-looking fragments (like partial aromatic rings) will not cause mol generation to fail
        mol = Chem.MolFromSmiles(smiles, sanitize=False)

    atoms = mol.GetAtoms()
    # Scale up image size up to a limit, so that structures aren't tiny
    if scaleImage:
        num_atoms = len(list(mol.GetAtoms()))
        if num_atoms > 10:
            scaling_factor = num_atoms / 10
            scaling_factor = 2 if scaling_factor > 2 else scaling_factor
            x = int(scaling_factor * x)
            y = int(scaling_factor * y)

    # Create a transparent background: no more white background in molecule images
    cairo_drawer = Draw.rdMolDraw2D.MolDraw2DCairo(x, y)
    drawing_options = cairo_drawer.drawOptions()
    drawing_options.clearBackground = False

    constant_indices = [atom.GetIdx() for atom in atoms if atom.GetAtomMapNum() == 7]
    remove_maps = [atom.SetAtomMapNum(0) for atom in atoms]
    constant_atom_colors = {i: (1, 0.7, 0.7) for i in constant_indices}
    atom_radii = {i: 0.8 for i in constant_indices}

    try:
        # This generates molecule images with transparent background
        cairo_drawer.DrawMolecule(mol, highlightAtoms=constant_indices, highlightAtomColors=constant_atom_colors, highlightAtomRadii=atom_radii)
        cairo_drawer.FinishDrawing()
        buffered = io.BytesIO(cairo_drawer.GetDrawingText())
    except:
        # Note: molecule images generated as below will have a white background
        img = Draw.MolToImage(mol, size=(x,y), highlightAtoms=constant_indices, highlightColor=(1, 0.7, 0.7))
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
    buffered.seek(0)
    return StreamingResponse(buffered, media_type="image/png")


def align_mols(mol, pattern, pattern_format='smiles', match_whole_pattern=True):
    # Modifies mol in place

    # If we are going to use atom highlights to show differences between mol and pattern in a depiction, then we set match_whole_pattern=True
    # But if we are just aligning mol to pattern, and don't need atom highlights, we don't care about the whole pattern match,
    #   just the match of mol to the largest fragment in pattern that we will align to
    if match_whole_pattern:
        assert pattern_format == 'smiles'
        # Remove dummy atoms, otherwise substructure matching fails
        pattern_mol = Chem.MolFromSmiles(pattern)
        for atom in pattern_mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(1)
        pattern_mol = Chem.RemoveHs(pattern_mol)

    # For template ("pattern") molecule with multiple fragments, they will mess up the alignment, so we need to pick only the largest fragment
    max_index = 0
    max_length = 0
    if pattern_format == 'smiles':
        pattern_frags = pattern.split(".")
        for i, frag in enumerate(pattern_frags):
            if len(frag) > max_length:
                max_length = len(frag)
                max_index = i
        largest_frag_mol = Chem.MolFromSmiles(pattern_frags[max_index])
        for atom in largest_frag_mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(1)
        largest_frag_mol = Chem.RemoveHs(largest_frag_mol)

    # Sometimes we want to align to a precise user-specified structure from a sketcher, for example aligning to a core mol provided by extract_variable_and_enumeration_core
    elif pattern_format == 'rdkit_mol':
        pattern_frags = Chem.GetMolFrags(pattern, asMols=True, sanitizeFrags=False)
        for i, frag in enumerate(pattern_frags):
            num_atoms = len(frag.GetAtoms())
            if num_atoms > max_length:
                max_length = num_atoms
                max_index = i
        largest_frag_mol = pattern_frags[max_index]

    for molecule in (mol, largest_frag_mol):
        if len(molecule.GetConformers()) == 0:
            rdDepictor.Compute2DCoords(molecule)

    if match_whole_pattern:
        match = mol.GetSubstructMatch(pattern_mol)
    largest_frag_match = mol.GetSubstructMatch(largest_frag_mol)
    if match_whole_pattern:
        # SSS isn't perfect. If no match, revert to unaligned, unhighlighted image
        if len(match) == 0 or len(largest_frag_match) == 0:
            return None
    elif len(largest_frag_match) == 0:
        return None
    TemplateAlign.AlignMolToTemplate2D(mol,largest_frag_mol,match=largest_frag_match,clearConfs=True)

    # The alignment often causes parts of the molecule to overlap/clash.
    # Address this by recomputing the 2d coords of all atoms, except the ones we wanted to align to
    # Might not be necessary when we use rdDepictor.SetPreferCoordGen(True)

    if not match_whole_pattern:
        return None
    else:
        # Return the substructure match which is used further in smilesToImage_aligned
        return match


@app.get("/smilesToImage_aligned")
async def smilesToImage_aligned(smiles: str, pattern, x: int, y: int):
    """
    Utility endpoint for generating images in the frontend.

    Parameters:

    - **smiles**: SMILES string of structure from which to generate an image
    - **pattern**: SMILES string of pattern to which we want to align the molecule image generated from **smiles** input
    - **x**: integer width of image
    - **y**: integer height of image
    """

    # This setting might help with preventing clashes/overlaps in the drawn structures
    rdDepictor.SetPreferCoordGen(True)
    # Create a transparent background: no more white background in molecule images
    mol = Chem.MolFromSmiles(smiles)

    # modifies mol in place to align 2D depiction to largest fragment in pattern, and returns substructure matching indices of pattern in mol
    match = align_mols(mol, pattern)
    if match is None:
        # SSS isn't perfect. If no match, revert to unaligned, unhighlighted image
        return _smilesToImage(smiles, x, y)

    fragment_indices = [x for x in range(len(mol.GetAtoms())) if x not in set(match)]

    fragment_atom_colors = {x: (0.7, 1, 0.7) for x in fragment_indices}

    ## Optionally increase highlight radii
    atom_radii = {x: 0.8 for x in fragment_indices}

    cairo_drawer = Draw.rdMolDraw2D.MolDraw2DCairo(x, y)
    drawing_options = cairo_drawer.drawOptions()
    drawing_options.clearBackground = False

    Draw.rdMolDraw2D.PrepareAndDrawMolecule(
        cairo_drawer, mol, highlightAtoms=fragment_indices, highlightAtomColors=fragment_atom_colors, highlightAtomRadii=atom_radii
    )

    buffered = io.BytesIO(cairo_drawer.GetDrawingText())
    buffered.seek(0)
    return StreamingResponse(buffered, media_type="image/png")

if __name__ == '__main__':
    uvicorn.run(app)
