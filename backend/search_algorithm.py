import pandas as pd
import time
from rdkit.Chem import AllChem as Chem
import copy
import logging
import ss_select
import smiles_syntax
import backend_api
logging.basicConfig(level=logging.INFO)


def query_v3_statement(query):

    req, opt = query["REQUIRED_properties"], query["OPTIONAL_properties"]
    num_props = len(req) + len(opt)
    query_type = query["query_type"]
    constant, variable1, variable2 = query["constant"], query["variable1"], query["variable2"]    

    assert variable1 or variable2

    query_id = query["query_id"]
    num_frags = 0

    num_frags_check_mol = variable1 if variable1 is not None else variable2
    for atom in num_frags_check_mol.GetAtoms():
        if atom.GetAtomicNum() == 37:
            # Rb atoms were inserted as dummies anywhere a cut was performed; we count # of Rb in a fragment to determine number of fragmentations: num_frags
            num_frags += 1
    if query_type == "exact":
        # Add explicit H so that SSS functions like an exact search
        # Conversion to binary and back seems necessary to prevent RDKit error when calling Chem.AddHs
        variable1 = Chem.AddHs(Chem.Mol(variable1.ToBinary())) if variable1 is not None else None
        variable2 = Chem.AddHs(Chem.Mol(variable2.ToBinary())) if variable2 is not None else None

        #variable1 = Chem.Mol(variable1.ToBinary())
        #variable1 = Chem.AddHs(variable1)

    variable1_smiles = Chem.MolToSmiles(variable1) if variable1 else None
    variable2_smiles = Chem.MolToSmiles(variable2) if variable2 else None

    #elif query_type == "many_to_many":
        #variable2 = variable1

    # A bunch of acrobatics is necessary in multicut queries
    # A multicut fragment is possibly only represented in DB in one of multiple possible ways, e.g. [*:1]CO[*:2] and/or [*:2]CO[*:1]
    # We need to query all possibilities, and make sure the cut points in variable map to the cut points in constant; that's what the below code does
    variable1_perms, variable2_perms, constant_perms = ss_select.permute_variables_constant(variable1=variable1, variable2=variable2, constant=constant)

    # variable_perms contains all possible variable permutations of variable1 and variable2, and constant_perms contains all possible constant permutations

    # variable_perms = variable1_perms + variable2_perms if variable2_perms else variable1_perms
    if variable1_perms and variable2_perms:
        variable_perms = variable1_perms + variable2_perms
    elif variable1_perms:
        variable_perms = variable1_perms
    elif variable2_perms:
        variable_perms = variable2_perms

    # Convert to binary for passing to postgres server
    variable_pkls = [mol.ToBinary() for mol in variable_perms]
    constant_pkls = [mol.ToBinary() for mol in constant_perms] if constant_perms else None

    variable_ids = ss_select.map_unique_permutations(variable_perms)
    if variable1 and variable2:
        variable1_ids = variable_ids[0 : int(len(variable_ids) / 2)]
        variable2_ids = variable_ids[int(len(variable_ids) / 2) :]
    elif variable1:
        variable1_ids = variable_ids
        variable2_ids = [None]*len(variable1_ids)
    elif variable2:
        variable2_ids = variable_ids
        variable1_ids = [None]*len(variable2_ids)

    constant_ids = ss_select.map_unique_permutations(constant_perms) if constant_perms else [None]*len(variable1_ids)

    # We need to record environment smiles for multicut transformations where a symmetrical variable fragment has an environment with less symmetry than the constant fragment
    # For now will just toggle this when there is any symmetry in the variable fragment, and when constant fragment is present; not sure if all of such cases require this feature
    # For example, ethylene as variable fragment (2-cut) being transformed to a methylated ethylene, with environment of O on one end and C on the other
    # Sometimes the same rule is used for 2 different MMPs where the O and C would be swapped on the ethylene, causing multiple rule/environment combinations to appear in the same cluster
    # If we track and include the user-defined environment_smarts in the GROUP BY along with rule_id, this should fix it, and enable accurate enumerations
    query['insert_environment_smarts'] = False
    if num_frags > 1 and constant_pkls is not None:
        query['insert_environment_smarts'] = True
        # We need smarts for proper representation of fragments, for use in substructure matching during compound enumeration
        constant_smarts_perms = [Chem.MolToSmarts(perm) for perm in constant_perms]
        # We need smiles for gluing to the variable fragments using smiles_syntax.convert_labeled_wildcards_to_closures
        constant_smiles_perms = get_constant_smiles_perms(constant_perms)
    """
    if num_frags == 2 and constant_pkls is not None:
        if len(set(variable1_ids)) < 2 or query_type == 'many_to_many':
            query['insert_environment_smarts'] = True
            constant_smiles_perms = get_constant_smiles_perms(constant_perms)
    if num_frags == 3 and constant_pkls is not None:
        if len(set(variable1_ids)) < 6 or query_type == 'many_to_many':
            query['insert_environment_smarts'] = True
            constant_smiles_perms = get_constant_smiles_perms(constant_perms)
    """

    # v1_c and v2_c have mappings between variable and constant, for all unique combos of variable1/constant and variable2/constant
    v1_c = set(zip(variable1_ids, constant_ids))
    v2_c = set(zip(variable2_ids, constant_ids))
    # For each pairing, we filter from_construct by v1/c, and to_construct by v2/c (and vice/versa in a unioned query with use_original_direction = False)
    pairings = []
    for v1, c1 in v1_c:
        for v2, c2 in v2_c:
            if c1 == c2:
                pairings.append({
                    'v1_id': v1,
                    'v2_id': v2,
                    'c_id': c1
                })

    ############################################################################################################################
    # Use common table expressions for filtering rule_smiles and constant_smiles
    # This helps with organization, but also performance, as the CTEs can be reused
    ############################################################################################################################

    CTE_expression = "WITH \n"
    CTE_definitions = []
    # We will insert molecule binaries into mol_pkls in the same order we add items to CTE_definitions
    mol_pkls = []
    CTE_template = """{}
AS (
    SELECT A.id AS id{}
    FROM {} A
    WHERE {}
    AND A.num_frags = {} {} {}
)"""

    # Here we are creating CTEs (1 per unique variable permutation) with rule_smiles filtered by variable
    for i in set(variable_ids):
        # We need to use a plpgsql procedure (below), to force the query planner to expect a certain number of rows to result from the substructure search
        # Currently the row estimate is set to 100 for rule_smiles (see last row of the procedure)
        # Otherwise, the planner tends to vastly underestimate how many rows will be returned, causing bad plans such as nested loop joins between large relations
        #  which can, for example, cause a 15 second query (with below code) to take 10+ minutes
        """
CREATE OR REPLACE FUNCTION filter_rule_smiles ( pkl IN BYTEA, nf IN INTEGER, min_heavies IN INTEGER DEFAULT 0, max_heavies IN INTEGER DEFAULT 0) RETURNS TABLE(id INTEGER)
AS $$
DECLARE
BEGIN
    IF (min_heavies > 0 AND max_heavies = 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf
            AND rule_smiles.num_heavies >= min_heavies;
    ELSIF (min_heavies > 0 AND max_heavies > 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf
            AND rule_smiles.num_heavies >= min_heavies
            AND rule_smiles.num_heavies <= max_heavies;
    ELSIF (min_heavies = 0 AND max_heavies > 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf
            AND rule_smiles.num_heavies <= max_heavies;
    ELSIF (min_heavies = 0 AND max_heavies = 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf;
    END IF;
END;
$$ LANGUAGE plpgsql
   ROWS 100
        """
        CTE_definitions.append(f"rule_smiles_{i} AS (SELECT id FROM filter_rule_smiles(%s, {num_frags}, {query['variable_min_heavies']}, {query['variable_max_heavies']}))")
        mol_pkls.append(variable_pkls[i])

    newline = "\n"
    construct_conditions = {}
    for direction in ('from', 'to'):
        construct_conditions[direction] = [
            f"{direction}_construct.num_frags = {num_frags}",
            f"{direction}_construct.rule_smiles_num_heavies >= {query['variable_min_heavies']}" if query['variable_min_heavies'] > 0 else "",
            f"{direction}_construct.rule_smiles_num_heavies <= {query['variable_max_heavies']}" if query['variable_max_heavies'] > 0 else "",
            f"{direction}_construct.compound_num_heavies >= {query['compound_min_heavies']}" if query['compound_min_heavies'] > 0 else "",
            f"{direction}_construct.compound_num_heavies <= {query['compound_max_heavies']}" if query['compound_max_heavies'] > 0 else "",
        ]
        construct_conditions[direction] = [c for c in construct_conditions[direction] if c != ""]

    if constant:
        for i in set(constant_ids):
            # We need to use a plpgsql procedure (below), to force the query planner to expect a certain number of rows to result from the substructure search
            # Currently the row estimate is set to 100,000 for constant_smiles (see last row of the procedure)
            # Otherwise, the planner tends to vastly underestimate how many rows will be returned, causing bad plans such as nested loop joins between large relations
            #  which can, for example, cause a 15 second query (with below code) to take 10+ minutes
            """
CREATE OR REPLACE FUNCTION filter_constant_smiles ( pkl IN BYTEA, nf IN INTEGER ) RETURNS TABLE(id integer)
AS $$
DECLARE
BEGIN
  RETURN QUERY
  SELECT constant_smiles.id AS id
  FROM constant_smiles
  WHERE constant_smiles.smiles_mol@>mol_from_pkl(pkl)
  AND constant_smiles.num_frags = nf;
END;
$$ LANGUAGE plpgsql
   ROWS 100000
            """
            CTE_definitions.append(f"constant_smiles_{i} AS (SELECT id FROM filter_constant_smiles(%s, {num_frags}))")
            mol_pkls.append(constant_pkls[i])
            for direction in ('from', 'to'):
                CTE_definitions.append(f"""{direction}_construct_{i}
AS (
    SELECT {direction}_construct.id AS id, {direction}_construct.rule_smiles_id AS rule_smiles_id, {direction}_construct.constant_id AS constant_id, {direction}_construct.compound_id AS compound_id
    FROM {direction}_construct 
    INNER JOIN constant_smiles_{i} ON constant_smiles_{i}.id = {direction}_construct.constant_id
    WHERE {(newline + '    AND ').join([c for c in construct_conditions[direction]])}
)""")

    ############################################################################################################################
    # Define the main query body
    ############################################################################################################################

    not_yet_unioned_statements = []

    for p in pairings:

        v1_id, v2_id, c_id = p['v1_id'], p['v2_id'], p['c_id']
        from_params = ("from", "to", "TRUE")
        to_params = ("to", "from", "FALSE")
        for direction, opposite_direction, use_original_direction in (from_params, to_params):
        
            if direction == "to" and variable1 and variable2 and (variable1_smiles == variable2_smiles) and query['aggregation_type'] == 'group_by_fragment':
                # No need to scan for pairs in the "TO" direction, because the "FROM" query is symmetrical
                # Only implement this de-duplication for 'group_by_fragment' queries, because the downstream data processing queries assume we only scanned the 'from' direction
                # Otherwise, the user may actually want to see both directions, even if redundant
                continue

            construct_1_rs_join, construct_2_rs_join = "", ""
            construct_1, construct_2 = f"{direction}_construct", f"{opposite_direction}_construct"
            construct_join_conditions = f"from_construct.constant_id = to_construct.constant_id AND from_construct.rule_smiles_id != to_construct.rule_smiles_id"
            construct_1_frags_condition = f"\nWHERE {direction}_construct.num_frags = {num_frags}"
            construct_2_frags_condition = f"\nWHERE {opposite_direction}_construct.num_frags = {num_frags}"
            construct_1_conditions, construct_2_conditions = construct_conditions[direction], construct_conditions[opposite_direction]
            if c_id is not None:
                construct_1 = f"{direction}_construct_{c_id} " + construct_1
                construct_2 = f"{opposite_direction}_construct_{c_id} " + construct_2
                # We already filtered constructs by these conditions in CTEs (see above), when joining to constant
                construct_1_conditions, construct_2_conditions = [], []
            if v1_id is not None:
                construct_1_rs_join = f"\n        INNER JOIN rule_smiles_{v1_id} ON rule_smiles_{v1_id}.id = {direction}_construct.rule_smiles_id"
                construct_1 = f"""(
        SELECT {direction}_construct.id AS id, {direction}_construct.rule_smiles_id AS rule_smiles_id, {direction}_construct.constant_id AS constant_id, {direction}_construct.compound_id AS compound_id
        FROM {construct_1} {construct_1_rs_join} {newline+'WHERE' if len(construct_1_conditions) > 0 else ''} {(newline + '    AND ').join([condition for condition in construct_1_conditions])}
    ) {direction}_construct"""

            if v2_id is not None:
                construct_2_rs_join = f"\n        INNER JOIN rule_smiles_{v2_id} ON rule_smiles_{v2_id}.id = {opposite_direction}_construct.rule_smiles_id"
                construct_2 = f"""(
        SELECT {opposite_direction}_construct.id AS id, {opposite_direction}_construct.rule_smiles_id AS rule_smiles_id, {opposite_direction}_construct.constant_id AS constant_id, {opposite_direction}_construct.compound_id AS compound_id
        FROM {construct_2} {construct_2_rs_join} {newline+'WHERE' if len(construct_2_conditions) > 0 else ''} {(newline + '    AND ').join([condition for condition in construct_2_conditions])}
    ) {opposite_direction}_construct"""

            env_smarts_perm = f", '{constant_smarts_perms[c_id]}' AS environment_smarts" if query['insert_environment_smarts'] == True else ""
            env_smiles_perm = f", '{constant_smiles_perms[c_id]}' AS environment_smiles" if query['insert_environment_smarts'] == True else ""
            block = f"""
    SELECT {query_id} AS query_id, {use_original_direction} AS use_original_direction, from_construct.id AS from_construct_id, to_construct.id AS to_construct_id,
        from_construct.compound_id AS from_compound_id, to_construct.compound_id AS to_compound_id,
        from_construct.rule_smiles_id AS from_smiles_id, to_construct.rule_smiles_id AS to_smiles_id{env_smarts_perm}{env_smiles_perm}
    FROM {construct_1}
    INNER JOIN {construct_2} ON {construct_join_conditions} \n"""

            if c_id is None:
                if v1_id is None or v2_id is None:
                    block += "\nWHERE "
                    if v1_id is None:
                        block += '\n    AND '.join([condition for condition in construct_1_conditions])
                    if v2_id is None:
                        block += '\n    AND '.join([condition for condition in construct_2_conditions])

            not_yet_unioned_statements.append(block)

    no_props = "\n\n    UNION ALL\n".join(not_yet_unioned_statements)
    CTE_definitions.append(f"no_props AS ({no_props}\n)")
    CTE_expression = "WITH \n" + ",\n".join(CTE_definitions) + "\n\n"
    env_columns = ", environment_smarts, from_smiles_env" if query['insert_environment_smarts'] == True else ""
    insertion = f"INSERT INTO query_result (query_id, rule_id, use_original_direction, from_construct_id, to_construct_id{env_columns})"
    env_columns_aliases = ", no_props.environment_smarts, no_props.environment_smiles" if query['insert_environment_smarts'] == True else ""

    not_yet_unioned_prop_selects = []
    opt_iterator = list(zip(opt, ["INNER JOIN"]*len(opt)))
    for prop, join in opt_iterator:
        prop_select = f"""
SELECT query_id, rule.id, use_original_direction, from_construct_id, to_construct_id{env_columns_aliases}
FROM no_props
{join} compound_property A_{prop} ON A_{prop}.compound_id = from_compound_id
{join} compound_property B_{prop} ON B_{prop}.compound_id = to_compound_id
{join} property_name {prop} ON A_{prop}.property_name_id = {prop}.id AND B_{prop}.property_name_id = {prop}.id AND {prop}.name = '{prop}'
INNER JOIN rule ON rule.from_smiles_id = no_props.from_smiles_id AND rule.to_smiles_id = no_props.to_smiles_id
"""
        not_yet_unioned_prop_selects.append(prop_select)

    statement = CTE_expression + insertion + "\n\nUNION\n".join(not_yet_unioned_prop_selects)

    # asyncpg requires $1, $2, etc. as placeholders for query variables, but we used %s in the above construction
    for i in range(len(mol_pkls)):
        statement = statement.replace('%s', '$' + str(i + 1), 1)

    return statement, mol_pkls, query

def aggregate_transforms_v2(filters, use_environment, property_metadata={}):

    # We will keep track of all the properties we see, so that we can do inner joins to these property tables later in the FROM block
    unique_prop_values = {}
    # We also need to keep track of unique combinations of properties and change types, because for each such unique combo we could calculate multiple stats
    unique_prop_changes = {}

    statistics_expressions = ""
    statistics_names = ""
    inverse_statistics_names = ""
    stats_nulls = ""
    m2m_aggs = ""

    for stat in filters.statistics:

        A_value, B_value = f"A_{stat.property_name}.value", f"B_{stat.property_name}.value"

        if stat.base == 'raw':
            if property_metadata[stat.property_name]['base'] == 'log':
                A_value, B_value = f"10 ^ ({A_value})", f"10 ^ ({B_value})"
            elif property_metadata[stat.property_name]['base'] == 'negative_log':
                A_value, B_value = f"10 ^ (-1 * {A_value})", f"10 ^ (-1 * {B_value})"

        if 'unit' in property_metadata[stat.property_name]:
            if stat.units == 'uM' and property_metadata[stat.property_name]['unit'] == 'M':
                A_value, B_value = f"(1000000 * ({A_value}))", f"(1000000 * ({B_value}))"
            elif stat.units == 'M' and property_metadata[stat.property_name]['unit'] == 'uM':
                A_value, B_value = f"(({A_value}) / 1000000)", f"(({B_value}) / 1000000)"

        # At present, formatting avg_value and avg_name is an "arbitrary decision" not controllable with API parameters
        if property_metadata[stat.property_name]['base'] == 'log':
            avg_value = f"({A_value} + {B_value}) / 2.0"
            avg_name = f"Average_log_{stat.property_name}"
        if property_metadata[stat.property_name]['base'] == 'negative_log':
            avg_value = f"(-1*{A_value} + -1*{B_value}) / 2.0"
            avg_name = f"Average_log_{stat.property_name}"
        if property_metadata[stat.property_name]['base'] == 'raw':
            avg_value = f"({A_value} + {B_value}) / 2.0"
            avg_name = f"Average_{stat.property_name}"

        unique_prop_changes[(stat.property_name, stat.change_type, stat.base, stat.units)] = {'A_value': A_value, 'B_value': B_value, 'name': f"{stat.property_name}_{stat.change_type}"}
        unique_prop_values[(stat.property_name, stat.base, stat.units)] = {'A_value': A_value, 'B_value': B_value, 'avg_name': avg_name, 'avg_value': avg_value}

        if stat.statistic == "median":
            name = f"{stat.statistic}_{stat.property_name}_{stat.change_type}"
            statistics_names += f", {name} AS {name}"
            if stat.change_type == 'delta':
                inverse_statistics_names += f", -1 * {name} AS {name}"
            elif stat.change_type == 'fold_change':
                inverse_statistics_names += f", 1 / {name} AS {name}"
            statistics_expressions += f""",
        PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY {stat.property_name}_{stat.change_type}) AS {name}"""

        if filters.aggregation_type == 'individual_transforms':
            stats_nulls += f", NULL AS {name}"
        elif filters.aggregation_type == 'group_by_fragment':
            stats_nulls += f", NULL AS percent_increased_{stat.property_name}"
            if stat.change_type == 'delta':
                m2m_aggs += f", AVG(({name} > 0)::int) * 100 AS percent_increased_{stat.property_name}"
            elif stat.change_type == 'fold_change':
                m2m_aggs += f", AVG(({name} > 1)::int) * 100 AS percent_increased_{stat.property_name}"

    # Apply filters to compounds A and/or B, for various properties, with various mathematical conditions
    range_conditions = ""
    inverse_range_conditions = ""
    unique_rf_props = set()
    for rf in filters.range_filters:

        if rf.compound == "A":
            other_compound = "B"
        elif rf.compound == "B":
            other_compound = "A"

        rf_value = f"{rf.compound}_{rf.property_name}.value"
        inverse_rf_value = f"{other_compound}_{rf.property_name}.value"

        if rf.base == 'raw':
            if property_metadata[rf.property_name]['base'] == 'log':
                rf_value, inverse_rf_value = f"10 ^ ({rf_value})", f"10 ^ ({inverse_rf_value})"
            elif property_metadata[rf.property_name]['base'] == 'negative_log':
                rf_value, inverse_rf_value = f"10 ^ (-1 * {rf_value})", f"10 ^ (-1 * {inverse_rf_value})"

        if rf.units == 'uM' and property_metadata[rf.property_name].get('unit') == 'M':
            rf_value, inverse_rf_value = f"1000000 * ({rf_value})", f"1000000 * ({inverse_rf_value})"
        elif rf.units == 'M' and property_metadata[rf.property_name].get('unit') == 'uM':
            rf_value, inverse_rf_value = f"({rf_value}) / 1000000", f"({inverse_rf_value}) / 1000000"
        
        unique_rf_props.add(rf.property_name)

        range_conditions += f"""
    AND {rf_value} {rf.operator} {rf.value}"""
        inverse_range_conditions += f"""
    AND {inverse_rf_value} {rf.operator} {rf.value}"""

    # Join all the necessary compound property tables needed, for both getting values (unique_prop_values), and filtering based on values (unique_rf_props)
    property_joins = ""
    property_joins_props = set(prop for prop,base,units in unique_prop_values.keys()).union(unique_rf_props)
    for prop in property_joins_props:
        property_joins += f"""
    INNER JOIN compound_property A_{prop} ON A_{prop}.compound_id = from_construct.compound_id
    INNER JOIN compound_property B_{prop} ON B_{prop}.compound_id = to_construct.compound_id
    INNER JOIN property_name {prop} ON A_{prop}.property_name_id = {prop}.id AND B_{prop}.property_name_id = {prop}.id AND {prop}.name = '{prop}'"""

    # Obtain values for the properties of interest, for each compound, A and B
    prop_values = ""
    inverse_prop_values = ""
    minmax_prop_values = ""
    minmax_prop_values_nulls = ""
    for (prop, base, units), vdict in unique_prop_values.items():
        ## Here we could format units / base for inserting into the property names in the SQL, so we can distinguish different combinations of units + base + prop values
        # But ultimately decided to simplify the code by not including base/units in the field names
        # This will probably lead to bugs if the user tries to simultaneously query different base/unit combinations of the same property, for instance log Soly_pH7 (M), and Soly_pH7 (uM)
        if base == 'raw':
            base = ""
        else:
            base = "_" + base
        if units is None:
            units = ""
        else:
            units = "_" + units

        prop_values += f", {vdict['A_value']} AS A_{prop}, {vdict['B_value']} AS B_{prop}, {vdict['avg_value']} AS {vdict['avg_name']}"
        inverse_prop_values += f", {vdict['B_value']} AS A_{prop}, {vdict['A_value']} AS B_{prop}, {vdict['avg_value']} AS {vdict['avg_name']}"
        minmax_prop_values += f""", MIN(A_{prop}) AS MIN_A_{prop}, MAX(A_{prop}) AS MAX_A_{prop}, MIN(B_{prop}) AS MIN_B_{prop}, 
    MAX(B_{prop}) AS MAX_B_{prop}, MIN({vdict['avg_name']}) AS MIN_{vdict['avg_name']}, MAX({vdict['avg_name']}) AS MAX_{vdict['avg_name']}"""
        minmax_prop_values_nulls += ", NULL, NULL, NULL, NULL, NULL, NULL"

    changes = ""
    inverse_changes = ""
    minmax_changes = ""
    minmax_changes_nulls = ""
    for (prop, change, base, units), vdict in unique_prop_changes.items():
        # Format units / base for inserting into the property names in the SQL, so we can distinguish different combinations of units + base + prop values
        if base == 'raw':
            base = ""
        else:
            base = "_" + base
        if units is None:
            units = ""
        else:
            units = "_" + units

        if change == "fold_change":
            changes += f", {vdict['B_value']} / {vdict['A_value']} AS {vdict['name']}"
            inverse_changes += f", {vdict['A_value']} / {vdict['B_value']} AS {vdict['name']}"
        elif change == "delta":
            changes += f", {vdict['B_value']} - {vdict['A_value']} AS {vdict['name']}"
            inverse_changes += f", {vdict['A_value']} - {vdict['B_value']} AS {vdict['name']}"
        minmax_changes += f", MIN({vdict['name']}) AS MIN_{vdict['name']}, MAX({vdict['name']}) AS MAX_{vdict['name']}"

        minmax_changes_nulls += ", NULL, NULL"

    # Here we control how the pairs stored in query_result are flipped
    # For one_to_all results, we want to flip the pairs where use_original_direction = FALSE, because they are stored in the opposite direction relative to the user's query
    # However, for many_to_many queries where the from_smiles = to_smiles, then we only ever searched in the original direction; we did not query in the opposite direction
    # This is because the opposite direction, when from_smiles = to_smiles, is by definition the inverse of the forward direction, and would be a waste of space to store in the query_result table
    # Therefore all many_to_many results, when from_smiles = to_smiles, will have use_original_direction = TRUE, and we need to flip all of these, which we accomplish by setting flip_boolean = 'TRUE'
    ## TODO if we want to do many_to_many queries with from_smiles != to_smiles, we need to reexamine this code and probably update
    if filters.aggregation_type == 'individual_transforms':
        flip_boolean = 'FALSE'
    elif filters.aggregation_type == 'group_by_fragment':
        flip_boolean = 'TRUE'

    if use_environment:
        ungrouped_distinct_CTE = """
ungrouped_distinct AS (
    SELECT DISTINCT ON (ungrouped.from_smiles_env, ungrouped.to_smiles_env, ungrouped.from_compound_id, ungrouped.to_compound_id) *
    FROM ungrouped
),
"""
        final_ungrouped_selection = "ungrouped_distinct"
        grouped_CTE_select = f"from_smiles_env, to_smiles_env"
        one_to_all_NULLs = "NULL AS from_smiles_env, NULL AS to_smiles_env"
        many_to_many_select = "to_smiles_env, count(to_smiles_env) AS transform_count, sum(pair_count) AS pair_count, ARRAY_AGG(from_smiles_env) AS from_smiles_env_array"
        many_to_many_NULLs = "NULL AS to_smiles_env, NULL AS transform_count, NULL AS pair_count, NULL AS from_smiles_env_array"
        many_to_many_groupby = "to_smiles_env"
    else:
        ungrouped_distinct_CTE = ""
        final_ungrouped_selection = "ungrouped"
        grouped_CTE_select = f"rule_id, from_smiles, to_smiles"
        one_to_all_NULLs = "NULL AS rule_id, NULL AS from_smiles, NULL AS to_smiles"
        many_to_many_select = "ARRAY_AGG(rule_id) AS rule_id_array, to_smiles, count(rule_id) AS transform_count, sum(pair_count) AS pair_count, ARRAY_AGG(from_smiles) AS from_smiles_array"
        many_to_many_NULLs = "NULL AS rule_id_array, NULL AS to_smiles, NULL AS transform_count, NULL AS pair_count, NULL AS from_smiles_array"
        many_to_many_groupby = "to_smiles"

    CTE = f"""
WITH ungrouped AS (
    SELECT query_result.rule_id AS rule_id, from_smiles.smiles AS from_smiles, to_smiles.smiles AS to_smiles, query_result.from_smiles_env AS from_smiles_env,
        query_result.to_smiles_env AS to_smiles_env, from_construct.compound_id AS from_compound_id, to_construct.compound_id AS to_compound_id {prop_values} {changes}

    FROM query_result
    INNER JOIN from_construct ON query_result.from_construct_id = from_construct.id
    INNER JOIN to_construct ON query_result.to_construct_id = to_construct.id
    INNER JOIN rule_smiles from_smiles ON from_construct.rule_smiles_id = from_smiles.id
    INNER JOIN rule_smiles to_smiles ON to_construct.rule_smiles_id = to_smiles.id {property_joins}
    WHERE query_result.query_id = {filters.query_id}
    AND query_result.use_original_direction = TRUE {range_conditions}
    
    UNION ALL
    
    SELECT -1 * query_result.rule_id AS rule_id, to_smiles.smiles AS from_smiles, from_smiles.smiles AS to_smiles, query_result.to_smiles_env AS from_smiles_env, 
        query_result.from_smiles_env AS to_smiles_env, to_construct.compound_id AS from_compound_id, from_construct.compound_id AS to_compound_id {inverse_prop_values} {inverse_changes}

    FROM query_result
    INNER JOIN from_construct ON query_result.from_construct_id = from_construct.id
    INNER JOIN to_construct ON query_result.to_construct_id = to_construct.id
    INNER JOIN rule_smiles from_smiles ON from_construct.rule_smiles_id = from_smiles.id
    INNER JOIN rule_smiles to_smiles ON to_construct.rule_smiles_id = to_smiles.id {property_joins}
    WHERE query_result.query_id = {filters.query_id}
    AND query_result.use_original_direction = {flip_boolean} {inverse_range_conditions}
), {ungrouped_distinct_CTE}
grouped AS (
    SELECT {grouped_CTE_select}, count(rule_id) AS pair_count {statistics_expressions}
    FROM {final_ungrouped_selection}
    GROUP BY {grouped_CTE_select}
    ORDER BY pair_count DESC
)"""

    if filters.aggregation_type == 'individual_transforms':

        statement = CTE + f"""
SELECT {one_to_all_NULLs}, NULL AS pair_count {stats_nulls} {minmax_prop_values} {minmax_changes}
FROM {final_ungrouped_selection}

UNION ALL

SELECT * {minmax_prop_values_nulls} {minmax_changes_nulls} FROM grouped
"""
        return statement

    elif filters.aggregation_type == 'group_by_fragment':

        statement = CTE + f"""
SELECT {many_to_many_NULLs} {stats_nulls} {minmax_prop_values} {minmax_changes}
FROM {final_ungrouped_selection}

UNION ALL

SELECT {many_to_many_select} {m2m_aggs} {minmax_prop_values_nulls} {minmax_changes_nulls}
FROM grouped
GROUP BY {many_to_many_groupby}
"""
        return statement

    # Currently not used, but could be used for a more general many_to_many search where from_smiles != to_smiles
    elif filters.aggregation_type == 'old_group_by_fragment':

        # Here we combine the original transforms with flipped versions thereof.
        # We want original + flipped in the same list, so we can aggregate by unique SMILES on one side of the equation, and catch everything.
        CTE += f""",
one_to_all AS (
    SELECT rule_id AS rule_id, from_smiles AS from_smiles, to_smiles AS to_smiles, pair_count AS pair_count {statistics_names}
    FROM grouped

    UNION ALL

    SELECT -1 * rule_id AS rule_id, to_smiles AS from_smiles, from_smiles AS to_smiles, pair_count AS pair_count {inverse_statistics_names}
    FROM grouped
)
"""        
        statement = CTE + f"""
SELECT NULL AS rule_id_array, NULL AS to_smiles, NULL AS transform_count, NULL AS pair_count, NULL AS from_smiles_array {stats_nulls} {minmax_prop_values} {minmax_changes}
FROM ungrouped

UNION ALL

SELECT ARRAY_AGG(rule_id) AS rule_id_array, to_smiles, count(rule_id) AS transform_count, sum(pair_count) AS pair_count, ARRAY_AGG(from_smiles) AS from_smiles_array
    {m2m_aggs} {minmax_prop_values_nulls} {minmax_changes_nulls}
FROM one_to_all
GROUP BY to_smiles
"""
        return statement

def get_plot_data_statement(filters, property_metadata={}):

    # We will keep track of all the properties we see, so that we can do inner joins to these property tables later in the FROM block
    props_for_joins = set()
    axis_expressions = set()
    axis_names, inverse_axis_names = "", ""

    for axis in (filters.x_data, filters.y_data):
        props_for_joins.add(axis.property_name)

        A_value, B_value = f"A_{axis.property_name}.value", f"B_{axis.property_name}.value"

        if axis.base == 'raw':
            if property_metadata[axis.property_name]['base'] == 'log':
                A_value, B_value = f"(10 ^ {A_value})", f"(10 ^ {B_value})"
            elif property_metadata[axis.property_name]['base'] == 'negative_log':
                A_value, B_value = f"(10 ^ (-1 * {A_value}))", f"(10 ^ (-1 * {B_value}))"

        if 'unit' in property_metadata[axis.property_name]:
            if axis.units == 'uM' and property_metadata[axis.property_name]['unit'] == 'M':
                A_value, B_value = f"(1000000 * {A_value})", f"(1000000 * {B_value})"
            elif axis.units == 'M' and property_metadata[axis.property_name]['unit'] == 'uM':
                A_value, B_value = f"({A_value} / 1000000)", f"({B_value} / 1000000)"

        # We always need the A and B fields for each unique property, but use set.add to avoid duplication if we select two different change_types for the same property
        axis_expressions.add(f", {A_value} AS A_{axis.property_name}")
        axis_expressions.add(f", {B_value} AS B_{axis.property_name}")

        if axis.change_type == 'A':
            axis_names += f", A_{axis.property_name}"
            inverse_axis_names += f", B_{axis.property_name} AS A_{axis.property_name}"
        elif axis.change_type == 'B':
            axis_names += f", B_{axis.property_name}"
            inverse_axis_names += f", A_{axis.property_name} AS B_{axis.property_name}"
        elif axis.change_type == "average":
            axis_expressions.add(f", ({A_value} + {B_value}) / 2.0 AS Average_{axis.property_name}")
            axis_names += f", Average_{axis.property_name}"
            inverse_axis_names += f", Average_{axis.property_name}"
        elif axis.change_type == "fold_change": 
            axis_expressions.add(f", {B_value} / {A_value} AS {axis.property_name}_fold_change")
            axis_names += f", {axis.property_name}_fold_change"
            inverse_axis_names += f", 1 / {axis.property_name}_fold_change AS {axis.property_name}_fold_change"
        elif axis.change_type == "delta":
            axis_expressions.add(f", {B_value} - {A_value} AS {axis.property_name}_delta")
            axis_names += f", {axis.property_name}_delta"
            inverse_axis_names += f", -1 * {axis.property_name}_delta AS {axis.property_name}_delta"
            
    axis_expressions = ''.join(axis_expressions)

    range_conditions, inverse_range_conditions = "", ""
    unique_rf_props = set()
    for rf in filters.range_filters:

        if rf.compound == "A":
            other_compound = "B"
        elif rf.compound == "B":
            other_compound = "A"

        rf_value = f"{rf.compound}_{rf.property_name}"
        inverse_rf_value = f"{other_compound}_{rf.property_name}"

        """
        if rf.base == 'raw':
            if default_bases_units[rf.property_name]['base'] == 'log':
                rf_value, inverse_rf_value = f"10 ^ ({rf_value})", f"10 ^ ({inverse_rf_value})"
            elif default_bases_units[rf.property_name]['base'] == 'negative_log':
                rf_value, inverse_rf_value = f"10 ^ (-1 * {rf_value})", f"10 ^ (-1 * {inverse_rf_value})"

        if rf.units == 'uM' and default_bases_units[rf.property_name].get('unit') == 'M':
            rf_value, inverse_rf_value = f"1000000 * {rf_value}", f"1000000 * {inverse_rf_value}"
        elif rf.units == 'M' and default_bases_units[rf.property_name].get('unit') == 'uM':
            rf_value, inverse_rf_value = f"{rf_value} / 1000000", f"{inverse_rf_value} / 1000000"
        """

        unique_rf_props.add(rf.property_name)

        range_conditions += f"\nAND {rf_value} {rf.operator} {rf.value}"
        inverse_range_conditions += f"\nAND {inverse_rf_value} {rf.operator} {rf.value}"

    # Obtain values for the properties of interest, for each compound, A and B
    property_joins = ""
    for prop in props_for_joins.union(unique_rf_props):
        property_joins += f"""
    INNER JOIN compound_property A_{prop} ON A_{prop}.compound_id = from_construct.compound_id
    INNER JOIN compound_property B_{prop} ON B_{prop}.compound_id = to_construct.compound_id
    INNER JOIN property_name {prop} ON A_{prop}.property_name_id = {prop}.id AND B_{prop}.property_name_id = {prop}.id AND {prop}.name = '{prop}'"""

    m2m_selects, m2m_joins, m2m_names, inverse_m2m_names = '', '', '', ''
    if filters.grouped_by_environment == False:
        if filters.aggregation_type == 'individual_transforms':
            m2m_selects = ', query_result.rule_id AS rule_id'
            m2m_names = ', rule_id'
            inverse_m2m_names = ', -1 * rule_id AS rule_id'
        elif filters.aggregation_type == 'group_by_fragment':
            m2m_selects = ', query_result.rule_id AS rule_id, from_rule_smiles.smiles AS from_smiles, to_rule_smiles.smiles AS to_smiles'
            m2m_joins = """\nINNER JOIN rule_smiles from_rule_smiles ON from_rule_smiles.id = from_construct.rule_smiles_id
    INNER JOIN rule_smiles to_rule_smiles ON to_rule_smiles.id = to_construct.rule_smiles_id"""
            m2m_names = ', rule_id, to_smiles'
            inverse_m2m_names = ', -1 * rule_id AS rule_id, from_smiles AS to_smiles'
    elif filters.grouped_by_environment == True:
        m2m_selects = ', query_result.from_smiles_env AS from_smiles_env, query_result.to_smiles_env AS to_smiles_env,\n'
        m2m_selects += 'from_construct.compound_id AS from_compound_id, to_construct.compound_id AS to_compound_id'  
        if filters.aggregation_type == 'individual_transforms':
            # We need both of the below fields to encode "rule" groupings of one-to-all multicut queries
            m2m_names = ', from_smiles_env, to_smiles_env, from_compound_id, to_compound_id'
            inverse_m2m_names = ', to_smiles_env AS from_smiles_env, from_smiles_env AS to_smiles_env, to_compound_id AS from_compound_id, from_compound_id AS to_compound_id'
        elif filters.aggregation_type == 'group_by_fragment':
            # We need only to_smiles_env for many_to_many since we group/color by 'TO' fragment
            m2m_names = ', to_smiles_env, from_compound_id, to_compound_id'
            inverse_m2m_names = ', from_smiles_env AS to_smiles_env, to_compound_id AS from_compound_id, from_compound_id AS to_compound_id'

    final_select, forward_conditions, reverse_conditions = [], "", ""
    if filters.ids:
        if filters.grouped_by_environment == False:
            # Then ids are integers (rule_ids), and it's a 1-cut query
            forward_conditions, reverse_conditions = rule_ids_to_conditions(filters.ids, filters.aggregation_type)
        elif filters.grouped_by_environment == True:
            # Then ids are strings, and it's a 2-cut or 3-cut query
            if filters.aggregation_type == 'individual_transforms':
                # Then ids are strings: fromSmilesEnv_toSmilesEnv
                forward_conditions, reverse_conditions = smiles_envs_to_conditions(filters.ids)

            elif filters.aggregation_type == 'group_by_fragment':
                # Then ids are strings: toSmilesEnv
                forward_conditions = f"to_smiles_env IN ({str(filters.ids)[1:-1]})"
                reverse_conditions = f"from_smiles_env IN ({str(filters.ids)[1:-1]})"

    if forward_conditions:
        final_select.append(f"""
SELECT from_construct_id, to_construct_id {m2m_names} {axis_names}
FROM mixed_directions
WHERE {forward_conditions} {range_conditions}\n""")

    if reverse_conditions:
        # Make construct_id values negative to indicate that they are flipped. We need to encode this info for when we fetch individual pair data from points in the scatterplot
        final_select.append(f"""
SELECT -1 * to_construct_id AS from_construct_id, -1 * from_construct_id AS to_construct_id {inverse_m2m_names} {inverse_axis_names} 
FROM mixed_directions
WHERE {reverse_conditions} {inverse_range_conditions}\n""")

    if len(final_select) > 1:
        final_select = "\nUNION ALL\n".join(final_select)
    else:
        final_select = "".join(final_select)

    if filters.grouped_by_environment == True:
        final_select = f""",
nondistinct AS ({final_select})
SELECT DISTINCT ON (from_compound_id, to_compound_id) from_construct_id, to_construct_id {m2m_names.replace(', from_compound_id', '').replace(', to_compound_id', '')} {axis_names}
FROM nondistinct
"""

    statement = f"""
WITH mixed_directions AS (
    SELECT query_result.use_original_direction AS use_original_direction,
    from_construct.id AS from_construct_id, to_construct.id AS to_construct_id {m2m_selects}
    {axis_expressions}
    FROM query_result
    INNER JOIN from_construct ON query_result.from_construct_id = from_construct.id
    INNER JOIN to_construct ON query_result.to_construct_id = to_construct.id {property_joins} {m2m_joins}
    WHERE query_result.query_id = {filters.query_id} 
){final_select}
"""
    return statement

def get_pair_data_statement(pairs, property_metadata={}):

    # We will keep track of all the properties we see, so that we can do inner joins to these property tables later in the FROM block
    props_for_joins = set()
    axis_expressions = {}
    inverse_axis_expressions = {}

    for axis in pairs.prop_changes:
        props_for_joins.add(axis.property_name)

        A_value, B_value = f"A_{axis.property_name}.value", f"B_{axis.property_name}.value"

        if axis.base == 'raw':
            if property_metadata[axis.property_name]['base'] == 'log':
                A_value, B_value = f"(10 ^ {A_value})", f"(10 ^ {B_value})"
            elif property_metadata[axis.property_name]['base'] == 'negative_log':
                A_value, B_value = f"(10 ^ (-1 * {A_value}))", f"(10 ^ (-1 * {B_value}))"

        if 'unit' in property_metadata[axis.property_name]:
            if axis.units == 'uM' and property_metadata[axis.property_name]['unit'] == 'M':
                A_value, B_value = f"(1000000 * {A_value})", f"(1000000 * {B_value})"
            elif axis.units == 'M' and property_metadata[axis.property_name]['unit'] == 'uM':
                A_value, B_value = f"({A_value} / 1000000)", f"({B_value} / 1000000)"

        # We always need the A and B fields for each unique property, but use dict to avoid duplication if we select two different change_types for the same property
        axis_expressions[f"A_{axis.property_name}"] = f", {A_value} AS A_{axis.property_name}"
        axis_expressions[f"B_{axis.property_name}"] = f", {B_value} AS B_{axis.property_name}"

        inverse_axis_expressions[f"A_{axis.property_name}"] = f", {B_value} AS A_{axis.property_name}"
        inverse_axis_expressions[f"B_{axis.property_name}"] = f", {A_value} AS B_{axis.property_name}"

        if axis.change_type == "fold_change": 
            axis_expressions[f"{axis.property_name}_fold_change"] = f", {B_value} / {A_value} AS {axis.property_name}_fold_change"
            inverse_axis_expressions[f"{axis.property_name}_fold_change"] = f", {A_value} / {B_value} AS {axis.property_name}_fold_change"
        elif axis.change_type == "delta":
            axis_expressions[f"{axis.property_name}_delta"] = f", {B_value} - {A_value} AS {axis.property_name}_delta"
            inverse_axis_expressions[f"{axis.property_name}_delta"] = f", {A_value} - {B_value} AS {axis.property_name}_delta"
            
    # Need the ordering of the axis expressions to match, otherwise the wrong columns can be glued together in the UNION
    axis_expressions = ''.join([b for a,b in sorted(axis_expressions.items())])
    inverse_axis_expressions = ''.join([b for a,b in sorted(inverse_axis_expressions.items())])

    # Obtain values for the properties of interest, for each compound, A and B
    property_joins = ""
    for prop in props_for_joins:
        property_joins += f"""
    INNER JOIN compound_property A_{prop} ON A_{prop}.compound_id = from_construct.compound_id
    INNER JOIN compound_property B_{prop} ON B_{prop}.compound_id = to_construct.compound_id
    INNER JOIN property_name {prop} ON A_{prop}.property_name_id = {prop}.id AND B_{prop}.property_name_id = {prop}.id AND {prop}.name = '{prop}'"""

    forward_pairs_A, forward_pairs_B, reversed_pairs_A, reversed_pairs_B = [], [], [], []
    for a,b in pairs.construct_id_pairs:
        # a and b are negative if the order is flipped compared to the originally stored direction in the DB
        # For example, if the MMP is stored as b -> a in the DB, but we wanted to display the MMP in the a -> b direction, we 'remember' we flipped it by making them negative
        # if a is negative, then b should always be negative (a and b always have the same sign)
        if a < 0:
            reversed_pairs_A.append(-1 * b)
            reversed_pairs_B.append(-1 * a)
        else:
            forward_pairs_A.append(a)
            forward_pairs_B.append(b)

    pre_union = []
    if len(forward_pairs_A) > 0:
        pre_union.append(f"""
SELECT A_compound.public_id AS A_ID, A_compound.clean_smiles AS A, B_compound.public_id AS B_ID, B_compound.clean_smiles AS B, constant_smiles.smiles AS constant
{axis_expressions}
FROM (SELECT unnest(ARRAY[{str(forward_pairs_A)[1:-1]}]) AS from, unnest(ARRAY[{str(forward_pairs_B)[1:-1]}]) AS to) input_constructs
INNER JOIN from_construct ON input_constructs.from = from_construct.id
INNER JOIN to_construct ON input_constructs.to = to_construct.id 
INNER JOIN compound A_compound ON A_compound.id = from_construct.compound_id
INNER JOIN compound B_compound ON B_compound.id = to_construct.compound_id 
INNER JOIN constant_smiles ON from_construct.constant_id = constant_smiles.id {property_joins}\n""")

    if len(reversed_pairs_A) > 0:
        pre_union.append(f"""
SELECT B_compound.public_id AS A_ID, B_compound.clean_smiles AS A, A_compound.public_id AS B_ID, A_compound.clean_smiles AS B, constant_smiles.smiles AS constant
{inverse_axis_expressions}
FROM (SELECT unnest(ARRAY[{str(reversed_pairs_A)[1:-1]}]) AS from, unnest(ARRAY[{str(reversed_pairs_B)[1:-1]}]) AS to) input_constructs
INNER JOIN from_construct ON input_constructs.from = from_construct.id
INNER JOIN to_construct ON input_constructs.to = to_construct.id 
INNER JOIN compound A_compound ON A_compound.id = from_construct.compound_id
INNER JOIN compound B_compound ON B_compound.id = to_construct.compound_id 
INNER JOIN constant_smiles ON from_construct.constant_id = constant_smiles.id {property_joins}\n""")

    statement = "\nUNION ALL\n".join(pre_union)

    return statement


def get_all_raw_data_statement(inputs, req, opt):

    prop_values, inverse_prop_values = "", ""
    joins = f"""
    INNER JOIN from_construct ON query_result.from_construct_id = from_construct.id
    INNER JOIN to_construct ON query_result.to_construct_id = to_construct.id
    INNER JOIN compound compound_A ON compound_A.id = from_construct.compound_id
    INNER JOIN compound compound_B ON compound_B.id = to_construct.compound_id
    INNER JOIN constant_smiles ON from_construct.constant_id = constant_smiles.id
    INNER JOIN rule_smiles from_smiles ON from_construct.rule_smiles_id = from_smiles.id
    INNER JOIN rule_smiles to_smiles ON to_construct.rule_smiles_id = to_smiles.id"""

    req_iterator = list(zip(req, ["INNER JOIN"]*len(req)))
    opt_iterator = list(zip(opt, ["LEFT OUTER JOIN"]*len(opt)))
    for prop, join in req_iterator + opt_iterator:
        prop_values += f", A_{prop}.value AS A_{prop}, B_{prop}.value AS B_{prop}"
        inverse_prop_values += f", B_{prop}.value AS A_{prop}, A_{prop}.value AS B_{prop}"
        joins += f"""
    INNER JOIN property_name {prop} ON {prop}.name = '{prop}'
    {join} compound_property A_{prop} ON A_{prop}.compound_id = from_construct.compound_id AND A_{prop}.property_name_id = {prop}.id
    {join} compound_property B_{prop} ON B_{prop}.compound_id = to_construct.compound_id AND B_{prop}.property_name_id = {prop}.id
"""
    joins += f"WHERE query_result.query_id = {inputs.query_id}\n"

    statement = f"""
    SELECT compound_A.public_id AS A_ID, compound_A.clean_smiles AS A_smiles, compound_B.public_id AS B_ID, compound_B.clean_smiles AS B_smiles,
    query_result.rule_id AS rule_id, from_smiles.smiles AS from_smiles, to_smiles.smiles AS to_smiles, constant_smiles.smiles AS constant_smiles 
    {prop_values}

    FROM query_result {joins}
    AND query_result.use_original_direction = TRUE
"""

    if inputs.query.advanced_options.aggregation_type == 'individual_transforms':

        statement += f"""
    UNION ALL
    
    SELECT compound_B.public_id AS A_ID, compound_B.clean_smiles AS A_smiles, compound_A.public_id AS B_ID, compound_A.clean_smiles AS B_smiles,
    -1 * query_result.rule_id AS rule_id, to_smiles.smiles AS from_smiles, from_smiles.smiles AS to_smiles, constant_smiles.smiles AS constant_smiles 
    {inverse_prop_values}

    FROM query_result {joins}
    AND query_result.use_original_direction = FALSE
"""
    return statement

def get_constant_smiles_perms(constant_perms):

    # We need to use constant (environment) smiles that can be converted with the smiles_syntax dependency into fragments with closures,
    # because later we will need to glue the constant smiles to the variable smiles
    constant_smiles_list = []
    He_Ne_Ar_set = set((2, 10, 18))
    for perm in constant_perms:
        # Don't modify in place
        perm = copy.deepcopy(perm)
        for atom in perm.GetAtoms():
            if atom.GetAtomicNum() not in He_Ne_Ar_set:
                # Later we will glue the perm onto the variable, but we want to 'remember' which atoms came from the constant. We do this with an atom map number = 7
                atom.SetAtomMapNum(7)
        constant_smiles_list.append(ss_select.get_multifragment_gluable_smiles(perm))

    return constant_smiles_list

def select_variables_and_env(query_id):

# We start off with the smiles representation of our environment in the from_smiles_env column
# After we retrieve this environment_smiles here, we no longer need it, and will overwrite the value there with the glued-together from_smiles_env
    statement = f"""
SELECT query_result.id, from_smiles.smiles AS from_smiles, to_smiles.smiles AS to_smiles, query_result.from_smiles_env AS environment_smiles
FROM query_result
INNER JOIN from_construct ON from_construct.id = query_result.from_construct_id
INNER JOIN to_construct ON to_construct.id = query_result.to_construct_id
INNER JOIN rule_smiles from_smiles ON from_smiles.id = from_construct.rule_smiles_id
INNER JOIN rule_smiles to_smiles ON to_smiles.id = to_construct.rule_smiles_id
WHERE query_result.query_id = {query_id}
"""
    return statement

def glue_variables_and_env(data):
    glued_smiles = {}
    for row in data:
        glued_smiles[row['id']] = {}
        env_smiles = row['environment_smiles']
        for direction in ('from', 'to'):
            variable_smiles = row[direction + '_smiles']
            for patt in ('[*:1]', '[*:2]', '[*:3]'):
                if patt in variable_smiles and patt not in env_smiles:
                    variable_smiles = variable_smiles.replace(patt, '[At:7]')
            variable_to_weld = smiles_syntax.convert_labeled_wildcards_to_closures(variable_smiles)
            env_to_weld = smiles_syntax.convert_labeled_wildcards_to_closures(env_smiles)
            glued_smiles[row['id']][direction] = Chem.MolToSmiles(Chem.MolFromSmiles(f"{variable_to_weld}.{env_to_weld}", sanitize=False))

    data_tuples = [(row_id, glued_smiles[row_id]['from'], glued_smiles[row_id]['to']) for row_id in glued_smiles]
    data_tuples = str(data_tuples)[1:-1]
    statement = f"""
UPDATE query_result
SET from_smiles_env = data.from_glued, to_smiles_env = data.to_glued
FROM (VALUES {data_tuples}) AS data(id, from_glued, to_glued)
WHERE query_result.id = data.id
"""
    return statement

def rule_ids_to_conditions(ids, aggregation_type):
    forward_conditions, reverse_conditions = "", ""
    assert ids
    positive_rule_ids = [a for a in ids if a >= 0]
    if len(positive_rule_ids) > 0:
        forward_conditions = f"rule_id IN ({str(positive_rule_ids)[1:-1]})"
        if aggregation_type == 'individual_transforms':
            forward_conditions = f"use_original_direction = TRUE\nAND {forward_conditions}"
    # We make the rule_ids positive to match the originally stored rule_ids (based on foreign key), but "remember" that these rule_ids were flipped by having them in a separate list
    negative_rule_ids = [-1 * a for a in ids if a < 0]
    if len(negative_rule_ids) > 0:
        reverse_conditions = f"rule_id IN ({str(negative_rule_ids)[1:-1]})"
        if aggregation_type == 'individual_transforms':
            reverse_conditions = f"use_original_direction = FALSE\nAND {reverse_conditions}"

    return forward_conditions, reverse_conditions

def smiles_envs_to_conditions(ids):
    forward_prejoin, reverse_prejoin = [], []
    for from_to in ids:
        from_smiles_env, to_smiles_env = from_to.split('_')
        forward_prejoin.append(f"""(from_smiles_env = '{from_smiles_env}' AND to_smiles_env = '{to_smiles_env}')""")
        reverse_prejoin.append(f"""(from_smiles_env = '{to_smiles_env}' AND to_smiles_env = '{from_smiles_env}')""")
    forward_conditions = "(" + f"\nOR ".join(forward_prejoin) + ")"
    reverse_conditions = "(" + f"\nOR ".join(reverse_prejoin) + ")"
    return forward_conditions, reverse_conditions

def ids_to_transforms(ids, query_id, aggregation_type, use_environment):

    ids = tuple(ids.keys())
    select_forward = "SELECT DISTINCT to_smiles.smiles AS to_smiles, query_result.environment_smarts AS environment_smarts"
    select_reverse = "SELECT DISTINCT from_smiles.smiles AS to_smiles, query_result.environment_smarts AS environment_smarts"

    if aggregation_type == "individual_transforms":
        if use_environment == False:
            select_forward += ", from_smiles.smiles AS from_smiles, query_result.rule_id AS transform_id"
            select_reverse += ", to_smiles.smiles AS from_smiles, -1 * query_result.rule_id AS transform_id"
            forward_conditions, reverse_conditions = rule_ids_to_conditions(ids, aggregation_type)
        else:
            select_forward += ", from_smiles.smiles AS from_smiles, (query_result.from_smiles_env || '_' || query_result.to_smiles_env) AS transform_id"
            select_reverse += ", to_smiles.smiles AS from_smiles, (query_result.to_smiles_env || '_' || query_result.from_smiles_env) AS transform_id"
            forward_conditions, reverse_conditions = smiles_envs_to_conditions(ids)

        from_block_forward = f"""
FROM query_result
INNER JOIN from_construct ON from_construct.id = query_result.from_construct_id
INNER JOIN to_construct ON to_construct.id = query_result.to_construct_id
INNER JOIN rule_smiles from_smiles ON from_smiles.id = from_construct.rule_smiles_id
INNER JOIN rule_smiles to_smiles ON to_smiles.id = to_construct.rule_smiles_id
WHERE query_result.query_id = {query_id}
"""
        from_block_reverse = from_block_forward

    elif aggregation_type == "group_by_fragment":
        assert use_environment == True
        select_forward += ", query_result.to_smiles_env AS transform_id"
        select_reverse += ", query_result.from_smiles_env AS transform_id"
        forward_conditions = f"to_smiles_env IN ({str(ids)[1:-1]})"
        reverse_conditions = f"from_smiles_env IN ({str(ids)[1:-1]})"

        from_block_forward = f"""
FROM query_result
INNER JOIN to_construct ON to_construct.id = query_result.to_construct_id
INNER JOIN rule_smiles to_smiles ON to_smiles.id = to_construct.rule_smiles_id
WHERE query_result.query_id = {query_id}
"""
        from_block_reverse = f"""
FROM query_result
INNER JOIN from_construct ON from_construct.id = query_result.from_construct_id
INNER JOIN rule_smiles from_smiles ON from_smiles.id = from_construct.rule_smiles_id
WHERE query_result.query_id = {query_id}
"""
    preunion = []
    if forward_conditions:
        preunion.append(f"{select_forward} {from_block_forward} AND {forward_conditions}")
    if reverse_conditions:
        preunion.append(f"{select_reverse} {from_block_reverse} AND {reverse_conditions}")
    statement = "\n UNION \n".join(preunion) if preunion else ""

    return statement

def flip_pairs(pairs):
    flipped_pairs = []

    for pair in pairs:
            prop_values = pair[8:]
            switched_prop_values = []
            for i in range(0, len(prop_values), 2):
                switched_prop_values += [prop_values[i+1], prop_values[i]]
            flipped_pairs.append(pair[2:4] + pair[0:2] + [-1 * pair[4]] + [pair[6]] + [pair[5]] + [pair[7]] + switched_prop_values)

    return flipped_pairs

def frag_to_q_mol(frag_smiles : str):
    return Chem.MolToMolBlock(Chem.AddHs(Chem.MolFromSmiles(frag_smiles.replace('[*:1]', '[Rb]')))).replace('0  0  0  0  0  1  0  0  0  0  0  0', '0  0  0  0  0  0  0  0  0  0  0  0')

def flip_inferred_pairs(pairs, num_props):
    flipped_pairs = []
    for pdata in pairs:
            flipped_pairs.append(tuple([d for d in pdata[(2 + num_props)*3 : -5]] + [c for c in pdata[(2 + num_props)*2 : (2 + num_props)*3]] + \
                [b for b in pdata[(2 + num_props) : (2 + num_props)*2]] + [a for a in pdata[0 : (2 + num_props)]] + [-1*pdata[-2], pdata[-1], pdata[-3]] + \
                [-1*pdata[-5], pdata[-4]]))
    return flipped_pairs

# Trying out a SQL version of inferring pairs
# Currently will return inferred pairs in the form of xA -> xB, yB -> yC, where lowercase = constant, uppercase = variable, x != y, A != C
# Currently only applying variable constraints to A and C. If we apply the variable constraint also to B, we will get fewer results

# The final desired function which accomplishes all of the following in one query:
# from to from to
# from to to from
# to from from to
# to from to from

# flipping should already be taken care of when variable1 != variable2 and we want everything to be in the direction of variable1 -> variable2
# Specify get_flipped = True for variable1 = variable2 type comparisons, when you want to be able to see everything going to one fragment

def infer_pairs(query, cursor, print_statement=False, num_frags=1, get_flipped = False):
    
    pairs = {}
    req, opt = query["REQUIRED_properties"], query["OPTIONAL_properties"]
    num_props = len(req + opt)
    
    # A and D are the compounds to which we apply variable constraints
    # B and C are the "connection points" for inferring changes from A to D
    # We apply the constant constraint to both pairs A->B and C->D
    # These two pairs must have different constants (so that A->D is not a real matched pair)

    selectA = "X.A_public_id, X.A_clean_smiles"
    XsubA = "A.public_id AS A_public_id, A.clean_smiles AS A_clean_smiles"
    selectB = "X.B_public_id, X.B_clean_smiles"
    XsubB = "B.public_id AS B_public_id, B.clean_smiles AS B_clean_smiles"

    selectC = "Y.C_public_id, Y.C_clean_smiles"
    YsubC = "C.public_id AS C_public_id, C.clean_smiles AS C_clean_smiles"
    selectD = "Y.D_public_id, Y.D_clean_smiles"
    YsubD = "D.public_id AS D_public_id, D.clean_smiles AS D_clean_smiles"

    A_cp_labels, B_cp_labels, C_cp_labels, D_cp_labels, pname_labels = [], [], [], [], []
    X_prop_joins_A, X_prop_joins_B, Y_prop_joins_C, Y_prop_joins_D = '', '', '', ''
    label = "Z"
    
    # INNER JOINS for required properties - all compounds in all pairs must have these properties
    # LEFT OUTER JOINS for optional properties - will retrieve these properties if available
    req_iterator = list(zip(req, ["INNER JOIN"]*len(req)))
    opt_iterator = list(zip(opt, ["LEFT OUTER JOIN"]*len(opt)))

    for prop, join in req_iterator + opt_iterator:

        A_cp_label = get_next_label(label)
        A_cp_labels.append(A_cp_label)
        B_cp_label = get_next_label(A_cp_label)
        B_cp_labels.append(B_cp_label)
        C_cp_label = get_next_label(B_cp_label)
        C_cp_labels.append(C_cp_label)
        D_cp_label = get_next_label(C_cp_label)
        D_cp_labels.append(D_cp_label)
        pname_label = get_next_label(D_cp_label)
        pname_labels.append(pname_label)
        label = pname_label

        selectA += ", X." + A_cp_label + "_value"
        XsubA += ", " + A_cp_label + ".value AS " + A_cp_label + "_value"
        selectB += ", X." + B_cp_label + "_value"
        XsubB += ", " + B_cp_label + ".value AS " + B_cp_label + "_value"
        selectC += ", Y." + C_cp_label + "_value"
        YsubC += ", " + C_cp_label + ".value AS " + C_cp_label + "_value"
        selectD += ", Y." + D_cp_label + "_value"
        YsubD += ", " + D_cp_label + ".value AS " + D_cp_label + "_value"

        X_prop_joins_A += join + " (compound_property " + A_cp_label + " INNER JOIN property_name " + pname_label + " ON " + \
            A_cp_label + ".property_name_id = " + pname_label + ".id AND " + pname_label + ".name = '" + prop + "') ON A.id = " + \
            A_cp_label + ".compound_id\n"
        X_prop_joins_B += join + " (compound_property " + B_cp_label + " INNER JOIN property_name " + pname_label + " ON " + \
            B_cp_label + ".property_name_id = " + pname_label + ".id AND " + pname_label + ".name = '" + prop + "') ON B.id = " + \
            B_cp_label + ".compound_id\n"
        Y_prop_joins_C += join + " (compound_property " + C_cp_label + " INNER JOIN property_name " + pname_label + " ON " + \
            C_cp_label + ".property_name_id = " + pname_label + ".id AND " + pname_label + ".name = '" + prop + "') ON C.id = " + \
            C_cp_label + ".compound_id\n"
        Y_prop_joins_D += join + " (compound_property " + D_cp_label + " INNER JOIN property_name " + pname_label + " ON " + \
            D_cp_label + ".property_name_id = " + pname_label + ".id AND " + pname_label + ".name = '" + prop + "') ON D.id = " + \
            D_cp_label + ".compound_id\n"

    constant, variable1, variable2 = query["constant"], query["variable1"], query["variable2"]    
    
    X_constant, X_variable, Y_constant, Y_variable = '', '', '', ''

    if constant:
        X_constant = "AND sss(M.smiles_blob, '{}') = 1".format(constant)            
        Y_constant = "AND sss(N.smiles_blob, '{}') = 1".format(constant)

    if variable1:
        X_variable = "AND sss(I.smiles_blob, '{}') = 1".format(variable1)
    if variable2:
        Y_variable = "AND sss(J.smiles_blob, '{}') = 1".format(variable2)    

    # Need to enforce the desired join order with subqueries (views) and optimizer hints, for estimated 1000x+ improvement in query speed

    from_X = """
FROM (
    SELECT /*+ no_merge */ F.rule_smiles_id AS F_rule_smiles_id, F.constant_id AS F_constant_id, K.id AS K_id, I.smiles AS I_smiles, O.smiles AS O_smiles, {}, {}
    FROM
    rule_smiles I
    INNER JOIN from_construct E ON E.rule_smiles_id = I.id {}
    INNER JOIN constant_smiles M ON M.id = E.constant_id {}
    INNER JOIN compound A ON A.id = E.compound_id
    {}
    INNER JOIN to_construct F ON F.constant_id = E.constant_id
    INNER JOIN compound B ON B.id = F.compound_id
    {}
    INNER JOIN rule K ON K.from_smiles_id = E.rule_smiles_id AND K.to_smiles_id = F.rule_smiles_id
    INNER JOIN rule_smiles O ON O.id = F.rule_smiles_id
    WHERE E.num_frags = {}
    AND F.num_frags = {}

    UNION ALL

    SELECT /*+ no_merge */ F.rule_smiles_id AS F_rule_smiles_id, F.constant_id AS F_constant_id, -1*K.id AS K_id, I.smiles AS I_smiles, O.smiles AS O_smiles, {}, {}
    FROM
    rule_smiles I
    INNER JOIN to_construct E ON E.rule_smiles_id = I.id {}
    INNER JOIN constant_smiles M ON M.id = E.constant_id {}
    INNER JOIN compound A ON A.id = E.compound_id
    {}
    INNER JOIN from_construct F ON F.constant_id = E.constant_id
    INNER JOIN compound B ON B.id = F.compound_id
    {}
    INNER JOIN rule K ON K.to_smiles_id = E.rule_smiles_id AND K.from_smiles_id = F.rule_smiles_id
    INNER JOIN rule_smiles O ON O.id = F.rule_smiles_id
    WHERE E.num_frags = {}
    AND F.num_frags = {} 
) X
    """.format(
        XsubA, XsubB, X_variable, X_constant, X_prop_joins_A, X_prop_joins_B, str(num_frags), str(num_frags),
        XsubA, XsubB, X_variable, X_constant, X_prop_joins_A, X_prop_joins_B, str(num_frags), str(num_frags)
    )

    # Below joins are in reverse order compared to above (to -> from, rather than from -> to),
    # because we're coming at the inferred connections from the other direction

    from_Y = """
INNER JOIN (
    SELECT /*+ no_merge */ G.rule_smiles_id AS G_rule_smiles_id, G.constant_id AS G_constant_id, L.id AS L_id, J.smiles AS J_smiles, {}, {}
    FROM
    rule_smiles J
    INNER JOIN to_construct H ON H.rule_smiles_id = J.id {}
    INNER JOIN constant_smiles N ON N.id = H.constant_id {}
    INNER JOIN compound D ON D.id = H.compound_id
    {}
    INNER JOIN from_construct G ON G.constant_id = H.constant_id
    INNER JOIN compound C ON C.id = G.compound_id
    {}
    INNER JOIN rule L ON L.from_smiles_id = G.rule_smiles_id AND L.to_smiles_id = H.rule_smiles_id
    WHERE H.num_frags = {}
    AND G.num_frags = {}

    UNION ALL

    SELECT /*+ no_merge */ G.rule_smiles_id AS G_rule_smiles_id, G.constant_id AS G_constant_id, -1*L.id AS L_id, J.smiles AS J_smiles, {}, {}
    FROM
    rule_smiles J
    INNER JOIN from_construct H ON H.rule_smiles_id = J.id {}
    INNER JOIN constant_smiles N ON N.id = H.constant_id {}
    INNER JOIN compound D ON D.id = H.compound_id
    {}
    INNER JOIN to_construct G ON G.constant_id = H.constant_id
    INNER JOIN compound C ON C.id = G.compound_id
    {}
    INNER JOIN rule L ON L.to_smiles_id = G.rule_smiles_id AND L.from_smiles_id = H.rule_smiles_id
    WHERE H.num_frags = {}
    AND G.num_frags = {}    
) Y
ON Y.G_rule_smiles_id = X.F_rule_smiles_id AND X.F_constant_id != Y.G_constant_id
    """.format(YsubC, YsubD, Y_variable, Y_constant, Y_prop_joins_D, Y_prop_joins_C, str(num_frags), str(num_frags),
               YsubC, YsubD, Y_variable, Y_constant, Y_prop_joins_D, Y_prop_joins_C, str(num_frags), str(num_frags)
    )

    #select = "explain plan for SELECT " + selectA + ", " + selectB + ", " + selectC + ", " + selectD + ", X.K_id, X.I_smiles, X.O_smiles, Y.L_id, Y.J_smiles"
    select = "SELECT " + selectA + ", " + selectB + ", " + selectC + ", " + selectD + ", X.K_id, X.I_smiles, X.O_smiles, Y.L_id, Y.J_smiles"
    statement = select + from_X + from_Y

    if print_statement:
        print(statement)
        return
        
    t5 = time.time()
    cursor.execute(statement)
    t6 = time.time()
    print("Time in seconds for search: " + str(round(t6-t5, 2)))
    #return
    t7 = time.time()
    pairs = [a for a in cursor]
    t8 = time.time()
    print("Time in seconds for download: " + str(round(t8-t7, 2)))
    print("Number of pairs: " + str(len(pairs)))

    if get_flipped == True:
        flipped_pairs = flip_inferred_pairs(pairs, num_props)
        return pairs + flipped_pairs

    return pairs

# for analyzing inferred pairs raw data
def measure_raw(raw):
    print("total pairs = " + str(len(raw)))
    raw_u = [a[:-5] for a in raw]
    print("unique pairs = " + str(len(raw_u)))
    return
    
# for analyzing inferred pairs raw data
def A_minus_B(raw1, raw2):
    raw1_u = [a[:-5] for a in raw1]
    raw2_u = [a[:-5] for a in raw2]
    diff = len(set(raw1_u) - set(raw2_u))
    print('A_uniques - B_uniques = ' + str(diff))

def inferred_raw_to_csv(query : dict, raw : list, csv_path : str):
    # unpack data from cursor into a dataframe, save as csv

    props = query['REQUIRED_properties'] + query['OPTIONAL_properties']

    columns = []

    for cmpd in ['A', 'B', 'C', 'D']:
        columns += [cmpd + '_ID', cmpd]
        for prop in props:
            columns.append(cmpd + '_' + prop)

    columns += ['AB_rule_id', 'FROM', 'CONNECT', 'CD_rule_id', 'TO']

    df = pd.DataFrame(raw, columns = columns)
    df.to_csv(csv_path)

    return
