"""Test that database built successfully, 
also test ability to query database from backend."""

from backend.backend_api import get_matcher_conn
import pytest

# Success of this test depends on using the "test" dataset in backend/initialize_db

@pytest.mark.asyncio
async def test_db_build_complete():

    conn = await get_matcher_conn(schema='public')

    result = await conn.fetch('SELECT num_compounds, num_rules, num_pairs, num_rule_environments, num_rule_environment_stats FROM dataset')
    result = result[0]
    expected = {
        'num_compounds': 16,
        'num_rules': 14,
        'num_pairs': 108,
        'num_rule_environments': 100,
        'num_rule_environment_stats': 106,
    }
    for key in result.keys():
        assert result[key] == expected[key]

    expected = {
        'from_construct': 18,
        'to_construct': 18,
        'snapshot': 8,
    }
    for table in expected.keys():
        count = await conn.fetchval(f"SELECT count(id) FROM {table}")
        assert expected[table] == count

    await conn.close()
