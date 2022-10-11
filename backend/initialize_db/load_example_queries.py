import sys, os

# backend_api.py is one level up
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from backend_api import snap_write

import json
import asyncio
from models import ExampleQuery

async def load_examples(json_file, schema='public'):
    with open(json_file) as f:
        examples = json.load(f)
        for ex in examples:
            ex = ExampleQuery.parse_obj(ex)
            result = await snap_write(ex, schema=schema)

if __name__ == "__main__":
    assert len(sys.argv) > 1 and len(sys.argv) < 4, '\n\nload_example_queries example usage:\n\npython load_example_queries.py example_queries.json public\n'

    loop = asyncio.get_event_loop()
    work = [loop.create_task(load_examples(sys.argv[1], schema=sys.argv[2]))]
    loop.run_until_complete(asyncio.wait(work))
    loop.close()