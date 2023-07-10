DASH_PAGE_PREFIX = 'rep'


def create_id(name: str) -> str:
    return f'{DASH_PAGE_PREFIX}_{name}'
