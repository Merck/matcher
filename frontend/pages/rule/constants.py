DASH_PAGE_PREFIX = 'rule'


def create_id(name: str) -> str:
    return f'{DASH_PAGE_PREFIX}_{name}'
