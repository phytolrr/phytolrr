
import json
from bottle import response
import settings


VALID_VERSIONS = {'1': 1, '2': 2, 1: 1, 2: 2, '3': 3, 3: 3}


def get_version_arg(version):
    return VALID_VERSIONS.get(version, 3)


def response_error(message):
    response.headers['Content-Type'] = 'application/json'
    response.status = 400
    return json.dumps({
        'message': message
    })


def _set_cache_head(with_cache):
    if settings.mode == settings.MODE_DEV:
        response.headers['Cache-Control'] = 'no-store'
        return
    if with_cache:
        response.headers['Cache-Control'] = 'public, max-age=300'
        return


def response_ok(body, with_cache=False):
    response.headers['Content-Type'] = 'application/json'
    _set_cache_head(with_cache)

    response.status = 200
    if isinstance(body, dict):
        return json.dumps(body)
    else:
        return body
