import redis
from settings import cache
from web_service import service_utils

pool = redis.ConnectionPool(host=cache.host, port=cache.port)
r = redis.Redis(connection_pool=pool)

def traffic_controller(callback):
    def wrapper(*args, **kwargs):
        body = callback(*args, **kwargs)
        return body
    return wrapper
