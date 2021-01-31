import logging
import os
from logging.handlers import RotatingFileHandler
import bottle
import settings
from web_service.motif_service import *
from web_service.sequence_service import *
from web_service.lrr_service import *


@bottle.get('/<filepath:path>')
def sequence_index_default(filepath):
    return bottle.static_file(filepath, root='../frontend/dist')


@bottle.get('/')
def get_index():
    return sequence_index_default('lrr_db.html')


@bottle.get('/findlrr')
def get_index():
    return sequence_index_default('find_lrr.html')


@bottle.get('/about')
def get_index():
    return sequence_index_default('about.html')


def bottle_logger_plugin(fn):
    def _log_to_logger(*args, **kwargs):
        actual_response = fn(*args, **kwargs)
        logging.info('%s %s %s %s' % (bottle.request.remote_addr,
                                        bottle.request.method,
                                        bottle.request.url,
                                        bottle.response.status))
        return actual_response
    return _log_to_logger


def _init_log():
    if settings.log.path is not None:
        log_path = os.path.abspath(settings.log.path)
        if not os.path.isdir(log_path):
            logging.error(str.format("The logging path {} does not exists", log_path))
            return
        print("Begin to init logging path, path: " + log_path)
        rotating_handler = RotatingFileHandler(os.path.join(log_path, 'nsites.log'),
                                               maxBytes=settings.log.rotate_size,
                                               backupCount=50,
                                               encoding=settings.log.encoding)
        rotating_handler.setLevel(settings.log.level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        rotating_handler.setFormatter(formatter)

        logger = logging.getLogger()
        logger.setLevel(settings.log.level)
        logger.handlers = []
        logger.addHandler(rotating_handler)

        bottle.install(bottle_logger_plugin)


if __name__ == '__main__':
    import sys
    print(str.format("Configs: {}", settings.repr_all_configs()))
    _init_log()
    logging.warning("Begin to start the bottle server")
    bottle.run(host=settings.server.host, port=settings.server.port)
