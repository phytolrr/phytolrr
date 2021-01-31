# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019-2021 phytolrr.com. All rights reserved.

from settings import db as settings
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from contextlib import contextmanager


__all__ = ['Base', 'Session', 'session_scope', 'query_session', 'engine']

Base = None
engine = None
Session = None


def reconnect():
    global Base, engine, Session

    if engine is not None:
        engine.dispose()

    Base = declarative_base()
    if settings.db_type == 'sqlite3':
        engine = sqlalchemy.create_engine("sqlite://")
    else:
        url = sqlalchemy.engine.url.URL(settings.db_type, settings.user, settings.password,
                                        settings.host, settings.port, settings.database)
        engine = sqlalchemy.create_engine(url, pool_size=50, pool_recycle=300)
    Session = sessionmaker(bind=engine)


reconnect()

@contextmanager
def session_scope():
    """Provide a transactional scope around a series of operations."""
    session = Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.expunge_all()
        session.close()

@contextmanager
def query_session():
    """Provide a transactional scope around a series of operations."""
    session = Session()
    try:
        yield session
    except:
        raise
    finally:
        session.expunge_all()
        session.close()
