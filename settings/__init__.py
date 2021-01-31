# -*- coding: utf-8 -*
# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019-2021 phytolrr.com. All rights reserved.

import configparser
from settings import db
from settings import log
from settings import server

MODE_DEV = 'development'
MODE_PRODUCT = 'production'

mode = MODE_DEV

def read_in_configs(file_path):
    parser = configparser.ConfigParser()
    parser.read([file_path])
    if parser.has_section('database'):
        for key, value in parser.items('database'):
            setattr(db, key, value)
        if db.db_type == 'mysql':
            db.db_type = 'mysql+pymysql'
    if parser.has_section('general'):
        if parser.has_option('general', 'mode'):
            global mode
            mode = parser.get('general', 'mode')
    if parser.has_section('log'):
        for key, value in parser.items('log'):
            setattr(log, key, value)
    if parser.has_section('server'):
        for key, value in parser.items('server'):
            setattr(server, key, value)


def check_configs():
    pass


def repr_config(module):
    des = str.format("\n{}:\n", module.__name__)
    for item in dir(module):
        if item.startswith("__"):
            continue
        if item == 'password':
            des += str.format("    {}={}\n", item, '****')
        else:
            des += str.format("    {}={}\n", item, getattr(module, item))
    return des


def repr_all_configs():
    config_description = str.format("\nmode: {}\n", mode)
    config_description += repr_config(server)
    config_description += repr_config(db)
    config_description += repr_config(log)

    return config_description


read_in_configs("settings.ini")
check_configs()
