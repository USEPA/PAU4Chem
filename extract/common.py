# common.py (PAU4Chem)
# !/usr/bin/env python3
# coding=utf-8


"""TODO describe this file/module."""

# Importing libraries

import yaml
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

__config = None


def config():
    """Open, load, and return the configuration yaml file."""
    global __config
    if not __config:
        with open(dir_path + '/config.yaml', mode='r') as f:
            __config = yaml.load(f, Loader=yaml.FullLoader)
    return __config
