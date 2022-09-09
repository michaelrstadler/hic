#!/usr/bin/env python

"""name.py: Declarative description of what program does

More detailed description.

TO DO:
-
"""

__author__      = "Michael Stadler"
__copyright__   = "Copyright 2022, California, USA"
__version__ = "1.0.0"

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument("-f", "--training_data_folder", type=str,  required=True,
                help="Folder containing training data in folders labeled left and right.")
    
    args = parser.parse_args()
    return args

args = parse_args()
