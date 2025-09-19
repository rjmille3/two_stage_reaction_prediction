"""Various constants for use by the kernel modules"""

import logging

"""Delimiter for feature and matrix files"""
TEXT_DELIM = "\t"

"""Delimiter for atom pair weight specification strings and feature index:count mappings"""
KEY_DELIM = ":"
ITEM_DELIM = ","

"""Prefix to identify feature:index mapping rows of feature dictionary text files"""
FEATURE_PREFIX = "#"

"""Define a really small number to check convergence and to avoid division by zero."""
EPSILON = 1E-12;
