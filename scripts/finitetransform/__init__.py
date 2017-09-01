# -*- coding: utf-8 -*-
"""
Finite Transform Library (FTL) is a collection of algorithms for image processing and reconstruction.

There are a number of sub-modules:

radon - this houses the finite Radon transform algorithms
mojette - this houses the aperiodic Radon transform algorithms
tomo - this houses the traditional reconstruction algorithms

@author: shakes
"""
import os.path as osp

pkg_dir = osp.abspath(osp.dirname(__file__))
data_dir = osp.join(pkg_dir, '../data')

__version__ = '0.1'

del osp
