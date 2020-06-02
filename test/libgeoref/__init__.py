#!/usr/bin/env python3

import ctypes as ct
import os

X = 'libgeoref_0.1.0-' + os.environ['COMP_ARCH'] + '_' + os.environ['BASE_ARCH']

libgeoref = ct.cdll.LoadLibrary(
    os.path.join(os.environ['SSM_DEV'], 'workspace', X,
                 os.environ['SSM_DEV'][1:], 'workspace', X, 'lib',
                 'libgeoref.so'))
