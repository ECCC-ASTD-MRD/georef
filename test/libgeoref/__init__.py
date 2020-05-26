#!/usr/bin/env python3

import ctypes as ct
import os

libgeoref = ct. cdll.LoadLibrary(os.path.join(os.environ['SSM_DEV'], 'workspace',
                 '_-intel-19.0.3.199_ubuntu-18.04-skylake-64', 'lib',
                 'libgeoref.so'))
                 