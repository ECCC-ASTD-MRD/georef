#!/usr/bin/env python3

import ctypes as ct
import os

# libgeoref = ct. cdll.LoadLibrary(os.path.join(os.environ['SSM_DEV'], 'workspace',
#                  'libgeoref_0.1.0-intel-19.0.3.199_ubuntu-18.04-skylake-64', 
#                  os.environ['SSM_DEV'], 'workspace',
#                  'libgeoref_0.1.0-intel-19.0.3.199_ubuntu-18.04-skylake-64', 'lib',
#                  'libgeoref.so'))
libgeoref = ct. cdll.LoadLibrary('/tmp/map007/188159/tmp.k5ydkEdmBp/workspace/libgeoref_0.1.0-intel-19.0.3.199_ubuntu-18.04-skylake-64/tmp/map007/188159/tmp.k5ydkEdmBp/workspace/libgeoref_0.1.0-intel-19.0.3.199_ubuntu-18.04-skylake-64/lib/libgeoref.so')
                 