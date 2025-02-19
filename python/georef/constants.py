"""Constants used by the georef package."""

# Error codes
GEOREF_SUCCESS = 1
GEOREF_ERROR = 0

# Grid types
GRID_UNDEF = 0
GRID_LATLON = 1
GRID_PS = 2
GRID_YINYANG = 3
GRID_MASKED = 4

# Interpolation types
IR_UNDEF = 0
IR_NEAREST = 1
IR_LINEAR = 2
IR_CUBIC = 3
IR_AVERAGE = 4

# Extrapolation types
ER_UNDEF = 0
ER_NONE = 1
ER_VALUE = 2
ER_NEAREST = 3
ER_MAXIMUM = 4
ER_MINIMUM = 5

# Vector interpolation methods
IV_UNDEF = 0
IV_FAST = 1
IV_WITHIN = 2
IV_INTERSECT = 3
IV_CENTROID = 4 