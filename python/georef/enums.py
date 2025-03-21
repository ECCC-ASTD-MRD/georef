import enum

# Using IntEnum is important in the context of ctypes because it allows us to
# to pass instances of this enum to functions marked as taking an int argument.
class RefInterpR(enum.IntEnum):
   IR_UNDEF                          = 0
   IR_NEAREST                        = 1
   IR_LINEAR                         = 2
   IR_CUBIC                          = 3
   IR_NORMALIZED_CONSERVATIVE        = 4
   IR_CONSERVATIVE                   = 5
   IR_MAXIMUM                        = 6
   IR_MINIMUM                        = 7
   IR_SUM                            = 8
   IR_AVERAGE                        = 9
   IR_VARIANCE                       = 10
   IR_SQUARE                         = 11
   IR_NORMALIZED_COUNT               = 12
   IR_COUNT                          = 13
   IR_VECTOR_AVERAGE                 = 14
   IR_NOP                            = 15
   IR_ACCUM                          = 16
   IR_BUFFER                         = 17
   IR_SUBNEAREST                     = 18
   IR_SUBLINEAR                      = 19
   IV_FAST                           = 20
   IV_WITHIN                         = 21
   IV_INTERSECT                      = 22
   IV_CENTROID                       = 23
   IV_ALIASED                        = 24
   IV_CONSERVATIVE                   = 25
   IV_NORMALIZED_CONSERVATIVE        = 26
   IV_POINT_CONSERVATIVE             = 27
   IV_LENGTH_CONSERVATIVE            = 28
   IV_LENGTH_NORMALIZED_CONSERVATIVE = 29
   IV_LENGTH_ALIASED                 = 30

class RefExtrapR(enum.IntEnum):
    ER_UNDEF   = 0
    ER_MAXIMUM = 1
    ER_MINIMUM = 2
    ER_VALUE   = 3
    ER_ABORT   = 4

class RefInterpVR(enum.IntEnum):
    IV_UNDEF                          = 0
    IV_FAST                           = 1
    IV_WITHIN                         = 2
    IV_INTERSECT                      = 3
    IV_CENTROID                       = 4
    IV_ALIASED                        = 5
    IV_CONSERVATIVE                   = 6
    IV_NORMALIZED_CONSERVATIVE        = 7
    IV_POINT_CONSERVATIVE            = 8
    IV_LENGTH_CONSERVATIVE           = 9
    IV_LENGTH_NORMALIZED_CONSERVATIVE = 10
    IV_LENGTH_ALIASED                = 11

class RefCombineR(enum.IntEnum):
    CB_REPLACE   = 0
    CB_MIN       = 1
    CB_MAX       = 2
    CB_SUM       = 3
    CB_AVERAGE   = 4


class TDefType(enum.IntEnum):
    TD_UNKNOWN = 0
    TD_BINARY  = 1
    TD_UBYTE   = 2
    TD_BYTE    = 3
    TD_UINT16  = 4
    TD_INT16   = 5
    TD_UINT32  = 6
    TD_INT32   = 7
    TD_UINT64  = 8
    TD_INT64   = 9
    TD_FLOAT32 = 10
    TD_FLOAT64 = 11
