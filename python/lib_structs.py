from ctypes import *

# From RPN.h
class TRPNFile(Structure):
    _fields_ = [
        ("CId", c_char_p),
        ("Name", c_char_p),
        ("Open", c_int),
        ("Id", c_uint),
        ("NRef", c_int),
        ("Mode", c_char)
    ]

# From RPN.h
class TRPNHeader(Structure):
    _fields_ = [
        ("File", POINTER(TRPNFile)),
        ("FID", c_int),
        ("KEY", c_int),
        ("DATEO", c_int),
        ("DATEV", c_int),
        ("DEET", c_int),
        ("NPAS", c_int),
        ("NBITS", c_int),
        ("DATYP", c_int),
        ("IP1", c_int),
        ("IP2", c_int),
        ("IP3", c_int),
        ("NI", c_int),
        ("NJ", c_int),
        ("NK", c_int),
        ("NIJ", c_int),
        ("IG", c_int * 4),
        ("IGREF", c_int * 4),
        ("XG", c_float * 4),
        ("XGREF", c_float * 4),
        ("SWA", c_int),
        ("LNG", c_int),
        ("DLTF", c_int),
        ("UBC", c_int),
        ("EX1", c_int),
        ("EX2", c_int),
        ("EX3", c_int),
        ("TYPVAR", c_char * 3),
        ("NOMVAR", c_char * 5),
        ("ETIKET", c_char * 13),
        ("GRTYP", c_char * 2),
        ("GRREF", c_char * 2)
    ]

class TGeoRef(Structure):
    pass
# # From GeoRef.h
# class TGeoRef(Structure):
#     _fields_ = [
#         ("Name", c_char_p),
#         ("NRef", c_int),
#         ("RefFrom", POINTER(TGeoRef)),
#         ("Subs", POINTER(POINTER(TGeoRef))),
#         ("NbSub", c_int),
#         ("Type", c_int),
#         ("BD", c_int),
#         ("NX", c_int),
#         ("NY", c_int),
#         ("X0", c_int),
#         ("Y0", c_int),
#         ("X1", c_int),
#         ("Y1", c_int),
#         ("Extension", c_int),
#         ("GRTYP", c_char * 3),
#         ("Hemi", c_int),
#         ("NbSet", c_int),
#         ("Sets", POINTER(TGridSet)),
#         ("LastSet", POINTER(TGridSet)),
#         ("NIdx", c_uint),
#         ("Idx", POINTER(c_uint)),
#         ("Lat", POINTER(c_double)),
#         ("Lon", POINTER(c_double)),
#         ("AX", POINTER(c_float)),
#         ("AY", POINTER(c_float)),
#         ("NCX", POINTER(c_float)),
#         ("NCY", POINTER(c_float)),
#         ("Hgt", POINTER(c_float)),
#         ("Wght", POINTER(c_double)),
#         ("String", c_char_p),
#         ("Transform", POINTER(c_double)),
#         ("InvTransform", POINTER(c_double)),
#         ("RotTransform", POINTER(TRotationTransform)),
#         ("Options", TGeoOptions),
#         ("Project", CFUNCTYPE(c_int, POINTER(TGeoRef), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int)),
#         ("UnProject", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), c_double, c_double, c_int, c_int)),
#         ("XY2LL", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int)),
#         ("LL2XY", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int)),
#         ("RPNHead", TRPNHeader),
#         ("i1", c_int),
#         ("i2", c_int),
#         ("j1", c_int),
#         ("j2", c_int),
#         ("mask", POINTER(c_int)),
#         ("mymaskgrid", POINTER(TGeoRef)),
#         ("mymaskgridi0", c_int),
#         ("mymaskgridi1", c_int),
#         ("mymaskgridj0", c_int),
#         ("mymaskgridj1", c_int)
#     ]

# From GeoRef.h
class TGeoOptions(Structure):
    _fields_ = [
        ('Interp', c_int),          # TRef_InterpR: Interpolation degree
        ('Extrap', c_int),          # TRef_ExtrapR: Extrapolation method
        ('InterpVector', c_int),     # TRef_InterpV: Vector interpolation method
        ('Combine', c_int),          # TRef_Combine: Aggregation type
        ('Transform', c_int),        # Apply transformation or stay within master referential
        ('CIndex', c_int),          # C Indexing (starts at 0)
        ('Symmetric', c_int),        # Symmetric
        ('WeightNum', c_int),       # Weight number
        ('Segment', c_int),         # How much segmentation (Conservatives/Geometric modes)
        ('Sampling', c_int),        # Sampling interval
        ('PolarCorrect', c_char),   # Apply polar corrections
        ('VectorMode', c_char),     # Process data as vector
        ('DistTreshold', c_float),  # Distance treshold for point clouds
        ('LonRef', c_float),        # Longitude referential (-180.0,0.0)
        ('NoData', c_float),        # NoData Value (Default: NaN)
        ('Table', POINTER(c_double)), # Data table to check of values to check for
        ('lutDef', POINTER(POINTER(c_double))), # Lookup table
        ('lutSize', c_int),         # Number of lookup elements
        ('lutDim', c_int),          # Dimension of the lookup elements
        ('Ancilliary', POINTER(c_double)) # Pre calculated field (ex: variance, average,...)
    ]

# From GeoRef.h
class TRotationTransform(Structure):
    _fields_ = [
        ("Lat", c_double),
        ("Lon", c_double),
        ("Angle", c_double),
        ("SinTheta", c_double),
        ("CosTheta", c_double),
        ("SinPhi", c_double),
        ("CosPhi", c_double)
    ]

# From GeoRef.h
class TGeoZone(Structure):
    _fields_ = [
        ("npts", c_int),
        ("x", POINTER(c_double)),
        ("y", POINTER(c_double)),
        (("idx", POINTER(c_int)))
    ]

# From Vect.h
class Vect2d(Structure):
    _fields_ = [
        ("X", c_double),
        ("Y", c_double)
    ]

# From Vect.h
class Vect3d(Structure):
    _fields_ = [
        ("X", c_double),
        ("Y", c_double),
        ("Z", c_double)
    ]

# From Array.h
class T3DArray(Structure):
    _fields_ = [
        ("Data", POINTER(c_double * 3)),  # Vect3d is a double[3]
        ("Value", c_double),
        ("Size", c_ulong)
    ]

# From ZRefInterp.h
class TZRefInterp(Structure):
    _fields_ = [
        ("ZRefSrc", POINTER(TZRef)),      # Source vertical reference
        ("ZRefDest", POINTER(TZRef)),     # Destination vertical reference
        ("Indexes", POINTER(c_int32)),    # Interpolation array indices
        ("NIJ", c_int32),                 # 2D dimensions
        ("Same", c_int32)                 # Flag indicating source and destination are the same
    ]
    
# From Def.h
class TDef(Structure):
    _fields_ = [
        ("Buffer", POINTER(c_double)),
        ("Aux", POINTER(c_double)),
        ("Accum", POINTER(c_int32)),
        ("Mask", c_char_p),
        ("Data", c_char_p * 4),
        ("Mode", c_char_p),
        ("Dir", c_char_p),
        ("Pres", POINTER(c_float)),
        ("PresLS", POINTER(c_float)),
        ("Height", POINTER(c_float)),
        ("Sub", POINTER(c_float)),
        ("Pick", POINTER(c_void_p)),  # OGRGeometryH
        ("Poly", POINTER(c_void_p)),  # OGRGeometryH
        ("Segments", c_void_p),       # TList*
        ("NoData", c_double),
        ("Type", c_int32),            # TDef_Type
        ("NI", c_int32),
        ("NJ", c_int32),
        ("NK", c_int32),
        ("NC", c_int32),
        ("NIJ", c_int32),
        ("Idx", c_int32),
        ("CellDim", c_int32),
        ("CoordLimits", (c_double * 2) * 2),
        ("Limits", (c_int32 * 2) * 3),
        ("Level", c_int32),
        ("Sample", c_int32),
        ("SubSample", c_int32),
        ("Alias", c_char)
    ]

# From GeoScan.h
class TGeoScan(Structure):
    _fields_ = [
        ("X", c_double),          # X coordinate
        ("Y", c_double),          # Y coordinate
        ("I", c_int32),           # Grid I index
        ("J", c_int32),           # Grid J index
        ("D", c_double),          # Distance
        ("S", c_double),          # Segment length
        ("N", c_int32),           # Number of points
        ("C", c_int32),           # Current point
        ("DI", c_int32),          # I increment
        ("DJ", c_int32),          # J increment
        ("Side", c_int32),        # Side number
        ("Done", c_int32),        # Scan done flag
        ("In", c_int32),          # Inside flag
        ("NI", c_int32),          # Grid dimension I
        ("NJ", c_int32),          # Grid dimension J
        ("Wrap", c_int32)         # Wrapping flag
    ]

# From GeoSet.h
class TGeoSet(Structure):
    _fields_ = [
        ("NbZone", c_int32),                # Number of zones
        ("Zone", POINTER(c_int32)),         # Zone array
        ("NbIdx", c_int32),                 # Number of indices
        ("Idx", POINTER(c_int32)),          # Index array
        ("NbYY", c_int32),                  # Number of YY points
        ("YY", POINTER(c_int32)),           # YY array
        ("YYLat", POINTER(c_double)),       # YY latitudes
        ("YYLon", POINTER(c_double)),       # YY longitudes
        ("YYDist", POINTER(c_double)),      # YY distances
        ("YYW", POINTER(c_double)),         # YY weights
        ("YYIdxS", POINTER(c_int32)),       # YY source indices
        ("YYIdxD", POINTER(c_int32)),       # YY destination indices
        ("YYNb", c_int32),                  # Number of YY points
        ("YYNbT", c_int32),                 # Total number of YY points
        ("YYNbIdx", c_int32),               # Number of YY indices
        ("YYType", c_int32),                # YY type
        ("YYDone", c_int32),                # YY done flag
        ("NbXY", c_int32),                  # Number of XY points
        ("XY", POINTER(c_int32)),           # XY array
        ("XYLat", POINTER(c_double)),       # XY latitudes
        ("XYLon", POINTER(c_double)),       # XY longitudes
        ("XYDist", POINTER(c_double)),      # XY distances
        ("XYW", POINTER(c_double)),         # XY weights
        ("XYIdxS", POINTER(c_int32)),       # XY source indices
        ("XYIdxD", POINTER(c_int32)),       # XY destination indices
        ("XYNb", c_int32),                  # Number of XY points
        ("XYNbT", c_int32),                 # Total number of XY points
        ("XYNbIdx", c_int32),               # Number of XY indices
        ("XYType", c_int32),                # XY type
        ("XYDone", c_int32)                 # XY done flag
    ]

# From ZRef.h
class TZRef(Structure):
    _fields_ = [
        ("Name", c_char_p),
        ("VGD", c_char_p),
        ("Levels", POINTER(c_float)),
        ("A", POINTER(c_float)),
        ("B", POINTER(c_float)),
        ("P0", POINTER(c_float)),
        ("P0LS", POINTER(c_float)),
        ("PCube", POINTER(c_float)),
        ("LevelNb", c_int32),
        ("NRef", c_int32),
        ("Version", c_int32),
        ("Type", c_int32),
        ("SLEVE", c_int32),
        ("POff", c_float),
        ("PTop", c_float),
        ("PRef", c_float),
        ("RCoef", c_float * 2),
        ("ETop", c_float),
        ("Style", c_int32)  # TZRef_IP1Mode
    ]

# From Coord.h
class TCoord(Structure):
    _fields_ = [
        ("Lat", c_double),          # Latitude
        ("Lon", c_double),          # Longitude
        ("Elev", c_double),         # Elevation
        ("X", c_double),            # X coordinate
        ("Y", c_double),            # Y coordinate
        ("Z", c_double),            # Z coordinate
        ("I", c_int32),             # Grid I index
        ("J", c_int32),             # Grid J index
        ("K", c_int32),             # Grid K index
        ("T", c_int32),             # Time index
        ("H", c_int32),             # Height index
        ("Level", c_float),         # Level value
        ("Time", c_double),         # Time value
        ("Height", c_float),        # Height value
        ("Pressure", c_float),      # Pressure value
        ("Distance", c_double),     # Distance value
        ("Value", c_float),         # Data value
        ("Valid", c_int32)          # Validity flag
    ]

# From OGR_Layer.h
class OGR_Layer(Structure):
    _fields_ = [
        ("Name", c_char * 128),             # Layer name
        ("Type", c_int32),                  # Layer type
        ("Format", c_char * 32),            # Format name
        ("NFeature", c_int32),              # Number of features
        ("Feature", POINTER(c_void_p)),     # Array of feature pointers
        ("FID", POINTER(c_int32)),          # Feature IDs
        ("Selection", POINTER(c_int32)),    # Selection array
        ("NbSelected", c_int32),            # Number of selected features
        ("Proj", c_char * 2048),            # Projection string
        ("Extent", c_double * 4),           # Spatial extent [minx,miny,maxx,maxy]
        ("Modified", c_int32),              # Modified flag
        ("Mode", c_int32),                  # Layer mode
        ("Active", c_int32),                # Active flag
        ("Transparent", c_int32),           # Transparency flag
        ("Width", c_int32),                 # Line width
        ("Color", c_uint32),                # Color value
        ("Alpha", c_float),                 # Alpha transparency
        ("Icon", c_int32),                  # Icon type
        ("Size", c_int32),                  # Icon size
        ("Outline", c_int32),               # Outline flag
        ("Fill", c_int32),                  # Fill flag
        ("Dash", c_int32),                  # Dash pattern
        ("Stipple", c_int32),               # Stipple pattern
        ("Bkg", c_uint32),                  # Background color
        ("BkgAlpha", c_float),              # Background alpha
        ("Desc", c_char * 2048)             # Description
    ]

# From OGR_Core.h
class OGREnvelope(Structure):
    _fields_ = [
        ("MinX", c_double),
        ("MaxX", c_double),
        ("MinY", c_double),
        ("MaxY", c_double)
    ]

# From GridSet.h
class TGridSet(Structure):
    _fields_ = [
        ("RefFrom", POINTER(TGeoRef)),
        ("zones", TGeoZone * 5), #TODO: define SET_NZONES 5
        ("flags", c_int),
        ("x", POINTER(c_double)), ("y", POINTER(c_double)),
        (("mask_in", POINTER(c_int))), (("mask_out", POINTER(c_int))),
        ("yin_maskout", POINTER(c_float)), ("yan_maskout", POINTER(c_float)),
        (("yinlat", POINTER(c_int))), (("yinlon", POINTER(c_int))), (("yanlat", POINTER(c_int))), (("yanlon", POINTER(c_int))),
        (("yin2yin_lat", POINTER(c_int))), (("yin2yin_lon", POINTER(c_int))), (("yan2yin_lat", POINTER(c_int))), (("yan2yin_lon", POINTER(c_int))),
        (("yin2yan_lat", POINTER(c_int))), (("yin2yan_lon", POINTER(c_int))), (("yan2yan_lat", POINTER(c_int))), (("yan2yan_lon", POINTER(c_int))),
        (("yin2yin_x", POINTER(c_int))), (("yin2yin_y", POINTER(c_int))), (("yan2yin_x", POINTER(c_int))), (("yan2yin_y", POINTER(c_int))),
        (("yin2yan_x", POINTER(c_int))), (("yin2yan_y", POINTER(c_int))), (("yan2yan_x", POINTER(c_int))), (("yan2yan_y", POINTER(c_int))),
        ("yincount_yin", c_int), ("yancount_yin", c_int), ("yincount_yan", c_int), ("yancount_yan", c_int),
        ("wts", c_int),
        ("y", POINTER(c_double)),
        (("mask", POINTER(c_int))), (("idx", POINTER(c_int)))
    ]


TGeoRef._fields_ = [
    ("Name", c_char_p),
    ("NRef", c_int),
    ("RefFrom", POINTER(TGeoRef)),
    ("Subs", POINTER(POINTER(TGeoRef))), #TODO: pointer of pointer
    ("NbSub", c_int),
    ("Type", c_int),
    ("BD", c_int),
    ("NX", c_int),("NY", c_int),("X0", c_int),("Y0", c_int),("X1", c_int),("Y1", c_int),
    ("Extension", c_int),
    ("GRTYP", c_char * 3),
    ("Hemi", c_int),
    ("NbSet", c_int),
    ("Sets", POINTER(TGridSet)),("LastSet", POINTER(TGridSet)),
    ("NIdx", c_uint),("Idx", POINTER(c_uint)),
    ("Lat", POINTER(c_double)),("Lon", POINTER(c_double)),
    ("AX", POINTER(c_float)),("AY", POINTER(c_float)),("NCX", POINTER(c_float)),("NCY", POINTER(c_float)),("Hgt", POINTER(c_float)),
    ("Wght", POINTER(c_double)),
    ("String", c_char_p),
    ("Transform", POINTER(c_double)),("InvTransform", POINTER(c_double)),
    ("RotTransform", POINTER(TRotationTransform)),
    ("GCPTransform", c_void_p),  # void pointer for now
    ("TPSTransform", c_void_p),
    ("RPCTransform", c_void_p),
    ("QTree", c_void_p),  # or create TQTree structure if available
    ("LLExtent", OGREnvelope),  # already defined in your structures
    ("Function", c_void_p),  # or appropriate function pointer type
    ("InvFunction", c_void_p),
    ("Spatial", c_void_p),
    ("Loc", TCoord),  # already defined in your structures
    ("CTH", c_double),
    ("STH", c_double),
    ("ResR", c_double),
    ("ResA", c_double),
    ("R", c_int),
    ("Options", TGeoOptions),

    # TODO: Is it the correct way ?
    # typedef int    (TGeoRef_Project)   (struct TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform);
    ("Project", CFUNCTYPE(c_int, POINTER(TGeoRef), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int)),

    # typedef int    (TGeoRef_UnProject) (struct TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform);
    ("UnProject", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), c_double, c_double, c_int, c_int)),

    # typedef int    (TGeoRef_LL2XY)     (struct TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb);
    ("XY2LL", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int)),

    # typedef int    (TGeoRef_XY2LL)     (struct TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
    ("LL2XY", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int)),
                         
#    TGeoRef_Value     *Value;                              ///< Valeur a un point xy
#    TGeoRef_Height    *Height;

    ("RPNHead", TRPNHeader),
    ("i1", c_int),("i2", c_int),("j1", c_int),("j2", c_int),
    ("mask", POINTER(c_int)),
    ("mymaskgrid", POINTER(TGeoRef)),
    ("mymaskgridi0", c_int),("mymaskgridi1", c_int),
    ("mymaskgridj0", c_int),("mymaskgridj1", c_int),
    ("NC_Id", c_int),
    ("NC_NXDimId", c_int),
    ("NC_NYDimId", c_int),
    ]
