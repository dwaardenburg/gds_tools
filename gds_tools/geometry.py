import numpy as np
import gdspy as gp
import gds_tools as gdst

#=======================
# Generate a hollow box \\
#=========================================================================
# Arguments:    size            :   (x, y) specifiying box outer size   ||
#               border_width    :   width of the border                 ||
#=========================================================================
def hollow_box(size, border_width, layer = 0):

    outer = gp.Rectangle((0, 0), size)
    inner = gp.Rectangle((border_width, border_width), (size[0] - border_width, size[1] - border_width))

    hbox = gp.boolean(outer, inner, 'not', layer = layer)
    ends = {'A': (0, size[1] / 2), 'B': (size[0] / 2, 0), 'C': (size[0] / 2, size[1]), 'D': (size[0], size[1] / 2), 'CENTER': (size[0] / 2, size[1] / 2)}
    epsz = {'A': size[1], 'CENTER': None}

    return gdst.classes.GDStructure(hbox, ends, epsz)

#========================
# Generate a regular box \\
#=========================================================================
# Arguments:    size            :   (x, y) specifiying box outer size   ||
#=========================================================================
def box(size, layer = 0, datatype = 0):

    box = gp.Rectangle((0, 0), size, layer = layer, datatype = datatype)
    ends = {'A': (0, size[1] / 2), 'B': (size[0] / 2, 0), 'C': (size[0] / 2, size[1]), 'D': (size[0], size[1] / 2), 'CENTER': (size[0] / 2, size[1] / 2),
            'W': (0, 0), 'X': (size[0], 0), 'Y': (size[0], size[1]), 'Z': (0, size[1])}
    epsz = {'A': size[1], 'D': size[1], 'B': size[0], 'C': size[0], 'CENTER': None, 'W': size[0], 'X': size[0], 'Y': size[0], 'Z': size[0]}

    return gdst.classes.GDStructure(box, ends, epsz)

#===================
# Generate a circle \\
#=========================================================================
# Arguments:    r           :   radius                                  ||
#    (optional) resolution  :   number of points in [0..2pi]            ||
#=========================================================================
def circle(r, resolution = 100, layer = 0):

    circle = gp.Round((0, 0), r, number_of_points = resolution, layer = layer)
    ends = {'CENTER': (0, 0), 'BOTTOM': (0, -r)}
    epsz = {'CENTER': r, 'BOTTOM': None}

    return gdst.classes.GDStructure(circle, ends, epsz)

#================================================
# Define healing function for transmission lines \\
#=========================================================================
# Initialization:   radius      :   radius of healing circle            ||
#                   endpoint    :   (x, y) of circle center             ||
#       (optional)  npoints     :   number of points to use for circum. ||
#       (optional)  layer       :   layer to put healer on              ||
#=========================================================================
def circle(radius, endpoint, npoints = 100, layer = 0, datatype = 0):

    c = gp.Round(endpoint, radius, number_of_points = npoints, layer = layer, datatype = datatype)
    c = gdst.classes.GDStructure(c, {'A': endpoint}, {'A': 2*radius})

    return c

def triangle(side_1, side_2, angle, layer = 0):

    p = [(0, 0), (side_1, 0), gdst.funcs.VecRot(angle, (side_2, 0))]

    poly = gp.Polygon(p, layer = layer)
    ends = {gdst.alphabet[i]: p[i] for i in range(0, len(p))}
    epsz = {gdst.alphabet[i]: None for i in range(0, len(p))}

    return gdst.classes.GDStructure(poly, ends, epsz)
