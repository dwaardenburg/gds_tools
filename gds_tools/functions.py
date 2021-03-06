
import sys, os, copy, cProfile, pstats, io

import numpy as np
import gdspy as gp
import gds_tools as gdst
from recordclass import recordclass

def profile(fnc):
    
    """A decorator that uses cProfile to profile a function"""
    
    def inner(*args, **kwargs):
        
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner

def RotMat(rad):
    #==========================
    # Generate rotation matrix \\
    #=========================================================================
    # Arguments:    rad     :   radians to rotate about origin              ||
    #=========================================================================
    return np.matrix([[np.cos(rad), -np.sin(rad)], [np.sin(rad), np.cos(rad)]])

def VecRot(rad, vec, origin = (0, 0)):
    #=========================
    # Perform vector rotation \\
    #=========================================================================
    # Arguments:    rad     :   radians to rotate about origin              ||
    #               vec     :   input vector (2-x-n list)                   ||
    #=========================================================================
    return (RotMat(rad).dot(np.array(vec) - np.array(origin)) + np.array(origin)).tolist()[0]

def instruction_parse(s, args = None):
    #============================
    # Simple instructions parser \\
    #=========================================================================
    # Parses a string and converts it to a dictionary.                      ||
    #                                                                       ||
    # Arguments:    s       :   input strung                                ||
    #               args    :   dictionary with keys for variable placement ||
    #=========================================================================
    if args:
        for a in args:
            s = s.replace('{'+a+'}', str(args[a]))

    dic = {}
    key = ''
    val = ''
    pos = 'key'
    for i, c in enumerate(s):
        if c == ' ' or c == '\n' or c == '\t': # ignore whitespace
            continue
        elif c == ':':
            pos = 'val'
        elif c != ':' and c != ',' and pos == 'key':
            key += c
        elif c != ':' and c != ',' and pos == 'val':
            val += c
        elif c == ',':
            True # do nothing
        else:
            print('Error: unknown parameter, could not parse.')
            return False

        if c == ',' or (i + 1) == len(s):
            val = eval(val.replace('pi', str(np.pi)))
            dic[key] = float(val)
            key = ''
            val = ''
            pos = 'key'
            if (i + 1) == len(s):
                break

    return dic

def flatten(objectlist, endpoints, endpoint_dims, layer = 0):
    #===========================
    # FLatten a list of objects \\
    #=========================================================================
    # Flattening will cause all objects in the objectlist to be placed in   ||
    # one single layer and remove boundaries between them if there are any. ||
    # All layer information will become lost! If you just want to combine   ||
    # structures while keeping layer information, use cluster()             ||
    #                                                                       ||
    # Arguments:    objectlist      :   list of objects (GDStructure)       ||
    #               endpoints       :   dictionary of new endpoints         ||
    #               endpoint_dims   :   dictionary of new endpoint sizes    ||
    #=========================================================================
    # Define function to allow for recursive walk through list and pick out all
    # compound structures
    def stacker(inlist):

        outlist = []
        for i in inlist:
            if i.compound:
                outlist += [i] + stacker(i.compound)
            else:
                outlist += [i]

        return outlist

    objectlist = stacker(objectlist)

    ends = copy.deepcopy(endpoints)
    epsz = copy.deepcopy(endpoint_dims)

    objs = []
    for i in objectlist:
        objs.append(i.structure)

    return gdst.classes.GDStructure(gp.boolean(objs, None, 'or', layer = layer), ends, epsz)

def lattice(cell, repeat, spacing):
    #============================
    # Generate a crystal lattice \\
    #=========================================================================
    # Arguments:    cell        :   unit cell as gdspy Cell object          ||
    #               repeat      :   (n_x, n_y) vector with amount of cells  ||
    #               spacing     :   sapce between unit cells                ||
    #=========================================================================
    array = gp.CellArray(cell, repeat[0], repeat[1], spacing)
    ends = {'A': (0, spacing[1] * repeat[1] / 2), 'B': (spacing[0] * (repeat[0] - 1) / 2, spacing[1] * (repeat[1] - 1/2))}
    epsz = {'A': 0, 'B': 0}

    return gdst.classes.GDStructure(array, ends, epsz)

def lattice_cutter(lattice, objectlist, mode = 'and', layer = 0):
    #=====================================
    # Cut a lattice up using boolean \\
    #=========================================================================
    # Arguments:    lattice     :   output of lattice() function            ||
    #               objectlist  :   list of objects that intersect lattice  ||
    #    (optional) mode        :   what boolean operation to apply         ||
    #    (optional) layer       :   layer to put resulting structure on     ||
    #=========================================================================
    if type(objectlist) is not type([]):
        objectlist = [objectlist]

    for i in objectlist:
        if i.compound:
            lattice = lattice_cutter(lattice, i.compound)
        lattice.structure = gp.boolean(lattice.structure, i.structure, mode, layer = layer)

    return lattice

# def cluster(objectlist, endpoints, endpoint_dims, endpoint_directions = None, method = None, args = {}, ignore_conflict = False):
#     #============================
#     # Cluster list of structures \\
#     #=========================================================================
#     # Will combine GDStructure objects in a list so they will act as one    ||
#     # single structure. They will keep their individual properties.         ||
#     #                                                                       ||
#     # Arguments:    objectlist      :   list of GDStructure objects         ||
#     #               endpoints       :   dict of endpoints                   ||
#     #               endpoint_dims   :   dict of endpoint sizes              ||
#     #=========================================================================
#     if not ignore_conflict:
#         for i in endpoints:
#             if i in objectlist[0].endpoints:
#                 print('Error: this endpoint already exists in first element of given list of structures to combine.\nCannot continue, as overwriting existing endpoints will result in broken compound structure.\nTo fix this, give different names to your new endpoints as input to this functions.\n')
#                 return False

#     objectlist[0].endpoints = endpoints
#     objectlist[0].endpoint_dims = endpoint_dims
#     if endpoint_directions != None:
#         objectlist[0].endpoint_directions = endpoint_directions
#     else:
#         objectlist[0].endpoint_directions = {}
#     objectlist[0].compound += objectlist[1:]

#     objectlist[0].method = method
#     objectlist[0].__args_tuple__ = recordclass('args', args.keys())
#     objectlist[0].args = objectlist[0].__args_tuple__(**args)

#     for i, object in enumerate(objectlist[1:]):
#         objectlist[0].next['COMPOUND_'+str(i)] = object
#         object.prev['COMPOUND_0'] = objectlist[0]

#     return objectlist[0]

def add(cell, elements, recursive=False):
    #================================
    # Add structures to a gdspy cell \\
    #=========================================================================
    # Arguments:    cell        :   gdspy cell object                       ||
    #               elements  :   list of GDStructure objects             ||
    #=========================================================================
    if not isinstance(elements, list):
        elements = [elements]

    for element in elements:
        if isinstance(element, list):
            gdst.add(cell, element)
        else:
            if isinstance(element, gdst.classes.GDSComponent):
                for polygon in element.polygons:
                    cell.add(polygon)
            else:
                cell.add(element)

def mirror(p):
    #============================
    # Mirror points about y-axis \\
    #=========================================================================
    # Arguments:    p   :   list of (x, y) points                           ||
    #=========================================================================
    for i, val in enumerate(p):
        p[i] = (-val[0], val[1])

    return p

def symm_coords(points, mirror_x = True, mirror_y = True):
    
    if not isinstance(points, list):
        points = [points]
    
    output_points = copy.deepcopy(points)

    if mirror_y:
        for i, val in enumerate(points):
            output_points.append((-val[0], val[1]))
    if mirror_x:
        for i, val in enumerate(points):
            output_points.append((val[0], -val[1]))
    if mirror_x and mirror_y:
        for i, val in enumerate(points):
            output_points.append((-val[0], -val[1]))
    return output_points

def save(cell, filename, unit = 1e-6, precision = 1e-9):
    #=====================
    # Save cell to a file \\
    #=========================================================================
    # Arguments:    cell        :   gdspy cell object or a list of cells    ||
    #               filename    :   filename to write to (relative path)    ||
    #=========================================================================
    writer = gp.GdsWriter(filename, unit = unit, precision = precision)
    if type(cell) == type([]):
        for cell_item in cell:
            writer.write_cell(cell_item)
    else:
        writer.write_cell(cell)

    return writer.close()

def biquadratic_func(x):
    return x ** 2 * (2 - x ** 2)

def rotate_reference_cell(reference, angle, center = (0, 0)):
    dx = np.cos(angle) * (reference.origin[0] - center[0]) - np.sin(angle) * (reference.origin[1] - center[1]) + center[0]
    dy = np.sin(angle) * (reference.origin[0] - center[0]) + np.cos(angle) * (reference.origin[1] - center[1]) + center[1]
    angle_deg = np.degrees(angle)
    reference.rotation += angle_deg
    reference.translate(dx - reference.origin[0], dy - reference.origin[1])

def inside(points, cellref, dist, nop = 3, precision = 0.001):
    #=====================
    # Save cell to a file \\
    #=========================================================================
    # Arguments:    points      :   list of points to check                 ||
    #               cellref     :   gdspy cell reference object             ||
    #               dist        :   distance from points to search          ||
    #               nop         :   number of probe points within dist      ||
    #               precision   :   gdspy.inside precision parameter        ||
    #=========================================================================
    # Force uneven
    if nop % 2 == 0:
        nop += 1

    search_ps = []
    for p in points:
        px = np.linspace(p[0] - dist/2, p[0] + dist/2, nop)
        py = np.linspace(p[1] - dist/2, p[1] + dist/2, nop)
        search_ps.append([[i, j] for i in px for j in py])

    return gp.inside(search_ps, cellref, precision = precision)

def convert_to_dxf(filename):

    print("-- Converting to DXF --")
    # Convert GDS to DXF with Klayout
    os.system('/Applications/klayout.app/Contents/MacOS/klayout -zz -rd input="{}.gds" -rd output="{}.dxf" -r convert.rb'.format(filename, filename))

def bounding_box_center(object):
    bounding_box = object.get_bounding_box()
    bounding_box_x = (bounding_box[1][0] + bounding_box[0][0]) / 2
    bounding_box_y = (bounding_box[1][1] + bounding_box[0][1]) / 2
    return (bounding_box_x, bounding_box_y)