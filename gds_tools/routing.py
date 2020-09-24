import time, heapq, copy

import numpy as np
import numpy.ma as ma
import scipy.spatial as sp
import gdspy as gp
import gds_tools as gdst

#===========================
# Make a transmission line  \\
#=========================================================================
# Arguments:    dxy         :   (dx, dy), lengths of the line           ||
#               width       :   width of the line                       ||
#    (optional) width_end   :   different width at the end              ||
#    (optional) layer       :   choose layer to put spiral onto         ||
#=========================================================================
def line(dxy, width, width_end = False, layer = 0, datatype = 0):

    width_end = width if (not width_end) else width_end

    p = [(0, 0), dxy]
    w = [width, width_end]

    poly = gp.PolyPath(p, w, layer = layer, datatype = datatype)
    ends = {'A': p[0], 'B': p[-1]}
    epsz = {'A': w[0], 'B': w[-1]}

    return gdst.classes.GDStructure(poly, ends, epsz)

#================================================
# Generate a meander \\ Author: Marta Pita Vidal \\
#=========================================================================
# Arguments:    n               :   number of meanders                  ||
#               width           :   width of the meander path           ||
#               radius          :   radius of the corner circles        ||
#               meander_length  :   length from corner to corner        ||
#               connect_length  :   length of initial transmission line ||
#    (optional) layer           :   layer to put meander on             ||
#=========================================================================
def meander(n, width, radius, meander_length, connect_length, layer = 0):

    # Aliases
    w = width
    d = radius
    N = n
    position = (0, 0)
    xMeander = meander_length / 2
    yBase = connect_length

    # Start path
    meander_path = gp.Path(w, (0, 0))
    meander_path.segment(yBase, np.pi/2, layer = layer).turn(w/2, 'r', layer = layer).segment(xMeander - w/2, layer = layer).turn(d, 'll', layer = layer).segment(xMeander, layer = layer)

    # Generate periodic meanders
    for i in range(N - 1):
        meander_path.segment(xMeander, layer = layer).turn(d, 'rr', layer = layer).segment(2 * xMeander, layer = layer).turn(d, 'll', layer = layer).segment(xMeander, layer = layer)

    # Finalize
    meander_path.segment(xMeander, layer = layer).turn(d, 'rr', layer = layer).segment(xMeander - w/2, layer = layer).turn(w/2, 'l', layer = layer).segment(yBase, '+y', layer = layer)

    # Create object
    struct = gdst.classes.GDStructure(meander_path, {'A': (0, 0), 'B': (meander_path.x, meander_path.y)}, {'A': w, 'B': w})

    return struct

#===================================
# Make a transmission line splitter \\
#=========================================================================
# Arguments:    n       :   number of splits                            ||
#               width   :   width of the paths                          ||
#               length  :   total length from start to end              ||
#               space   :   spacing between splits                      ||
#    (optional) layer   :   choose layer to put split onto              ||
#=========================================================================
def splitter(n, width, length, space, layer = 0):

    p = [(0, 0), (length / 2 - width / 2, 0), (length / 2 - width / 2, (n*width + (n-1)*space)/2 - width / 2)]

    ends = {'A': (0, -width / 2)}
    epsz = {'A': width}

    for i in range(0, n):
        p += [(length, (n*width + (n-1)*space)/2 - i*(space + width) - width / 2)]

        ends[gdst.alphabet[i + 1]] = (length, p[-1][1] - width / 2)
        epsz[gdst.alphabet[i + 1]] = width

        p += [(length, (n*width + (n-1)*space)/2 - i*(space + width) - width - width / 2)]
        p += [(length / 2 + width / 2, (n*width + (n-1)*space)/2 - i*(space + width) - width - width / 2)]
        p += [(length / 2 + width / 2, (n*width + (n-1)*space)/2 - i*(space + width) - width - space - width / 2)]

    p = p[:-1]

    p += [(length / 2 - width / 2, (n*width + (n-1)*space)/2 - (n-1)*(space + width) - width - width / 2)]
    p += [(length / 2 - width / 2, -width)]
    p += [(0, -width)]

    poly = gp.Polygon(p, layer = layer)

    return gdst.classes.GDStructure(poly, ends, epsz)

# Helper function, returns indeces of neighbouring points (u, d, l, r) in a grid
def lookaround(index, row_len, pref = 'y'):

    n = []
    n.append(index - row_len)
    n.append(index + row_len)
    if index % row_len != 0:
        n.append(index - 1)
    n.append(index + 1)

    if pref == 'x':
        n = n[::-1]

    return n

# Autoroute using Lee's routing algorithm
def router(cell, fr_str, fr_ep, to_str, to_ep, width, bmul, grid_s = 1, xr = False, yr = False, uniform_width = False, precision = 0.001, pref = 'y', switch_pref = False, layer = 0, debug = False, nop = 21, dist_multi = 2, pathmethod = 'poly', detect = 'native'):

    fr = fr_str.endpoints[fr_ep]
    to = to_str.endpoints[to_ep]

    border_s = bmul * grid_s
    box = [fr, to]

    # Make sure first box coord is always the top-left corner and add additional border points
    xs = [box[0][0], box[1][0]]
    ys = [box[0][1], box[1][1]]
    box = [[min(xs) - border_s, max(ys) + border_s], [max(xs) + border_s, min(ys) - border_s]]

    # Build list of gridpoints that are outside all structures in cell
    lxr = int((box[1][0] - box[0][0]) / grid_s) + 1
    lyr = int((box[0][1] - box[1][1]) / grid_s) + 1

    if type(xr) not in [list, np.ndarray]:
        xr = np.linspace(box[0][0], box[1][0], lxr)
    else:
        lxr = len(xr)

    if type(yr) not in [list, np.ndarray]:
        yr = np.linspace(box[0][1], box[1][1], lyr)
    else:
        lyr = len(yr)

    p = []
    p_d_to = []
    p_d_fr = []
    for y in yr:
        for x in xr:

            c = (x, y)
            p.append(c)

            # Compute squared Euclidean distance from to and fr coords for each gridpoint
            # For optimization we don't need to sqrt() since minimal in squared is minimal in sqrt
            dist_fr = (x - fr[0])**2 + (y - fr[1])**2
            dist_to = (x - to[0])**2 + (y - to[1])**2

            p_d_fr.append(dist_fr)
            p_d_to.append(dist_to)

    p_i = np.array(p)
    p_d_fr = np.array(p_d_fr)
    p_d_to = np.array(p_d_to)

    # Build list of points that are inside a structure
    cell_ref = gp.CellReference(cell)
    if detect == 'native':
        inside = np.array(gp.inside(p, cell_ref, precision = precision))
    elif detect == 'custom':
        inside = np.array(gdst.funcs.inside(p, cell_ref, dist = dist_multi*grid_s, nop = nop, precision = precision))
    else:
        raise ValueError('Parameter \'detect\' is only allowed to have values [\'native\', \'custom\'], cannot continue')

    p_d_fr_min = np.min(p_d_fr[np.argwhere(inside == False)])
    p_d_to_min = np.min(p_d_to[np.argwhere(inside == False)])

    # Get p_i index of starting values
    start_i = np.argwhere(p_d_fr == p_d_fr_min).tolist()
    end_i = np.argwhere(p_d_to == p_d_to_min).tolist()

    start_i = [item for sublist in start_i for item in sublist]
    end_i = [item for sublist in end_i for item in sublist]

    start = p_i[start_i]
    end = p_i[end_i]

    # Now start stepping from start to end, labelling all gridpoints accordingly by the number of steps required from starting point to reach it
    n = [0] * 4
    lp = len(p)
    p_g = [0] * lp
    path_found = False
    k = 0

    while not path_found and start_i:

        k += 1
        next_start_i = []

        if debug:
            print(start_i)

        for i in start_i:

            # Look up, down, left and right, store the index
            n = lookaround(i, lxr, pref = pref)

            # Check if any of the neighbouring points are not in a structure and not in p_g
            for nb in n:

                if nb in end_i:
                    path_found = True
                    p_g[nb] = k
                    final_index = nb

                    if debug:
                        # Visualize
                        circ = gp.Round(p[nb], 0.1, layer = 10)
                        cell.add(circ)
                        txt = gp.Text(str(k), 0.5, p[nb], layer = 11)
                        cell.add(txt)

                    break

		         # Point is out of bounds, marked as structure (< 0) or already has a step value (> 0)
                if nb < 0 or nb >= lp or p_g[nb] != 0 or (i % lxr == 0 and nb % lxr == 1) or (i % lxr == 1 and nb % lxr == 0):
                    continue # Skip this iteration

                if inside[nb]:
                    p_g[nb] = -1
                else:
                    p_g[nb] = k
                    next_start_i.append(nb)

                    if debug:
                        # Visualize
                        circ = gp.Round(p[nb], 0.1, layer = 1)
                        cell.add(circ)
                        txt = gp.Text(str(k), 0.5, p[nb], layer = 2)
                        cell.add(txt)

        start_i = copy.copy(next_start_i)

    # Routing ended, checking whether we succeeded
    if not path_found:
        print('>> ERROR: No existing route was found.')
        return False

    print('>> Found a route in ' + str(k) + ' steps.')

    # Backtrace path
    this_index = final_index
    backtraced = [to, p[final_index]]
    switched = False
    for this_k in range(k, -1, -1):

        # Change move preference after switch_pref moves
        if switch_pref and not switched and this_k < switch_pref*k:
            pref = 'x' if pref == 'y' else 'y'
            switched = True

        n = lookaround(this_index, lxr, pref = pref)
        for nb in n:
            if nb < lp and p_g[nb] == this_k:
                this_index = nb
                backtraced.append(p[nb])
                break

    backtraced.append(fr)

    if debug:
        print('>> Points of found route:')
        print(backtraced)

    # Generate list of widths for the route
    if not uniform_width:
        to_w = to_str.endpoint_dims[to_ep]
        fr_w = fr_str.endpoint_dims[fr_ep]
        ws = [to_w if to_w != None else width]*2 + [width]*(len(backtraced)-4) + [fr_w if fr_w != None else width]*2
    else:
        ws = width

    # Create backtraced path
    if pathmethod == 'poly':
        r = gp.PolyPath(backtraced, ws, layer = layer)
    elif pathmethod == 'flex':
        r = gp.FlexPath(backtraced, ws, corners = 'smooth', layer = layer)
    else:
        raise ValueError('Parameter \'pathmethod\' only has allowed values [\'poly\', \'flex\']')

    ends = {'A': backtraced[0], 'B': backtraced[-1]}
    epsz = {'A': ws[0] if not uniform_width else width, 'B': ws[-1] if not uniform_width else width}

    structure = gdst.classes.GDStructure(r, ends, epsz)

    # Connect to 'to' and 'from' structures
    fr_str.next['AUTOROUTE_A'] = structure
    to_str.next['AUTOROUTE_B'] = structure
    structure.prev['AUTOROUTE_A'] = fr_str
    structure.prev['AUTOROUTE_B'] = to_str

    return structure

class GdsMap:
    def __init__(self, cell, grid_size, buffer=None, layers=None, precision=0.001):

        self.grid_size = grid_size

        start_time = time.time()
        print('--- Constructing map ---')

        if layers is None:
            layers = list(cell.get_layers())
        if not isinstance(layers, list):
            layers = [layers]
        
        if buffer is None:
            buffer = grid_size

        bounding_box = cell.get_bounding_box()
        x_length = int((bounding_box[1][0] - bounding_box[0][0] + 2 * grid_size) / grid_size) + 1
        y_length = int((bounding_box[1][1] - bounding_box[0][1] + 2 * grid_size) / grid_size) + 1
        x_linspace = np.linspace(bounding_box[0][0] - grid_size,
                                 bounding_box[1][0] + grid_size,
                                 x_length)
        y_linspace = np.flip(np.linspace(bounding_box[0][1] - grid_size,
                                         bounding_box[1][1] + grid_size,
                                         y_length))
        x_array, y_array = np.meshgrid(x_linspace, y_linspace)
        xy_array = np.array(list(zip(x_array.ravel(),
                                     y_array.ravel()))).reshape(*x_array.shape, 2)
        self.mask = np.zeros((y_length, x_length))

        cell_copy = cell.copy(cell.name + "_copy",
                              deep_copy=True)
        cell_copy.remove_polygons(lambda pts, layer, datatype: layer not in layers)
        map_gdspy = np.array(gp.inside(np.array(xy_array).reshape(-1, 2),
                                       gp.CellReference(cell_copy),
                                       precision=precision)).reshape((y_length, x_length))

        neighbours, grid_distance = get_neighbours(self, buffer)

        for i in range(y_length):
            for j in range(x_length):
                if not map_gdspy[i][j]:
                    if self.mask[i][j] != 1:
                        self.mask[i][j] = 0
                else:
                    if grid_distance == 1:
                        self.mask[i][j] = self.mask[i - 1][j] = self.mask[i + 1][j] = self.mask[i][j - 1] = self.mask[i][j + 1] = 1
                    else:
                        for neighbour in neighbours:
                            try:
                                self.mask[i + neighbour[0]][j + neighbour[1]] = 1
                            except:
                                pass

        elapsed_time = round(time.time() - start_time, 3)
        print('--- Constructed map in ' + str(elapsed_time) + 's ---')

        self.x_length = x_length
        self.y_length = y_length
        self.x_linspace = x_linspace
        self.y_linspace = y_linspace
        self.x_array = x_array
        self.y_array = y_array
        self.xy_array = xy_array

def heuristic(point_a, point_b, grid_size = 1, sqrt_two=1.4):
    delta_x = abs(point_a[0] - point_b[0])
    delta_y = abs(point_a[1] - point_b[1])
    edge_dist = grid_size
    diag_dist = sqrt_two * grid_size
    distance = edge_dist * max(delta_x, delta_y) + (diag_dist - edge_dist) * min(delta_x, delta_y)
    return distance

def array_heuristic(x_array, y_array, point, grid_size, sqrt_two=1.4):
    delta_x = np.absolute(x_array - point[0])
    delta_y = np.absolute(y_array - point[1])
    edge_dist = grid_size
    diag_dist = sqrt_two * grid_size
    distance = edge_dist * np.maximum(delta_x, delta_y) + (diag_dist - edge_dist) * np.minimum(delta_x, delta_y)
    return distance

def bipoint_angle(point_a, point_b):
    delta_x = point_b[0] - point_a[0]
    delta_y = point_b[1] - point_a[1]
    theta = np.arctan2(delta_y, delta_x)
    return theta

def blocked(point, delta_x, delta_y, mask):
    if point[0] + delta_x < 0 or point[0] + delta_x >= mask.shape[1]:
        return True
    if point[1] + delta_y < 0 or point[1] + delta_y >= mask.shape[0]:
        return True
    if delta_x != 0 and delta_y != 0:
        if mask[point[1]][point[0] + delta_x] == 1 and mask[point[1] + delta_y][point[0]] == 1:
            return True
        if mask[point[1] + delta_y][point[0] + delta_x] == 1:
            return True
    else:
        if delta_x != 0:
            if mask[point[1]][point[0] + delta_x] == 1:
                return True
        else:
            if mask[point[1] + delta_y][point[0]] == 1:
                return True
    return False

def jump_point_astar(mask, start_point, end_point, grid_size=100, sqrt_two=1.4):
    closed_set = set()
    came_from = {}
    g_score = {start_point: 0}
    f_score = {start_point: heuristic(start_point, end_point, grid_size, sqrt_two)}

    point_queue = []

    heapq.heappush(point_queue, (f_score[start_point], start_point))

    while point_queue:

        current_pos = heapq.heappop(point_queue)[1]
        if current_pos == end_point:
            route = []
            while current_pos in came_from:
                route.append(current_pos)
                current_pos = came_from[current_pos]
            route.append(start_point)
            route = list(reversed(route[::]))
            return route

        closed_set.add(current_pos)
        for delta_x, delta_y in [(0, 1),
                                 (0, -1),
                                 (1, 0),
                                 (-1, 0),
                                 (1, 1),
                                 (1, -1),
                                 (-1, 1),
                                 (-1, -1)]:

            if blocked(current_pos, delta_x, delta_y, mask):
                continue

            neighbour = current_pos[0] + delta_x, current_pos[1] + delta_y

            if delta_x != 0 and delta_y != 0:
                tentative_g_score = g_score[current_pos] + sqrt_two * grid_size
            else:
                tentative_g_score = g_score[current_pos] + grid_size

            if neighbour in closed_set:
                continue

            if tentative_g_score < g_score.get(
                    neighbour, 0) or neighbour not in [i[1] for i in point_queue]:
                came_from[neighbour] = current_pos
                g_score[neighbour] = tentative_g_score
                f_score[neighbour] = tentative_g_score + heuristic(
                    neighbour, end_point, grid_size, sqrt_two)
                heapq.heappush(point_queue, (f_score[neighbour], neighbour))
    return False

def get_neighbours(GdsMap, distance):
    grid_distance = int((distance - distance % GdsMap.grid_size) / (2 * GdsMap.grid_size) + 1)
    if grid_distance == 1:
        neighbours = [(0, 0), (0, 1), (0, -1), (1, 0), (-1, 0)]
    else:
        neighbours = [(i, j)
                        for i in range(-grid_distance, grid_distance + 1)
                        for j in range(-grid_distance, grid_distance + 1)]
    return neighbours, grid_distance

def pathfinder(GdsMap,
               end_points,
               directions = None,
               path_width = None,
               segment_distances = 0,
               reference_point = None,
               path_buffer_distance = 0):

    if not isinstance(path_width, list):
        path_width = [path_width]
    path_width = np.array(path_width)

    if not isinstance(path_buffer_distance, list):
        path_buffer_distance = [path_buffer_distance] * len(path_width)
    path_buffer_distance = np.array(path_buffer_distance)

    if not isinstance(segment_distances, list):
        segment_distances = [segment_distances]

    if directions != None:
        from_point = (end_points[0][0] + path_width[0] * np.cos(directions[0]), end_points[0][1] + path_width[0] * np.sin(directions[0]))
        to_point = (end_points[1][0] + path_width[-1] * np.cos(directions[1]), end_points[1][1] + path_width[-1] * np.sin(directions[1]))
    else:
        from_point = end_points[0]
        to_point = end_points[1]

    reference_point = reference_point if reference_point is not None else to_point

    start_point_index = ma.MaskedArray.argmin(
        ma.array(
            array_heuristic(
                GdsMap.x_array,
                GdsMap.y_array, 
                from_point,
                GdsMap.grid_size),
            mask=GdsMap.mask))

    start_point = (start_point_index % GdsMap.x_length,
                   start_point_index // GdsMap.x_length)
    
    end_point_index = ma.MaskedArray.argmin(
        ma.array(
            array_heuristic(
                GdsMap.x_array,
                GdsMap.y_array,
                to_point,
                GdsMap.grid_size),
            mask=GdsMap.mask))
    
    end_point = (end_point_index % GdsMap.x_length,
                 end_point_index // GdsMap.x_length)

    start_time = time.time()
    route = jump_point_astar(GdsMap.mask, start_point, end_point)
    elapsed_time = round(time.time() - start_time, 3)

    if not route:
        print('--- Path not found after ' + str(elapsed_time) + 's ---')
        return False
    print('--- Found route in ' + str(elapsed_time) + 's ---')

    routed_path = []
    routed_path.append(end_points[0])
    for point in route:
        routed_path.append(GdsMap.xy_array[point[1]][point[0]])
    routed_path.append(end_points[1])

    if directions != None:
        routed_path[1] = from_point
        routed_path[-2] = to_point

    index = 0
    neighbours, _ = get_neighbours(GdsMap, path_width[index] + 2 * path_buffer_distance[index])

    # Remove trivial points
    remove_index_list = []
    new_path_point_list = []
    path_length = len(routed_path)
    for i in range(1, path_length - 2):
        for neighbour in neighbours:
            GdsMap.mask[route[i - 1][1] + neighbour[1]][route[i - 1][0] + neighbour[0]] = 1
        if bipoint_angle(routed_path[i - 1], routed_path[i]) == bipoint_angle(routed_path[i], routed_path[i + 1]):
            if heuristic(routed_path[i], reference_point) < segment_distances[index]:
                new_path_point_list.append(routed_path[i])
                index += 1
                neighbours, _ = get_neighbours(GdsMap, path_width[index] + 2 * path_buffer_distance[index])
            else:
                remove_index_list.append(i)

    for index in reversed(remove_index_list):
        del routed_path[index]
    path_length = len(routed_path)
    print('--- Removed ' + str(len(remove_index_list)) +
          ' trivial points from route ---')

    route_segments = []
    for segment_point in new_path_point_list:
        for i, route_point in enumerate(routed_path):
            if np.array_equal(route_point, segment_point):
                route_segments.append(routed_path[:i + 1])
                routed_path = routed_path[i:]
                break
    route_segments.append(routed_path)

    return route_segments
