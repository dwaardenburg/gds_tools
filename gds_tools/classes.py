import copy

import numpy as np
import gdspy as gp
import gds_tools as gdst
import design_elements as elmt
from operator import add, sub
from recordclass import recordclass

class GDStructure:

    def __init__(self, structure, endpoints = {}, endpoint_dims = {}, endpoint_directions = None, method = None, args = {}):
        """Define class to store structures in

        Args:
            structure (gp.Path()): gdspy description of object (i.e. ouput of gp.PolyPath()).
            endpoints (dict): positions of the endpoints to which you connect other structures.
            endpoint_dims (dict): dimensions of the endpoints (width/height/etc.), to be used for connecting e.g. transmission lines.
            endpoint_directions ([type], optional): [description]. Defaults to None.
            method ([type], optional): [description]. Defaults to None.
            args (dict, optional): [description]. Defaults to {}.
        """

        self.type = 'GDStructure'
        self.structure = structure
        self.endpoints = endpoints
        self.endpoint_dims = endpoint_dims
        self.endpoint_directions = endpoint_directions
        self.compound = []
        self.prev = {}
        self.next = {}
        self.method = method

        self.__args_tuple__ = recordclass('args', args.keys())
        self.args = self.__args_tuple__(**args)

        self.gen = self.generate

    def generate(self):
        if self.method != None:
            return self.method(**self.args.__dict__)
        else:
            raise TypeError('This GDStructure does not support the .generate() method (yet)')

    def copy(self):
        """Makes a new copy of the structure,
        i.e. allocates new memory and copies the data contents to it.

        Returns:
        gds_tools.GDStructure: pointer of copy of original
        """

        # Create deepcopy of self and return the copy
        # Do not copy.deepcopy(self) directly, as this will also copy the connections and screw things up
        # Instead, just fill a new instance of the class with copies of relevant internal structure
        new_obj = gdst.classes.GDStructure(copy.deepcopy(self.structure), copy.deepcopy(self.endpoints), copy.deepcopy(self.endpoint_dims), copy.deepcopy(self.endpoint_directions))

        compound = []
        for i, c in enumerate(self.compound):
            c_copy = c.copy()
            compound.append(c_copy)
            new_obj.next[i] = c_copy

        new_obj.compound = compound

        return new_obj

    #==================
    # Rotate structure \\
    #=========================================================================
    # When rotating, all endpoints need to move wih the rotation,           ||
    # for this reason we need our own rotation function.                    ||
    #                                                                       ||
    # Arguments:    rad         :   amount of radians to rotate             ||
    #               signal_from :   tells linked structure who sent command ||
    #=========================================================================
    def rotate(self, rad, signal_from = None):

        if not isinstance(signal_from, list):
            signal_from = [signal_from]

        # Rotate the structure with standard gdspy .rotation() function
        if self not in signal_from:
            signal_from.append(self)
        self.structure.rotate(rad)

        # Rotate the endpoints
        for i in self.endpoints:
            self.endpoints[i] = tuple(gdst.functions.VecRot(rad, self.endpoints[i]))
            if self.endpoint_directions:
                self.endpoint_directions[i] += rad

        # Rotate all connected structures
        for i in self.prev:
            if self.prev[i] not in signal_from and self.prev:
                signal_from.append(self.prev[i])
                self.prev[i].rotate(rad, signal_from = signal_from)

        for i in self.next:
            if self.next[i] not in signal_from and self.next:
                signal_from.append(self.next[i])
                self.next[i].rotate(rad, signal_from = signal_from)

        return self

    #=====================
    # Get structure layer \\
    #=========================================================================
    # Gdspy has some inconsistencies in how it stores layers, this getter   ||
    # function should make life easier to obtain the layer of the object    ||
    # stored in self.structure                                              ||
    #=========================================================================
    def getlayer(self):
        return self.structure.layers[0] if hasattr(self.structure, 'layers') else (self.structure.layer if hasattr(self.structure, 'layer') else 0)

    #=====================
    # Translate structure \\
    #=========================================================================
    # Arguments:    delta       :   vector (x, y) how much you want to move ||
    #               signal_from :   tells linked structure who sent command ||
    #=========================================================================
    def translate(self, delta, signal_from = None):

        if not isinstance(signal_from, list):
            signal_from = [signal_from]

        # Translate structure with standard gdspy .translate() function
        if self not in signal_from:
            signal_from.append(self)
        self.structure.translate(delta[0], delta[1])

        # Translate all the endpoints
        for i in self.endpoints:
            self.endpoints[i] = tuple(map(add, self.endpoints[i], delta))

        # Translate all connected structures
        for i in self.prev:
            if self.prev[i] not in signal_from and self.prev:
                signal_from.append(self.prev[i])
                self.prev[i].translate(delta, signal_from = signal_from)

        for i in self.next:
            if self.next[i] not in signal_from and self.next:
                signal_from.append(self.next[i])
                self.next[i].translate(delta, signal_from = signal_from)

        return self

    #=====================
    # Mirror structure \\
    #=========================================================================
    # Arguments:    p1, p2  :   p1 and p2 are (x, y) coords forming         ||
    #                           the mirror line                             ||
    #=========================================================================
    def mirror(self, p1, p2, signal_from = None):

        if not isinstance(signal_from, list):
            signal_from = [signal_from]

        # Mirror the gdspy shape
        if self not in signal_from:
            signal_from.append(self)
        self.structure.mirror(p1, p2 = p2)

        # Process the endpoints
        for k, v in self.endpoints.items():

            p1 = list(p1)
            p2 = list(p2)

            if p1[0] == p2[0]:
                p1[0] += 1E-20

            if p1[1] == p2[1]:
                p1[1] += 1E-20

            # y = ax + c : mirror line
            a = (p2[1] - p1[1]) / (p2[0] - p1[0])
            c = p1[1] - a * p1[0]
            d = (v[0] + (v[1] - c)*a) / (1 + a**2)

            v2x = 2*d - v[0]
            v2y = 2*d*a - v[1] + 2*c

            self.endpoints[k] = (v2x, v2y)

        # Ripple through all connected shapes
        for i in self.prev:
            if self.prev[i] not in signal_from and self.prev:
                signal_from.append(self.prev[i])
                self.prev[i].mirror(p1, p2, signal_from = signal_from)

        for i in self.next:
            if self.next[i] not in signal_from and self.next:
                signal_from.append(self.next[i])
                self.next[i].mirror(p1, p2, signal_from = signal_from)

        return self

    #=================
    # Heal connection \\
    #=========================================================================
    # When connection endpoints, there is always a little space that is not ||
    # overlapped. This function will fill this overlap.                     ||
    #                                                                       ||
    # Arguments:    endpoint    :   (x, y) center of healer                 ||
    #    (optional) npoints     :   number of points to use for boundary    ||
    #=========================================================================
    def heal(self, endpoint, npoints = 100, r = 'auto', layer = None, datatype = None):

        if type(r) == str and r == 'auto':
            r = self.endpoint_dims[endpoint] / 2

        healer = gdst.heal.circle(r, self.endpoints[endpoint], npoints = npoints, layer = layer if layer != None else self.getlayer(), datatype = datatype if datatype != None else self.structure.datatypes[0])

        self.prev['HEAL_' + endpoint] = healer
        self.compound += [healer]

        return self

    #===================
    # Connect structure \\
    #=========================================================================
    # Arguments:    ep_self :   reference to A, B, C, etc.                  ||
    #                           name of endpoint of this structure          ||
    #               to      :   GDStructure class input to which to connect ||
    #               ep_to   :   endpoint to which to connect to             ||
    #               offset  :   (x, y) offset from target endpoints         ||
    #               obj_link:   use Python object reference for linked list ||
    #                           entry rather than endpoint key              ||
    #=========================================================================
    def connect(self, ep_self, to, ep_to, offset = (0, 0), obj_link = False, rotate = True):

        # Rotate to correct orientation
        if self.endpoint_directions != None and to.endpoint_directions != None and self.endpoint_directions and to.endpoint_directions and rotate:
            self.rotate((to.endpoint_directions[ep_to] - self.endpoint_directions[ep_self] - np.pi) % (2 * np.pi))

        # Move to connect endpoints
        delta = tuple(map(sub, to.endpoints[ep_to], self.endpoints[ep_self]))
        delta = tuple(map(add, delta, offset))
        self.mov(delta)

        # Update linked list
        if not obj_link:
            to.next[ep_to] = self
            self.prev[ep_self] = to
        else:
            to.next[self] = self
            self.prev[to] = to

        return self

    #======================
    # Disconnect structure \\
    #=========================================================================
    # Only removes the references in the linked lists                       ||
    #=========================================================================
    def disconnect(self):

        for e in self.endpoints:
            if e in self.next:
                for ne in self.next[e].endpoints:
                    if ne in self.next[e].prev and self.next[e].prev[ne] == self:
                        del self.next[e].prev[ne]
                del self.next[e]
            if e in self.prev:
                for pe in self.prev[e].endpoints:
                    if pe in self.prev[e].next and self.prev[e].next[pe] == self:
                        del self.prev[e].next[pe]
                del self.prev[e]

        return self

    # Aliases
    rot = rotate
    mov = translate
    con = connect
    dis = disconnect

class Waveguide:

    def __init__(self,
                 path_list,
                 path_types,
                 layer_specs,
                 hole_radius = None,
                 hole_distance = None,
                 path_holey_ground = None,
                 hole_references = None,
                 length = 0,
                 endpoints = None,
                 endpoint_directions = None):

        """Create a waveguide object that allows an array of path elements to be treated as an single path.

        Args:
            path_list (list): the gp.path classes to be bundled.
            layer_list (list): the layers of each of paths.
            holey_ground_path ([type], optional): [description]. Defaults to None.
            hole_references ([type], optional): [description]. Defaults to None.
            hole_distance (float, optional): distance between hole centres. Defaults to None.
            hole_radius (float, optional): radius of the holes. Defaults to None.
            length (float, optional): length of the path. Defaults to 0.
            endpoints (dict, optional): dictionairy defining the endpoints (x, y). Defaults to None.
            endpoint_directions (dict, optional): dictionairy defining the endpoint directions. Defaults to None.
        """

        self.type = 'Waveguide'
        self.paths = path_list
        self.path_types = path_types
        self.layer_specs = layer_specs
        self.structure = path_list[self.path_types['stripline'][0]]
        self.compound = [gdst.classes.GDStructure(path, {}, {}) for path in self.paths]
        self.holey_unit_cells = []
        self.path_holey_ground = [] if path_holey_ground == None else path_holey_ground
        self.hole_references = [] if hole_references == None else hole_references
        self.compound.append(gdst.classes.GDStructure(self.path_holey_ground, {}, {}))
        self.compound.append(gdst.classes.GDStructure(self.hole_references, {}, {}))

        if 'holey_unit_cell_index' in globals():
            global holey_unit_cell_index
        else:
            global holey_unit_cell_index
            holey_unit_cell_index = 0

        if endpoints == None:
            self.endpoints = {}
            self.endpoints['A'] = (self.structure.x, self.structure.y)
            self.endpoints['B'] = (self.structure.x, self.structure.y)
        else:
            self.endpoints = endpoints

        self.endpoint_directions = {}
        if endpoint_directions == None:
            if isinstance(self.structure.direction, str):
                if self.structure.direction == '+x':
                    self.endpoint_directions['A'] = np.pi
                elif self.structure.direction == '-x':
                    self.endpoint_directions['A'] = 0
                elif self.structure.direction == '+y':
                    self.endpoint_directions['A'] = - np.pi / 2
                elif self.structure.direction == '-y':
                    self.endpoint_directions['A'] = np.pi / 2
            else:
                self.endpoint_directions['A'] = self.structure.direction
            self.endpoint_directions['B'] = self.endpoint_directions['A'] + np.pi
        else:
            self.endpoint_directions = endpoint_directions

        self.endpoint_widths = [0] * len(self.paths)
        self.endpoint_distances = [0] * len(self.paths)
        self.endpoint_path_numbers = [0] * len(self.paths)
        self.holey_ground_paths = []

        self.requested_hole_distance = hole_distance
        self.hole_radius = hole_radius

        for i, path in enumerate(self.paths):
            self.endpoint_widths[i] = 2 * path.w
            self.endpoint_path_numbers[i] = path.n
            self.endpoint_distances[i] = path.distance

        if self.requested_hole_distance != None and self.hole_radius != None:
            self.holey_ground_paths = self.path_types['holey_grounds']
            holey_path_index = self.path_types['holey_grounds'][0]
            self.holey_ground_layer = self.layer_specs[holey_path_index]['layer']
            self.holey_ground_datatype = self.layer_specs[holey_path_index]['datatype']
            self.holey_ground_spec = self.layer_specs[holey_path_index]

            self.hole_distances = [0] * len(self.paths)
            self.round_verticies = 6
            hole_list = []
            for i in self.holey_ground_paths:
                effective_width = self.endpoint_widths[i] - 2 * self.hole_radius
                total_vertical_hole_distances = int((effective_width - effective_width % self.requested_hole_distance) / self.requested_hole_distance)
                self.hole_distances[i] = effective_width / total_vertical_hole_distances
                if self.endpoint_path_numbers[i] == 1:
                    for j in range(total_vertical_hole_distances + 1):
                        hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_width / 2), self.hole_radius, number_of_points = self.round_verticies))
                    for j in range(total_vertical_hole_distances):
                        hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_width / 2), self.hole_radius, number_of_points = self.round_verticies))
                if self.endpoint_path_numbers[i] == 2:
                    for j in range(total_vertical_hole_distances + 1):
                        hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_width / 2 - self.endpoint_distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                        hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_width / 2 + self.endpoint_distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                    for j in range(total_vertical_hole_distances):
                        hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_width / 2 - self.endpoint_distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                        hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_width / 2 + self.endpoint_distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
            self.holey_unit_cell_width = (np.sqrt(3) / 2) * self.hole_distances[self.holey_ground_paths[0]]
            self.holey_unit_cell_polygon = gp.boolean(hole_list, [], 'or', **self.holey_ground_spec).translate(- self.holey_unit_cell_width / 2, 0)

            self.active_unit_cell = 0
            try:
                self.holey_unit_cells.append(gp.Cell('Holey_unit_cell_' + str(holey_unit_cell_index) + '_' + str(self.active_unit_cell)))
                holey_unit_cell_index += 1
            except:
                holey_unit_cell_index += 1
                self.holey_unit_cells.append(gp.Cell('Holey_unit_cell_' + str(holey_unit_cell_index) + '_' + str(self.active_unit_cell)))

            self.holey_unit_cells[self.active_unit_cell].add(gp.boolean(hole_list, [], 'or', **self.holey_ground_spec).translate(- self.holey_unit_cell_width / 2, 0))
        else:
            self.hole_distances = self.requested_hole_distance
            self.hole_radius = self.hole_radius

        self.length = length

        self.prev = {}
        self.next = {}

    def construct_holey_unit_cell(self, widths, distances):
        self.hole_distances = [0] * len(self.paths)
        hole_list = []
        for i in self.holey_ground_paths:
            effective_path_width = widths[i] - 2 * self.hole_radius
            total_vertical_hole_distances = int((effective_path_width - effective_path_width % self.requested_hole_distance) / self.requested_hole_distance)
            self.hole_distances[i] = effective_path_width / total_vertical_hole_distances
            if self.endpoint_path_numbers[i] == 1:
                for j in range(total_vertical_hole_distances + 1):
                    hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_path_width / 2), self.hole_radius, number_of_points = self.round_verticies))
                for j in range(total_vertical_hole_distances):
                    hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_path_width / 2), self.hole_radius, number_of_points = self.round_verticies))
            if self.endpoint_path_numbers[i] == 2:
                for j in range(total_vertical_hole_distances + 1):
                    hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_path_width / 2 - distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                    hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_path_width / 2 + distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                for j in range(total_vertical_hole_distances):
                    hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_path_width / 2 - distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                    hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_path_width / 2 + distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
        self.holey_unit_cell_width = (np.sqrt(3) / 2) * self.hole_distances[self.holey_ground_paths[0]]

        global holey_unit_cell_index
        self.active_unit_cell += 1
        self.holey_unit_cells.append(gp.Cell('Holey_unit_cell_' + str(holey_unit_cell_index) + '_' + str(self.active_unit_cell)))
        self.holey_unit_cells[self.active_unit_cell].add(gp.boolean(hole_list, [], 'or', **self.holey_ground_spec).translate(- self.holey_unit_cell_width / 2, 0))

    def segment(self, length, direction = None, final_widths = None, final_distances = None, axis_offset = 0):

        if isinstance(self.structure.direction, str) or isinstance(direction, str):
            if self.structure.direction == '+x':
                self.endpoint_directions['B'] = 0
            elif self.structure.direction == '-x':
                self.endpoint_directions['B'] = np.pi
            elif self.structure.direction == '+y':
                self.endpoint_directions['B'] = np.pi / 2
            elif self.structure.direction == '-y':
                self.endpoint_directions['B'] = - np.pi / 2
        else:
            self.endpoint_directions['B'] = self.structure.direction

        if final_widths == None and final_distances == None:
            if self.hole_distances != None:
                effective_length = length - self.holey_unit_cell_width
                total_unit_cells = int((effective_length - effective_length % 2 * self.holey_unit_cell_width) / (2 * self.holey_unit_cell_width))
                unit_cell_distance = effective_length / total_unit_cells

                for i in range(total_unit_cells):
                    referenced_unit_cell = gp.CellReference(self.holey_unit_cells[self.active_unit_cell], origin = (self.structure.x, self.structure.y), rotation = np.degrees(self.endpoint_directions['B']))

                    dx = np.cos(self.endpoint_directions['B']) * (self.holey_unit_cell_width + unit_cell_distance) / 2
                    dy = np.sin(self.endpoint_directions['B']) * (self.holey_unit_cell_width + unit_cell_distance) / 2
                    referenced_unit_cell.translate(dx, dy)

                    dx = np.cos(self.endpoint_directions['B']) * i * unit_cell_distance
                    dy = np.sin(self.endpoint_directions['B']) * i * unit_cell_distance

                    referenced_unit_cell.translate(dx, dy)
                    self.hole_references.append(referenced_unit_cell)

            for i, path in enumerate(self.paths):
                if i not in self.holey_ground_paths:
                    path.segment(length, direction, None, None, axis_offset, **self.layer_specs[i])
                else:
                    path.x = self.structure.x
                    path.y = self.structure.y
                    path.direction = self.structure.direction
        else:
            holey_ground = []
            holey_index = 0
            for i, path in enumerate(self.paths):
                if i not in self.holey_ground_paths:
                    path.segment(length, direction, final_widths[i], final_distances[i], axis_offset, **self.layer_specs[i])
                    path.w = final_widths[i] / 2
                    path.distance = final_distances[i]
                else:
                    if self.requested_hole_distance != None and self.hole_radius != None:
                        holey_ground.append(gp.Path(2 * path.w, initial_point=(path.x, path.y), number_of_paths = path.n, distance = path.distance))
                        holey_ground[holey_index].direction = self.structure.direction
                        holey_ground[holey_index].segment(length, direction, final_widths[i], final_distances[i], axis_offset)
                        holey_index += 1
                    path.x = self.structure.x
                    path.y = self.structure.y
                    path.direction = self.structure.direction

            self.endpoint_widths = final_widths
            self.endpoint_distances = final_distances
            if self.requested_hole_distance != None and self.hole_radius != None:
                self.path_holey_ground.append(elmt.lattices.generate_holey_ground(gp.boolean(holey_ground, [], 'or', **self.holey_ground_spec),
                                                                                     self.hole_radius,
                                                                                     self.requested_hole_distance,
                                                                                     self.round_verticies,
                                                                                     **self.holey_ground_spec))
                self.construct_holey_unit_cell(widths = final_widths, distances = final_distances)
        self.length += length
        if length < 0:
            print("Warning: Negative segment length detected, waveguide lengths will no longer be correct.")

        self.endpoints['B'] = (self.structure.x, self.structure.y)

    def klopfenstein_taper(self, length, final_widths, final_distances, point_distance):
        holey_ground = []
        holey_index = 0
        for i, path in enumerate(self.paths):
            if i not in self.holey_ground_paths:
                direction = path.direction
                d = {'+x': 0, '-x': np.pi, '+y': np.pi / 2, '-y': np.pi / 2}
                direction = direction if isinstance(direction, float) else (d[direction])
                x, y = path.x, path.y
                path.rotate(- direction, center=(x, y))
                if self.endpoint_distances[i] == 0:
                    width = lambda x: self.endpoint_widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.endpoint_widths[i])
                    path.parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                else:
                    distance = lambda x: self.endpoint_distances[i] + gdst.functions.biquadratic_func(x) * (final_distances[i] - self.endpoint_distances[i])
                    width = lambda x: self.endpoint_widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.endpoint_widths[i])
                    path.parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, final_distance = distance, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                path.rotate(direction, center=(x, y))
            else:
                if self.requested_hole_distance != None and self.hole_radius != None:
                    holey_ground.append(gp.Path(2 * path.w, initial_point=(path.x, path.y), number_of_paths = path.n, distance = path.distance))
                    holey_ground[holey_index].direction = self.structure.direction
                    if self.endpoint_distances[i] == 0:
                        width = lambda x: self.endpoint_widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.endpoint_widths[i])
                        holey_ground[holey_index].parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                    else:
                        distance = lambda x: self.endpoint_distances[i] + gdst.functions.biquadratic_func(x) * (final_distances[i] - self.endpoint_distances[i])
                        width = lambda x: self.endpoint_widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.endpoint_widths[i])
                        holey_ground[holey_index].parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, final_distance = distance, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                    holey_index += 1
                path.x = self.structure.x
                path.y = self.structure.y
                path.direction = self.structure.direction

        if self.requested_hole_distance != None and self.hole_radius != None:
                self.path_holey_ground.append(elmt.lattices.generate_holey_ground(gp.boolean(holey_ground, [], 'or', **self.holey_ground_spec),
                                                                                     self.hole_radius,
                                                                                     self.requested_hole_distance,
                                                                                     self.round_verticies,
                                                                                     **self.holey_ground_spec))
                self.construct_holey_unit_cell(widths = final_widths, distances = final_distances)

        self.endpoints['B'] = (self.structure.x, self.structure.y)

        for i, path in enumerate(self.paths):
            path.w = final_widths[i] / 2
            path.distance = final_distances[i]
        self.endpoint_widths = final_widths
        self.endpoint_distances = final_distances
        # if self.requested_hole_distance != None and self.hole_radius != None:
        #     self.construct_holey_unit_cell(widths = final_widths, distances = final_distances)
        self.length += length
        if length < 0:
            print("Warning: Negative tapering length length detected, waveguide lengths will no longer be correct.")

    def turn(self, radius, angle, tolerance = 0.01, number_of_points = None, max_points = 199, final_width = None, final_distance = None):

        angles = {'r': - np.pi / 2, 'l': np.pi / 2, 'rr': - np.pi, 'll': np.pi}
        angle = angle if not isinstance(angle, str) else angles[angle]

        turn_length = radius * abs(angle)
        self.length += turn_length
        if turn_length < 0:
            print("Warning: Negative turn length detected, waveguide lengths will no longer be correct.")

        if self.hole_distances != None:
            effective_length = turn_length - self.holey_unit_cell_width
            total_unit_cells = int((effective_length - effective_length % 2 * self.holey_unit_cell_width) / (2 * self.holey_unit_cell_width))
            unit_cell_distance = effective_length / total_unit_cells

            for i in range(total_unit_cells):
                referenced_unit_cell = referenced_unit_cell = gp.CellReference(self.holey_unit_cells[self.active_unit_cell], origin = (self.structure.x, self.structure.y), rotation = np.degrees(self.endpoint_directions['B']))

                relative_radius = radius * np.sign(angle)
                dx_0 = np.cos(self.endpoint_directions['B'] + np.pi / 2) * relative_radius
                dy_0 = np.sin(self.endpoint_directions['B'] + np.pi / 2) * relative_radius
                turn_center = (self.structure.x + dx_0, self.structure.y + dy_0)
                dtheta = (self.holey_unit_cell_width + unit_cell_distance) / (2 * relative_radius)
                gdst.functions.rotate_reference_cell(referenced_unit_cell, dtheta, turn_center)

                dtheta = unit_cell_distance / relative_radius

                gdst.functions.rotate_reference_cell(referenced_unit_cell, i * dtheta, turn_center)

                self.hole_references.append(referenced_unit_cell)


        for i, path in enumerate(self.paths):
            if i not in self.holey_ground_paths:
                path.turn(radius, angle, tolerance = tolerance, number_of_points = number_of_points, max_points = max_points, final_width = final_width, final_distance = final_distance, **self.layer_specs[i])
            else:
                path.x = self.structure.x
                path.y = self.structure.y
                path.direction = self.structure.direction

        self.endpoints['B'] = (self.structure.x, self.structure.y)
        self.endpoint_directions['B'] = self.structure.direction

    def arc(self, radius, initial_angle, final_angle, number_of_points = 0.01, max_points = 199, final_width = None, final_distance = None):
        for i, path in enumerate(self.paths):
            path.arc(radius, initial_angle, final_angle, number_of_points, max_points, final_width, final_distance, **self.layer_specs[i])
        self.endpoints['B'] = (self.structure.x, self.structure.y)
        self.endpoint_directions['B'] = self.structure.direction
        arc_length = radius * abs(final_angle - initial_angle)
        self.length += arc_length
        if arc_length < 0:
            print("Warning: Negative arc length detected, waveguide lengths will no longer be correct.")

    def fillet(self, radius, points_per_2pi = 128, max_points = 199, precision = 0.001):
        for path in self.paths:
            path.fillet(radius, points_per_2pi, max_points, precision)

    def fracture(self, max_points = 199, precision = 0.001):
        for path in self.paths:
            path.fracture(max_points, precision)

    def parametric(self, curve_function, curve_derivative = None, number_of_evaluations = 99, max_points = 199, final_width = None, final_distance = None):
        for i, path in enumerate(self.paths):
            path.parametric(curve_function, curve_derivative, number_of_evaluations, max_points, final_width, final_distance, **self.layer_specs[i])
        self.endpoints['B'] = (self.structure.x, self.structure.y)
        self.endpoint_directions['B'] = self.structure.direction

    def scale(self, scalex, scaley = None, center = (0, 0)):
        for path in self.paths:
            path.scale(scalex, scaley, center)
        scaley = 1 if scaley == None else scaley
        self.endpoints['A'] = ((self.endpoints['A'][0] - center[0]) * scalex * scaley + center[0], (self.endpoints['A'][1] - center[1]) * scalex * scaley + center[1])
        self.endpoints['B'] = ((self.endpoints['B'][0] - center[0]) * scalex * scaley + center[0], (self.endpoints['B'][1] - center[1]) * scalex * scaley + center[1])

    def meander(self, meander_direction, meanders, meander_length = None, total_length = None, inintial_length = None, io_difference = 0, turn_radius = None, number_of_points = 0.01, max_points = 199):
        direction = {'r': - np.pi, 'l': np.pi}
        if turn_radius == None:
            turn_radius = 14 * (self.paths[0].w + 2 * self.paths[1].w)
        length_turns = (meanders - 1) * turn_radius * np.pi
        if inintial_length == None:
            inintial_length = meander_length
        if io_difference == 0:
            final_length = inintial_length
        if meander_length == None and inintial_length == None:
            length_turns = (meanders - 1) * turn_radius * np.pi
            meander_length = (total_length - length_turns) / meanders
            if meander_length < 0:
                print("Warning: Enforced turn radius results in a negative meander length.")
            inintial_length = meander_length
            final_length = meander_length
        if meander_length == None and inintial_length != None:
            final_length = inintial_length + io_difference
            meander_length = (total_length - inintial_length - length_turns - final_length) / (meanders - 2)
            if meander_length < 0:
                print("Warning: Enforced turn radius results in a negative meander length.")
        j = 0
        self.segment(length = inintial_length)
        self.turn(radius = turn_radius, angle = direction[meander_direction], number_of_points = number_of_points, max_points = max_points)
        while j != meanders - 2:
            if j % 2 == 0:
                self.segment(length = meander_length)
                self.turn(radius = turn_radius, angle = - direction[meander_direction], number_of_points = number_of_points, max_points = max_points)
            else:
                self.segment(length = meander_length)
                self.turn(radius = turn_radius, angle = direction[meander_direction], number_of_points = number_of_points, max_points = max_points)
            j += 1
        self.segment(length = final_length)

        total_length = length_turns + (meanders - 2) * meander_length + inintial_length + final_length
        if total_length < 0:
            print("Warning: Negative total meander length detected, waveguide lengths will no longer be correct.")

        self.endpoints['B'] = (self.structure.x, self.structure.y)
        self.endpoint_directions['B'] = self.structure.direction

    def translate(self, dx, dy):
        for path in self.paths:
            path.translate(dx, dy)
        for endpoint in self.endpoints:
            self.endpoints[endpoint] = (self.endpoints[endpoint][0] + dx, self.endpoints[endpoint][1] + dy)
        for unit_cell in self.path_holey_ground:
            unit_cell.translate(dx, dy)
        for unit_cell in self.hole_references:
            unit_cell.translate(dx, dy)

    def rotate(self, angle, center = (0, 0)):
        for path in self.paths:
            path.rotate(angle, center)
        for endpoint in self.endpoints:
            centered_endpoint = np.array(self.endpoints[endpoint]) - np.array(center)
            centered_endpoint = np.array(gdst.functions.VecRot(angle, centered_endpoint))
            self.endpoints[endpoint] = tuple(centered_endpoint + center)
            self.endpoint_directions[endpoint] += angle
        for unit_cell in self.path_holey_ground:
            unit_cell.rotate(angle, center)
        for unit_cell in self.hole_references:
            gdst.functions.rotate_reference_cell(unit_cell, angle, center)

    def connect(self, endpoint_self, destination_structure, destination_endpoint, rotate = True):

        # Rotate to correct orientation
        if rotate:
            self.rotate((destination_structure.endpoint_directions[destination_endpoint] - self.endpoint_directions[endpoint_self] - np.pi) % (2 * np.pi))

        # Move to connect endpoints
        delta = tuple(map(sub, destination_structure.endpoints[destination_endpoint], self.endpoints[endpoint_self]))
        self.translate(delta[0], delta[1])

        return self

    def copy(self, dx = 0, dy = 0):

        # Create deepcopy of self and return the copy
        # Do not copy.deepcopy(self) directly, as this will also copy the connections and screw things up
        # Instead, just fill a new instance of the class with copies of relevant internal structure

        new_waveguide = Waveguide(path_list = copy.deepcopy(self.paths),
                                  path_types = copy.deepcopy(self.path_types),
                                  layer_specs = copy.deepcopy(self.layer_specs),
                                  hole_radius = copy.deepcopy(self.hole_radius),
                                  hole_distance = copy.deepcopy(self.requested_hole_distance),
                                  path_holey_ground = copy.deepcopy(self.path_holey_ground),
                                  hole_references = copy.deepcopy(self.hole_references),
                                  length = copy.deepcopy(self.length),
                                  endpoints = copy.deepcopy(self.endpoints),
                                  endpoint_directions = copy.deepcopy(self.endpoint_directions))

        for reference in new_waveguide.hole_references:
            reference.ref_cell = new_waveguide.holey_unit_cells[0]

        new_waveguide.translate(dx, dy)

        return new_waveguide
