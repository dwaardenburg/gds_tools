import copy

import numpy as np
import gdspy as gp
import gds_tools as gdst
from operator import add, sub

class GDSComponent:

    def __init__(
        self,
        elements,
        align_points = None,
        align_angles = None,
        align_dimensions = None):

        self.align_points = {} if align_points is None else align_points
        self.align_angles = {} if align_angles is None else align_angles
        self.align_dimensions = {} if align_dimensions is None else align_dimensions

        self.polygons = []

        self.previous = []
        self.next = []

        if not isinstance(elements, list):
            elements = [elements]
        
        for element in elements:
            if isinstance(element, GDSComponent):
                self.next.append(element)
                element.previous.append(self)
            else:
                self.polygons.append(element)
    
    #===================
    # Connect structure \\
    #=========================================================================
    # Arguments:    align_self :   reference to A, B, C, etc.                  ||
    #                           name of endpoint of this structure          ||
    #               to      :   GDSComponent class input to which to connect ||
    #               align_to   :   endpoint to which to connect to             ||
    #               offset  :   (x, y) offset from target align_points         ||
    #               obj_link:   use Python object reference for linked list ||
    #                           entry rather than endpoint key              ||
    #=========================================================================
    def connect(
        self,
        align_self,
        to_component,
        align_to,
        offset = (0, 0),
        rotate = True):

        # Rotate to correct orientation
        if self.align_angles != None and to_component.align_angles != None and self.align_angles and to_component.align_angles and rotate:
            self.rotate((to_component.align_angles[align_to] - self.align_angles[align_self] - np.pi) % (2 * np.pi))

        # Move to connect align_points
        delta = tuple(map(sub, to_component.align_points[align_to], self.align_points[align_self]))
        delta = tuple(map(add, delta, offset))
        self.translate(delta[0], delta[1])

        # Update linked list
        to_component.next.append(self)
        self.previous.append(to_component)

    #======================
    # Disconnect structure \\
    #=========================================================================
    # Only removes the references in the linked lists                       ||
    #=========================================================================
    def disconnect(self):

        for e in self.align_points:
            if e in self.next:
                for ne in self.next[e].align_points:
                    if ne in self.next[e].prev and self.next[e].prev[ne] == self:
                        del self.next[e].prev[ne]
                del self.next[e]
            if e in self.previous:
                for pe in self.previous[e].align_points:
                    if pe in self.previous[e].next and self.previous[e].next[pe] == self:
                        del self.previous[e].next[pe]
                del self.previous[e]

        return self

    def copy(
        self,
        deep_copy = False,
        signal_from = None):
        """Makes a new copy of the structure,
        i.e. allocates new memory and copies the data contents to it.

        Returns:
        gds_tools.GDSComponent: pointer of copy of original
        """

        if not isinstance(signal_from, list):
            signal_from = [signal_from]

        # Create deepcopy of self and return the copy
        # Do not copy.deepcopy(self) directly, as this will also copy the connections and screw things up
        # Instead, just fill a new instance of the class with copies of relevant internal structure
        new_component = copy.deepcopy(self)
        new_component.previous = []
        new_component.next = []

        for i, component in enumerate(self.next):
            if self.next[i] not in signal_from and self.next:
                signal_from.append(self.next[i])
                new_component.next.append(component.copy(signal_from = signal_from))

        return new_component

    #==================
    # Rotate structure \\
    #=========================================================================
    # When rotating, all align_points need to move wih the rotation,           ||
    # for this reason we need our own rotation function.                    ||
    #                                                                       ||
    # Arguments:    rad         :   amount of radians to rotate             ||
    #               signal_from :   tells linked structure who sent command ||
    #=========================================================================
    def rotate(
        self,
        angle,
        signal_from = None):

        if not isinstance(signal_from, list):
            signal_from = [signal_from]
            
        # Rotate the structure with standard gdspy .rotation() function
        if self not in signal_from:
            signal_from.append(self)

        for polygon in self.polygons:
            polygon.rotate(angle)

        # Rotate the align_points
        for key in self.align_points:
            self.align_points[key] = tuple(gdst.functions.VecRot(angle, self.align_points[key]))
            if self.align_angles:
                self.align_angles[key] += angle

        # Rotate all connected structures
        for previous_component in self.previous:
            if previous_component not in signal_from and self.previous:
                signal_from.append(previous_component)
                previous_component.rotate(angle, signal_from = signal_from)

        for next_component in self.next:
            if next_component not in signal_from and self.next:
                signal_from.append(next_component)
                next_component.rotate(angle, signal_from = signal_from)

        return self

    #=====================
    # Translate structure \\
    #=========================================================================
    # Arguments:    delta       :   vector (x, y) how much you want to move ||
    #               signal_from :   tells linked structure who sent command ||
    #=========================================================================
    def translate(
        self,
        dx,
        dy,
        signal_from = None):

        if not isinstance(signal_from, list):
            signal_from = [signal_from]

        # Translate structure with standard gdspy .translate() function
        if self not in signal_from:
            signal_from.append(self)
        
        for polygon in self.polygons:
            polygon.translate(dx, dy)

        # Translate all the align_points
        for i in self.align_points:
            self.align_points[i] = tuple(map(add, self.align_points[i], (dx, dy)))

        # Translate all connected structures
        for previous_component in self.previous:
            if previous_component not in signal_from and self.previous:
                signal_from.append(previous_component)
                previous_component.translate(dx, dy, signal_from = signal_from)

        for next_component in self.next:
            if next_component not in signal_from and self.next:
                signal_from.append(next_component)
                next_component.translate(dx, dy, signal_from = signal_from)

        return self

class PathComponent(GDSComponent):

    def __init__(
        self,
        widths,
        initial_points = None,
        path_counts = [1],
        distances = [0],
        types = None,
        layer_specs = []
        ):

        """Create a waveguide object that allows an array of path elements to be treated as an single path.

        Args:
            path_list (list): the gp.path classes to be bundled.
            layer_list (list): the layers of each of paths.
            holey_ground_path ([type], optional): [description]. Defaults to None.
            hole_references ([type], optional): [description]. Defaults to None.
            hole_distance (float, optional): distance between hole centres. Defaults to None.
            hole_radius (float, optional): radius of the holes. Defaults to None.
            length (float, optional): length of the path. Defaults to 0.
            align_points (dict, optional): dictionairy defining the align_points (x, y). Defaults to None.
            align_angles (dict, optional): dictionairy defining the endpoint directions. Defaults to None.
        """

        self.widths = widths
        self.initial_points = [(0, 0)] * len(widths) if initial_points is None else initial_points
        self.path_counts = path_counts
        self.distances = distances
        self.types = ['stripline'] if types is None else types
        self.layer_specs = layer_specs
        self.length = 0
        self.paths = []
        self.type_index = {}
        for i, width in enumerate(self.widths):
            self.paths.append(
                gp.Path(
                    width,
                    initial_point = self.initial_points[i],
                    number_of_paths = self.path_counts[i],
                    distance = self.distances[i]
                    )
                )
            self.type_index[self.types[i]] = i
        
        self.main_path = self.paths[0]
        self.n = self.main_path.n

        super().__init__(self.paths)

        self.initial_keys = []
        self.final_keys = []
        if self.n == 1:
            self.initial_keys = ['A']
            self.final_keys = ['B']
        else:
            for n in range(self.n):
                self.initial_keys.append('A' + str(n))
                self.final_keys.append('B' + str(n))

        align_keys = self.initial_keys + self.final_keys
        self.align_points = dict.fromkeys(align_keys)
        self.align_angles = dict.fromkeys(align_keys)
        for key in align_keys:
            self.align_dimensions[key] = {}
            for i, type in enumerate(self.types):
                self.align_dimensions[key][type] = {}
                self.align_dimensions[key][type]['number'] = self.paths[i].n
                self.align_dimensions[key][type]['hole_radius'] = False
                self.align_dimensions[key][type]['hole_distance'] = False
        self.update_endpoints(self.initial_keys)
        for key in self.initial_keys:
            self.align_angles[key] = self.direction + np.pi
        self.update_endpoints(self.final_keys)
        
        self.requested_hole_distance = None
        self.hole_radius = None
        self.path_holey_ground = []
        self.hole_references = []
        self.holey_ground_paths = []
        self.holey_unit_cells = []

    def update_endpoints(
        self,
        align_keys 
        ):
        self.x = self.main_path.x
        self.y = self.main_path.y
        self.w = self.main_path.w
        self.distance = self.main_path.distance

        dir_str = {'+x': 0, '-x': np.pi, '+y': np.pi / 2, '-y': np.pi / 2}
        if isinstance(self.main_path.direction, float):
            self.direction = self.main_path.direction
        else:
            self.direction = dir_str[self.main_path.direction]

        cos_dir = np.cos(self.direction)
        sin_dir = np.sin(self.direction)

        for n in range(self.n):
            d0 = n * self.distance - (self.n - 1) * self.distance * 0.5
            self.align_points[align_keys[n]] = (self.x + d0 * sin_dir, self.y - d0 * cos_dir)
            self.align_angles[align_keys[n]] = self.direction

        for key in align_keys:
            for i, type in enumerate(self.types):
                self.align_dimensions[key][type]['width'] = 2 * self.paths[i].w
                self.align_dimensions[key][type]['distance'] = self.paths[i].distance

    def segment(
        self,
        length,
        direction = None,
        final_widths = None,
        final_distances = None,
        axis_offset = 0):

        if direction != None:
            dir_str = {'+x': 0, '-x': np.pi, '+y': np.pi / 2, '-y': np.pi / 2}
            self.direction = direction if isinstance(direction, float) else dir_str[direction]
        if final_widths == None:
            final_widths = [None] * len(self.paths)
        if final_distances == None:
            final_distances = [None] * len(self.paths)

        for i, path in enumerate(self.paths):
            path.segment(length, direction, final_widths[i], final_distances[i], axis_offset, **self.layer_specs[i])

        self.length += length
        if length < 0:
            print("Warning: Negative segment length detected, waveguide lengths will no longer be correct.")
        self.update_endpoints(self.final_keys)

    def klopfenstein(
        self,
        length,
        final_widths,
        final_distances,
        point_distance):

        for i, path in enumerate(self.paths):
            direction = path.direction
            dir_str = {'+x': 0, '-x': np.pi, '+y': np.pi / 2, '-y': np.pi / 2}
            direction = direction if isinstance(direction, float) else dir_str[direction]
            x, y = path.x, path.y
            path.rotate(- direction, center=(x, y))
            width = lambda x: self.widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.widths[i])
            if self.distances[i] == 0:
                distance = None
            else:
                distance = lambda x: self.distances[i] + gdst.functions.biquadratic_func(x) * (final_distances[i] - self.distances[i])
            path.parametric(
                    lambda x: (x * length, 0),
                    lambda x: (length, 0),
                    final_width = width,
                    final_distance = distance,
                    max_points = 4094,
                    number_of_evaluations = int(2 * length / point_distance),
                    **self.layer_specs[i]
                    )
            path.rotate(direction, center=(x, y))

        self.length += length
        if length < 0:
            print("Warning: Negative tapering length detected, waveguide lengths will no longer be correct.")
        self.update_endpoints(self.final_keys)

    def turn(
        self,
        radius,
        angle,
        tolerance = 0.01,
        number_of_points = None,
        max_points = 199,
        final_width = None,
        final_distance = None):

        angles = {'r': - np.pi / 2, 'l': np.pi / 2, 'rr': - np.pi, 'll': np.pi}
        angle = angle if not isinstance(angle, str) else angles[angle]

        for i, path in enumerate(self.paths):
            path.turn(radius, angle, tolerance = tolerance, number_of_points = number_of_points, max_points = max_points, final_width = final_width, final_distance = final_distance, **self.layer_specs[i])

        turn_length = radius * abs(angle)
        self.length += turn_length
        if turn_length < 0:
            print("Warning: Negative turn length detected, waveguide lengths will no longer be correct.")
        self.update_endpoints(self.final_keys)

    def arc(
        self,
        radius,
        initial_angle,
        final_angle,
        number_of_points = 0.01,
        max_points = 199,
        final_width = None,
        final_distance = None):

        for i, path in enumerate(self.paths):
            path.arc(radius, initial_angle, final_angle, number_of_points, max_points, final_width, final_distance, **self.layer_specs[i])
        
        arc_length = radius * abs(final_angle - initial_angle)
        self.length += arc_length
        if arc_length < 0:
            print("Warning: Negative arc length detected, waveguide lengths will no longer be correct.")
        self.update_endpoints(self.final_keys)

    def fillet(
        self,
        radius,
        points_per_2pi = 128,
        max_points = 199,
        precision = 0.001):

        for path in self.paths:
            path.fillet(radius, points_per_2pi, max_points, precision)

    def fracture(
        self,
        max_points = 199,
        precision = 0.001):

        for path in self.paths:
            path.fracture(max_points, precision)

    def parametric(
        self,
        curve_function,
        curve_derivative = None,
        number_of_evaluations = 99,
        max_points = 199,
        final_width = None,
        final_distance = None):

        for i, path in enumerate(self.paths):
            path.parametric(curve_function, curve_derivative, number_of_evaluations, max_points, final_width, final_distance, **self.layer_specs[i])

        self.update_endpoints(self.final_keys)

    def scale(
        self,
        scalex,
        scaley = None,
        center = (0, 0)):

        for path in self.paths:
            path.scale(scalex, scaley, center)

        print("Warning: Unfinished function, path endpoints are no longer correct.")

        scaley = 1 if scaley == None else scaley
        for key in self.align_points:
            self.align_points[key] = ((self.align_points[key][0] - center[0]) * scalex * scaley + center[0], (self.align_points[key][1] - center[1]) * scalex * scaley + center[1])

    def meander(
        self,
        meander_direction,
        meanders,
        meander_length = None,
        total_length = None,
        inintial_length = None,
        io_difference = 0,
        turn_radius = None,
        number_of_points = 0.01,
        max_points = 199):

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

class HoleyPath:

    def __init__(
        self,
        width,
        initial_point = (0, 0),
        number_of_paths = 1,
        distance = 0,
        type = 'holey_ground',
        hole_radius = None,
        hole_distance = None,
        layer_specs = []
        ):

        self.width = width
        self.number_of_paths = number_of_paths
        self.distance = distance
        self.type = type
        self.layer_specs = layer_specs

        self.x = []
        self.y = []
        self.w = []
        self.n = []
        self.direction = []
        self.distance = []
        self.length = []
        self.properties = []

        self.holey_unit_cells = []
        self.holey_ground_paths = []
        self.path_holey_ground = []
        self.hole_references = []

        self.requested_hole_distance = hole_distance
        self.hole_radius = hole_radius

        if 'holey_unit_cell_index' in globals():
            global holey_unit_cell_index
        else:
            global holey_unit_cell_index
            holey_unit_cell_index = 0

        self.holey_ground_paths = self.path_types['holey_grounds']
        holey_path_index = self.path_types['holey_grounds'][0]
        self.holey_ground_layer = self.layer_specs[holey_path_index]['layer']
        self.holey_ground_datatype = self.layer_specs[holey_path_index]['datatype']
        self.holey_ground_spec = self.layer_specs[holey_path_index]

        self.hole_distances = [0] * len(self.paths)
        self.round_verticies = 6

        hole_list = []
        for i in self.holey_ground_paths:
            effective_width = self.widths[i] - 2 * self.hole_radius
            total_vertical_hole_distances = int((effective_width - effective_width % self.requested_hole_distance) / self.requested_hole_distance)
            self.hole_distances[i] = effective_width / total_vertical_hole_distances
            if self.endpoint_path_counts[i] == 1:
                for j in range(total_vertical_hole_distances + 1):
                    hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_width / 2), self.hole_radius, number_of_points = self.round_verticies))
                for j in range(total_vertical_hole_distances):
                    hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_width / 2), self.hole_radius, number_of_points = self.round_verticies))
            if self.endpoint_path_counts[i] == 2:
                for j in range(total_vertical_hole_distances + 1):
                    hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_width / 2 - self.distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                    hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_width / 2 + self.distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                for j in range(total_vertical_hole_distances):
                    hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_width / 2 - self.distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
                    hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_width / 2 + self.distances[i] / 2), self.hole_radius, number_of_points = self.round_verticies))
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


        self.length = length

    def construct_holey_unit_cell(
        self,
        widths,
        distances):

        self.hole_distances = [0] * len(self.paths)
        hole_list = []
        for i in self.holey_ground_paths:
            effective_path_width = widths[i] - 2 * self.hole_radius
            total_vertical_hole_distances = int((effective_path_width - effective_path_width % self.requested_hole_distance) / self.requested_hole_distance)
            self.hole_distances[i] = effective_path_width / total_vertical_hole_distances
            if self.endpoint_path_counts[i] == 1:
                for j in range(total_vertical_hole_distances + 1):
                    hole_list.append(gp.Round((0, j * self.hole_distances[i] - effective_path_width / 2), self.hole_radius, number_of_points = self.round_verticies))
                for j in range(total_vertical_hole_distances):
                    hole_list.append(gp.Round(((np.sqrt(3) / 2) * self.hole_distances[i], (j + 1 / 2) * self.hole_distances[i] - effective_path_width / 2), self.hole_radius, number_of_points = self.round_verticies))
            if self.endpoint_path_counts[i] == 2:
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

    def segment(
        self,
        length,
        direction = None,
        final_widths = None,
        final_distances = None,
        axis_offset = 0):

        if final_widths == None and final_distances == None:
            effective_length = length - self.holey_unit_cell_width
            total_unit_cells = int((effective_length - effective_length % 2 * self.holey_unit_cell_width) / (2 * self.holey_unit_cell_width))
            unit_cell_distance = effective_length / total_unit_cells

            for i in range(total_unit_cells):
                referenced_unit_cell = gp.CellReference(self.holey_unit_cells[self.active_unit_cell], origin = (self.main_path.x, self.main_path.y), rotation = np.degrees(self.align_angles['B']))

                dx = np.cos(self.align_angles['B']) * (self.holey_unit_cell_width + unit_cell_distance) / 2
                dy = np.sin(self.align_angles['B']) * (self.holey_unit_cell_width + unit_cell_distance) / 2
                referenced_unit_cell.translate(dx, dy)

                dx = np.cos(self.align_angles['B']) * i * unit_cell_distance
                dy = np.sin(self.align_angles['B']) * i * unit_cell_distance

                referenced_unit_cell.translate(dx, dy)
                self.hole_references.append(referenced_unit_cell)

            for i, path in enumerate(self.paths):
                if i not in self.holey_ground_paths:
                    path.segment(length, direction, None, None, axis_offset, **self.layer_specs[i])
                else:
                    path.x = self.main_path.x
                    path.y = self.main_path.y
                    path.direction = self.direction
        else:
            holey_ground = []
            holey_index = 0
            
            if self.requested_hole_distance != None and self.hole_radius != None:
                holey_ground.append(gp.Path(2 * path.w, initial_point=(path.x, path.y), number_of_paths = path.n, distance = path.distance))
                holey_ground[holey_index].direction = self.direction
                holey_ground[holey_index].segment(length, direction, final_widths[i], final_distances[i], axis_offset)
                holey_index += 1
            path.x = self.main_path.x
            path.y = self.main_path.y
            path.direction = self.direction

            self.widths = final_widths
            self.distances = final_distances
            if self.requested_hole_distance != None and self.hole_radius != None:
                self.path_holey_ground.append(generate_holey_ground(gp.boolean(holey_ground, [], 'or', **self.holey_ground_spec),
                                                                                     self.hole_radius,
                                                                                     self.requested_hole_distance,
                                                                                     self.round_verticies,
                                                                                     **self.holey_ground_spec))
                self.construct_holey_unit_cell(widths = final_widths, distances = final_distances)
        self.length += length
        if length < 0:
            print("Warning: Negative segment length detected, waveguide lengths will no longer be correct.")

        self.align_points['B'] = (self.main_path.x, self.main_path.y)

    def klopfenstein(
        self,
        length,
        final_widths,
        final_distances,
        point_distance):

        holey_ground = []
        holey_index = 0
        for i, path in enumerate(self.paths):
            if i not in self.holey_ground_paths:
                direction = path.direction
                d = {'+x': 0, '-x': np.pi, '+y': np.pi / 2, '-y': np.pi / 2}
                direction = direction if isinstance(direction, float) else (d[direction])
                x, y = path.x, path.y
                path.rotate(- direction, center=(x, y))
                if self.distances[i] == 0:
                    width = lambda x: self.widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.widths[i])
                    path.parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                else:
                    distance = lambda x: self.distances[i] + gdst.functions.biquadratic_func(x) * (final_distances[i] - self.distances[i])
                    width = lambda x: self.widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.widths[i])
                    path.parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, final_distance = distance, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                path.rotate(direction, center=(x, y))
            else:
                if self.requested_hole_distance != None and self.hole_radius != None:
                    holey_ground.append(gp.Path(2 * path.w, initial_point=(path.x, path.y), number_of_paths = path.n, distance = path.distance))
                    holey_ground[holey_index].direction = self.direction
                    if self.distances[i] == 0:
                        width = lambda x: self.widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.widths[i])
                        holey_ground[holey_index].parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                    else:
                        distance = lambda x: self.distances[i] + gdst.functions.biquadratic_func(x) * (final_distances[i] - self.distances[i])
                        width = lambda x: self.widths[i] + gdst.functions.biquadratic_func(x) * (final_widths[i] - self.widths[i])
                        holey_ground[holey_index].parametric(lambda x: (x * length, 0), lambda x: (length, 0), final_width = width, final_distance = distance, max_points = 4094, number_of_evaluations = int(2 * length / point_distance), **self.layer_specs[i])
                    holey_index += 1
                path.x = self.main_path.x
                path.y = self.main_path.y
                path.direction = self.direction

        if self.requested_hole_distance != None and self.hole_radius != None:
                self.path_holey_ground.append(generate_holey_ground(gp.boolean(holey_ground, [], 'or', **self.holey_ground_spec),
                                                                                     self.hole_radius,
                                                                                     self.requested_hole_distance,
                                                                                     self.round_verticies,
                                                                                     **self.holey_ground_spec))
                self.construct_holey_unit_cell(widths = final_widths, distances = final_distances)

        self.align_points['B'] = (self.main_path.x, self.main_path.y)

        for i, path in enumerate(self.paths):
            path.w = final_widths[i] / 2
            path.distance = final_distances[i]
        self.widths = final_widths
        self.distances = final_distances
        # if self.requested_hole_distance != None and self.hole_radius != None:
        #     self.construct_holey_unit_cell(widths = final_widths, distances = final_distances)
        self.length += length
        if length < 0:
            print("Warning: Negative tapering length length detected, waveguide lengths will no longer be correct.")

    def turn(
        self,
        radius,
        angle,
        tolerance = 0.01,
        number_of_points = None,
        max_points = 199,
        final_width = None,
        final_distance = None):

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
                referenced_unit_cell = referenced_unit_cell = gp.CellReference(self.holey_unit_cells[self.active_unit_cell], origin = (self.main_path.x, self.main_path.y), rotation = np.degrees(self.align_angles['B']))

                relative_radius = radius * np.sign(angle)
                dx_0 = np.cos(self.align_angles['B'] + np.pi / 2) * relative_radius
                dy_0 = np.sin(self.align_angles['B'] + np.pi / 2) * relative_radius
                turn_center = (self.main_path.x + dx_0, self.main_path.y + dy_0)
                dtheta = (self.holey_unit_cell_width + unit_cell_distance) / (2 * relative_radius)
                gdst.functions.rotate_reference_cell(referenced_unit_cell, dtheta, turn_center)

                dtheta = unit_cell_distance / relative_radius

                gdst.functions.rotate_reference_cell(referenced_unit_cell, i * dtheta, turn_center)

                self.hole_references.append(referenced_unit_cell)


        for i, path in enumerate(self.paths):
            if i not in self.holey_ground_paths:
                path.turn(radius, angle, tolerance = tolerance, number_of_points = number_of_points, max_points = max_points, final_width = final_width, final_distance = final_distance, **self.layer_specs[i])
            else:
                path.x = self.main_path.x
                path.y = self.main_path.y
                path.direction = self.direction

        self.align_points['B'] = (self.main_path.x, self.main_path.y)
        self.align_angles['B'] = self.direction

    def arc(
        self,
        radius,
        initial_angle,
        final_angle,
        number_of_points = 0.01,
        max_points = 199,
        final_width = None,
        final_distance = None):

        for i, path in enumerate(self.paths):
            path.arc(radius, initial_angle, final_angle, number_of_points, max_points, final_width, final_distance, **self.layer_specs[i])
        self.align_points['B'] = (self.main_path.x, self.main_path.y)
        self.align_angles['B'] = self.direction
        arc_length = radius * abs(final_angle - initial_angle)
        self.length += arc_length
        if arc_length < 0:
            print("Warning: Negative arc length detected, waveguide lengths will no longer be correct.")

    def fillet(
        self,
        radius,
        points_per_2pi = 128,
        max_points = 199,
        precision = 0.001):

        for path in self.paths:
            path.fillet(radius, points_per_2pi, max_points, precision)

    def fracture(
        self,
        max_points = 199,
        precision = 0.001):

        for path in self.paths:
            path.fracture(max_points, precision)

    def parametric(
        self,
        curve_function,
        curve_derivative = None,
        number_of_evaluations = 99,
        max_points = 199,
        final_width = None,
        final_distance = None):

        for i, path in enumerate(self.paths):
            path.parametric(curve_function, curve_derivative, number_of_evaluations, max_points, final_width, final_distance, **self.layer_specs[i])
        self.align_points['B'] = (self.main_path.x, self.main_path.y)
        self.align_angles['B'] = self.direction

    def scale(
        self,
        scalex,
        scaley = None,
        center = (0, 0)):

        for path in self.paths:
            path.scale(scalex, scaley, center)
        scaley = 1 if scaley == None else scaley
        self.align_points['A'] = ((self.align_points['A'][0] - center[0]) * scalex * scaley + center[0], (self.align_points['A'][1] - center[1]) * scalex * scaley + center[1])
        self.align_points['B'] = ((self.align_points['B'][0] - center[0]) * scalex * scaley + center[0], (self.align_points['B'][1] - center[1]) * scalex * scaley + center[1])

    def meander(
        self,
        meander_direction,
        meanders,
        meander_length = None,
        total_length = None,
        inintial_length = None,
        io_difference = 0,
        turn_radius = None,
        number_of_points = 0.01,
        max_points = 199):

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

        self.align_points['B'] = (self.main_path.x, self.main_path.y)
        self.align_angles['B'] = self.direction

    def translate(
        self,
        dx,
        dy):

        for path in self.paths:
            path.translate(dx, dy)
        for endpoint in self.align_points:
            self.align_points[endpoint] = (self.align_points[endpoint][0] + dx, self.align_points[endpoint][1] + dy)
        for unit_cell in self.path_holey_ground:
            unit_cell.translate(dx, dy)
        for unit_cell in self.hole_references:
            unit_cell.translate(dx, dy)

    def rotate(
        self,
        angle,
        center = (0, 0)):
        for path in self.paths:
            path.rotate(angle, center)
        for endpoint in self.align_points:
            centered_endpoint = np.array(self.align_points[endpoint]) - np.array(center)
            centered_endpoint = np.array(gdst.functions.VecRot(angle, centered_endpoint))
            self.align_points[endpoint] = tuple(centered_endpoint + center)
            self.align_angles[endpoint] += angle
        for unit_cell in self.path_holey_ground:
            unit_cell.rotate(angle, center)
        for unit_cell in self.hole_references:
            gdst.functions.rotate_reference_cell(unit_cell, angle, center)

    def copy(
        self,
        dx = 0,
        dy = 0):

        # Create deepcopy of self and return the copy
        # Do not copy.deepcopy(self) directly, as this will also copy the connections and screw things up
        # Instead, just fill a new instance of the class with copies of relevant internal structure

        new_waveguide = copy.deepcopy(self)
  
        self.previous = []
        self.next = []

        new_waveguide.translate(dx, dy)

        return new_waveguide

def generate_holey_ground(holey_ground, hole_radius, hole_distance, round_verticies, layer, datatype):
    bb = holey_ground.get_bounding_box()
    bb_width = abs(bb[1][0] - bb[0][0]) - 2 * hole_radius
    bb_height = abs(bb[1][1] - bb[0][1]) - 2 * hole_radius
    bb_x = (bb[0][0] + bb[1][0]) / 2
    bb_y = (bb[0][1] + bb[1][1]) / 2
    bb_min_x = bb_x - bb_width / 2
    bb_max_x = bb_x + bb_width / 2
    bb_min_y = bb_y - bb_height / 2
    bb_max_y = bb_y + bb_height / 2

    tot_h_double_distances = int((bb_width - bb_width % (np.sqrt(3) * hole_distance)) / (np.sqrt(3) * hole_distance))
    h_double_distance = bb_width / (tot_h_double_distances - 1 / 2)
    tot_v_distances = int((bb_height - bb_height % hole_distance) / hole_distance)
    v_distance = bb_height / (tot_v_distances - 1 / 2)

    circle = gp.Round((0, 0), hole_radius, number_of_points = round_verticies)
    circle_points = np.array(circle.polygons[0])
    circle_points_x = np.tile(circle_points[:,0], tot_h_double_distances)
    circle_points_y = np.tile(circle_points[:,1], tot_v_distances)

    circle_centers_x_0 = np.repeat(np.linspace(bb_min_x, bb_max_x, tot_h_double_distances), round_verticies)
    circle_centers_y_0 = np.repeat(np.linspace(bb_min_y, bb_max_y, tot_v_distances), round_verticies)
    coords_x_0 = np.tile(circle_points_x + circle_centers_x_0, tot_v_distances)
    coords_y_0 = np.tile(circle_points_y, tot_h_double_distances) + np.repeat(circle_centers_y_0, tot_h_double_distances)
    circle_polygons_0 = np.array([coords_x_0, coords_y_0]).T.reshape((tot_h_double_distances * tot_v_distances, round_verticies, 2))

    circle_centers_x_1 = np.repeat(np.linspace(bb_min_x + h_double_distance / 2, bb_max_x + h_double_distance / 2, tot_h_double_distances), round_verticies)
    circle_centers_y_1 = np.repeat(np.linspace(bb_min_y + v_distance / 2, bb_max_y + v_distance / 2, tot_v_distances), round_verticies)
    coords_x_1 = np.tile(circle_points_x + circle_centers_x_1, tot_v_distances)
    coords_y_1 = np.tile(circle_points_y, tot_h_double_distances) + np.repeat(circle_centers_y_1, tot_h_double_distances)
    circle_polygons_1 = np.array([coords_x_1, coords_y_1]).T.reshape((tot_h_double_distances * tot_v_distances, round_verticies, 2))

    circle_polygons = np.concatenate((circle_polygons_0, circle_polygons_1))
    polygonset = gp.PolygonSet(circle_polygons)
    return gp.fast_boolean(holey_ground, polygonset, 'and', layer = layer, datatype = datatype)
