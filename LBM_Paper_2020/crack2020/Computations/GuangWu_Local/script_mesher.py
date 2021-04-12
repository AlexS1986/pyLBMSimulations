#!/usr/bin/env python
# -*- coding: utf-8 -*-
import SolidLBM.mesher.mesh.Mesh as Mesh_Module
import SolidLBM.mesher.mesh.Point as Point_Module

cell_size = 1.0/42.0
#cell_size = 1.0/7.0

point1 = Point_Module.Point(cell_size/2.0, cell_size/2.0, 0.0)
mesh1 = Mesh_Module.Mesh(name='disc_with_crack', working_directory='./', cell_size=cell_size,
                 seed_point=point1)
mesh1.create_mesh_neighbor_points()
mesh1.plot_mesh()
mesh1.print_to_file()
mesh1.tmp_print_initial_conditions_to_file()
mesh1.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh1.print_alt_bc_to_file()
