import SolidLBM.mesher.mesh.Point as Point_Module
import SolidLBM.mesher.mesh.Mesh as Mesh_Module
import SolidLBM.solver.solution.Computation as Computation_Module
import SolidLBM.mesher.mesh.InitialCondition as Condition_Module

#h = 0.25
h = 0.05 # 0.5
seed_point = Point_Module.Point(-0.25-h/2, 0.0, 0.0)
#seed_point = Point_Module.Point(0.0, 0.0, 0.0)

mesh1 = Mesh_Module.Mesh(name='block', working_directory='./', cell_size=h, seed_point=seed_point)
mesh1.create_mesh_neighbor_points()
mesh1.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh1.plot_mesh()
mesh1.print_to_file()
mesh1.print_alt_bc_to_file()
mesh1.tmp_print_initial_conditions_to_file()


computation1 = Computation_Module.Computation(working_directory='./', name='block')
computation1.execute_command_sequence()
