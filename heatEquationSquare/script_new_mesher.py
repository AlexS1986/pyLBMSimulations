import SolidLBM.mesher.mesh.Point as Point_Module
import SolidLBM.mesher.mesh.Mesh as Mesh_Module
import SolidLBM.mesher.mesh.InitialCondition as Condition_Module
import SolidLBM.solver.solution.Computation as Computation_Module

# this is to generate the init file
#cond1 = Condition_Module.Condition(name='disc_with_hole', working_directory='/Users/Sikang Yan/Documents/pyLBM/examples/circle/')
#cond1.init_condition()
#cond1.tmp_print_initial_conditions_to_file()

seed_point = Point_Module.Point(0.25, 0.0,0.0)
mesh1 = Mesh_Module.Mesh(name='heatEquation', working_directory='/Users/alex/Work/LBM/PythonImplementierung/pyLBM/examples/heatEquationSquare/', cell_size=0.05, seed_point=seed_point)


mesh1.create_mesh_neighbor_points()
#mesh1.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh1.plot_mesh()
mesh1.print_to_file()
#mesh1.tmp_print_initial_conditions_to_file()

# Points = mesh1.Points
# Config = mesh1.Config


computation1 = Computation_Module.Computation(working_directory='/Users/alex/Work/LBM/PythonImplementierung/pyLBM/examples/heatEquationSquare/', name='heatEquation')
computation1.execute_command_sequence()
