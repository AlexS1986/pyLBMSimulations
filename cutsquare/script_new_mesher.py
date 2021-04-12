import SolidLBM.mesher.mesh.Point as Point_Module
import SolidLBM.mesher.mesh.Mesh as Mesh_Module
import SolidLBM.solver.solution.Computation as Computation_Module


seed_point = Point_Module.Point(0.25, 0.0,0.0)
mesh1 = Mesh_Module.Mesh(name='cutsquare', working_directory='/Users/alex/Work/LBM/PythonImplementierung/pyLBM/examples/cutsquare/', cell_size=0.025, seed_point=seed_point)
mesh1.create_mesh_neighbor_points()
mesh1.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh1.plot_mesh()
mesh1.print_to_file()
mesh1.tmp_print_initial_conditions_to_file()


computation1 = Computation_Module.Computation(working_directory='/Users/alex/Work/LBM/PythonImplementierung/pyLBM/examples/cutsquare/', name='cutsquare')
computation1.execute_command_sequence()
