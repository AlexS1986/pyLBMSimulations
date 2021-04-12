import MeshGeneration as mg
import Solver as sol


point1 = mg.Point(0.0,0.0,0.0)
point1.PositionIndex = list([0, 0])
mesh1 = mg.Mesh(name='cutsquare', working_directory='./', cell_size=0.125, seed_point=point1)
mesh1.create_mesh_neighbor_points()
mesh1.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh1.plot_mesh()
mesh1.print_to_file()
mesh1.tmp_print_initial_conditions_to_file()


computation1 = sol.Computation(working_directory='./', name='cutsquare')
computation1.execute_command_sequence()
