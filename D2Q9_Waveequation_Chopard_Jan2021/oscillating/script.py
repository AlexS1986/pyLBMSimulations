import MeshGeneration as mg
import Solver as sol


point1 = mg.Point(0.25,0.0,0.0)
#point1.PositionIndex = list([0, 0])
mesh1 = mg.Mesh(name='disc_with_hole', working_directory='/Users/alex/Desktop/circle_GAMM/', cell_size=0.025, seed_point=point1)
mesh1.create_mesh_neighbor_points()
mesh1.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh1.plot_mesh()
mesh1.print_to_file()
mesh1.tmp_print_initial_conditions_to_file()


computation1 = sol.Computation(working_directory='/Users/alex/Desktop/circle_GAMM/', name='disc_with_hole')
computation1.execute_command_sequence()