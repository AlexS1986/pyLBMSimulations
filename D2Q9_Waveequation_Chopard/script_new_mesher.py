import SolidLBM.mesher.mesh.Point as Point_Module
import SolidLBM.mesher.mesh.Mesh as Mesh_Module
import SolidLBM.mesher.mesh.InitialCondition as Condition_Module
import SolidLBM.solver.solution.Computation as Computation_Module



seed_point = Point_Module.Point(0.25, 0.0,0.0)
mesh1 = Mesh_Module.Mesh(name='disc_with_hole', working_directory='Pfad setzen', cell_size=0.025, seed_point=seed_point,mesh_type="D2Q9")


mesh1.create_mesh_neighbor_points()
mesh1.plot_mesh()
mesh1.print_to_file()



computation1 = Computation_Module.Computation(working_directory='Pfad setzen', name='disc_with_hole')
computation1.execute_command_sequence()
