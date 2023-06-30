import SolidLBM.util.geometry.Point as Point_Module
import SolidLBM.mesher.mesh.Mesh as Mesh_Module


h = 2**-7
name='dynamic'
seed_point = Point_Module.Point(h/2, h/2, 0.0)

mesh = Mesh_Module.Mesh(
    name=name,
    working_directory='./',
    cell_size=h,
    seed_point=seed_point,
    mesh_type="D2Q9",
)
mesh.create_mesh_neighbor_points(verbose=True)
mesh.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh.plot_mesh()
mesh.print_to_file()
mesh.print_alt_bc_to_file()
mesh.tmp_print_initial_conditions_to_file()

