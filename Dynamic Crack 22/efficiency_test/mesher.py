import logging

from SolidLBM.util.geometry.Point import Point
from SolidLBM.mesher.mesh.Mesh import Mesh

name = 'strip'
h = 2**-4

logging.basicConfig(
    filename = name + '.log',
    encoding = 'utf-8',
    format = '%(levelname)s -- %(filename)s (%(lineno)s) -- %(message)s',
    filemode = 'w',
    level = logging.INFO,
)

seed_point = Point(2 + h/2.0, h/2.0, 0.0)

mesh = Mesh(
    name=name,
    working_directory='./',
    cell_size=h,
    seed_point=seed_point,
)

mesh.create_mesh_neighbor_points(verbose=True)
mesh.print_to_file()
mesh.tmp_print_initial_conditions_to_file()
mesh.compute_cell_volumes_areas_boundary_names_at_boundary_points()
mesh.print_alt_bc_to_file()
mesh.plot_mesh(dpi=300)

