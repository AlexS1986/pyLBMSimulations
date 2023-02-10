import logging

from SolidLBM.util.geometry.Point import Point
from SolidLBM.mesher.mesh.Mesh import Mesh
from SolidLBM.solver.solution.Computation import Computation


def mesher():
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

    return mesh


def solver():
    comp = Computation(working_directory='./', name=name)
    comp.execute_command_sequence()

    return comp


name = 'strip'

logging.basicConfig(
        filename = name + '.log',
        encoding = 'utf-8',
        format = '%(levelname)s -- %(filename)s (%(lineno)s) -- %(message)s',
        filemode = 'w',
        level = logging.DEBUG,
    )

h = 2**-5
seed_point = Point(1 + h/2.0, h/2.0, 0.0)

mesh = mesher()
comp = solver()

from post_eval import *
name = 'strip'

kc = 0.042
fig_sif = plot_sif(name, norm=kc, tip=1)
add_k_crit(fig_sif, 1)
add_velocity(fig_sif, name, color='tab:green', tip=1)
add_length(fig_sif, name)
fig_sif = add_legend(fig_sif)
save_plot(fig_sif, './sif_k', dpi=300,
)

