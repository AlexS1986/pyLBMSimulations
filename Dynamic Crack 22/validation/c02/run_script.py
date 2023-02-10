import logging

from SolidLBM.util.geometry.Point import Point
from SolidLBM.mesher.mesh.Mesh import Mesh
from SolidLBM.solver.solution.Computation import Computation

from post_eval import *

name = 'strip'
h = 2**-4

seed_point = Point(1 + h/2.0, h/2.0, 0.0)

logging.basicConfig(
        filename = name + '.log',
        encoding = 'utf-8',
        format = '%(levelname)s -- %(filename)s (%(lineno)s) -- %(message)s',
        filemode = 'w',
        level = logging.DEBUG,
    )

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

mesh = mesher()
comp = solver()

l = 1.0
w0 = 0.1
fig_sif = plot_sif(name)
add_mandal_orig(fig_sif, name, w0*2, l)
#add_regression(fig_sif, name, guess=2*w0, model='c')
#add_lowpass(fig_sif, name)
add_stat_mean(fig_sif, name, tf=15)
add_stat_median(fig_sif, name, tf=15)
save_plot(fig_sif, './sif_h4', dpi=300)

