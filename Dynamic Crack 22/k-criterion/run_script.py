import logging

from SolidLBM.util.geometry.Point import Point
from SolidLBM.mesher.mesh.Mesh import Mesh
from SolidLBM.solver.solution.Computation import Computation


name = 'strip'
h = 2**-6

logging.basicConfig(
        filename = name + '.log',
        encoding = 'utf-8',
        format = '%(levelname)s -- %(filename)s (%(lineno)s) -- %(message)s',
        filemode = 'w',
        level = logging.WARNING,
    )

seed_point = Point(1 + h/2.0, h/2.0, 0.0)

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

comp = Computation(working_directory='./', name=name)

import time
from resource import *

process_start = time.process_time()
user_start = getrusage(RUSAGE_SELF).ru_utime
sys_start = getrusage(RUSAGE_SELF).ru_stime

comp.execute_command_sequence()

process_time = time.process_time() - process_start
user_time = getrusage(RUSAGE_SELF).ru_utime - user_start
sys_time = getrusage(RUSAGE_SELF).ru_stime - sys_start

print("propagation")
print("process time: ", process_time)
print("resource user:", user_time)
print("resource sys: ", sys_time)

from post_eval import *

kc = 0.006
fig_sif = plot_sif(name, norm=kc, tip=1)
add_k_crit(fig_sif, 1)
add_velocity(
    fig_sif, name, tip=1,
    color='tab:green',
)
add_length(fig_sif, name)
fig_sif = add_legend(fig_sif)
save_plot(
    fig_sif, './sif_k', dpi=450, 
    tex=True,
)

