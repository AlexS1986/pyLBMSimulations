import logging
from resource import *
import numpy as np

from SolidLBM.solver.solution.Computation import Computation

name = 'strip'

logging.basicConfig(
    filename = name + '.log',
    encoding = 'utf-8',
    format = '%(levelname)s -- %(filename)s (%(lineno)s) -- %(message)s',
    filemode = 'w',
    level = logging.INFO,
)

comp = Computation(working_directory='./', name=name)
# comp.execute_command_sequence()

""" COMMAND SEQUENCE """
i_max = 2000

total_timings = np.empty(i_max)
crack_timings = np.empty(i_max)

comp.command_init()                                 # INITIAL

for i in range(i_max):
    print('iteration {} / {}'.format(i+1, i_max), end='\r')

    total_time = getrusage(RUSAGE_SELF).ru_utime
    comp.command_equi()                             # EQUILIBRIUM
    comp.command_colli()                            # COLLISION
    comp.command_stre()                             # STREAM
    comp.command_update_bc()                        # UPDATE_BOUNDARY_CONDITION
    comp.command_boun()                             # APPLY_BOUNDARY_CONDITION
    comp.command_update_distribution_functions()    # UPDATE_DISTRIBUTION_FUNCTIONS
    comp.command_integrate()                        # INTEGRATE

    crack_time = getrusage(RUSAGE_SELF).ru_utime
    comp.command_dyn_crack()                        # DYNAMIC_CRACK
    crack_timings[i] = getrusage(RUSAGE_SELF).ru_utime - crack_time

    comp.command_time()                             # TIME
    total_timings[i] = getrusage(RUSAGE_SELF).ru_utime - total_time

print('finished ')

# analysis of timings
lbm_timings = total_timings - crack_timings
mean = {
    'total': total_timings.mean(),
    'lbm': lbm_timings.mean(),
    'crack': crack_timings.mean(),
}

median = {
    'total': np.median(total_timings),
    'lbm': np.median(lbm_timings),
    'crack': np.median(crack_timings),
}

summed = {
    'total': total_timings.sum(),
    'lbm': lbm_timings.sum(),
    'crack': crack_timings.sum(),
}

relative = {
    'lbm': summed['lbm'] / summed['total'],
    'crack': summed['crack'] / summed['total'],
}

# output of results
print('mean:', mean)
print('sum: ', summed)
print('rel: ', relative)
print('summary:')
print('mean: {:1.2e} s, std: {:1.2e} s, min: {:1.4e} s, max: {:1.4e} s'.format(
    mean['crack'], np.std(crack_timings),
    crack_timings.min(), crack_timings.max()))
print('statistics:')
print('time spent on crack: {:.3f} s of {:.3f} s; {:.2%} of time'.format(
    summed['crack'], summed['total'], relative['crack']))


