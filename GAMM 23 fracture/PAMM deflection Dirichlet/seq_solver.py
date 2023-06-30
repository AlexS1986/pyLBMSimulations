import SolidLBM.solver.solution.Computation as Computation_Module


write_interval = 1

''' === INITIALIZATION === '''

name='deflection'
comp = Computation_Module.Computation(
    working_directory='./',
    name=name,
    verbose=True,
)
comp.command_init()
latt = comp.Lattice

''' === SIMULATION === '''

i = 0
while comp.current_time <= comp.Parameters['max_time']:
    print(f"    time {comp.current_time: >4f} / {comp.Parameters['max_time']: >4f},                iteration {i: >3d}", end='\r')

    comp.command_equi()
    comp.command_colli()
    comp.command_stre()
    comp.command_update_bc()
    comp.command_boun()
    comp.command_update_distribution_functions()
    comp.command_integrate()
    comp.command_dyn_crack()
    comp.command_time()
    if i % write_interval == 0:
        comp.command_vtk()
        comp.command_output_crack()
    
    i += 1
    
print('\n... done')

''' === POST-PROCESSING === '''

import crackeval.combi.crackpath as cpath
import crackeval.combi.jvec as jvec
import crackeval.utils as cutils


fig_c = cpath.plot(angle=-70)
cutils.save_plot(fig_c, 'cpath')

fig_j = jvec.plot_g(steps=True)
cutils.save_plot(fig_j, 'g-v-a')

