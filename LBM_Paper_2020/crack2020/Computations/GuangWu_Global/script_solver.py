import SolidLBM.solver.solution.Computation as Computation_Module

computation1 = Computation_Module.Computation(working_directory='./', name='disc_with_crack')
computation1.execute_command_sequence()
