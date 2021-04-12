import SolidLBM.mesher.mesh.Boundary as Boundary_Module
<<<<<<< HEAD

boun1 = Boundary_Module.Boundary('/Users/alex/Work/LBM/PythonImplementierung/pyLBM/examples/circle/block.geo')
#boun1.read_input_file()
boun1.plot_boundary()

=======
import SolidLBM.mesher.mesh.Mesh as Mesh_Module

boun1 = Boundary_Module.Boundary('/Users/Sikang Yan/Documents/pyLBM/examples/circle/block.geo')
boun1.read_input_file()
boun1.plot_boundary()


# mesh size is here 0.025
#boun1.SetMeshableBoundary(0.025, '/Users/Sikang Yan/Documents/pyLBM/examples/circle/')
#boun1.correctBoundary(0.025, '/Users/Sikang Yan/Documents/pyLBM/examples/circle/')


#boun2 = Boundary_Module.Boundary('/Users/Sikang Yan/Documents/pyLBM/examples/circle/heatEquation.geo')
#boun2.plot_boundary()

>>>>>>> Sikang
