#!/usr/bin/env python
# coding: utf-8

# Bilaplaciano desacoplado
import mfem.par as mfem
from mpi4py import MPI
import numpy as np
from mfem.common.arg_parser import ArgParser
# Paso de argumentos en la llamada
parser = ArgParser()
parser.add_argument('-r', '--refine',
                    action='store', default=0, type=int,
                    help="Número de refinamientos de la malla")
parser.add_argument('-v', '--visualization',
                    action='store_true',
                    help='Activar visualización GLVis')

args = parser.parse_args()
ref_levels = args.refine
visual = args.visualization

inicio = MPI.Wtime()

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

serial_mesh = mfem.Mesh(100)
dim = serial_mesh.Dimension()

for x in range(ref_levels):
    serial_mesh.UniformRefinement()


mesh = mfem.ParMesh(MPI.COMM_WORLD, serial_mesh)
del serial_mesh

order = 1
h2 = mfem.H1_FECollection(order, dim)
U_space = mfem.ParFiniteElementSpace(mesh, h2)


bd = mfem.intArray([1,1])
# Asignamos etiquetas al espacio de elementos 
boundary_dofs = mfem.intArray()
U_space.GetBoundaryTrueDofs(boundary_dofs)


ux = mfem.ParGridFunction(U_space)
vx = mfem.ParGridFunction(U_space)

f1 = mfem.PWConstCoefficient(mfem.Vector([0.,0.]))
f2 = mfem.PWConstCoefficient(mfem.Vector([2*np.pi,-2.*np.pi]))
ux.ProjectBdrCoefficient(f1,bd)
vx.ProjectBdrCoefficient(f2,bd)


# Segundo miembro 
class ffunc(mfem.PyCoefficient):
    def EvalValue(self,x):
        return -np.pi**3*(np.pi*x[0]*np.sin(np.pi*x[0]) - 4*np.cos(np.pi*x[0]))
fcoef = ffunc()



# Forma lineal 
gform = mfem.ParLinearForm(U_space)
gform.AddDomainIntegrator(mfem.DomainLFIntegrator(fcoef))
gform.Assemble()

cVarf = mfem.ParBilinearForm(U_space)
cVarf.AddDomainIntegrator(mfem.DiffusionIntegrator())
cVarf.Assemble()


# ### Formulación del sistema para v

A = mfem.HypreParMatrix()
B = mfem.Vector()
V = mfem.Vector()
cVarf.FormLinearSystem(boundary_dofs, vx, gform, A, V, B)

# Precondicionador tipo Gauss-Seidel
M = mfem.HypreBoomerAMG(A)
M.SetPrintLevel(0)
cg = mfem.CGSolver(MPI.COMM_WORLD)
cg.SetRelTol(1e-12)
cg.SetMaxIter(2000)
cg.SetPrintLevel(0)
cg.SetPreconditioner(M)
cg.SetOperator(A)
cg.Mult(B, V)

# Asignamos solución a la función grid
cVarf.RecoverFEMSolution(V, gform, vx)

# Problema para u
vcoef = mfem.GridFunctionCoefficient(vx)

fform = mfem.ParLinearForm(U_space)
fform.AddDomainIntegrator(mfem.DomainLFIntegrator(mfem.ProductCoefficient(-1.,vcoef)))
fform.Assemble()

U = mfem.Vector()
cVarf.FormLinearSystem(boundary_dofs,ux,fform,A,U,B)
cg.Mult(B, U)


# Asignamos solución a la función grid
cVarf.RecoverFEMSolution(U, fform, ux)

# #### Cálculo del error
# Solución exacta
class usol(mfem.PyCoefficient):
    def EvalValue(self,x):
        return x[0]*np.sin(np.pi*x[0])
        
class vsol(mfem.PyCoefficient):
    def EvalValue(self,x):
        return np.pi*(-np.pi*x[0]*np.sin(np.pi*x[0]) + 2*np.cos(np.pi*x[0]))

uc = usol()
vc = vsol()

irs = [mfem.IntRules.Get(i,2) for i in range(mfem.Geometry.NumGeom)]

err_u = ux.ComputeL2Error(uc, irs)
err_v = vx.ComputeL2Error(vc, irs)
if not myid:
    print(err_u)
    print(err_v)

