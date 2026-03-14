#!/usr/bin/env python
# coding: utf-8

import mfem.par as mfem
from mpi4py import MPI
from mfem.common.arg_parser import ArgParser
# Paso de argumentos en la llamada
parser = ArgParser()
parser.add_argument('-r', '--refine',
                    action='store', default=2, type=int,
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

# mallado
serial_mesh = mfem.Mesh('mallas/conductor.mesh',1 ,1)

for i in range(ref_levels):
    serial_mesh.UniformRefinement()

mesh= mfem.ParMesh(MPI.COMM_WORLD,serial_mesh)

# Espacio de elementos finitos
fec = mfem.H1_FECollection(1,  mesh.Dimension())
fespace = mfem.ParFiniteElementSpace(mesh, fec)


# #### Coeficiente $\kappa$ definido mediante `PWConstCoefficient`
kappa = mfem.PWConstCoefficient(mfem.Vector([5.,0.,0.,0.,1.]))

# Formas bilineal y lineal
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(kappa))
a.Assemble()

# Segundo miembro (f(x)=0)
b = mfem.ParLinearForm(fespace)
b.Assign(0.)


# ### Condiciones en la frontera

ess_list = [0]*mesh.bdr_attributes.Max()
ess_list[1] = 20.
ess_list[99] = 60.
ess = mfem.Vector(ess_list)

border = mfem.PWConstCoefficient(ess)

ess_lt = [0]*mesh.bdr_attributes.Max()
ess_lt[1] = 1
ess_lt[99] = 1
essl = mfem.intArray(ess_lt)

x = mfem.ParGridFunction(fespace)
x.Assign(0.)
x.ProjectBdrCoefficient(border,essl)

ess_tdof_list = mfem.intArray()
fespace.GetEssentialTrueDofs(essl, ess_tdof_list)


# ### Formulación del sistema

A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
print(f"Proc: {myid:2d}. Tamaño del sistema: {A.Height():8d}")

# Resolución del sistema
# Precondicionador
M = mfem.HypreBoomerAMG(A)
M.SetPrintLevel(0)

# Solver
cg = mfem.CGSolver(MPI.COMM_WORLD)
cg.SetRelTol(1e-12)
cg.SetMaxIter(2000)
cg.SetPrintLevel(0)
cg.SetPreconditioner(M)
cg.SetOperator(A)
cg.Mult(B, X)

# Asignamos solución a la función grid
a.RecoverFEMSolution(X, b, x)


if not myid:
    print("Tiempo:",MPI.Wtime()-inicio)

# Visualización
if visual:
    u_sock = mfem.socketstream("localhost", 19916)
    u_sock << "parallel " << num_procs << " " << myid << "\n"
    u_sock.precision(8)
    u_sock << "solution\n" << mesh << x
