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


meshfile = 'mallas/beam-tet.mesh'

serial_mesh = mfem.Mesh(meshfile, 1, 1)

for i in range(ref_levels):
    serial_mesh.UniformRefinement()

if not myid:
    print(f'Número de elementos: {serial_mesh.GetNE():4d}')
    print(f'Número de vértices: {serial_mesh.GetNV():4d}')

mesh = mfem.ParMesh(MPI.COMM_WORLD, serial_mesh)

del serial_mesh
dim = mesh.Dimension()

# ### Espacio de elementos finitos vectorial

fec = mfem.H1_FECollection(2, dim)
fespace = mfem.ParFiniteElementSpace(mesh, fec, dim)


# ### Condiciones frontera Dirichlet en etiq 1

ess_tdof_list = mfem.intArray()
ess_bdr = mfem.intArray([1,0,0])
fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

# ### Fuerza aplicada
g = mfem.VectorConstantCoefficient(mfem.Vector([0., -1.e-2, 0.]))

# ### Formulación variacional 
bN = mfem.intArray([0,1,0])
b = mfem.ParLinearForm(fespace)
b.AddBoundaryIntegrator(mfem.VectorBoundaryLFIntegrator(g), bN)
b.Assemble()

lamb = mfem.Vector(mesh.attributes.Max())
lamb.Assign(1.0)
lamb[0] = lamb[1] * 50
lambda_c = mfem.PWConstCoefficient(lamb)
mu = mfem.Vector(mesh.attributes.Max())
mu.Assign(1.0)
mu[0] = mu[1] * 50
mu_c = mfem.PWConstCoefficient(mu)

a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_c, mu_c))
a.Assemble()


# ### Resolución
u = mfem.ParGridFunction(fespace)
u.Assign(0.0)

A = mfem.HypreParMatrix()
B = mfem.Vector()
U = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, u, b, A, U, B)


M = mfem.HypreBoomerAMG(A)
M.SetPrintLevel(0)

# Solver
cg = mfem.CGSolver(MPI.COMM_WORLD)
cg.SetRelTol(1e-12)
cg.SetMaxIter(2000)
cg.SetPrintLevel(0)
cg.SetPreconditioner(M)
cg.SetOperator(A)
cg.Mult(B, U)

a.RecoverFEMSolution(U, b, u)


# ### Visualización de resultados

deform_mesh = mfem.Mesh(mesh)
deform_mesh.SetNodalFESpace(fespace)
nodes = deform_mesh.GetNodes()
nodes += u

if visual:
    u_sock = mfem.socketstream("localhost", 19916)
    u_sock << "parallel " << num_procs << " " << myid << "\n"
    u_sock.precision(8)
    u_sock << "solution\n" << deform_mesh << u
    u_sock << "zoom 2\n"
    u_sock << "keys m\n"

# # ### Obtención de valores de la solución en puntos concretos
count , elem, ip  = mesh.FindPoints([[2,0,0],[4,0,0],[8,0,0]])

if not myid:
    print("\nDesplazamientos en [2,0,0], [4,0,0], [8,0,0]")

# Atención a los procesadores donde no se encuentra el punto (elem[j]<0)
for j in range(count):
    if not elem[j]<0: 
#        for i in range(dim):
#            print(u.GetValue(elem[j],ip[j],i),end=' ')
#        print()
#    else:
        print("",end='\r')
#for j in range(count):
#    print(myid, count, j, elem[j])

# # ### Evaluación de integrales
energy = B*U
tot_energy = MPI.COMM_WORLD.allreduce(energy, op = MPI.SUM)
print(myid, energy)
if not myid:
    print(tot_energy)

if not myid:
    print("\nTiempo: ",MPI.Wtime()-inicio)
