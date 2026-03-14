#!/usr/bin/env python
# coding: utf-8


import mfem.par as mfem
from mpi4py import MPI
from mfem.common.arg_parser import ArgParser

# Paso de argumentos en la llamada
parser = ArgParser()
parser.add_argument('-r', '--refine',
                    action='store', default=0, type=int,
                    help="Número de refinamientos de la malla")
args = parser.parse_args()
ref_levels = args.refine

inicio = MPI.Wtime()

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

# malla serial
serial_mesh = mfem.Mesh.MakeCartesian2D(100,10,mfem.Geometry.TRIANGLE,False, 10.,1.)

for i in range(ref_levels):
    serial_mesh.UniformRefinement()

if not myid:
    print(f'Número de elementos: {serial_mesh.GetNE():4d}')
    print(f'Número de vértices: {serial_mesh.GetNV():4d}')

# malla paralela
mesh = mfem.ParMesh(MPI.COMM_WORLD,serial_mesh)

del serial_mesh
dim = mesh.Dimension()


# ### Espacio de elementos finitos vectorial
fec = mfem.H1_FECollection(2, dim)
fespace = mfem.ParFiniteElementSpace(mesh, fec, dim)


# Condición Dirichlet ${\bf u}={\bf 0}$ en la frontera 4
ess_tdof_list = mfem.intArray()
border = [0]*mesh.bdr_attributes.Max()
# Condición Dirichlet en frontera 4
border[3] = 1

# Extraemos los nodos del espacio asociados a la condición frontera en `ess_tdof_list`
ess_bdr = mfem.intArray(border)
fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

# Formulación Variacional
E = 21.e5
nu = 0.28
mu = 1/(2*(1.+nu))
lamb = nu/ ( (1 + nu)*(1-2*nu) )

Ecoeff = mfem.ConstantCoefficient(E)

a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.ElasticityIntegrator(Ecoeff,lamb,mu))
a.Assemble()

# fuerza RHS
f = mfem.VectorConstantCoefficient(mfem.Vector([0.,-1.]))

b = mfem.ParLinearForm(fespace)
b.AddDomainIntegrator(mfem.VectorBoundaryLFIntegrator(f))
b.Assemble()


# ### Resolución

# Función para valores en la frontera
u = mfem.ParGridFunction(fespace)
u.Assign(0.0)

# Preparamos resolución del sistema
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


# Asignamos solución a la función grid
a.RecoverFEMSolution(U, b, u)

# ### Visualización de resultados
deform_mesh = mfem.ParMesh(mesh)
deform_mesh.SetNodalFESpace(fespace)
nodes = deform_mesh.GetNodes()

# Copia separada de u
x = mfem.ParGridFunction(u) 
x *= 100
nodes += x


# Visualización
u_sock = mfem.socketstream("localhost", 19916)
u_sock << "parallel " << num_procs << " " << myid << "\n"
u_sock.precision(8)
u_sock << "solution\n" << deform_mesh << x
u_sock << "view 0 0\n"
u_sock << "zoom 2\n"
u_sock << "keys m\n"


# ### Extracción de componentes
fesp = mfem.ParFiniteElementSpace(mesh, fec) 
u1 = mfem.ParGridFunction(fesp,mfem.Vector(u.GetDataArray()))
u2 = mfem.ParGridFunction(fesp,mfem.Vector(u.GetDataArray()),u.Size()//2)

# # ##### Valores en los nodos
y1 = mfem.Vector()
y2 = mfem.Vector()
u.GetNodalValues(y1,1)
u.GetNodalValues(y2,1)


# # ### Interpolación en otros espacios
fec0 = mfem.H1_FECollection(1, dim)
fespa = mfem.ParFiniteElementSpace(mesh, fec0)

# Función de H1, con valores en los nodos de u1
z1 = mfem.ParGridFunction(fespa,y1)

# # Proyectamos u1 (de H2) en una función de H1, con `ProjectGridFunction`
w1 = mfem.ParGridFunction(fespa)
w1.ProjectGridFunction(u1) # proyecta los valores de u1 en w1

# Calculamos max(abs(w1-z1))
w1 -= z1
w1.Abs()
maxw1 = w1.Max()
maxw= MPI.COMM_WORLD.allreduce(maxw1, op=MPI.MAX)
if not myid:
    print(maxw)

# Con ComputeMaxError
zz1 = mfem.GridFunctionCoefficient(u1)
reso = z1.ComputeMaxError(zz1)
if not myid:
    print(reso)

    
# # También lo podemos hacer projectando desde P1 a P2:
z2 = mfem.ParGridFunction(fesp) # función en P2
w2 = mfem.ParGridFunction(fespa,y2) # define un función en P1 con los valores de y2
z2.ProjectGridFunction(w2) # proyecta los valores de w2 en z2 (de P1 a P2)

# Calculamos max(abs(w2-z2))
z2 -= u2
z2.Abs()
maxz2 = z2.Max()
maxz= MPI.COMM_WORLD.allreduce(maxz2, op=MPI.MAX)
if not myid:
    print(maxz)

if not myid:
    print("Tiempo: ",MPI.Wtime()-inicio)


