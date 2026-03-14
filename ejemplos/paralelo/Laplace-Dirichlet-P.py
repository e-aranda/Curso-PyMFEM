#!/usr/bin/env python
# coding: utf-8


import mfem.par as mfem
from mpi4py import MPI

inicio = MPI.Wtime()

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

# malla original
serial_mesh = mfem.Mesh('mallas/doblecirc.mesh')

if myid==0:
    print("Número de elementos:", serial_mesh.GetNE())
    print("Número de vértices:",serial_mesh.GetNV())

mesh= mfem.ParMesh(MPI.COMM_WORLD,serial_mesh)

dim = mesh.Dimension()
order = 1
fec = mfem.H1_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(mesh, fec)

print(f'Soy el proceso {myid:2d}, con número de incógnitas: {fespace.GetTrueVSize():8d}')
print(f'Soy el proceso {myid:2d}, con número de incógnitas global: {fespace.GlobalTrueVSize():8d}')


# forma bilineal
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator())
a.Assemble()

# Segundo miembro (f=1)
one = mfem.ConstantCoefficient(1.0)

b = mfem.ParLinearForm(fespace)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
b.Assemble()

# Condiciones frontera
boundary_dofs = mfem.intArray()
# Obtenemos etiquetas del espacio de elementos finitos: todas las fronteras son Dirichlet
fespace.GetBoundaryTrueDofs(boundary_dofs)

# Función para definir condiciones frontera
x = mfem.ParGridFunction(fespace)
x.Assign(0.)
inner = [0,1]
inner_border = mfem.intArray(inner)

coef = mfem.ConstantCoefficient(3.)
x.ProjectBdrCoefficient(coef,inner_border)

class Outer_Fun(mfem.PyCoefficient):
    def EvalValue(self,x):
        return x[0] + 2*x[1]        
u0 = Outer_Fun()

outer = [1,0]
outer_border = mfem.intArray(outer)
x.ProjectBdrCoefficient(u0,outer_border)

# Formulación del sistema
A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()

a.FormLinearSystem(boundary_dofs, x, b, A, X, B)

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
    print("Tiempo: ",MPI.Wtime()-inicio)

# Visualización
u_sock = mfem.socketstream("localhost", 19916)
u_sock << "parallel " << num_procs << " " << myid << "\n"
u_sock.precision(8)
u_sock << "solution\n" << mesh << x

