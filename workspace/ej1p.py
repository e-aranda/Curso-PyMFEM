#!/usr/bin/env python
# coding: utf-8

import mfem.par as mfem
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

serial_mesh = mfem.Mesh.MakeCartesian2D(20,20,mfem.Geometry.SQUARE)
mesh = mfem.ParMesh(MPI.COMM_WORLD,serial_mesh)
del serial_mesh

fec = mfem.H1_FECollection(1,  mesh.Dimension())
fespace = mfem.ParFiniteElementSpace(mesh, fec)
print('Número de incógnitas: de ',myid, fespace.GetTrueVSize())

# forma bilineal
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator())
a.Assemble()

# Segundo miembro (f=1)
one = mfem.ConstantCoefficient(1.0)

b = mfem.ParLinearForm(fespace)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
b.Assemble()

# ### Condiciones frontera
boundary_dofs = mfem.intArray()
# Asignamos etiquetas al espacio de elementos finitos
fespace.GetBoundaryTrueDofs(boundary_dofs)


# #### Definimos función para asignar el valor en la frontera
x = mfem.ParGridFunction(fespace)
x.Assign(0.0)


# ### Formulación del sistema

A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(boundary_dofs, x, b, A, X, B)

# Precondicionador tipo Gauss-Seidel
M = mfem.HypreBoomerAMG(A)
M.SetPrintLevel(0)
cg = mfem.CGSolver(MPI.COMM_WORLD)
cg.SetRelTol(1e-12)
cg.SetMaxIter(2000)
cg.SetPrintLevel(0)
cg.SetPreconditioner(M)
cg.SetOperator(A)
cg.Mult(B, X)

# Asignamos solución a la función grid
a.RecoverFEMSolution(X, b, x)

# visualización con GLVis
u_sock = mfem.socketstream("localhost", 19916)
u_sock << "parallel " << num_procs << " " << myid << "\n"
u_sock.precision(8)
u_sock << "solution\n" << mesh << x


