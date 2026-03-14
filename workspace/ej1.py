#!/usr/bin/env python
# coding: utf-8

import mfem.ser as mfem

mesh = mfem.Mesh.MakeCartesian2D(20,20,mfem.Geometry.SQUARE)

fec = mfem.H1_FECollection(1,  mesh.Dimension())
fespace = mfem.FiniteElementSpace(mesh, fec)
print('Número de incógnitas:',fespace.GetTrueVSize())

# forma bilineal
a = mfem.BilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator())
a.Assemble()

# Segundo miembro (f=1)
one = mfem.ConstantCoefficient(1.0)

b = mfem.LinearForm(fespace)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))
b.Assemble()

# ### Condiciones frontera
boundary_dofs = mfem.intArray()
# Extracción de los grados de libertad asociados a la frontera (al completo)
fespace.GetBoundaryTrueDofs(boundary_dofs)

# #### Definimos función para asignar el valor en la frontera
x = mfem.GridFunction(fespace)
x.Assign(0.0)

# ### Formulación del sistema
A = mfem.SparseMatrix()
B = mfem.Vector()
X = mfem.Vector()

a.FormLinearSystem(boundary_dofs, x, b, A, X, B)

mfem.CG(A, B, X, 0, 200, 1e-12, 0.0)

# Asignamos solución a la función grid
a.RecoverFEMSolution(X, b, x)

# visualización con GLVis
u_sock = mfem.socketstream("localhost", 19916)
u_sock.precision(8)
u_sock << "solution\n" << mesh << x

