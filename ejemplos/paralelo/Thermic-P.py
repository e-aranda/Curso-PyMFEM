#!/usr/bin/env python
# coding: utf-8

import mfem.par as mfem
from mpi4py import MPI
import numpy as np

inicio = MPI.Wtime()

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

# mallado
serial_mesh = mfem.Mesh('mallas/termico.mesh')

mesh = mfem.ParMesh(MPI.COMM_WORLD,serial_mesh)
del serial_mesh

fec = mfem.H1_FECollection(1,  mesh.Dimension())
fespace = mfem.ParFiniteElementSpace(mesh, fec)


# #### Etiquetas frontera

# para las condiciones Robin
ess_robin = mfem.intArray([1,0,1,0])

# para la condiciones Dirichlet
ess_dirich = mfem.intArray([0,1,0,1])

# arrays para definir todas las etiquetas frontera
ess_tdof_list = mfem.intArray()

# recopilamos etiquetas para pasarlas al solver
fespace.GetEssentialTrueDofs(ess_dirich, ess_tdof_list)


# Coeficientes
alfa = 0.25
alfa_coeff = mfem.ConstantCoefficient(alfa)

ue_v = 25.
ue_coeff = mfem.ConstantCoefficient(ue_v*alfa)

class u0Coeff(mfem.PyCoefficient):
    def EvalValue(self,x):
        return 10.+90.*x[0]/6.
u0 = u0Coeff()    

kappa = mfem.PWConstCoefficient(mfem.Vector([0.2,2]))

# Coeficiente constante 1/dt
dt = 0.1
delta_t = mfem.ConstantCoefficient(1./dt)

# Coeficiente uold variable: representará u_{n-1} en la F.V. (habrá que dividirlo por dt)
uold = mfem.ParGridFunction(fespace)
uold_coef = mfem.ProductCoefficient(delta_t,mfem.GridFunctionCoefficient(uold))

# ### Formulación variacional

# forma bilineal
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(kappa))
a.AddDomainIntegrator(mfem.MassIntegrator(delta_t))
a.AddBoundaryIntegrator(mfem.MassIntegrator(alfa_coeff),ess_robin)
a.Assemble()

# forma lineal
b = mfem.ParLinearForm(fespace)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(uold_coef))
b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(ue_coeff),ess_robin)

# #### Dato frontera Dirichlet
x = mfem.ParGridFunction(fespace)
x.ProjectBdrCoefficient(u0,ess_dirich)    


# Dato inicial
uold.ProjectCoefficient(u0)
b.Assemble()
# Tiempo inicial y final    
t=0.
T = 5

# visualización inicial
u_sock = mfem.socketstream("localhost", 19916)
u_sock << "parallel " << num_procs << " " << myid << "\n"
u_sock.precision(8)
u_sock << "solution\n" << mesh << uold
u_sock << "view 0 0\n" 
u_sock << "pause\n"

# Preparamos el sistema  el solver
A = mfem.HypreParMatrix()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list,x,b,A,X,B)

M = mfem.HypreSmoother()
cg = mfem.CGSolver(MPI.COMM_WORLD)
cg.SetRelTol(1e-12)
cg.SetMaxIter(2000)
cg.SetPrintLevel(0)
cg.SetPreconditioner(M)
cg.SetOperator(A)

# Bucle: Euler implícito
while (t<T):    
    b.Update()
    b.Assemble()
    a.FormLinearSystem(ess_tdof_list,x,b,A,X,B)
    cg.Mult(B,X)
    a.RecoverFEMSolution(X, b, x)
    
    uold.ProjectGridFunction(x)
    
    cad = '{0:1.2f}'.format(t)
    u_sock << "parallel " << num_procs << " " << myid << "\n"
    u_sock << "solution\n" << mesh << x
    u_sock << "plot_caption '" << cad << "'"

    t += dt    
