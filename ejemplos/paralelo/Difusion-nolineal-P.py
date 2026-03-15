#!/usr/bin/env python
# coding: utf-8

import mfem.par as mfem
from mpi4py import MPI
import numpy as np
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
meshfile = 'mallas/circulo.mesh'
serial_mesh = mfem.Mesh(meshfile)

for i in range(ref_levels):
    serial_mesh.UniformRefinement()

dim = serial_mesh.Dimension()

mesh = mfem.ParMesh(MPI.COMM_WORLD, serial_mesh)
del serial_mesh

# FE
order = 1
fec = mfem.H1_FECollection(order, dim)

fespace = mfem.ParFiniteElementSpace(mesh, fec)

n_incog = fespace.GlobalTrueVSize()
if not myid:
    print("Número de incógnitas:",n_incog)


# Coeficiente difusión no lineal k(u) y su derivada
def pFunc(x):
    return 1.+x*x

def dpFunc(x):
    return 2*x


# #### Coeficientes para el integrador no lineal

# Coeficiente k(u)\nabla u
class NonlinearCoefficient(mfem.VectorPyCoefficientBase):
    def __init__(self, dim, gf):
        super(NonlinearCoefficient, self).__init__(dim,0)
        self.gf = gf
    def Eval(self, elvect, T, ip):        
        self.gf.GetGradient(T,elvect)
        val = self.gf.GetValue(T,ip)
        elvect *= pFunc(val)

# Coeficiente k(u)
class NonlinearDerivativeCoefficient(mfem.PyCoefficientBase):
    def __init__(self, gf):
        super(NonlinearDerivativeCoefficient, self).__init__(0) # inicialización en la clase padre
        self.gf = gf
    def Eval(self, T, ip):
        val = self.gf.GetValue(T,ip)
        return pFunc(val)

# Coeficiente k'(u)\nabla u
class NonlinearDerivativeCoefficient2(mfem.VectorPyCoefficientBase):
    def __init__(self, dim, gf):
        super(NonlinearDerivativeCoefficient2, self).__init__(dim,0) # inicialización en la clase padre
        self.gf = gf
    def Eval(self, elvect, T, ip):
        self.gf.GetGradient(T, elvect)
        val = self.gf.GetValue(T,ip)
        elvect *= dpFunc(val)

# segundo miembro
class rhsCoef(mfem.PyCoefficient):
    def EvalValue(self, x):
         return 12*x[0]**4 + 24*x[0]**2*x[1]**2 - 16*x[0]**2 +12*x[1]**4 - 16*x[1]**2 +8;
     

brhs = rhsCoef()

# Integrador para la forma no linael y su derivada
class NonlinearDiffusionIntegrator(mfem.PyNonlinearFormIntegrator):
    def __init__(self,fes, d):
        super().__init__()
        self.fes = fes
        self.gf = mfem.GridFunction(fes)
        self.dim = d

    # Corresponde a k(u)\nabla u \nabla v
    def AssembleElementVector(self, el, T, elfun, elvect):
        dofs = self.fes.GetElementDofs(T.ElementNo)        
        ardofs = mfem.intArray(dofs)
        self.gf.SetSubVector(ardofs,elfun)
        coeff = NonlinearCoefficient(self.dim, self.gf)
        integ = mfem.DomainLFGradIntegrator(coeff)        
        integ.AssembleRHSElementVect(el, T, elvect)
        
        integ2 = mfem.DomainLFIntegrator(brhs)
        elvect2 = mfem.Vector()
        integ2.AssembleRHSElementVect(el, T, elvect2)
        elvect.Add(-1., elvect2)


    # Corresponde a k(u)\nabla w\nabla v + k'(u)\nabla u w \nabla v
    def AssembleElementGrad(self, el, T, elfun, elmat):
        dofs = self.fes.GetElementDofs(T.ElementNo)        
        ardofs = mfem.intArray(dofs)
        
        self.gf.SetSubVector(ardofs, elfun)
        coeff = NonlinearDerivativeCoefficient(self.gf)
        integ1 = mfem.DiffusionIntegrator(coeff)
        integ1.AssembleElementMatrix(el,T,elmat)

        coeff2 = NonlinearDerivativeCoefficient2(self.dim, self.gf)
        integ2 = mfem.MixedScalarWeakDivergenceIntegrator(coeff2)
        elmat2 = mfem.DenseMatrix()
        integ2.AssembleElementMatrix(el, T, elmat2)
        
        elmat.Add(1., elmat2)

# Asignación de valores en la frontera
x = mfem.ParGridFunction(fespace)
x.Assign(0.)

# 1. Crear la forma no lineal y añadir nuestro integrador
nform = mfem.ParNonlinearForm(fespace)
nonlin = NonlinearDiffusionIntegrator(fespace, dim)
nform.AddDomainIntegrator(nonlin)

b = mfem.ParLinearForm(fespace)
b.Assign(0.)

bdr = mfem.intArray([1]*mesh.bdr_attributes.Max())
nform.SetEssentialBC(bdr, b)


# Segundo miembro y vector solución
B  = b.ParallelAssemble()

X = mfem.Vector(fespace.GetTrueVSize())
x.ParallelProject(X)


# 2. Configurar el Solver Lineal (interno para el paso de Newton)
lsolver = mfem.CGSolver(MPI.COMM_WORLD)
lsolver.SetRelTol(1e-5)
lsolver.SetMaxIter(300)
lsolver.SetPrintLevel(0)

# Precondicionador
amg = mfem.HypreBoomerAMG()
amg.SetPrintLevel(0)
lsolver.SetPreconditioner(amg)

# 3. Configurar el NewtonSolver
newton_solver = mfem.NewtonSolver(MPI.COMM_WORLD)
newton_solver.SetOperator(nform)
newton_solver.SetSolver(lsolver)
newton_solver.SetMaxIter(100)
newton_solver.SetPrintLevel(1)
newton_solver.SetRelTol(1e-7)

newton_solver.Mult(B,X)    


# recuperamos la solución a la función grid
x.Distribute(X)


if not myid:
    print("Tiempo: ",MPI.Wtime()-inicio)


if visual:
    u_sock = mfem.socketstream("localhost", 19916)
    u_sock << "parallel " << num_procs << " " << myid << "\n"
    u_sock.precision(8)
    u_sock << "solution\n" << mesh << x 




