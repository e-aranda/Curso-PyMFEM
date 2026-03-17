#!/usr/bin/env python
# coding: utf-8

# Problema  no lineal: -\Delta u + G(u) = f

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
meshfile = 'mallas/star.mesh'
serial_mesh = mfem.Mesh(meshfile)

for i in range(ref_levels):
    serial_mesh.UniformRefinement()

dim = serial_mesh.Dimension()

mesh = mfem.ParMesh(MPI.COMM_WORLD, serial_mesh)
del serial_mesh

# espacio de elementos finitos
order = 1
fec = mfem.H1_FECollection(order, dim)

fespace = mfem.ParFiniteElementSpace(mesh, fec)
# Número de incógnitas
n_incog = fespace.GlobalTrueVSize()
if not myid:
    print("Número de incógnitas", n_incog)

# función no lineal G(u) y derivada G'(u)
alpha = -0.5
def pFunc(x):
    return np.exp(alpha*x)

def dpFunc(x):
    return alpha*np.exp(alpha*x)

# Coeficientes para las formulaciones variacionales con G(u) y G'(u)    
class NonlinearCoefficient(mfem.PyCoefficientBase):
    def __init__(self, gf):
        super(NonlinearCoefficient, self).__init__(0) # inicialización en la clase padre
        self.gf = gf
    def Eval(self, T, ip):
        val = self.gf.GetValue(T,ip)
        return pFunc(val)
        
class NonlinearDerivativeCoefficient(mfem.PyCoefficientBase):
    def __init__(self, gf):
        super(NonlinearDerivativeCoefficient, self).__init__(0) # inicialización en la clase padre
        self.gf = gf
    def Eval(self, T, ip):
        val = self.gf.GetValue(T,ip)
        return dpFunc(val)

# Integrador para la forma no lineal: precisa método para la integral y su derivada

class NonlinearMassIntegrator(mfem.PyNonlinearFormIntegrator):
    def __init__(self,fes):
        super().__init__()
        self.fes = fes
        self.gf = mfem.GridFunction(fes)

    def AssembleElementVector(self, el, T, elfun, elvect):
        dofs = self.fes.GetElementDofs(T.ElementNo)        
        ardofs = mfem.intArray(dofs)
        self.gf.SetSubVector(ardofs,elfun)
        coeff = NonlinearCoefficient(self.gf)
        integ = mfem.DomainLFIntegrator(coeff)
        integ.AssembleRHSElementVect(el,T,elvect)

    def AssembleElementGrad(self, el, T, elfun, elmat):
        dofs = self.fes.GetElementDofs(T.ElementNo)
        ardofs = mfem.intArray(dofs)
        self.gf.SetSubVector(ardofs, elfun)
        coeff = NonlinearDerivativeCoefficient(self.gf)
        integ = mfem.MassIntegrator(coeff)
        integ.AssembleElementMatrix(el,T,elmat)

        
# Fijamos condición frontera
x = mfem.ParGridFunction(fespace)
x.Assign(1.)

# Forma No lineal
nform = mfem.ParNonlinearForm(fespace)
nonlin = NonlinearMassIntegrator(fespace)
nform.AddDomainIntegrator(nonlin)
nform.AddDomainIntegrator(mfem.DiffusionIntegrator())

# segundo miembro
b = mfem.ParLinearForm(fespace)
b.Assign(0.)

# Condiciones frontera 
bdat = mfem.intArray([1]*mesh.bdr_attributes.Max())
nform.SetEssentialBC(bdat,b)



# Segundo miembro
B  = b.ParallelAssemble()

# Vector solución (proyectamos los valores en la frontera)
X = mfem.Vector(fespace.GetTrueVSize())
x.ParallelProject(X)

# Solver Lineal (interno para Newton)
lsolver = mfem.CGSolver(MPI.COMM_WORLD)
# Precondicionador
#amg = mfem.HypreBoomerAMG();
#amg.SetPrintLevel(0);
#lsolver.SetPreconditioner(amg);

# NewtonSolver
newton_solver = mfem.NewtonSolver(MPI.COMM_WORLD)
newton_solver.SetOperator(nform)
newton_solver.SetSolver(lsolver)
newton_solver.SetPrintLevel(1)
newton_solver.SetMaxIter(100)
newton_solver.SetRelTol(1e-5)

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

