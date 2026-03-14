#!/usr/bin/env python
# coding: utf-8

import mfem.par as mfem
from mpi4py import MPI
import numpy as np
from scipy.special import erfc
from mfem.common.arg_parser import ArgParser
# Paso de argumentos en la llamada
parser = ArgParser()
parser.add_argument('-r', '--refine',
                    action='store', default=2, type=int,
                    help="Número de refinamientos de la malla")
parser.add_argument('-v', '--visualization',
                    action='store_true',
                    help='Activar visualización GLVis')
parser.add_argument('-vs', '--visual_step',
                    action='store', default=100, type=int,
                    help='Define cada cuántas iteraciones se muestra la visualización')
parser.add_argument('-o', '--ode_solver',
                    action='store', default=2, type=int,
                    help="""Selección del ODE solver:
                    ForwardEulerSolver (1),
                    RK2Solver (2)
                    RK3SSolver (3)
                    RK4Solver (4)
                    RK6Solver (6)
                    BackwardEulerSolver (11)""")

vis_step = 100 # parámetro para visualización


args = parser.parse_args()
ref_levels = args.refine
visual = args.visualization
vis_step = args.visual_step
ode_solver_type = args.ode_solver

inicio = MPI.Wtime()
comm = MPI.COMM_WORLD
num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank

# Clase para el problema de evolución
class FE_Evolution(mfem.PyTimeDependentOperator):
    def __init__(self, M, K):
        mfem.PyTimeDependentOperator.__init__(self, M.Height())

        self.K = K
        self.M = M
        self.T = None
        
        self.z = mfem.Vector(M.Height())

        self.M_prec = mfem.HypreSmoother()
        self.M_prec.SetType(mfem.HypreSmoother.Jacobi)
        self.M_solver = mfem.CGSolver(comm)
        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.SetOperator(self.M)
        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(1.e-8)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(500)
        self.M_solver.SetPrintLevel(0)

        # Solver implícito
        self.T_prec = mfem.HypreSmoother()
        self.T_solver = mfem.CGSolver(comm)
        self.T_solver.iterative_mode = False
        self.T_solver.SetRelTol(1.e-8)
        self.T_solver.SetAbsTol(0.0)
        self.T_solver.SetMaxIter(500)
        self.T_solver.SetPrintLevel(0)
        self.T_solver.SetPreconditioner(self.T_prec)        
        self.Tmat = None
        self.Tmat_e = None


    def Mult(self, x, y):
        self.K.Mult(x, self.z)
        self.z.Neg()
        self.M_solver.Mult(self.z, y)

    def ImplicitSolve(self, dt, x, y):
        if self.T is None:
            self.T = mfem.Add(1.0, self.M, dt, self.K)
            self.T_solver.SetOperator(self.T)
        
        self.K.Mult(x, self.z)
        self.z.Neg()
        self.T_solver.Mult(self.z, y)       
    
# malla serie
meshfile = 'mallas/periodic-hexagon.mesh'

mesh = mfem.Mesh(meshfile)
dim = mesh.Dimension()

for lev in range(ref_levels):
    mesh.UniformRefinement()


# #### Velocidad y dato inicial
bb_min, bb_max = mesh.GetBoundingBox()
center = (bb_min + bb_max)/2.0
        
class velocity_coeff(mfem.VectorPyCoefficient):
    def EvalValue(self, x):
        X = 2 * (x - center) / (bb_max - bb_min)
        v = [np.pi/2*X[1],  - np.pi/2*X[0]]
        return v

class u0_coeff(mfem.PyCoefficient):
    def EvalValue(self, x):
        X = 2 * (x - center) / (bb_max - bb_min)
        rx = 0.45
        ry = 0.25
        cx = 0.
        cy = -0.2
        w = 10.
        return (erfc(w * (X[0]-rx)) * erfc(-w*(X[0]+rx)) *
                erfc(w * (X[1]-cy-ry)) * erfc(-w*(X[1]-cy+ry)))/16.

velocity = velocity_coeff(dim)
u0 = u0_coeff()

# Malla paralela
pmesh = mfem.ParMesh(comm, mesh)
del mesh

# Solver para la ecuación en tiempo
if ode_solver_type == 1:
    ode_solver = mfem.ForwardEulerSolver()
elif ode_solver_type == 2:
    ode_solver = mfem.RK2Solver(1.0)
elif ode_solver_type == 3:
    ode_solver = mfem.RK3SSolver()
elif ode_solver_type == 4:
    ode_solver = mfem.RK4Solver()
elif ode_solver_type == 6:
    ode_solver = mfem.RK6Solver()
elif ode_solver_type == 11:
    ode_solver = mfem.BackwardEulerSolver()

# Elementos finitos discontinuos 
order = 3
fec = mfem.L2_FECollection(order, dim, mfem.BasisType.GaussLobatto)
fes = mfem.ParFiniteElementSpace(pmesh, fec)

# Formulación Variacional
m = mfem.ParBilinearForm(fes)
m.AddDomainIntegrator(mfem.MassIntegrator())
k = mfem.ParBilinearForm(fes)
k.AddDomainIntegrator(mfem.ConvectionIntegrator(velocity, 1.0))
k.AddInteriorFaceIntegrator(
    mfem.NonconservativeDGTraceIntegrator(velocity, 1.0))
k.AddBdrFaceIntegrator(
        mfem.NonconservativeDGTraceIntegrator(velocity, 1.0))

m.Assemble()
m.Finalize()
k.Assemble()
k.Finalize()

# Matrices paralelas para pasar al operador TimeDependent
M = m.ParallelAssemble()
K = k.ParallelAssemble()

# Condición inicial
u = mfem.ParGridFunction(fes)
u.ProjectCoefficient(u0)

# Socket de visualización
if visual:
    sol_sock = mfem.socketstream("localhost", 19916)
    sol_sock << "parallel " << num_procs << " " << myid << "\n"
    sol_sock.precision(8)
    sol_sock.send_solution(pmesh,  u)

# creamos el TimeDependentOperator
adv = FE_Evolution(M, K)
ode_solver.Init(adv)


# Parámetros de tiempo
t_final = 5.
dt = 0.001
t = 0.0
ti = 0 # para visualización

while True:
    if t > t_final - dt/2:
        break
    t, dt = ode_solver.Step(u, t, dt)
    ti +=1

    if ti % vis_step == 0:
        cad = "Time: {0:2.2f}".format(t)
        if not myid:
            print("time step:",ti,"time:",np.round(t, 3))
        if visual:
            sol_sock << "parallel " << num_procs << " " << myid << "\n"
            sol_sock.send_solution(pmesh,  u)
            sol_sock.send_text("plot_caption '{0:s}'".format(cad))


if not myid:
    print("Tiempo:",MPI.Wtime()-inicio)
