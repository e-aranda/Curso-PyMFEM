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


# función auxiliar g
def pFunc_exact(x):
    return np.exp(x[0])*np.sin(x[1])*np.cos(x[2])

# Funciones exactas (para comparación)
class uFunc_ex(mfem.VectorPyCoefficient): 
    def EvalValue(self, x):
        ret = [- np.exp(x[0])*np.sin(x[1])*np.cos(x[2]),
               - np.exp(x[0])*np.cos(x[1])*np.cos(x[2]),
               np.exp(x[0])*np.sin(x[1])*np.sin(x[2])]
        return ret

class pFunc_ex(mfem.PyCoefficient):
    def EvalValue(self, x):
        return pFunc_exact(x)

# Segundos miembros y condición frontera
class fFunc(mfem.VectorPyCoefficient):
    def EvalValue(self, x):
        return [0., 0., 0.]

class gFunc(mfem.PyCoefficient):
    def EvalValue(self, x):
        return -pFunc_exact(x)
        
class f_natural(mfem.PyCoefficient):
    def EvalValue(self, x):
        return -pFunc_exact(x)


### Malla serial y paralela
meshfile = 'mallas/fichera.mesh'
serial_mesh = mfem.Mesh(meshfile, 1, 1)
for x in range(ref_levels):
    serial_mesh.UniformRefinement()

mesh = mfem.ParMesh(MPI.COMM_WORLD, serial_mesh)
del serial_mesh

dim = mesh.Dimension()

# #### Espacios (Raviart-Thomas y L2)
order = 1
hdiv_coll = mfem.RT_FECollection(order, dim)
l2_coll = mfem.L2_FECollection(order, dim)

R_space = mfem.ParFiniteElementSpace(mesh, hdiv_coll)
W_space = mfem.ParFiniteElementSpace(mesh, l2_coll)

# Dimensiones de bloques
# Hay que distinguir entre GetVSize (grados de libertad en cada procesador (están repetidos)
# de TrueVSize (grados de libertad sin repetir)
block_offsets = mfem.intArray([0, R_space.GetVSize(), W_space.GetVSize()])
block_offsets.PartialSum()
block_trueOffsets = mfem.intArray([0, R_space.TrueVSize(), W_space.TrueVSize()])
block_trueOffsets.PartialSum()

# Vectores por bloques
x = mfem.BlockVector(block_offsets)
rhs = mfem.BlockVector(block_offsets)
# Vectores de TrueDofS
trueX = mfem.BlockVector(block_trueOffsets)
trueRhs = mfem.BlockVector(block_trueOffsets)
trueRhs.Assign(0.0)


fcoeff = fFunc(dim)
fnatcoeff = f_natural()
gcoeff = gFunc()
# Segundos miembros 
# .ParallelAssemble tiene en cuenta solo los TrueDofs
fform = mfem.ParLinearForm()
fform.Update(R_space, rhs.GetBlock(0), 0)
fform.AddDomainIntegrator(mfem.VectorFEDomainLFIntegrator(fcoeff))
fform.AddBoundaryIntegrator(
    mfem.VectorFEBoundaryFluxLFIntegrator(fnatcoeff))
fform.Assemble()
fform.ParallelAssemble(trueRhs.GetBlock(0))

gform = mfem.ParLinearForm()
gform.Update(W_space, rhs.GetBlock(1), 0)
gform.AddDomainIntegrator(mfem.DomainLFIntegrator(gcoeff))
gform.Assemble()
gform.ParallelAssemble(trueRhs.GetBlock(1))

# Formas bilineales
mVarf = mfem.ParBilinearForm(R_space)
mVarf.AddDomainIntegrator(mfem.VectorFEMassIntegrator())
mVarf.Assemble()
mVarf.Finalize()

bVarf = mfem.ParMixedBilinearForm(R_space, W_space)
bVarf.AddDomainIntegrator(mfem.VectorFEDivergenceIntegrator(mfem.ConstantCoefficient(-1.)))
bVarf.Assemble()
bVarf.Finalize()

# ### Matriz del sistema

matriz = mfem.BlockOperator(block_trueOffsets)

M = mVarf.ParallelAssemble()
B = bVarf.ParallelAssemble()
Bt = mfem.TransposeOperator(B)

matriz.SetBlock(0, 0, M)
matriz.SetBlock(0, 1, Bt)
matriz.SetBlock(1, 0, B)

# Solver
solver = mfem.MINRESSolver(MPI.COMM_WORLD)
solver.SetAbsTol(1.e-10)
solver.SetRelTol(1.e-8)
solver.SetMaxIter(1000)
solver.SetOperator(matriz)
solver.SetPrintLevel(0)
trueX.Assign(0.0)
solver.Mult(trueRhs, trueX)

if not myid:
    if solver.GetConverged():
        print("MINRES convergió en " + str(solver.GetNumIterations()) +
              " iteraciones, con residuo " + "{:g}".format(solver.GetFinalNorm()))
    else:
        print("MINRES no ha convergido tras " + str(solver.GetNumIterations()) +
          " iteraciones. Residuo:" + "{:g}".format(solver.GetFinalNorm()))

fin = MPI.Wtime()-inicio
if not myid:
    print("Tiempo en resolver: ",fin)

# #### Visualización de resultados

u = mfem.ParGridFunction()
p = mfem.ParGridFunction()
u.MakeRef(R_space, x.GetBlock(0), 0)
p.MakeRef(W_space, x.GetBlock(1), 0)

# Recuperamos los valores de la solución (que se han resuelto sobre trueX)
u.Distribute(trueX.GetBlock(0))
p.Distribute(trueX.GetBlock(1))


# #### Comparación con solución exacta
ucoeff = uFunc_ex(dim)
pcoeff = pFunc_ex()
order_quad = max(2, 2*order+1)

irs = [mfem.IntRules.Get(i, order_quad)
       for i in range(mfem.Geometry.NumGeom)]

norm_p = mfem.ComputeLpNorm(2, pcoeff, mesh, irs)
norm_u = mfem.ComputeLpNorm(2, ucoeff, mesh, irs)
err_u = u.ComputeL2Error(ucoeff, irs)
err_p = p.ComputeL2Error(pcoeff, irs)

if not myid:
    print("|| u_h - u_ex || / || u_ex || = " + "{:g}".format(err_u / norm_u))
    print("|| p_h - p_ex || / || p_ex || = " + "{:g}".format(err_p / norm_p))


if not myid:
    print("\nTiempo en comparar: ",MPI.Wtime()-fin)


if visual:
    u_sock = mfem.socketstream("localhost", 19916)
    u_sock << "parallel " << num_procs << " " << myid << "\n"
    u_sock.precision(8)
    u_sock << "solution\n" << mesh << u << "window_title 'Velocidad'\n"
    MPI.COMM_WORLD.Barrier()
    p_sock = mfem.socketstream("localhost", 19916)
    p_sock << "parallel " << num_procs << " " << myid << "\n"
    p_sock.precision(8)
    p_sock << "solution\n" << mesh << p << "window_title 'Presión'\n"

