#!/usr/bin/env python
# coding: utf-8

# Problema de Stokes paralelo

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


serial_mesh = mfem.Mesh.MakeCartesian2D(30,30,mfem.Geometry.SQUARE)
dim = serial_mesh.Dimension()
for i in range(ref_levels):
    serial_mesh.UniformRefinement()

mesh= mfem.ParMesh(MPI.COMM_WORLD,serial_mesh)
del serial_mesh

order = 2
h2 = mfem.H1_FECollection(order, dim)
h1 = mfem.H1_FECollection(order-1, dim)

U_space = mfem.ParFiniteElementSpace(mesh, h2, dim)
P_space = mfem.ParFiniteElementSpace(mesh, h1)

# Tamaño total del problema

dimR = U_space.GlobalTrueVSize()
dimW = P_space.GlobalTrueVSize()

if not myid:
    print("Total de Variables: {0} + {1} = {2}".format(dimR,dimW,dimR+dimW))

# Vectores por Bloques
block_offsets = mfem.intArray([0, U_space.GetVSize(), P_space.GetVSize()])
block_offsets.PartialSum()
block_trueOffsets = mfem.intArray([0, U_space.TrueVSize(), P_space.TrueVSize()])
block_trueOffsets.PartialSum()


x = mfem.BlockVector(block_offsets)
rhs = mfem.BlockVector(block_offsets)

trueX = mfem.BlockVector(block_trueOffsets)
trueRhs = mfem.BlockVector(block_trueOffsets)
trueRhs.Assign(0.0)


# Todas las fronteras son Dirichlet
boundary_dofs = mfem.intArray([1,1,1,1])

# Condición Dirichlet no nula
ess_bdr = mfem.intArray([0,0,1,0])

x.Assign(0.)
ux = mfem.ParGridFunction(U_space)
ux.MakeRef(U_space,x.GetBlock(0),0)

f1 = mfem.VectorConstantCoefficient([-1.,0.])
ux.ProjectBdrCoefficient(f1,ess_bdr)

# Coeficiente segundo miembro
class fFunc(mfem.VectorPyCoefficient):
    def EvalValue(self, x):
        return [1 , 1.]
        
fcoeff = fFunc(dim)
# RHS
fform = mfem.ParLinearForm()
fform.Update(U_space, rhs.GetBlock(0), 0)
fform.AddDomainIntegrator(mfem.VectorDomainLFIntegrator(fcoeff))
fform.Assemble()
fform.ParallelAssemble(trueRhs.GetBlock(0))

gform = mfem.ParLinearForm()
gform.Update(P_space, rhs.GetBlock(1), 0)
gform.Assemble()
gform.ParallelAssemble(trueRhs.GetBlock(1))

# Formas bilineales
mVarf = mfem.ParBilinearForm(U_space)
mVarf.AddDomainIntegrator(mfem.VectorDiffusionIntegrator())
mVarf.Assemble()
mVarf.EliminateEssentialBC(boundary_dofs, ux, rhs.GetBlock(0))
mVarf.Finalize()

bVarf = mfem.ParMixedBilinearForm(U_space, P_space)
bVarf.AddDomainIntegrator(mfem.VectorDivergenceIntegrator(mfem.ConstantCoefficient(-1.)))
bVarf.Assemble()
bVarf.EliminateTrialEssentialBC(boundary_dofs, ux, rhs.GetBlock(1))
bVarf.Finalize()

M = mVarf.ParallelAssemble()
B = bVarf.ParallelAssemble()
Bt = mfem.TransposeOperator(B)


# Para pasar los valores frontera a trueX y trueRhs
ux.ParallelAssemble(trueX.GetBlock(0))

# Copiamos valores de rhs a r0
r0 = mfem.ParGridFunction(U_space)
r0.MakeRef(U_space,rhs.GetBlock(0),0)
# los pasamos a trueRhs
r0.ParallelAssemble(trueRhs.GetBlock(0))


# Matriz del sistema
stokesOp = mfem.BlockOperator(block_trueOffsets)

stokesOp.SetBlock(0, 0, M)
stokesOp.SetBlock(0, 1, Bt)
stokesOp.SetBlock(1, 0, B)


solver = mfem.MINRESSolver(MPI.COMM_WORLD)
solver.SetAbsTol(1.e-10)
solver.SetRelTol(1.e-4)
solver.SetMaxIter(5000)
solver.SetOperator(stokesOp)
solver.SetPrintLevel(0)
trueX.Assign(0.)
solver.Mult(trueRhs, trueX)

if not myid:
    if solver.GetConverged():
        print("MINRES convergió en " + str(solver.GetNumIterations()) +
              " iteraciones, con residuo " + "{:g}".format(solver.GetFinalNorm()))
    else:
        print("MINRES no ha convergido tras " + str(solver.GetNumIterations()) +
          " iteraciones. Residuo:" + "{:g}".format(solver.GetFinalNorm()))

if not myid:
    print("Tiempo: ",MPI.Wtime()-inicio)

# #### Visualización de resultados

if visual:
    u = mfem.ParGridFunction()
    p = mfem.ParGridFunction()
    u.MakeRef(U_space, x.GetBlock(0), 0)
    p.MakeRef(P_space, x.GetBlock(1), 0)
    u.Distribute(trueX.GetBlock(0))
    p.Distribute(trueX.GetBlock(1))

    u_sock << "parallel " << num_procs << " " << myid << "\n"
    u_sock << "solution\n" << mesh << u << "window_title 'Velocidad'\n"
    MPI.COMM_WORLD.Barrier()
    p_sock = mfem.socketstream("localhost", 19916)
    p_sock << "parallel " << num_procs << " " << myid << "\n"
    p_sock.precision(8)
    p_sock << "solution\n" << mesh << p << "window_title 'Presión'\n"


