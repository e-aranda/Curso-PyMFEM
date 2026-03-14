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

# Inicialización de MPI
comm = MPI.COMM_WORLD
num_procs = comm.Get_size()
myid = comm.Get_rank()

inicio = MPI.Wtime()


# Clase para el Operador dependiente del tiempo
class ConductionOperator(mfem.PyTimeDependentOperator):
    def __init__(self, ess_, global_, M, K, b_vector):
        mfem.PyTimeDependentOperator.__init__(self, fespace.GetTrueVSize())

        self.ess_tdof_list = ess_
        self.global_tdof_size = global_
        self.M = M
        self.K = K
        self.neumann_vector = b_vector
        
        self.current_dt = -1.0
        self.z = mfem.Vector(fespace.GetTrueVSize())
        self.Mmat = self.M.ParallelAssemble()
        # En caso de que hay condición Dirichlet, debemos eliminarla de la matriz Mmat
        self.Mmat_e = mfem.HypreParMatrix()
        # Mmat_e guarda los elementos eliminados de las columnas de Mmat (para recuperar luego en rhs)
        if self.global_tdof_size > 0:
            self.Mmat_e = self.M.ParallelEliminateTDofs(self.ess_tdof_list, self.Mmat)

        
        self.Kmat = self.K.ParallelAssemble()
        # Full_Kmat es una copia (diferente si hay cond. Dirichlet (es necesario usar ambas)
        # Si no hay condición Dirichelt, no es necesaria
        self.Kmat_e = mfem.HypreParMatrix()
        if self.global_tdof_size > 0:
            self.Full_Kmat = mfem.HypreParMatrix(self.Kmat) # deep copy
            self.Kmat_e = self.K.ParallelEliminateTDofs(self.ess_tdof_list, self.Kmat)
        else:
            self.Full_Kmat = self.Kmat

        
        # Solver explícito
        self.M_prec = mfem.HypreSmoother()
        self.M_prec.SetType(mfem.HypreSmoother.Jacobi)
        self.M_solver = mfem.CGSolver(fespace.GetComm())
        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.SetOperator(self.Mmat)
        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(1.e-8)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(500)
        self.M_solver.SetPrintLevel(0)

        # Solver implícito
        self.T_prec = mfem.HypreSmoother()
        self.T_solver = mfem.CGSolver(fespace.GetComm())
        self.T_solver.iterative_mode = False
        self.T_solver.SetRelTol(1.e-8)
        self.T_solver.SetAbsTol(0.0)
        self.T_solver.SetMaxIter(500)
        self.T_solver.SetPrintLevel(0)
        self.T_solver.SetPreconditioner(self.T_prec)        
        self.Tmat = None
        self.Tmat_e = None

        
    def Mult(self, u, du_dt):
        # du_dt = M^{-1} * (-K * u + b)
        self.Full_Kmat.Mult(u, self.z)
        self.z.Neg()
        self.z += self.neumann_vector
        
        if self.global_tdof_size > 0:
            du_dt.SetSubVector(self.ess_tdof_list, 0.) # ponemos a cero, pues es la derivada
            # Ajustamos rhs
            mfem.EliminateBC(self.Mmat, self.Mmat_e, self.ess_tdof_list, du_dt, self.z)
            
        self.M_solver.Mult(self.z, du_dt)

    def ImplicitSolve(self, dt, u, du_dt):
        # Resuelve: du_dt = M^{-1} * [-K(u + dt*du_dt) + b]
        if self.Tmat is None or dt != self.current_dt:
            # T = M + dt * K
            self.Tmat = mfem.Add(1.0, self.Mmat, dt, self.Kmat)
            if self.global_tdof_size > 0:
                self.Tmat_e = mfem.Add(1.0, self.Mmat_e, dt, self.Kmat_e)
            self.current_dt = dt
            self.T_solver.SetOperator(self.Tmat)
            
        self.Full_Kmat.Mult(u, self.z)
        self.z.Neg()
        self.z += self.neumann_vector
        
        if self.global_tdof_size > 0:
            du_dt.SetSubVector(self.ess_tdof_list, 0.0) # ponemos a cero, pues es la derivada 
            # Eliminate essential BC specified by ess_dof_list from the solution du_dt to the r.h.s. z
            mfem.EliminateBC(self.Tmat, self.Tmat_e, self.ess_tdof_list, du_dt, self.z)
            
        self.T_solver.Mult(self.z, du_dt)






        
# Lectura de malla y malla paralela
mesh = mfem.Mesh('mallas/termico.mesh',1 ,1)

for i in range(ref_levels):
    mesh.UniformRefinement()

pmesh = mfem.ParMesh(comm, mesh)
del mesh

# Espacio de elementos finitos
order = 2
fec = mfem.H1_FECollection(order, pmesh.Dimension())
fespace = mfem.ParFiniteElementSpace(pmesh, fec)


# Construcción del problema


# Datos y coeficientes del problema
dt = 0.1
T = 5.
alpha_v = 0.25
alpha_coeff = mfem.ConstantCoefficient(alpha_v)
ue_v = 25.
ue_coeff = mfem.ConstantCoefficient(ue_v*alpha_v)

# Dato frontera e inicial
class u0Coeff(mfem.PyCoefficient):
    def EvalValue(self,x):
        return 10.+90.*x[0]/6.
u0 = u0Coeff()    

# Coeficiente difusión
kappa = mfem.PWConstCoefficient(mfem.Vector([0.2,2]))






# Atributos de frontera y lista de DoFs esenciales
ess_tdof_list = mfem.intArray()

ess_robin = mfem.intArray([1,0,1,0])
ess_dirich = mfem.intArray([0,1,0,1])

fespace.GetEssentialTrueDofs(ess_dirich, ess_tdof_list)
# Número total de nodos en la frontera tipo Dirichelet
global_tdof_size = fespace.GetParMesh().ReduceInt(ess_tdof_list.Size());

# Inicialización
u_grid = mfem.ParGridFunction(fespace)
u_grid.ProjectCoefficient(u0)

# Es necesario fijar los valores en la frontera de los TrueDoF para imponer la condición Dirichlet
# Obtenemos los TrueDOFS de u_grid
u_true = mfem.Vector()
u_grid.GetTrueDofs(u_true)

# Valores frontera
u_bd = mfem.ParGridFunction(fespace)
u_bd.ProjectBdrCoefficient(u0, ess_dirich)
# extraemos a un vector
aux1 = mfem.Vector()
u_bd.GetTrueDofs(aux1)

# extraemos los valores de los nodos de la frontera y asignamos a u_true
aux = mfem.Vector()
aux1.GetSubVector(ess_tdof_list,aux)
u_true.SetSubVector(ess_tdof_list, aux)


# Formas variacionales
# Matriz de masa de la parte du/dt
M = mfem.ParBilinearForm(fespace)
M.AddDomainIntegrator(mfem.MassIntegrator())
M.Assemble()
M.Finalize(0)


# Matriz de rigidez
# Coefiente condición Robin (parte bilineal)
alpha = mfem.RestrictedCoefficient(alpha_coeff,ess_robin)

# Coeficiente condición Robin (parte lineal)
ue = mfem.RestrictedCoefficient(ue_coeff,ess_robin)
 
K = mfem.ParBilinearForm(fespace)
K.AddDomainIntegrator(mfem.DiffusionIntegrator(kappa))
K.AddBoundaryIntegrator(mfem.MassIntegrator(alpha))
K.Assemble(0)
K.Finalize()

# Segundo miembro
b = mfem.ParLinearForm(fespace)
b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(ue))
b.Assemble()
b_vector = mfem.Vector(fespace.TrueVSize())
b.ParallelAssemble(b_vector)

# Operador TimeDependent
oper = ConductionOperator(ess_tdof_list, global_tdof_size, M, K, b_vector)

# asignamos solución a toda la GridFunction    
u_grid.SetFromTrueDofs(u_true)

t = 0.0
dt = 0.1
t_final = 5.0

# ODE Solver: SDIRK23 (similar al ejemplo C++)
###    ode_solver = mfem.SDIRK23Solver()
ode_solver = mfem.BackwardEulerSolver()

ode_solver.Init(oper)

if visual:
    u_sock = mfem.socketstream("localhost", 19916)
    u_sock << "parallel " << num_procs << " " << myid << "\n"
    u_sock.precision(8)
    u_sock << "solution\n" << pmesh << u_grid
    u_sock << "view 0 0\n"
    u_sock << "pause\n"

while t < t_final:
    t, dt = ode_solver.Step(u_true, t, dt)

    u_grid.SetFromTrueDofs(u_true)
    if visual:
        u_sock << "parallel " << num_procs << " " << myid << "\n";
        u_sock << "solution\n" << pmesh << u_grid
        cad = '{0:1.2f}'.format(t)
        u_sock << "plot_caption '" << cad << "'"   
    

if not myid:
    print("\nTiempo: ",MPI.Wtime()-inicio)

