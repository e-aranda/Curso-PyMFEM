#!/usr/bin/env python
# coding: utf-8


import nlopt
import mfem.par as mfem
import numpy as np
from mpi4py import MPI

from mfem.common.arg_parser import ArgParser

# Paso de argumentos en la llamada
parser = ArgParser()

parser.add_argument('-v', '--visualization',
                    action='store_true',
                    help='Activar visualización GLVis')

args = parser.parse_args()
visual = args.visualization

comm = MPI.COMM_WORLD
num_procs = comm.size
myid = comm.rank
inicio = MPI.Wtime()


class Area(mfem.PyCoefficientBase):
    def __init__(self, u):
        super(Area, self).__init__(0)
        self.u = u
        self.grad = mfem.Vector()

    def Eval(self, T, ip):
        self.u.GetGradient(T, self.grad)
        sig = 1 + self.grad*self.grad
        return np.sqrt(sig)


class InverArea(mfem.VectorPyCoefficientBase):
    def __init__(self, dim, u):
        super(InverArea, self).__init__(dim, 0)
        self.u = u

    def Eval(self, elvect, T, ip):
        self.u.GetGradient(T, elvect)
        sig = 1 + elvect*elvect
        elvect *= 1./np.sqrt(sig)


class Initial(mfem.PyCoefficient):
    def EvalValue(self, x):
        return np.cos(np.sin(x[0]*x[1]))


class Boundary(mfem.PyCoefficient):
    def EvalValue(self, x):
        return 0.
#       return np.cos(np.sin(x[0]*x[1]))


class Surface:
    def __init__(self, u_, ess_, apar_, disp_, sout_):
        self.ess = ess_
        self.apar = apar_
        self.disp = disp_
        self.mu_gf = u_
        self.sout = sout_
        
        f_ = u_.ParFESpace()
        self.mesh = f_.GetMesh()
        dim = self.mesh.Dimension()
        self.iter = 0
        self.zero = mfem.ConstantCoefficient(0.0)
        self.vol = mfem.ParLinearForm(f_)
        pp_coeff = Area(self.mu_gf)
        self.vol.AddDomainIntegrator(mfem.DomainLFIntegrator(pp_coeff))
        ipp_coeff = InverArea(dim, self.mu_gf)
        self.a = mfem.ParLinearForm(f_)
        self.a.AddDomainIntegrator(mfem.DomainLFGradIntegrator(ipp_coeff))
        self.X = mfem.ParGridFunction(f_)
        self.iter = 0

    def myfunc(self, x, grad):
        flag = 0
        if x is None and grad is None: # cuando llamamos para parar
            flag = 1
            x = np.empty(self.apar.sum())
        # distribuimos x entre los procesos
        parx0 = np.empty(self.apar[myid])
        comm.Scatterv([x, self.apar, self.disp, MPI.DOUBLE], parx0, 0)
        # el proceso !0 pregunta si hay mensaje
        if myid:
            flag = comm.Iprobe(source=0, tag=MPI.ANY_TAG,
                               status=None)
        if flag:
            return None

        # Asignamos valores de x a mu_gf y calculamos la integral
        self.mu_gf.SetFromTrueDofs(mfem.Vector(parx0))
        self.vol.Assemble()
        C = self.vol.Sum()
        C_val = comm.allreduce(C, op=MPI.SUM)

        # construimos el gradiente 
        if grad.size > 0 or myid:
            self.a.Assemble()
            self.X.SetVector(self.a, 0)
            self.X.ProjectBdrCoefficient(self.zero, self.ess)
            pX = self.X.GetTrueDofs()
            parg = pX.GetDataArray()
            # agrupamos gradiente
            comm.Gatherv(parg, [grad, self.apar, self.disp, MPI.DOUBLE], 0)

        self.iter += 1
        if not myid:
            print("Iter: {0} -- {1}".format(self.iter, C_val))

        if visual and not self.iter % 10:
            self.sout << "parallel " << num_procs << " " << myid << "\n"
            self.sout << "solution\n" << self.mesh << self.mu_gf

        return C_val


# Definición de malla y espacio de elementos finitos

#  malla
serial_mesh = mfem.Mesh.MakeCartesian2D(100,100,mfem.Geometry.SQUARE,sx = np.pi, sy = np.pi)

mesh = mfem.ParMesh(comm, serial_mesh)

ess = mfem.intArray([1, 1, 1, 1])

# espacio de elementos finitos
fec = mfem.H1_FECollection(1,  mesh.Dimension())
fespace = mfem.ParFiniteElementSpace(mesh, fec)

u = mfem.ParGridFunction(fespace)
local_size = np.array([fespace.GetTrueVSize()], int)

# Este array permitirá distribuir los vectores entre los procesos
# MPI_Allgather: recolectar el tamaño de todos los procesos
# array_par_n contendrá el tamaño de los datos de cada proceso
array_par_n = np.empty(num_procs, int)
comm.Allgather(local_size, array_par_n)

# El acumulado por procesador
disp = np.zeros(num_procs, dtype='i')
disp[1:] = np.cumsum(array_par_n[:-1])

# Dato inicial y condición frontera
coeff = Initial()
boundary = Boundary()
u.ProjectCoefficient(coeff)
u.ProjectBdrCoefficient(boundary, ess)

# Obtenemos el vector de dofs reales y lo pasamos a un array de numpy
x = u.GetTrueDofs()
xx = x.GetDataArray()

sout = mfem.socketstream("localhost", 19916)
if visual:    
    sout.precision(8)
    sout << "parallel " << num_procs << " " << myid << "\n"
    sout << "solution\n" << mesh << u
    sout.send_text("valuerange 0 1\n")
    sout.send_text("autoscale off\n")
    sout.send_text("keys c\n")

# el proceso 0 es el que hará la optimización
if not myid:
    n = fespace.GlobalTrueVSize()
    iter = 200
    z = np.empty(n)
else:
    z = None

# distribuimos el dato inicial xx en z
comm.Gatherv(xx, [z, array_par_n, disp, MPI.DOUBLE], root=0)
obj = Surface(u, ess, array_par_n, disp, sout)

if not myid:
    opt = nlopt.opt(nlopt.LD_MMA, n)
    opt.set_min_objective(obj.myfunc)
    opt.set_xtol_rel(1.e-8)
    opt.set_maxeval(iter)
    opt.set_exceptions_enabled(True)
    opt.get_exceptions_enabled()
    try:
        x0 = opt.optimize(z)
        print("Optimal value", opt.last_optimum_value())
        print("\nTiempo: ", MPI.Wtime()-inicio)
    except Exception as e:
        print("Exception found:", e)
    # finalizada la optimización enviamos mensaje a los procesos para que paren
    i = 1
    for p in range(1,num_procs):
        comm.send(i, dest=p)
    # llamamos a la función para que llegue el mensaje        
    obj.myfunc(None, None)
else:
    # el resto de procesos ejecuta la función hasta que se devuelva None
    while True:
        res = obj.myfunc(None, np.array([]))
        if res is None:
            x0 = np.array([])
            break

# recuperamos la solución del proceso 0 al resto
parx0 = np.empty(array_par_n[myid])
comm.Scatterv([x0, array_par_n, disp, MPI.DOUBLE], parx0, 0)
u.SetFromTrueDofs(mfem.Vector(parx0))

if visual:
    sout << "parallel " << num_procs << " " << myid << "\n"
    sout << "solution\n" << mesh << u

MPI.Finalize()
