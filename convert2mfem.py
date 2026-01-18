#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 11:22:02 2018

@author: earanda

script to convert a GMF, Comsol or FreeFem++ format mesh in MFEM format mesh 

use: 
  $ python3 convert2mfem.py file.mesh [output.mesh]

"""

import os
import argparse
from itertools import islice

# check if linea is the first line of freefem format mesh
def check_msh(linea):
    try:
        l = linea.split();
        if len(l) == 3:
            for i in l:
                int(i)
        else:
            return False
    except:
        return False
    return True

class conversor(object):
    def __init__(self,fichero,salida=None):           
        self.archivo = fichero
        if salida:
            self.output = salida
        else:
            nombre = os.path.basename(fichero)
            name = os.path.splitext(nombre)
            self.output = name[0] + "-mfem.mesh"
        
    def convertidor(self,tipo,salida=None):
        if salida:
            self.salida = salida
        else:
            self.salida = self.output
        if tipo == "GMF":
            self.tratamientoGMF()
            self.escribirGMF()
            print("File created:",self.salida)
        elif tipo == "FreeFem++":
            self.tratamientoFF()
            self.escribirFF()
            print("File created:",self.salida)
        elif tipo == "Comsol":
            malla = self.tratamientoCOMSOL()
            self.escribirCOMSOL(malla)
            print("File created:",self.salida)
        elif tipo == "MFEM" or tipo == "GMSH":
            return 1
        else:
            print("Format not recognised")
            print(self.fullfile[0])
            return -1
        return 0
        
            
    def del_z(self,cadena):
        """
        deleting last component of vertices ()
        """
        lis = cadena.strip().split()
        lis.pop()
        return ' '.join(lis) + '\n'
        
     
    def del_zz(self,cadena):
        """
        deleting third coordinate of vertices (and the last component)
        """
        lis = cadena.strip().split()
        lis.pop()
        lis.pop()
        return ' '.join(lis) + '\n'
    
    def process(self,cadena,attr):
        """
        write the element
        """
        lis = cadena.strip().split()
        lista = [lis[-1],str(attr) + ' ']
        resto = map(str,map(lambda x: x-1 ,map(int,lis[:-1])))
        return ' '.join(lista) + ' '.join(resto) + '\n'


            
    def tipomalla(self):
        
        with open(self.archivo) as myFile:
            self.fullfile = myFile.readlines()
            
        if self.fullfile == None:
            print("Error reading the file. Aborting...")
            exit()
            
        # get the format type
        if "MeshVersionFormatted" in self.fullfile[0]:
            print("Format recognised: GMF")
            tipo = "GMF"
        elif "MFEM" in self.fullfile[0]:
            print("Format recognised: MFEM")
            tipo = "MFEM"
        elif "MeshFormat" in self.fullfile[0]:
            print("Format recognised: GMSH")
            tipo = "GMSH"
        elif check_msh(self.fullfile[0]):
            print("Format recognised: FreeFem++")
            tipo = "FreeFem++"
        elif "COMSOL" in self.fullfile[0]:
            print("Format recognised: Comsol")
            tipo = "Comsol"
        else:
            tipo = None
        return tipo
        
        
    def tratamientoGMF(self):        
        
        self.tri = False
        self.quad = False
        self.tetra = False
        self.hexa = False
        zero = first = second = third = -2
        for num, line in enumerate(self.fullfile):
            if "Tetrahedra" in line:
                zero = num
                self.tetra = True
            elif "Hexahedra" in line:
                zero = num       
                self.hexa = True
            elif "Vertices" in line:
                first = num
            elif "Triangles" in line:
                second = num
                self.tri = True
            elif "Quadrilaterals" in line:
                second = num
                self.quad = True
            elif "Edges" in line:
                third = num
                        
                
        if zero < 0:
            self.tetra = False
            self.hexa = False
            print("There is no volume elements.\n Assuming planar mesh: deleting 3rd component in vertices\n")
            self.planar = True
        else:
            self.planar = False
            self.n0 = zero + 1
            self.q = zero + 2
            self.fq = self.n0 + int(self.fullfile[self.n0].strip())
            
        if first < 0:
            print("There is no Vertices. Aborting...")
            exit(1)
        else:
            self.n1 = first + 1
            self.p = first + 2
            self.fp = self.n1 + int(self.fullfile[self.n1].strip())
            
        if second < 0:
            print("There is no Triangles or Quadrilaterals.")
            exit(1)
        elif self.tri and self.quad:
            print("Mesh with Triangles and Quadrilaterals. No contempled.")
            exit(1)
        else:
            self.n2 = second + 1
            self.s = second + 2
            self.fs = self.n2 + int(self.fullfile[self.n2].strip())
            
        if third < 0:
            print("There is no Edges.\n")
            self.edges = False
        else:
            self.edges = True
            self.n3 = third + 1
            self.t = third + 2
            self.ft = self.n3 + int(self.fullfile[self.n3].strip())


    def escribirGMF(self):
        """
        Writing the exit file in MFEM format
        """
        tratar = self.del_z
        boun = True
        dimv = '3'
        if self.tetra:
            dim = '3'
            el = self.n0
            inic = self.q
            finic = self.fq
            attr = 4  
            tr = self.n2
            inict = self.s
            finict = self.fs
            attrt = 2
        elif self.hexa:
            dim = '3'
            el = self.n0
            inic = self.q
            finic = self.fq
            attr = 5
            tr = self.n2
            inict = self.s
            finict = self.fs
            attrt = 3
        else:
            dim = '2'
            el = self.n2
            inic = self.s
            finic = self.fs
            if self.tri:
                attr = 2
            if self.quad:
                attr=3
            if self.planar:
                dimv = '2'
                tratar = self.del_zz
            if self.edges:
                tr = self.n3
                inict = self.t
                finict = self.ft
                attrt = 1
            else:
                boun = False
                if self.planar:
                    print("Planar mesh. No edges... Something wrong?")
            
                
        with open(self.salida,'w') as f:
            f.write('MFEM mesh v1.0\n\n')
            f.write('dimension\n' + dim + '\n\n')
            f.write('elements\n')
            f.write(self.fullfile[el].strip() + '\n')
            for line in self.fullfile[inic:finic+1]:
                f.write(self.process(line,attr))
            if boun:    
                f.write('\n\nboundary\n')
                f.write(self.fullfile[tr].strip() + '\n')
                for line in self.fullfile[inict:finict+1]:
                    f.write(self.process(line,attrt))
            else:
                f.write('\n\nboundary\n0\n\n')
            
                
            f.write('\n\nvertices\n')
            f.write(self.fullfile[self.n1].strip() + '\n')
            f.write(dimv+'\n')
            for line in self.fullfile[self.p:self.fp+1]:
                f.write(tratar(line))

                   
                
    def tratamientoFF(self):
                                
        self.nv,self.ne,self.nb = map(int,self.fullfile[0].strip().split())

        self.dim = len(self.fullfile[1].strip().split()) - 1
        if self.dim == 2:
            print("Converting a 2D mesh ...")
            self.eleme = '2'            
            self.borde = '1'
        elif self.dim == 3:
            print("Converting a 3D mesh ...")
            self.eleme = '4'
            self.borde = '2'
        else:
            print("Ooooops. Something wrong. Aborting ...")
            exit(1)
            
            
    def escribirFF(self):
        """
Writing the exit file in MFEM mesh format
"""

        cabecera = """\
MFEM mesh v1.0

#
# MFEM Geometry Types (see mesh/geom.hpp):
#
# POINT       = 0
# SEGMENT     = 1
# TRIANGLE    = 2
# SQUARE      = 3
# TETRAHEDRON = 4
# CUBE        = 5
#
        """
        
        with open(self.salida,'w') as f:
        
            f.write(cabecera+'\n')
            
            # dimension
            f.write("dimension\n")
            f.write(str(self.dim) + '\n')
            
            # elements
            f.write("\nelements\n")
            f.write(str(self.ne) + '\n')
            for x in self.fullfile[self.nv+1:self.nv+self.ne+1]:
                y = x.strip().split()
                zz = str(int(y[-1]) + 1)
                z = map(str,map(lambda x: x-1 ,map(int,y[0:self.dim+1])))
                f.write(zz + " " + self.eleme + " " + ' '.join(z) + '\n')
                
            # boundary elements
            f.write("\nboundary\n")
            f.write(str(self.nb) + '\n')
            for x in self.fullfile[self.nv+self.ne+1:self.nv+self.ne+self.nb+1]:
                y = x.strip().split()
                zz = str(int(y[-1])+1)
                z = map(str,map(lambda x: x-1 ,map(int,y[0:self.dim])))
                f.write(zz + " " + self.borde + " " + ' '.join(z) + '\n')

            # vertices
            f.write("\nvertices\n")
            f.write(str(self.nv) + '\n' + str(self.dim) + '\n')
            
            for x in self.fullfile[1:self.nv+1]:
                y = x.strip().split()
                f.write(' '.join(y[0:self.dim]) + '\n')
      
               
                
                
    def tratamientoCOMSOL(self):
        """
        Built a dictionary with all mesh information
        """
        malla = {}
        malla['puntos'] = []
        
        with open(self.archivo) as myFile:
            # Read dimension
            for line in myFile:
                if "# sdim" in line:
                    dim = int(line.split()[0])
                    break
            print("Dimension detected: {0:1d}".format(dim))
            malla['dimension'] = dim
            
            # Read mesh points
            myFile.seek(0)
            for num,line in enumerate(myFile):
                if "number of mesh" in line:
                    n = int(line.split()[0])
                if "Mesh" in line and "coordinates" in line:
                    nummesh = num
                    break                    

                
            myFile.seek(0)
            for line in islice(myFile,nummesh+1,nummesh+n+1):
                malla['puntos'].append(line[:-2])
                
            # Read valid elements
            cadenas = ["edg # type name","tri # type name","quad # type name",
                       "tet # type name","hex # type name"]
            names1 = ["edg","tri","quad","tet","hex"]
            names2 = ["edg_g","tri_g","quad_g","tet_g","hex_g"]
           
            for x,y,z in zip(cadenas,names1,names2):
                lin1,lin2 = self.reading_object(myFile,x)
                malla[y] = lin1
                malla[z] = lin2
    
            elements = []
            for x,y in zip(names1,names2):
                if malla[y]:
                    elements.append(x)
            sel = set(elements)        
            sets = set(names1)
            if not sel.issubset(sets):
                print("Mesh error: elements no recognised")
                exit(1)
            elem = list(sel)
            if (len(elem))>2:
                elem.remove('edg')
            if 'tet' in elem:
                malla['ele'] = 'tet'
                malla['bound'] = 'tri'
            elif 'hex' in elem:        
                malla['ele'] = 'hex'
                malla['bound'] = 'quad'
            elif 'tri' in elem:
                malla['ele'] = 'tri'
                malla['bound'] = 'edg'
            elif 'quad' in elem:
                malla['ele'] = 'quad'
                malla['bound'] = 'edg'

            
        return malla
        

    def tratarCOMSOL(self,x,tipo):
        cad = x.split()
        if tipo == 'edg':
            newcad = cad[-1] + ' ' + cad[0]
        elif tipo == 'quad':
            newcad = cad[0] + ' ' + cad[1] + ' ' + cad[3] + ' ' + cad[2]
        else:
            newcad = x
        return newcad


    def escribirCOMSOL(self,malla):
        """
Writing the exit file in MFEM mesh format
"""

        cabecera = """\
MFEM mesh v1.0

        """
        
        diccion = {'edg' : '1', 'tri' : '2', 'quad' : '3', 'tet' : '4',
                   'hex' : '5'}
        with open(self.salida,'w') as f:
        
            f.write(cabecera+'\n')
            
            # dimension
            f.write("dimension\n")
            f.write(str(malla['dimension']) + '\n')
            
            
            # elements
            
            f.write("\nelements\n")
            f.write(str(len(malla[malla['ele']])) + '\n')
            for x,y in zip(malla[malla['ele']],malla[malla['ele']+'_g']):
                f.write(str(int(y)+1) + ' ' + diccion[malla['ele']] + ' ' + self.tratarCOMSOL(x,malla['ele']) + '\n')
            

            f.write("\nboundary\n")
            f.write(str(len(malla[malla['bound']])) + '\n')
            for x,y in zip(malla[malla['bound']],malla[malla['bound']+'_g']):
                f.write(str(int(y)+1) + ' ' + diccion[malla['bound']] + ' ' + self.tratarCOMSOL(x,malla['bound'])  + '\n')

            # vertices
            f.write("\nvertices\n")
            f.write(str(len(malla['puntos'])) + '\n')
            f.write(str(malla['dimension']) + '\n')            
            for x in malla['puntos']:
                f.write(x + '\n')        



             
            
   

    def reading_object(self,myFile,cadena):
        """
        Read object from myFile and get information
        """
        lista_objeto = []
        lista_geom = []
        # read line where object start
        myFile.seek(0)
        for num,line in enumerate(myFile):
            if cadena in line:
                tri = num 
                break   
        else:
            return 0,0
        # read number of objects
        myFile.seek(0)
        num = 0
        for line in islice(myFile,tri+1,tri+10):
            num += 1
            if "number of nodes per element" in line or \
               "number of vertices per element" in line:
                num_inic = tri + 1 + num 
                break
        numtri = 0
        myFile.seek(0)
        for line in islice(myFile,num_inic,num_inic+1):
            numtri = int(line.split()[0])         
            
        # read elements
        myFile.seek(0)
        for line in islice(myFile,num_inic+2,num_inic+numtri+2):
            lista_objeto.append(line[:-2])
        # read geometry information of object
        
    
        myFile.seek(0)
        num = 0
        for line in islice(myFile,num_inic+numtri+3,num_inic+numtri+10):
            num += 1
            if "number of geometric entity" in line:
                break
        else:
            print("No geometric entity indices found")
            return 0,0
        # change numbering
        for line in islice(myFile,num,num+numtri):
            lista_geom.append(line[:-2])
        
        return lista_objeto,lista_geom
             

if __name__== "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('mesh', action="store")
    parser.add_argument('--output',action="store",
                    help='Name of output mesh file')

    args = parser.parse_args()
    
    nombre = os.path.basename(args.mesh)
    name = os.path.splitext(nombre)

    # check if output was set
    if not args.output:
        salida = name[0] + "-mfem.mesh"
    else:
        salida = args.output
        
    w = conversor(args.mesh,salida)
    tipo = w.tipomalla()
    w.convertidor(tipo)
    
#    w.escribir()
    
