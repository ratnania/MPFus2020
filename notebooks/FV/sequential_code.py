#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:48:49 2020

@author: kissami
"""

from collections import OrderedDict 
import numpy as np
import meshio
import timeit

class Cells:
    nodeid = OrderedDict()
    faceid = OrderedDict()
    center = OrderedDict()
    volume = OrderedDict()
    

class Nodes:
    vertex = OrderedDict()
    bound = OrderedDict()


class Faces:
    nodeid = OrderedDict()
    cellid = OrderedDict()
    bound = OrderedDict()
    normal = OrderedDict()


#compute ghost value (neumann)
def ghost_value(w, Faces):
    w_ghost = OrderedDict()#[-1]*len(face_cellid)
    
    for i in range(len(Faces.cellid)):
        if Faces.bound[i] == 2:
            w_ghost[i] = w[Faces.cellid[i][0]]
        else :
            w_ghost[i] = w[Faces.cellid[i][0]]
    return w_ghost

#compute de normal vector
def VecteurNormal(a,b,bary):

    n = [None]*2
    s = [None]*2
    m = [None]*2
    normal = [None]*2
     
   
    n[0] = a[1] - b[1]
    n[1] = b[0] - a[0];
    
    m[0] = 0.5 * (a[0] + b[0]);
    m[1] = 0.5 * (a[1] + b[1]);
    
    s[0] = bary[0] - m[0] ;
    s[1] = bary[1] - m[1] ;
    
    if ( (s[0] * n[0] + s[1] * n[1])> 0):
        normal[0] = -1*n[0];
        normal[1] = -1*n[1];
    else:
        normal[0] = n[0];
        normal[1] = n[1];
        
    normal[0] = normal[0]#/longueur
    normal[1] = normal[1]#/longueur 
    
    
    return normal#, longueur

#create structure of cells, nodes and faces
def create_local_mesh(file):
    
    #read the gmsh file
    mesh = meshio.gmsh.read(file)

    #create the 3 nodes for each cell
    for i,j in mesh.cells.items():
        if i == "triangle":
            for k in range(len(j)):
                Cells.nodeid[k] = list(j[k])
                Cells.nodeid[k].sort()
     
    #create the node vertexes           
    Nodes.vertex = mesh.points
    #create the node boudaries
    Nodes.bound = [0]*len(Nodes.vertex)
    for i,j in mesh.cell_data.items():
        if i=="line":
            x = j.get('gmsh:physical')
  
    for i,j in mesh.cells.items():
        if i == "line":
            for k in range(len(j)):
                for l in range(2):
                    if x[k] > 2:
                        Nodes.bound[j[k][l]] = x[k]
    for i,j in mesh.cells.items():
         if i == "line":
             for k in range(len(j)):
                 for l in range(2):
                     if x[k] <= 2:
                         Nodes.bound[j[k][l]] = x[k]
        
    
    #compute de center and volume of cells      
    for i in range(len(Cells.nodeid)):
        s1 = Cells.nodeid[i][0]
        s2 = Cells.nodeid[i][1]
        s3 = Cells.nodeid[i][2]
        
        x1 = Nodes.vertex[s1][0]
        y1 = Nodes.vertex[s1][1]
        x2 = Nodes.vertex[s2][0]
        y2 = Nodes.vertex[s2][1]
        x3 = Nodes.vertex[s3][0]
        y3 = Nodes.vertex[s3][1]
        
        
        Cells.center[i] = (1./3 * (x1 + x2 + x3), 1./3*(y1 + y2 + y3))
        Cells.volume[i] = (1./2) * abs((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))

    #create faces node_id  
    Cellule = Cells.nodeid
    
    cellF = []
    faces = []
    k = 0
    for i in range(len(Cellule)):
        faces.append([Cellule[i][0], Cellule[i][1]])
        faces.append([Cellule[i][1], Cellule[i][2]])
        faces.append([Cellule[i][0], Cellule[i][2]])
        cellF.append([faces[k], faces[k+1], faces[k+2]])
        k = k+3
    
    faces =   set( tuple(x) for x in faces)
    faces = list(faces)
    
    
    
    facesdict = OrderedDict()
    for i in range(len(faces)):
        facesdict[faces[i]] = i
        Faces.nodeid[i] = faces[i]
        Faces.cellid[i] = (-1,-1)
            
     
    #create the 3 faces of each cell            
    for i in range(len(cellF)):
        Cells.faceid[i]  = [facesdict.get(tuple(cellF[i][0])), facesdict.get(tuple(cellF[i][1])), 
                            facesdict.get(tuple(cellF[i][2]))]
        
    #Faces.cellid = [(-1, -1)]*len(faces)
    for i in range(len(Cells.faceid)):
        for j in range(3):
            if Faces.cellid[Cells.faceid[i][j]] == (-1, -1):
                Faces.cellid[Cells.faceid[i][j]]= (i, -1)
            if Faces.cellid[Cells.faceid[i][j]][0] != i:
                Faces.cellid[Cells.faceid[i][j]]= (Faces.cellid[Cells.faceid[i][j]][0], i) 
                
    #boundary faces (1 : inlet ,2: outlet ,3 : upper wall, 4: down wall)
    for i in range(len(faces)):
        Faces.bound[i] = 0
        if (Faces.cellid[i][1] == -1):
            if (Nodes.bound[Faces.nodeid[i][0]] == Nodes.bound[Faces.nodeid[i][1]] ):
                Faces.bound[i] = Nodes.bound[Faces.nodeid[i][0]]
        Faces.normal[i] = VecteurNormal(Nodes.vertex[Faces.nodeid[i][0]], Nodes.vertex[Faces.nodeid[i][1]], 
                    Cells.center[Faces.cellid[i][0]])
                
    return Cells, Nodes, Faces
        
      
#compute the upwind flux  
def compute_flux(wl, wr,u, n):
        c = 0
        q = np.dot(u,n) 
        if (q >= 0):
            c = wl
        else:
            c = wr   
        flux = q*c
        return flux

#compute the rezidus using explicit scheme
def ExplicitScheme(w, u, w_ghost, Faces):
    
    rezidus  = [0]*len(w)
    for i in range(len(Faces.cellid)):
        wl = w[Faces.cellid[i][0]]
        n = Faces.normal[i]
        
        if (Faces.bound[i] == 0):
            wr = w[Faces.cellid[i][1]]
            flx = compute_flux(wl, wr,u[i], n)
            rezidus[Faces.cellid[i][0]] -= flx 
            rezidus[Faces.cellid[i][1]] += flx
      
        else:
            wr = w_ghost[i]
            flx = compute_flux(wl, wr,u[i],  n)
            rezidus[Faces.cellid[i][0]] -= flx
        
  
    return rezidus    

#save results in vtk format
def save_paraview_results(w, n, dt, cells, nodes):
    elements = {"triangle": np.array(list(cells))}
    points = [[-1, -1, -1] for i in range(len(nodes))]
    for i in range(len(nodes)):
        for j in range(3):
            points[i][j]  = nodes[i][j]
        
    data = {"wcell" : np.array(w)}  
    data = {"wcell": data}
    maxw = np.zeros(1)
    maxw = max(w)
    
    
    if(n%50 == 0):
        print("saving paraview results, iteration number ", n)
        print("max w =", maxw, "dt =", dt)
            
        meshio.write_points_cells("results/visu-"+str(n)+".vtu", 
                                  points, elements, cell_data=data)
    


#Début du code
start = timeit.default_timer()

cells, nodes, faces = create_local_mesh("veryfine_mesh.msh")


w=wn=np.zeros(len(cells.center))
vit = [0] * len(faces.nodeid)

for i in range(len(vit)):
    vit[i] = [2.,0]

nbelements = len(cells.center)
w=wn=np.zeros(nbelements)



for i in range(nbelements):
    if(cells.center[i][0] > 20 and cells.center[i][0] < 30 ):
        w[i] = 10


t = 0
Tfinal = 20
dt = 0

n = 0    
#loop over time
while(t<Tfinal):
    dt = 0.01
    t = t + dt
    n += 1
    w_ghost = ghost_value(w, faces)      
    rezidus = ExplicitScheme(w, vit, w_ghost, faces)
    
    for i in range(nbelements):
        wn[i]= w[i] + dt * (rezidus[i]/cells.volume[i])
        
        
    save_paraview_results(w, n, dt, cells.nodeid.values(), nodes.vertex)

    
    w = wn
    
stop = timeit.default_timer()

print(stop - start)