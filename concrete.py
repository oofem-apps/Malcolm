import math
import numpy as np

# import globvars
import warnings

from cross_sections import cross_section
from nodes import node
from elements import element

class concrete:
    def __init__(self, globvar, cs, mat):

        self.cs = cs
        self.mat = mat
        concrete_cs = cross_section(cs, mat) 
        globvar.cross_sections.append(concrete_cs)
        
    def create_concrete_mesh(self, globvar):

        Bx_net = globvar.Bx
        By_net = globvar.By

        # subtract concrete cover to produce mesh which approximately corresponds to the confined region
        if (globvar.flag_ignore_cover):
            Bx_net -= 2. * globvar.cover
            By_net -= 2. * globvar.cover
            
        n_elem_X = math.ceil(Bx_net/globvar.elem_size_X);
        n_elem_Y = math.ceil(By_net/globvar.elem_size_Y);
        n_elem_Z = math.ceil(globvar.H/globvar.elem_size_Z);


        ### NODES
        # master_nodes[ix][iy][iz]
        # ix -> x coordinate
        # iy -> y coordinate
        # iz -> z coordinate
        
        # test:            
        # for ix in range (n_elem_X+1):
        #    print( master_nodes[ix][0][0].XYZ )                 

        globvar.master_nodes = [[[ [] for _ in range(n_elem_Z+1)] for _ in range(n_elem_Y+1)] for _ in range(n_elem_X+1)]
        
        master_nodes_count = 0
        for iz in range (n_elem_Z+1):
            for iy in range (n_elem_Y+1):
                for ix in range (n_elem_X+1):
                    master_nodes_count +=  1
                    coords = [ ix * globvar.elem_size_X - Bx_net/2., iy * globvar.elem_size_Y - By_net/2., iz * globvar.elem_size_Z ]

                    if (iz == 0):
                        if (ix == 0 and iy == 0):
                            tag = "bc 3 1 1 1"
                        elif (ix == n_elem_X and iy == 0):
                            tag = "bc 3 0 1 1"
                        else:
                            tag = "bc 3 0 0 1"
                            
                    elif (iz == n_elem_Z):
                        tag = "doftype 3 2 2 2"
                        tag += " masterDofMan 4 " + str( (globvar.master_nodes[ix][iy][0]).nr ) # bottom node
                        tag += " " + str( globvar.master_node_shortening )
                        tag += " " + str( globvar.master_node_bending_x )
                        tag += " " + str( globvar.master_node_bending_y )

                        tag += " weights 4 1. 1."
                        tag += " " + (f"{(coords[1] / (globvar.By/2.)):.6f}")
                        tag += " " + (f"{(coords[0] / (globvar.Bx/2.)):.6f}")
                            
                    else:
                        tag = None
                    
                    master_node = node( nr = master_nodes_count, XYZ = coords, tag = tag)
                    
                    globvar.master_nodes[ix][iy][iz] = master_node

        globvar.ndofman = master_nodes_count

        ### ELEMENTS
        globvar.brick_elements.clear()

        # indices differences to easily select particular node
        dx = [0, 0, 1, 1, 0, 0, 1, 1]
        dy = [0, 1, 1, 0, 0, 1, 1, 0]
        dz = [1, 1, 1, 1, 0, 0, 0, 0]

        elem_nr = 0
        for iz in range (n_elem_Z):
            for iy in range (n_elem_Y):
                for ix in range (n_elem_X):
            
                    elem_nr += 1

                    brick_nodes = []
                    for i in range(8):

                        brick_node = globvar.master_nodes[ix + dx[i] ][iy + dy[i] ][iz + dz[i]]
                        brick_nodes.append(brick_node.nr)
                    brick_elem = element (globvar = globvar, nr = elem_nr, nodes = brick_nodes, cs = 1, tag = "concrete" )
                    globvar.brick_elements.append(brick_elem)
                    
        globvar.nelem = elem_nr 
        
