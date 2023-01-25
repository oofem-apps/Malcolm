# import globvars
import warnings


def write_heading(f, s):
    f.write(f"#\n#\n#\n")
    f.write(f"### {s}: ###\n")
    f.write(f"#\n")


def write_oofem_input_header(globvar, task, f):

    f.write(f"{globvar.filename_out}\n")
    f.write(f"#\n")
    f.write(f"Simulation of a multispiral rectangular column with dimensions Bx x By x H = {(globvar.Bx):.3f} x {(globvar.By):.3f} x {(globvar.H):.3f} m\n")
    f.write(f"#\n")
        
    # definition of material properties
        
    f.write(f"# Concrete cover = {(globvar.cover):.3f} m\n")
        
    f.write(f"# Reinforcement ratio VERTICAL direction RHO_V {(globvar.rho_vert*100):.3f}%\n")

    f.write(f"# Reinforcement ratio LATERAL direction RHO_L {(globvar.rho_lat*100):.3f}%\n")

    f.write(f"#\n# Rebar specification:\n")
    # definition of spiral reinforcement
    f.write(f"#\n# large spirals:\n")
    for reb in globvar.rebars:
        if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

            f.write(f"# Spiral {reb.tag}, rebar D = {(reb.give_diameter(globvar))*1000.:.1f} mm, XYZ = [{(reb.XYZ[0]):.3f},{(reb.XYZ[1]):.3f},{(reb.XYZ[2]):.3f}], axial radius = {(reb.radius):.3f} m, pitch = {(reb.pitch):.3f} m, cross-section nr. {int(reb.cs_nr)}, material nr. {int(reb.give_mat_nr(globvar))}\n")

    f.write(f"#\n# small spirals:\n")
    for reb in globvar.rebars:
        if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):

            f.write(f"# Spiral {reb.tag}, rebar D = {(reb.give_diameter(globvar))*1000.:.1f} mm, XYZ = [{(reb.XYZ[0]):.3f},{(reb.XYZ[1]):.3f},{(reb.XYZ[2]):.3f}], axial radius = {(reb.radius):.3f} m, pitch = {(reb.pitch):.3f} m, cross-section nr. {int(reb.cs_nr)}, material nr. {int(reb.give_mat_nr(globvar))}\n")

    # definition of axial reinforcement
    f.write(f"#\n# vertical rebars:\n")
    for reb in globvar.rebars:
        if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):
            f.write(f"# Rebar {reb.tag}, rebar D = {(reb.give_diameter(globvar))*1000.:.1f} mm, XYZ = ")

            for xyz in reb.XYZ:
                f.write(f" [{(xyz[0]):.3f},{(xyz[1]):.3f},{(xyz[2]):.3f}],")
            f.write(f" cross-section nr. {int(reb.cs_nr)}, material nr. {int(reb.give_mat_nr(globvar))}\n")                

    f.write(f"#\n")
    
        
def write_oofem_engineering_model_modules_domain (globvar, task, f):
       
    f.write(f"#\n")

    f.write(task.engng_definition)
    
    f.write(f"#\n")
    f.write(f"#\n")
    f.write(f"# export modules:\n")
    f.write(f"# concrete\n")
    f.write(f"vtkxml tstep_step 5 domain_all vars 2 1 4 cellvars 2 27 13 primvars 1 1 stype 1 regionsets 1 1\n")
    f.write(f"# spirals\n")
    f.write(f"vtkxml tstep_step 5 domain_all vars 2 1 4 cellvars 1 27 primvars 1 1 stype 1 regionsets 2 2 3\n")
    f.write(f"# rebars\n")
    f.write(f"vtkxml tstep_step 5 domain_all vars 2 1 4 cellvars 1 27 primvars 1 1 stype 1 regionsets 1 4\n")
    f.write(f"#\n")
    f.write(f"#\n")
    f.write(f"domain 3D\n")
    f.write(f"#\n")
    f.write( "OutputManager tstep_all dofman_output {"  f"{int(globvar.master_node_force)} {int(globvar.master_node_shortening)} {int(globvar.master_node_bending_x)} {int(globvar.master_node_bending_y)}" + "}\n")
    f.write(f"#\n")

    # include artificially defined nodes for loading and master-slave conditions
    total_nodes = globvar.ndofman + 4
    
    f.write(f"ndofman {int(total_nodes)} nelem {int(globvar.nelem)} ncrosssect {int(len(globvar.cross_sections))} nmat {int(len(globvar.materials))} nbc {int(len(globvar.bcs))} nic 0 nltf {int(len(globvar.ltfs))} nset {int(globvar.sets)}\n")


    

def write_oofem_input_mesh_nodes(globvar, task, f):


    write_heading(f, "NODES")

    x_range = len(globvar.master_nodes)
    y_range = len(globvar.master_nodes[0])
    z_range = len(globvar.master_nodes[0][0])
        
    for iz in range(z_range):
        for iy in range(y_range):
            for ix in range(x_range):
                m_node = globvar.master_nodes[ix][iy][iz]

                if (m_node.tag):
                    if ( "doftype" in m_node.tag ): # slave node
                        f.write(f"slavenode")
                    else: # slave node
                        f.write(f"node")
                else:
                    f.write(f"node")

                f.write(f" {int(m_node.nr)} coords 3 {(m_node.XYZ[0]):.6f} {(m_node.XYZ[1]):.6f} {(m_node.XYZ[2]):.6f}")


                # write boundary conditions
                if (not m_node.tag == None ):
                    f.write(" " + m_node.tag)
                
                f.write(f"\n")

                   
        
def write_oofem_input_hanging_nodes(globvar, task, f):

    write_heading(f, "HANGING NODES")
    
    prev_tag = None
    for h_node in globvar.hanging_nodes:

        if prev_tag != h_node.tag:
            f.write(f"#\n# {h_node.tag} :\n#\n")
            prev_tag = h_node.tag
            
        f.write(f"hangingNode {int(h_node.nr)} coords 3 {(h_node.XYZ[0]):.6f} {(h_node.XYZ[1]):.6f} {(h_node.XYZ[2]):.6f} doftype 3 2 2 2\n")


def write_oofem_input_dummy_nodes(globvar, task, f):

    write_heading(f, "DUMMY NODES")

    f.write("# loading (indirect control)\n")
    f.write(f"hangingNode {int(globvar.master_node_force)} coords 3 {(task.eccentricity_actual[0]):.6f} {(task.eccentricity_actual[1]):.6f} {(globvar.H):.6f} doftype 3 2 2 2 load 1 2\n")

    f.write("# shortening\n")
    f.write(f"node {int(globvar.master_node_shortening)} coords 3 {(0.):.6f} {(0.):.6f} {(globvar.H):.6f} bc 3 1 1 0\n")
    
    f.write("# rotation about x-axis\n")
    f.write(f"node {int(globvar.master_node_bending_x)} coords 3 {(0.):.6f} {(globvar.By/2.):.6f} {(globvar.H):.6f} bc 3 1 1 0\n")
    
    f.write("# rotation about y-axis\n")
    f.write(f"node {int(globvar.master_node_bending_y)} coords 3 {(globvar.Bx/2.):.6f} {(0.):.6f} {(globvar.H):.6f} bc 3 1 1 1\n")    
            
def write_oofem_input_brick_elements(globvar, task, f):

    write_heading(f, "BRICK ELEMENTS")

    for b_elem in globvar.brick_elements:

        f.write(f"LSpaceBB {int(b_elem.nr)} nodes 8")
        for n in range(8):
            f.write(f" {int(b_elem.nodes[n])}")
        f.write(f"\n")                    

def write_oofem_input_truss_elements(globvar, task, f):

    write_heading(f, "TRUSS ELEMENTS")

    prev_tag = None
    for t_elem in globvar.truss_elements:

        if prev_tag != t_elem.tag:
            f.write(f"#\n# {t_elem.tag} :\n#\n")
            prev_tag = t_elem.tag
            
        f.write(f"Truss3d {int(t_elem.nr)} nodes 2 {int(t_elem.nodes[0])} {int(t_elem.nodes[1])}\n")
            
def write_oofem_input_cross_sections(globvar, task, f):
    write_heading(f, "CROSS-SECTIONS")

    for cs in globvar.cross_sections:      
        f.write(f"SimpleCS {int(cs.nr)}")
        if (cs.cs_type == "rebar"):
            f.write(f" area {(cs.area):.6f}")
        f.write(f" material {int(cs.mat)} set {int(cs.elem_set)}\n")
        
def write_oofem_input_materials(globvar, task, f):
    
    write_heading(f, "MATERIAL DEFINITIONS")

    for mat in globvar.materials:
        mat_type = mat.give_material_type() 
        if (mat_type == "concrete"):
            f.write(f"Con2DPM {int(mat.nr)} d 0. talpha 0. E {(mat.Eci):.3f} n {(mat.nu):.3f} fc {(mat.fcm):.3f} ft {(mat.fctm):.3f} wf {(mat.wf):.4e} ecc {(mat.ecc):.3f} kinit {(mat.kinit):.3f} Hp {(mat.Hp):.3f} dilation {(mat.dilation):.3f} Ahard {(mat.Ahard):.4e} Bhard {(mat.Bhard):.4e} Chard {(mat.Chard):.4e} Dhard {(mat.Dhard):.4e} Asoft {(mat.Asoft):.3f} efc {(mat.efc):.4e} stype 1 wf1 {(mat.wf1):.3f} ft1 {(mat.ft1):.3f}\n")
            
        elif (mat_type == "rebar"):
            f.write(f"MisesMat {int(mat.nr)} d 0. talpha 0. E {(mat.E):.3f} n 0.3 sig0 {(mat.sig_0):.3f} H {(mat.H):.3f} omega_crit {(mat.omega_c):.3f} a {(mat.a):.3f}\n")
            
        elif (mat_type == "tendon"):
            print("gogo")
        else:
            warnings.warn("unknown material type")


def write_oofem_input_boundary_conditions(globvar, task, f):
    
    write_heading(f, "BOUNDARY CONDITIONS")

    for bc in globvar.bcs:      
        f.write(f"{bc}\n")

def write_oofem_input_time_functions(globvar, task, f):

    write_heading(f, "TIME FUNCTIONS")

    for ltf in globvar.ltfs:      
        f.write(f"{ltf}\n")
    
            
def write_oofem_input_elements_sets(globvar, task, f):

    write_heading(f, "SETS")
    
    # write element sets
    for i_cs in globvar.cross_sections:
    
        lower = None
        upper = None

        set_ranges = []
    
        for i_elem in i_cs.elements:
            if (lower == None):
                lower = i_elem
                prev = i_elem
                continue

            # hit last entry
            if (i_elem == i_cs.elements[-1]):
                upper = i_elem

                i_range = (lower,upper)
                set_ranges.append(i_range)
                break

            # continuous numbering
            if (i_elem == prev + 1):
                prev = i_elem
                # continue

                # suddenly discontinuous
            else:
                upper = prev

                i_range = (lower,upper)
                set_ranges.append(i_range)

                lower = i_elem
                prev = i_elem
                upper = None

        f.write(f"# elements set for cross-section nr. {int(i_cs.nr)} :\n")            
                    
        f.write(f"set {int(i_cs.nr)} elementranges" + " { ")            
        for i_range in set_ranges:
            f.write(f"({int(i_range[0])} {int(i_range[1])}) ")
        f.write("}\n")

            
def write_oofem_input_extractor_record(globvar, task, f):

    write_heading(f, "EXTRACTOR")
    
    f.write("#%BEGIN_CHECK%\n")
    f.write("#TIME\n")
    f.write("#LOADLEVEL\n")
    f.write(f"#NODE number {int(globvar.master_node_force)} dof 3 unknown d\n")
    f.write(f"#NODE number {int(globvar.master_node_shortening)} dof 3 unknown d\n")
    f.write(f"#NODE number {int(globvar.master_node_bending_x)} dof 3 unknown d\n")
    f.write(f"#NODE number {int(globvar.master_node_bending_y)} dof 3 unknown d\n")
    f.write(f"#REACTION number {int(globvar.master_node_shortening)} dof 3\n")
    f.write(f"#REACTION number {int(globvar.master_node_bending_x)} dof 3\n")
    f.write(f"#REACTION number {int(globvar.master_node_bending_y)} dof 3\n")
    f.write("#%END_CHECK%\n")


def write_program_info(globvar, task, f):
    
    f.write(f"----------------------------------------------------------------------\n")
    f.write(f"-                                                                    -\n")
    f.write(f"-    MaLCoLM version 2.0, December 2022                              -\n")
    f.write(f"-                                                                    -\n")
    f.write(f"-     Multi-spiral column simulation module                          -\n")
    f.write(f"-     OOFEM extension module                                         -\n")
    f.write(f"-                                                                    -\n")
    f.write(f"----------------------------------------------------------------------\n")
    f.write(f"-     Author:      Petr Havlasek                                     -\n")
    f.write(f"-     Institution: Department of Mechanics                           -\n")
    f.write(f"-                  Faculty of Civil Engineering                      -\n")
    f.write(f"-                  Czech Technical University in Prague              -\n")
    f.write(f"----------------------------------------------------------------------\n")
    f.write(f"- Acknowledgment:                                                    -\n")
    f.write(f"- Financial support for this work was provided by the Technology     -\n")
    f.write(f"- Agency of the Czech Republic (TAČR), project TM01000059 (CeSTaR 2) -\n")
    f.write(f"-                                                                    -\n")
    f.write(f"----------------------------------------------------------------------\n")
