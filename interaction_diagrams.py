import math
import numpy as np

class interaction_diagram:
    def __init__(self, globvar):

        self.MNc_steel = []
        self.MNc_concrete = []

    def clear_results(self):
        self.MNc_steel = []
        self.MNc_concrete = []

    def distance_to_CG(self,globvar,dd):
        # d = distance of the stress resultant from the compression surface
        return (globvar.By / 2. - dd)
        
    def strain(self,globvar,dd,c):
        # c = position of neutral axis from the compression surface
        # d = distance of the stress resultant from the compression surface
        if ( c > 0. ):
            eps = self.eps_u * (c-dd) / c
        else: # uniform tension equal to strain at yielding
            eps = globvar.fy_vert/globvar.Es
        return eps
  
    
    def steel_stress(self,globvar,dd,c):

        eps = self.strain(globvar,dd,c)
        sigma = eps * globvar.Es
        if abs(sigma) > globvar.fy_vert:
            sigma = np.sign(sigma) * globvar.fy_vert
        return sigma
    
    
    def compute_steel_id(self, globvar):

        for c in np.linspace(self.c_min,self.c_max,self.c_div):

            M, N = self.compute_steel_MN(globvar, c)
            MNc = [M, N, c]
            self.MNc_steel.append(MNc)

            
    def compute_steel_MN(self, globvar, c):

        M = 0.
        N = 0.

        for reb in globvar.rebars:
            if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):

                dd = globvar.By /2. - reb.XYZ[0][1]
                area = reb.give_area(globvar)
                F_i = area * self.steel_stress(globvar,dd,c)
                N += F_i   
                M += -1. * F_i * self.distance_to_CG(globvar,dd)
                
        return [M, N]
            
        
    def compute_concrete_id(self, globvar):

        
        for c in np.linspace(self.c_min,self.c_max,self.c_div):

            M, N = self.compute_concrete_MN(globvar, c)
            MNc = [M, N, c]
            self.MNc_concrete.append(MNc)


    def compute_concrete_MN(self, globvar, c):
        
        fcm = globvar.concretes[globvar.fcm]  
        
        # depth of compressed concrete block
        if ( self.beta * c < globvar.By ):
            a = self.beta * c
        else:
            a = globvar.By

        # resultant of stress in compressed concrete    
        Fc = -1. *  self.eta * fcm * a * globvar.Bx
        
        N = Fc
        M = -1. * Fc * self.distance_to_CG(globvar,a/2.) 

        return [M, N]

    # getting combination MN for cases when the interaction diagram is not sampled uniformly (MSR reinforcement e.g.
    def interpolate_concrete_MN(self,c_target):

        Mc, Nc, c = zip(* self.MNc_concrete)

        if ( c_target <= min(c) or c_target >= max(c) ):
            return [math.nan, math.nan]

        i = 0
        while (c[i] < c_target):
            i += 1

        c_ii = c[i]
        c_i = c[i-1]

        aux = (c_target - c[i-1])/(c[i] - c[i-1])
        M = Mc[i-1] + aux*(Mc[i] - Mc[i-1])
        N = Nc[i-1] + aux*(Nc[i] - Nc[i-1])

        return [M,N]


    def find_load_for_eccentricity(self,ecc_fem):

        Mc, Nc, c = zip(* self.MNc_concrete)
        Ms, Ns, c = zip(* self.MNc_steel)

        M = [sum(i) for i in zip(Mc, Ms)]
        M.reverse()
        N = [sum(i) for i in zip(Nc, Ns)]
        N.reverse()

        
        if (ecc_fem == 0.):
            return abs(N[0])

        else:
            N_prev = N[0]
            ecc_id_prev = abs(M[0]/N[0])
        
        for Mi, Ni in zip(M, N):
            ecc_id = abs(Mi/Ni)
            
            if (ecc_id) > ecc_fem:
                aux = (ecc_fem - ecc_id_prev)/(ecc_id - ecc_id_prev)
                return abs(N_prev + aux * (Ni - N_prev))
            else:
                ecc_id_prev = ecc_id
                N_prev = Ni

        return math.nan
            
        
class id_ACI(interaction_diagram):
    def __init__(self, globvar):
        super().__init__(globvar)

        # ultimate concrete strain
        self.eps_u = -3.e-3

        # reduction of concrete block (lambda is taken by python), beta_1 is used in ACI
        self.beta = 0.85
        fcm = globvar.concretes[globvar.fcm]      
        if (fcm > 28.):
            self.beta -= (fcm-28.)/7. * 0.05
            self.beta = max(self.beta, 0.65)
            
        # reduction of concrete strength (parameter gamma in ACI)
        self.eta = 0.85

        # c is the width of the compressed zone (i.e. distance from the compressed fibers to neutral axis)
        self.c_min = 0.
        self.c_max = globvar.By * abs(self.eps_u)/( abs(self.eps_u) - globvar.fy_vert/globvar.Es )
        #print(self.c_max)
        self.c_div = 100


class id_MC2010(interaction_diagram):
    def __init__(self, globvar):
        super().__init__(globvar)

        # there is no reduction in all three parameters for fc <= C50
        # ultimate concrete strain
        self.eps_u = -3.5e-3

        # reduction of concrete block (lambda is taken by python)
        self.beta = 0.8
            
        # reduction of concrete strength (parameter gamma in ACI)
        self.eta = 1.0

        # c is the width of the compressed zone (i.e. distance from the compressed fibers to neutral axis)
        self.c_min = 0.
        self.c_max = globvar.By * abs(self.eps_u)/( abs(self.eps_u) - globvar.fy_vert/globvar.Es )
        #print(self.c_max)
        self.c_div = 100


class id_MSR(interaction_diagram):
    def __init__(self, globvar):
        super().__init__(globvar)

        # strain at which concrete reaches its strength
        self.eps_c = -2.e-3
        # ultimate concrete strain
        self.eps_u = -3.5e-3

        # reduction of concrete block (lambda is taken by python)
        # value 0.875 is the average for different concrete grades as well as for the expected confinement
        self.beta = 0.875
            
        # reduction of concrete strength (parameter gamma in ACI)
        self.eta = 1.0

        # c is the width of the compressed zone (i.e. distance from the compressed fibers to neutral axis)
        self.c_min = 0.
        self.c_max = globvar.By * abs(self.eps_u)/( abs(self.eps_u) - globvar.fy_vert/globvar.Es )
        #print(self.c_max)
        self.c_div = 100

        self.eps_uc = 0.
        self.eps_cc = 0.

        self.single_confined_areas = []
        self.double_confined_areas = []

      
    def update_MSR_intersections(self, globvar):

        # erase former entries
        self.single_confined_areas = []
        self.double_confined_areas = []

        # perhaps loop over spirals which do not intersect should be added
        
        for reb_L in globvar.rebars:
            # outer loop over all large spirals and find all potential intersecting small spirals
            if ( (reb_L.give_rebar_type() == 'spiral') and (reb_L.tag[0] == 'L') ):

                for reb_S in globvar.rebars:
                    # outer loop over all large spirals and find all potential intersecting small spirals
                    if ( (reb_S.give_rebar_type() == 'spiral') and (reb_S.tag[0] == 'S') ):

                        # determine overlap of small and large spirals if any
                        # CtC = center-to-center distance
                        CtC = math.sqrt( (reb_L.XYZ[0]-reb_S.XYZ[0])**2 + (reb_L.XYZ[1]-reb_S.XYZ[1])**2 )
                        if (CtC < reb_L.radius + reb_L.radius ):

                            ###
                            # 1) calculate areas and auxiliary variables
                            ###
                            
                            # https://www.xarg.org/2016/07/calculate-the-intersection-points-of-two-circles/
                            # x-coordinate of intersection in LCS from center of spiral L
                            L_to_I1I2 = ( reb_L.radius**2 - reb_S.radius**2 + CtC**2 ) / (2.*CtC)

                            # x-coord in LCS from center of spiral S
                            S_to_I1I2 = CtC - L_to_I1I2

                            # y-coordinate of intersection
                            y_inter = math.sqrt ( reb_L.radius**2 - L_to_I1I2**2 )
                                                      
                            # first base vectors (orientaction from large to small)
                            e1 = [ ( reb_S.XYZ[0] - reb_L.XYZ[0] ) / CtC,
                                   ( reb_S.XYZ[1] - reb_L.XYZ[1] ) / CtC]
                            #e2 = [ -e1[1], e1[0] ] .... opposite direction needed
                            e2 = [ e1[1], -e1[0] ]

                            rot_matrix = np.array( [e1, e2] )

                            vec_1 = [L_to_I1I2, y_inter]
                            vec_2 = [L_to_I1I2, -y_inter]

                            # center of large spiral [x,y]
                            L_CG = np.array(reb_L.XYZ[0:2])
                            # center of small spiral [x,y]
                            S_CG = np.array(reb_S.XYZ[0:2])
                            
                            # first intersection x, y
                            I1 = L_CG + np.dot( rot_matrix, vec_1)
                            I2 = L_CG + np.dot( rot_matrix, vec_2)
                            # middle of intersections
                            I1I2 = (I1 + I2)/2.
                            # length of this chord (formed by intersections)
                            I1I2_length = np.linalg.norm(I1 - I2)

                            # this angle is always smaller than pi, even for circles close to each other
                            # central angle given by the center of small spiral and intersection with large spiral
                            theta_S = 2. * math.asin ( (I1I2_length / 2.) / reb_S.radius )

                            # need to correct theta to have a proper meaning (we use it further for calculating area etc.)
                            L_to_I1I2 = np.linalg.norm(I1I2 - L_CG)
                            if (CtC < L_to_I1I2):
                                theta_S = 2*math.pi - theta_S
                               
                            # center angle from large spiral
                            theta_L = 2. * math.asin ( (I1I2_length / 2.) / reb_L.radius )

                            # area of the Large segment = (rather small area - large spiral vs. chord
                            # area of segment = sector + triangle
                            AL_sector = theta_L * reb_L.radius**2 / 2.
                            AL_triangle = I1I2_length * L_to_I1I2 /2.
                            AL_segment = AL_sector - AL_triangle

                            # area of the Small sector - inside large spiral
                            AS_sector_inter = theta_S * reb_S.radius**2 / 2.
                            AS_triangle_inter =  I1I2_length * S_to_I1I2 /2.
                            AS_segment_inter = AS_sector_inter - AS_triangle_inter

                            # area of double confined area
                            A_D = AL_segment + AS_segment_inter
  
                            # area of single confined area in small spiral
                            AS_tot = math.pi * reb_S.radius**2

                            # area of the Small segment - outer part of the intersection I1-I2
                            # = small spiral minus segment of small spiral closer to the large spiral
                            AS_segment_outer = AS_tot - AS_segment_inter

                            AS_S = AS_segment_outer - AL_segment

                            ###
                            # 2) compute centers of gravity
                            ###

                            # distance of center of small segment (internal part) wrt center of small spiral (works also for theta > pi)
                            S_to_S_segment_inter = ( 4. * reb_S.radius * (math.sin(theta_S/2.))**3 ) / ( 3.* (theta_S-math.sin(theta_S) ) )

                            I1I2_to_S_segment_inter =  S_to_S_segment_inter - S_to_I1I2
                            
                            
                            # distance of large segment to center of large spiral
                            L_to_L_segment = ( 4. * reb_L.radius * (math.sin(theta_L/2.))**3 ) / ( 3.* (theta_L-math.sin(theta_L) ) )
                            # distance of large segment to intersection
                            I1I2_to_L_segment = L_to_L_segment - L_to_I1I2

                            # distance of center of small segment (external part) wrt center of small spiral
                            S_to_S_segment_outer = AS_segment_inter * S_to_S_segment_inter / AS_segment_outer
                            # modify - offset wrt intersection 
                            I1I2_to_S_segment_outer = S_to_S_segment_outer + S_to_I1I2

                            S_single_CG_LCS = np.array( [ ( AS_segment_outer *  I1I2_to_S_segment_outer - AL_segment * I1I2_to_L_segment ) / ( AS_segment_outer - AL_segment ), 0.] )
                            S_single_CG = I1I2 + np.dot( rot_matrix, S_single_CG_LCS)

                            # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
                            # self.topology_canvas.axes.plot( S_single_CG[0], S_single_CG[1], color='red', marker='2', markeredgewidth=2, markersize=10)


                            # condition wrt chord (intersection small vs. large)
                            D_CG_LCS = np.array( [ (I1I2_to_L_segment * AL_segment - I1I2_to_S_segment_inter * AS_segment_inter ) / ( AL_segment + AS_segment_inter), 0.] )

                            # center of intersection large-vs-small + transformed distance form the chord
                            D_CG = I1I2 + np.dot( rot_matrix, D_CG_LCS)

                            # self.topology_canvas.axes.plot( D_CG[0], D_CG[1], color='magenta', marker='1', markeredgewidth=2, markersize=10)


                            single_conf_A = confined_area("single",
                                                          I1, I2,
                                                          reb_S.XYZ[0]-reb_S.radius, reb_S.XYZ[0]+reb_S.radius,
                                                          reb_S.XYZ[1]-reb_S.radius, reb_S.XYZ[1]+reb_S.radius,
                                                          S_single_CG, AS_S,
                                                          [reb_S.tag, reb_L.tag] )
                            self.single_confined_areas.append(single_conf_A)
                            
                            double_conf_A = confined_area("double",
                                                          I1, I2,
                                                          reb_S.XYZ[0]-reb_S.radius, reb_S.XYZ[0]+reb_S.radius,
                                                          reb_S.XYZ[1]-reb_S.radius, reb_S.XYZ[1]+reb_S.radius,
                                                          D_CG, A_D,
                                                          [reb_S.tag, reb_L.tag] )
                            self.double_confined_areas.append(double_conf_A)


    def compute_confined_strength_increase(self, globvar, sigL):
        fcm = globvar.concretes[globvar.fcm] 
        delta_fcc = 3.5 * (sigL)**(3./4.) * fcm**(1./4.)
        # delta_fcc = 4.0*sigL
        return delta_fcc
        
    def compute_confined_strength(self, globvar, sigL):
        fcm = globvar.concretes[globvar.fcm] 
        fcc = fcm + self.compute_confined_strength_increase(globvar,sigL)
        return fcc



    # distance from compression surface where strain = eps_uc
    # this strain is at the axis of the small spiral on the compression face
    def give_surface_offset(self,globvar):

        return globvar.cover + (globvar.rebars_CS[globvar.DS]).diam / 2.
        

    # returns y-coordinate of a concrete block for particular selection of neutral axis "c"
    def give_concrete_block_limit_y(self,globvar,c):

        offset = self.give_surface_offset(globvar)
        
        return globvar.By/2. - ( offset + (c-offset) * self.beta )

    def set_limits_for_NO(self,globvar):

        offset = self.give_surface_offset(globvar) + 2.e-3
        
        self.c_min = offset + globvar.dS / self.beta
        self.c_max = offset + (globvar.By - 2.*offset - globvar.dS) / self.beta


    def compute_segment_characteristics(self,globvar,tag,line_y):


        
        CG_segment = [math.nan, math.nan]
        A_segment = math.nan
                    
        for spiral in globvar.rebars:
            if ( (spiral.give_rebar_type() == 'spiral') and (spiral.tag == tag) ):

                
                # spiral center of gravity
                CG_spiral = np.array(spiral.XYZ[0:2])

                # no contribution
                if ( line_y >= CG_spiral[1] + spiral.radius ):
                    CG_segment = [math.nan, math.nan]
                    A_segment = 0.

                # full contribution
                elif (line_y <= CG_spiral[1] - spiral.radius):
                    CG_segment = CG_spiral
                    A_segment = math.pi * (spiral.radius)**2

                # partial contribution
                else:
                    line_vs_spiral_1 = np.array( [ CG_spiral[0] - math.sqrt(spiral.radius**2. - ( line_y-CG_spiral[1] )**2. ) ,
                                               line_y ] )

                    line_vs_spiral_2 = np.array( [ CG_spiral[0] + math.sqrt(spiral.radius**2. - ( line_y-CG_spiral[1] )**2. ) ,
                                               line_y ] )
                            
                    # Length of the intersection - horizontal line vs. spiral
                    length_line_vs_spiral = np.linalg.norm(line_vs_spiral_1 - line_vs_spiral_2)

                    # angle from the center of the large spiral wrt intersection with line
                    theta_line_vs_spiral = 2. * math.asin( length_line_vs_spiral / (2. * spiral.radius) )

                    # allow angles larger than 180 deg - neutral axis is below centroid of large spiral
                    if ( line_y < CG_spiral[1] ):
                        theta_line_vs_spiral = 2.*math.pi - theta_line_vs_spiral
                            
                    # parts of the total area of the large spiral above neutral axis
                    # sector given by the angle
                    A_sector = theta_line_vs_spiral * spiral.radius**2. / 2.
                    # remainder - triangle - either needs to be added or subtracted
                    A_triangle = (CG_spiral[1]-line_y) * length_line_vs_spiral / 2.
        
                    # total area of the large spiral above neutral axis
                    A_segment = A_sector + A_triangle

                    # center of gravity large spiral cut by the neutral axis in global CS (small spirals not considered)
                    CG_segment = np.array( [ CG_spiral[0],
                                             CG_spiral[1] + ( 4. * spiral.radius * (math.sin(theta_line_vs_spiral/2.))**3 ) / ( 3.* (theta_line_vs_spiral-math.sin(theta_line_vs_spiral) ) ) ] )

                # watch the correct indent here! 
                break
        
        return [A_segment, CG_segment]
    

    def strain(self,globvar,dd,c):
        # c = position of neutral axis from the compression surface
        # d = distance of the stress resultant from the compression surface
        # assumption - the ultimate strain is not reached at concrete surface but just beneath the concrete cover + DS/2 (it can be assumed that small spiral has smaller rebar diameter than the large one)

        # transform to reflect the above asumption
        offset = self.give_surface_offset(globvar)
        
        c = c - offset
        dd = dd - offset

        # neutral axis below concrete cover
        if ( c > 0. ): 
            eps = self.eps_uc * (c-dd) / c
        # neutral axis above (compressed) concrete cover
        else: # uniform tension equal to strain at yielding
            eps = globvar.fy_vert/globvar.Es
        return eps



    def compute_steel_id(self, globvar):

        
        MNc = super().compute_steel_id(globvar)

        c = self.give_surface_offset(globvar)
        # c += 5.e-3
        M, N = self.compute_steel_MN(globvar, c)
        MNc = [M, N, c]
        self.MNc_steel.insert(0,MNc)
        
        c = globvar.By / self.beta
        M, N = self.compute_steel_MN(globvar, c)
        MNc = [M, N, c]
        self.MNc_steel.append(MNc)
        
    
    def compute_concrete_id(self, globvar):

        self.set_limits_for_NO(globvar)
        
        fcm = globvar.concretes[globvar.fcm]

        sigL_L = 0.
        ke_L = 0.
        sigL_S = 0.
        ke_S = 0.

        # evaluate expected lateral confinement in large and small spirals
        for reb in globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):
                sigL_L = reb.compute_confinement(globvar, globvar.fy_lat)
                ke_L = reb.compute_confinement_effectiveness()
                break

        for reb in globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):
                sigL_S = reb.compute_confinement(globvar, globvar.fy_lat)
                ke_S = reb.compute_confinement_effectiveness()
                break

            
        min_sigL_eff = min(sigL_L*ke_L, sigL_S*ke_S)
        max_sigL_eff = max(sigL_L*ke_L, sigL_S*ke_S)

        # strain at which the confined strength is attained
        self.eps_uc = self.eps_u - 0.2 * min_sigL_eff / fcm 
        self.eps_cc = self.eps_c * ( 1. + 5. * ( ( max_sigL_eff )/fcm - 1.) )
        

               
        for c in np.linspace(self.c_min,self.c_max,self.c_div):
        
            M, N = self.compute_concrete_MN(globvar, c)
            MNc = [M, N, c]
            self.MNc_concrete.append(MNc)

        c = self.give_surface_offset(globvar)
        M, N = self.compute_concrete_MN(globvar, c)
        MNc = [M, N, c]
        self.MNc_concrete.insert(0,MNc)            

        c = globvar.By / self.beta
        M, N = self.compute_concrete_MN(globvar, c)
        MNc = [M, N, c]
        self.MNc_concrete.append(MNc)
           
    def compute_concrete_MN(self, globvar, c):

        #fcm = 0.
        fcm = globvar.concretes[globvar.fcm] 
        
        # c is the distance from the compression surface
        # neutral axis in global CS
        NO_y = self.distance_to_CG(globvar, c)
        # block line limit in global CS
        block_y = self.give_concrete_block_limit_y(globvar, c)

        # compute properties for the large spiral - will be needed anyway, no matter what approach is adopted
        N = 0.
        M = 0.

        
        # prepare vertical coordinates given by the topology of small spirals
        # spiral_S_y_coords = []

        '''
        for reb_S in globvar.rebars:
            # outer loop over all large spirals and find all potential intersecting small spirals
            if ( (reb_S.give_rebar_type() == 'spiral') and (reb_S.tag[0] == 'S') ):
                pass
        '''

        # single and double-confined areas are treated separately


        #####
        # SINGLE - small spirals
        #####        
        for single in self.single_confined_areas:
            # case #1
            # the entire spiral is in the compressed block
            # >>> the spiral parts should be treated individually - single and double-confined area
            if (single.y_min > block_y):
       
                for reb_S in globvar.rebars:
                    if ( (reb_S.give_rebar_type() == 'spiral') and (reb_S.tag == single.tags[0]) ):
                        break

                # add full contribution of single confined area
                AS_segment = single.A
                CGS_segment = single.CG
                    
                sigL_S = reb_S.compute_confinement(globvar, globvar.fy_lat)
                ke_S = reb_S.compute_confinement_effectiveness()

                # conventional strength + confined strength increase on effectively confined area
                F = AS_segment * fcm
                F += AS_segment * ke_S * self.compute_confined_strength_increase(globvar,sigL_S)

                if (F != 0.):
                    N -= F
                    M +=  F * CGS_segment[1]

            # case #2
            # only portion of the spiral is in the compressed concrete block
            # >>> contribution of single confined area, no subtraction, can be added rightaway
            elif (single.y_max > block_y):

                for reb_S in globvar.rebars:
                    if ( (reb_S.give_rebar_type() == 'spiral') and (reb_S.tag == single.tags[0]) ):
                        break
                
                AS_segment, CGS_segment = self.compute_segment_characteristics(globvar,reb_S.tag,block_y)
                
                sigL_S = reb_S.compute_confinement(globvar, globvar.fy_lat)
                ke_S = reb_S.compute_confinement_effectiveness()

                             # conventional strength + confined strength increase on effectively confined area
                F = AS_segment * fcm
                F += AS_segment * ke_S * self.compute_confined_strength_increase(globvar,sigL_S)

                if (F != 0.):
                    N -= F
                    M +=  F * CGS_segment[1]

            # the entire spiral is in tension
            # >>> no contribution
            else:
                N -= 0.
                M += 0.

        #####
        # DOUBLE
        #####        
        for double in self.double_confined_areas:
            # the entire spiral is in the compressed block, otherwise the double-confined area cannot develop
            if (double.y_min > block_y):
       
                for reb_S in globvar.rebars:
                    if ( (reb_S.give_rebar_type() == 'spiral') and (reb_S.tag == double.tags[0]) ):
                        break

                for reb_L in globvar.rebars:
                    if ( (reb_L.give_rebar_type() == 'spiral') and (reb_L.tag == double.tags[1]) ):
                        break                    

                # add full contribution of double confined area
                AD_segment = double.A
                CGD_segment = double.CG
                    
                sigL_S = reb_S.compute_confinement(globvar, globvar.fy_lat)
                ke_S = reb_S.compute_confinement_effectiveness()
                
                sigL_L = reb_L.compute_confinement(globvar, globvar.fy_lat)
                ke_L = reb_L.compute_confinement_effectiveness()

                ke_D = min(ke_S,ke_L)
                sigL_D = sigL_S+sigL_L

                # conventional strength + confined strength increase on effectively confined area
                F = AD_segment * fcm
                F += AD_segment * ke_D * self.compute_confined_strength_increase(globvar,sigL_D)

                if (F != 0.):
                    N -= F
                    M +=  F * CGD_segment[1]

            # the entire spiral is in tension
            # >>> no contribution
            else:
                N -= 0.
                M += 0.       
      
        for reb_L in globvar.rebars:


            
            # outer loop over all large spirals and find all potential intersecting small spirals
            if ( (reb_L.give_rebar_type() == 'spiral') and (reb_L.tag[0] == 'L') ):
              
                AL_segment, CGL_segment = self.compute_segment_characteristics(globvar,reb_L.tag,block_y)

                sigL_L = reb_L.compute_confinement(globvar, globvar.fy_lat)
                ke_L = reb_L.compute_confinement_effectiveness()

                # treat double confined areas which need to be subtracted

                # cumulative area and static moment of segments
                AD_segments = 0.
                SD_segments = 0.

                AL_segment_net = AL_segment
                CGL_segment_net = CGL_segment
                
                for double in self.double_confined_areas:
                    if ( (reb_L.tag == double.tags[1]) and (double.y_min > block_y) ) :

                        # add full contribution of double confined area
                        AD_segments += double.A
                        SD_segments += double.CG[1] * double.A

                AL_segment_net -= AD_segments

                if (AL_segment_net):
                    CGL_segment_net[1] = ( CGL_segment[1] * AL_segment - SD_segments ) / AL_segment_net

                # conventional strength + confined strength increase on effectively confined area
                F = AL_segment_net * fcm
                F += AL_segment_net * ke_L * self.compute_confined_strength_increase(globvar,sigL_L)

                # to get rid of nans which appears in CG when A == 0.
                if (F != 0.):
                    N -= F
                    M +=  F * CGL_segment_net[1]
                
                break
            
   
        return [M, N]


class id_MSR_simple(id_MSR):
    def __init__(self, globvar):
        super().__init__(globvar)


    #TODO
    def set_limits_for_NO(self,globvar):

        #offset = self.give_surface_offset(globvar) 
        offset = globvar.cover
        self.c_min = offset 
        self.c_max = offset + (globvar.By - 2.*offset) / self.beta

        
    ##
    # MSR simple
    ##
    def compute_concrete_id(self, globvar):

        self.set_limits_for_NO(globvar)
        
        fcm = globvar.concretes[globvar.fcm]

        sigL_L = 0.
        ke_L = 0.
        sigL_S = 0.
        ke_S = 0.

        # evaluate expected lateral confinement in large and small spirals
        for reb in globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):
                sigL_L = reb.compute_confinement(globvar, globvar.fy_lat)
                ke_L = reb.compute_confinement_effectiveness()
                break

        for reb in globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):
                sigL_S = reb.compute_confinement(globvar, globvar.fy_lat)
                ke_S = reb.compute_confinement_effectiveness()
                break

            
        min_sigL_eff = min(sigL_L*ke_L, sigL_S*ke_S)
        max_sigL_eff = max(sigL_L*ke_L, sigL_S*ke_S)

        # strain at which the confined strength is attained
        self.eps_uc = self.eps_u - 0.2 * min_sigL_eff / fcm 
        self.eps_cc = self.eps_c * ( 1. + 5. * ( ( max_sigL_eff )/fcm - 1.) )
        
               
        for c in np.linspace(self.c_min,self.c_max,self.c_div):
        
            M, N = self.compute_concrete_MN(globvar, c)
            MNc = [M, N, c]
            self.MNc_concrete.append(MNc)

    ##
    # MSR simple
    ##
    def compute_concrete_MN(self, globvar, c):

        fcm = globvar.concretes[globvar.fcm] 
        
        # c is the distance from the compression surface
        # neutral axis in global CS
        NO_y = self.distance_to_CG(globvar, c)
        # block line limit in global CS
        block_y = self.give_concrete_block_limit_y(globvar, c)

        # compute properties for the large spiral - will be needed anyway, no matter what approach is adopted
        N = 0.
        M = 0.


        #####
        # PLAIN CONCRETE STRENGTH ON A RECTANGLE SUPERSCRIBED AROUND MSRs
        #####
        #offset = self.give_surface_offset(globvar)
        # this approach vvv might give better agreement with FEM 
        offset = globvar.cover
        F = (globvar.Bx - 2.*offset ) * (globvar.By/2. - offset - block_y )  * fcm
        CGS_block = (globvar.By/2. - offset + block_y)/2.
        
        if (F > 0.):
             N -= F
             M +=  F * CGS_block
        
             
        for reb in globvar.rebars:

            # outer loop over all large spirals and find all potential intersecting small spirals
            if ( reb.give_rebar_type() == 'spiral' ):

                A_segment, CG_segment = self.compute_segment_characteristics(globvar,reb.tag,block_y)

                sigL = reb.compute_confinement(globvar, globvar.fy_lat)
                ke = reb.compute_confinement_effectiveness()

                # confined strength increase on effectively confined area
                F = A_segment * ke * self.compute_confined_strength_increase(globvar,sigL)

                # to get rid of nans which appears in CG when A == 0.
                if (F > 0.):
                    N -= F
                    M +=  F * CG_segment[1]
  
        return [M, N]            

        
        

class confined_area:
    def __init__(self, confinement, I1, I2, x_min, x_max, y_min, y_max, CG, A, tags):

        # single or double confinement
        self.confinement = confinement        
        # intersection of the two spirals
        self.I1 = I1
        self.I2 = I2
        # outer dimensions of the small spirals
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        # center of gravity of the confined area
        self.CG = CG
        # area of the confined area
        self.A = A
        # tags specifying the intersecting spirals
        self.tags = tags
