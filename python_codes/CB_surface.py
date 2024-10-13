import numpy as np
import math as math

# Here are all the vertices (formed by one 6MR window each, identified by the three existing
# surface ligands)
Ai = np.array([[5253, 5160, 5161]])
Bi = np.array([[5250, 5382, 5309]])
Ci = np.array([[5232, 5252, 5233]])
Di = np.array([[5103, 5101, 5105]])
Ei = np.array([[5169, 5269, 5312]])
Fi = np.array([[5319, 5391, 5399]])
Gi = np.array([[5242, 5383, 5272]])
Hi = np.array([[5120, 5121, 5255]])

# Please visualize in ovito which ligand IDs above are at which vertice and also
# which vertice forms an edge and faces.

# --------------------------------------------------------------------------------
NA = 69384

number_of_configurations = 50

trj = np.zeros((1,5))
box_size = np.zeros((1,2))
ofi = open("dump.pos", 'r')
for it_1 in range (0, number_of_configurations):
    add=[]
    tmp = []
    ofi.readline()
    ofi.readline()
    ofi.readline()
    ofi.readline()
    ofi.readline()
    for it_3 in range (0,3):
        line = ofi.readline()
        for f in zip(line.split(' ')):
            # e = list(e)*
            f = np.array(f)
            tmp = np.append(tmp,float(f))
    tmp = np.array(tmp,float)
    tmp = tmp.reshape(3,2)
    box_size = np.concatenate((box_size, tmp), axis = 0)
    ofi.readline()
    for it_2 in range(0, (NA)):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e in zip(dump.split(' ')):
            # e = list(e)*
            e = np.array(e)
            add = np.append(add,float(e))
    add = np.array(add,float)
    add = add.reshape(NA,8)
    # -------------------------------
    add = np.delete(add, 7, 1)
    add = np.delete(add, 6, 1)
    add = np.delete(add, 5, 1)
    # -------------------------------
    trj = np.concatenate((trj, add), axis = 0)
 
trj = np.delete(trj,0, axis=0)
box_size = np.delete(box_size, 0, 0)

# -------------------------------------------------------------------------------------
# I will need to read this for later. It was built from simply copying pasting the bond
# section of the initial file used to start the CG-MD simulation of the nanoparticle.
number_of_bonds = 6900

bonds_mof = []
ofi = open("bonds_mof.dat", 'r')
for it_1 in range(0, number_of_bonds):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(4)):
            bonds_mof.append(float(e))
bonds_mof = np.array(bonds_mof,float)
bonds_mof = bonds_mof.reshape(number_of_bonds,4)

# ------------------------------------------------------------------------------------
# I will define here the histograms. There will be three: one for the of verteces, 
# one to the environment of the edges and another for the environment of the facets. 
# All of them will have histograms 12 angs long. The bins here will have 0.25 angs.
bin_size = 0.25
number_of_bins = int(12/bin_size)

bin_size_angle = 0.25
number_of_bins_angle = int(180/bin_size_angle)

# I will also output the density profile for the polymer phase. I wll not segregate superatoms types.
rho_grp = np.zeros((number_of_bins, number_of_bins_angle))
rho_edges = np.zeros((number_of_bins, number_of_bins_angle))
rho_facets = np.zeros((number_of_bins, number_of_bins_angle))

# ------------------------------------------------------------------------------------
M1 = 65.38
M2 = 81.0

number_of_configurations_considered = 50

for loop_0 in range (0, number_of_configurations_considered):
    
    output = np.zeros((NA,5))
    output[:,:] = trj[int((loop_0)*NA):int((loop_0 + 1)*NA),:]
    
    i = np.zeros((1,3))
    i[0,0] = box_size[int(loop_0*3),1] - box_size[int(loop_0*3),0]
    mod_i = ((i[0,0])**2 + (i[0,1])**2 + (i[0,2])**2)**(1/2)
    j = np.zeros((1,3))
    j[0,1] = box_size[int((loop_0*3) + 1),1] - box_size[int((loop_0*3) + 1),0]
    mod_j = ((j[0,0])**2 + (j[0,1])**2 + (j[0,2])**2)**(1/2)
    k = np.zeros((1,3))
    k[0,2] = box_size[int((loop_0*3) + 2),1] - box_size[int((loop_0*3) + 2),0]
    mod_k = ((k[0,0])**2 + (k[0,1])**2 + (k[0,2])**2)**(1/2)

    output[:,2] = (output[:,2]*i[0,0]) + box_size[int(loop_0*3),0]
    output[:,3] = (output[:,3]*j[0,1]) + box_size[int(loop_0*3) + 1,0]
    output[:,4] = (output[:,4]*k[0,2]) + box_size[int(loop_0*3) + 2,0]
    
    x_min = box_size[int(loop_0*3),0]
    y_min = box_size[int(loop_0*3) + 1,0]
    z_min = box_size[int(loop_0*3) + 2,0]
    
    # I think the first thing I will do is to shift the (xlo, ylo, zlo) of the simulation box
    # to the origin to facilitate further reasoning!
    translator = np.array([[x_min, y_min, z_min]])
    output[:,2] = output[:,2] - translator[0,0]
    output[:,3] = output[:,3] - translator[0,1]
    output[:,4] = output[:,4] - translator[0,2]
    # Naturally I need also to redefine the values of x_min, y_min and z_min defined previ-
    # ously.
    x_min = 0
    y_min = 0
    z_min = 0
    
    # -----------------------------------------------------------------------------------
    counter_xi = 0
    counter_yi = 0
    counter_zi = 0
    counter_xf = 0
    counter_yf = 0
    counter_zf = 0
    
    for it_1 in range (0, len(output)):
        # I will only run all of this stuff if I found a superatom of the MOF phase
        if (output[it_1,1] == 1) or (output[it_1,1] == 2):
            position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            # ---------------------------------------
            if position[0,0] > i[0,0]/2:
                counter_xf = counter_xf + 1
            if position[0,0] < i[0,0]/2:
                counter_xi = counter_xi + 1
            # ---------------------------------------
            if position[0,1] > j[0,1]/2:
                counter_yf = counter_yf + 1
            if position[0,1] < j[0,1]/2:
                counter_yi = counter_yi + 1
            # ---------------------------------------
            if position[0,2] > k[0,2]/2:
                counter_zf = counter_zf + 1
            if position[0,2] < k[0,2]/2:
                counter_zi = counter_zi + 1
        
    # By the end of this, I will be able to identify, for each of the three directions,
    # in which part of the simulation box most of the atoms are. 
    # Note that here I am considering a simulation box centered in the origin.
    
    # I will then have to shift the atoms that sit at the portion containing the least a-
    # mount of atoms AND which has been disrupted by periodic boundaries (I mean, it is 
    # possible to be sitting at the less populated side of a direction and not need to be 
    # moved).
    flag = 1
    while flag == 1:
        flag = 0
        # flag will only remain as 0 if, by the end of the it_1 for loop below, I never happen
        # to need to shift atoms due to crossing of the periodic boundary. Note that I am al-
        # ways shifting the superatoms in the direction of the more populated region, so that
        # I spend less time here.
        for it_1 in range (0, len(bonds_mof)):
            given_CG1 = bonds_mof[it_1,2]
            given_CG2 = bonds_mof[it_1,3]
            counter = 0
            # -----------------------------------------------------------------
            # Finding the positions. The variable counter is meant to speed up things
            # in the sense that I am creating a condition to break this for loop so  
            # that I dont need to scan it up to the end to find the position of the
            # two given superatoms.
            for it_2 in range (0, len(output)):
                if output[it_2,0] == given_CG1:
                    given_CG1_pos = np.array([[output[it_2,2], output[it_2,3], output[it_2,4]]])
                    line_CG1 = it_2
                    counter = counter + 1
                if output[it_2,0] == given_CG2:
                    given_CG2_pos = np.array([[output[it_2,2], output[it_2,3], output[it_2,4]]])
                    line_CG2 = it_2
                    counter = counter + 1
                if counter == 2:
                    break
            # -----------------------------------------------------------------
            dist_x = given_CG1_pos[0,0] - given_CG2_pos[0,0]
            dist_y = given_CG1_pos[0,1] - given_CG2_pos[0,1]
            dist_z = given_CG1_pos[0,2] - given_CG2_pos[0,2]
            
            # if the condition below is met, this means that I should be worrying about
            # translating superatoms downwards since the xi side is more populated. The
            # equal condition is included here but could be included in the other also.
            if int(counter_xi/counter_xf) >= 1:
                if (abs(dist_x) > mod_i/2) & (dist_x > 0):
                    output[line_CG1,2] = output[line_CG1,2] - i[0,0]
                    flag = 1
                if (abs(dist_x) > mod_i/2) & (dist_x < 0):
                    output[line_CG2,2] = output[line_CG2,2] - i[0,0]
                    flag = 1
            if int(counter_xi/counter_xf) < 1:
                if (abs(dist_x) > mod_i/2) & (dist_x > 0):
                    output[line_CG2,2] = output[line_CG2,2] + i[0,0]
                    flag = 1
                if (abs(dist_x) > mod_i/2) & (dist_x < 0):
                    output[line_CG1,2] = output[line_CG1,2] + i[0,0]
                    flag = 1
            # ---------------------------------------------------
            # Similarly, to the y direction:
            if int(counter_yi/counter_yf) >= 1:
                if (abs(dist_y) > mod_j/2) & (dist_y > 0):
                    output[line_CG1,3] = output[line_CG1,3] - j[0,1]
                    flag = 1
                if (abs(dist_y) > mod_j/2) & (dist_y < 0):
                    output[line_CG2,3] = output[line_CG2,3] - j[0,1]
                    flag = 1
            if int(counter_yi/counter_yf) < 1:
                if (abs(dist_y) > mod_j/2) & (dist_y > 0):
                    output[line_CG2,3] = output[line_CG2,3] + j[0,1]
                    flag = 1
                if (abs(dist_y) > mod_j/2) & (dist_y < 0):
                    output[line_CG1,3] = output[line_CG1,3] + j[0,1]
                    flag = 1
            # ---------------------------------------------------
            # Similarly, to the z direction:
            if int(counter_zi/counter_zf) >= 1:
                if (abs(dist_z) > mod_k/2) & (dist_z > 0):
                    output[line_CG1,4] = output[line_CG1,4] - k[0,2]
                    flag = 1
                if (abs(dist_z) > mod_k/2) & (dist_z < 0):
                    output[line_CG2,4] = output[line_CG2,4] - k[0,2]
                    flag = 1
            if int(counter_zi/counter_zf) < 1:
                if (abs(dist_z) > mod_k/2) & (dist_z > 0):
                    output[line_CG2,4] = output[line_CG2,4] + k[0,2]
                    flag = 1
                if (abs(dist_z) > mod_k/2) & (dist_z < 0):
                    output[line_CG1,4] = output[line_CG1,4] + k[0,2]
                    flag = 1
                    
        # The only way flag will be equal to 0 by the end of this for loop is if I ran it
        # entirely from the start without ever having to shift the nanoparticle, meaning I
        # have finished unwrapping it so that it is "inteira".
    # -----------------------------------------------------------------------------------
    # Once I finish the part below I can finally calculate the com of the mof without worrying
    # that it is broken across boundaries.
    numerator = np.array([[0.0, 0.0, 0.0]])
    denominator = 0
    for it_1 in range (0, len(output)):
        if (output[it_1,1] == 1):
            position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            numerator = numerator + M1*position
            denominator = denominator + M1
        if (output[it_1,1] == 2):
            position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            numerator = numerator + M2*position
            denominator = denominator + M2
    
    com_mof = numerator/denominator
    
    # -------------------------------------------------------------------------------
    # I think I will also already re-position the polymer superatoms so that they are
    # in the simulation box the closest possible to the center of mass of the MOF.
    for it_1 in range (0, NA):
        if output[it_1,1] > 2:
            atom_position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            dist_x = atom_position[0,0] - com_mof[0,0]
            dist_y = atom_position[0,1] - com_mof[0,1]
            dist_z = atom_position[0,2] - com_mof[0,2]
            # ------------------------------------------------------------------
            # Shifting to get the minimal distance (i.e., to put the polymer in the simulation
            # box domain that is the closest to the com of the MOF)
            if (abs(dist_x) > mod_i/2) & (dist_x > 0):
                atom_position[0,0] = atom_position[0,0] - i[0,0]
                output[it_1,2] = output[it_1,2] - i[0,0]
            if (abs(dist_x) > mod_i/2) & (dist_x < 0):
                atom_position[0,0] = atom_position[0,0] + i[0,0]
                output[it_1,2] = output[it_1,2] + i[0,0]
            if (abs(dist_y) > mod_j/2) & (dist_y > 0):
                atom_position[0,1] = atom_position[0,1] - j[0,1]
                output[it_1,3] = output[it_1,3] - j[0,1]
            if (abs(dist_y) > mod_j/2) & (dist_y < 0):
                atom_position[0,1] = atom_position[0,1] + j[0,1]
                output[it_1,3] = output[it_1,3] + j[0,1]
            if (abs(dist_z) > mod_k/2) & (dist_z > 0):
                atom_position[0,2] = atom_position[0,2] - k[0,2]
                output[it_1,4] = output[it_1,4] - k[0,2]
            if (abs(dist_z) > mod_k/2) & (dist_z < 0):
                atom_position[0,2] = atom_position[0,2] + k[0,2]
                output[it_1,4] = output[it_1,4] + k[0,2]
    
    # ---------------------------------------------------------------------------------
    # Now that I unwrapped the MOF nanoparticle, let me compute the position that is the
    # "center of analysis" for each vertex, edge and facet. I can do it based on the IDs
    # of the superatoms composing each vertex, as identified earlier (see the READ.me no-
    # te).
    group = np.zeros((8,3))
    for it_1 in range (0, len(output)):
        # ----------------------------------------------
        if output[it_1,0] == Ai[0,0]:
            group[0,:] = group[0,:] + output[it_1,2:5]
        if output[it_1,0] == Ai[0,1]:
            group[0,:] = group[0,:] + output[it_1,2:5]
        if output[it_1,0] == Ai[0,2]:
            group[0,:] = group[0,:] + output[it_1,2:5]
        # ----------------------------------------------
        if output[it_1,0] == Bi[0,0]:
            group[1,:] = group[1,:] + output[it_1,2:5]
        if output[it_1,0] == Bi[0,1]:
            group[1,:] = group[1,:] + output[it_1,2:5]
        if output[it_1,0] == Bi[0,2]:
            group[1,:] = group[1,:] + output[it_1,2:5]
        # ----------------------------------------------
        if output[it_1,0] == Ci[0,0]:
            group[2,:] = group[2,:] + output[it_1,2:5]
        if output[it_1,0] == Ci[0,1]:
            group[2,:] = group[2,:] + output[it_1,2:5]
        if output[it_1,0] == Ci[0,2]:
            group[2,:] = group[2,:] + output[it_1,2:5]
        # ----------------------------------------------
        if output[it_1,0] == Di[0,0]:
            group[3,:] = group[3,:] + output[it_1,2:5]
        if output[it_1,0] == Di[0,1]:
            group[3,:] = group[3,:] + output[it_1,2:5]
        if output[it_1,0] == Di[0,2]:
            group[3,:] = group[3,:] + output[it_1,2:5]
        # ----------------------------------------------
        if output[it_1,0] == Ei[0,0]:
            group[4,:] = group[4,:] + output[it_1,2:5]
        if output[it_1,0] == Ei[0,1]:
            group[4,:] = group[4,:] + output[it_1,2:5]
        if output[it_1,0] == Ei[0,2]:
            group[4,:] = group[4,:] + output[it_1,2:5]
        # ----------------------------------------------
        if output[it_1,0] == Fi[0,0]:
            group[5,:] = group[5,:] + output[it_1,2:5]
        if output[it_1,0] == Fi[0,1]:
            group[5,:] = group[5,:] + output[it_1,2:5]
        if output[it_1,0] == Fi[0,2]:
            group[5,:] = group[5,:] + output[it_1,2:5]
        # -----------------------------------------------
        if output[it_1,0] == Gi[0,0]:
            group[6,:] = group[6,:] + output[it_1,2:5]
        if output[it_1,0] == Gi[0,1]:
            group[6,:] = group[6,:] + output[it_1,2:5]
        if output[it_1,0] == Gi[0,2]:
            group[6,:] = group[6,:] + output[it_1,2:5]
        # -----------------------------------------------
        if output[it_1,0] == Hi[0,0]:
            group[7,:] = group[7,:] + output[it_1,2:5]
        if output[it_1,0] == Hi[0,1]:
            group[7,:] = group[7,:] + output[it_1,2:5]
        if output[it_1,0] == Hi[0,2]:
            group[7,:] = group[7,:] + output[it_1,2:5]
            
    # Now it suffices to calculate the position of the verteces in each group.
    group[:,:] = group[:,:]/3

    # -----------------------------------------------------------------
    # Let's now redefine the vertices based on their actual position.
    A = np.array([[group[0,0], group[0,1], group[0,2]]])
    B = np.array([[group[1,0], group[1,1], group[1,2]]])
    C = np.array([[group[2,0], group[2,1], group[2,2]]])
    D = np.array([[group[3,0], group[3,1], group[3,2]]])
    E = np.array([[group[4,0], group[4,1], group[4,2]]])
    F = np.array([[group[5,0], group[5,1], group[5,2]]])
    G = np.array([[group[6,0], group[6,1], group[6,2]]])
    H = np.array([[group[7,0], group[7,1], group[7,2]]])
    
    # ----------------------------------------------------------------
    # Now let me "point of analysis" for the edges (they are all equal), so
    # no need to split into groups
    AB = (B - A)/2 + A
    BC = (C - B)/2 + B
    CD = (D - C)/2 + C
    DA = (A - D)/2 + D
    EF = (F - E)/2 + E
    FG = (G - F)/2 + F
    GH = (H - G)/2 + G
    HE = (E - H)/2 + H
    DH = (H - D)/2 + D
    EA = (A - E)/2 + E
    FB = (B - F)/2 + F
    GC = (C - G)/2 + G
    
    grp_edges = np.concatenate((AB, BC, CD, DA, EF, FG, GH, HE, DH, EA, FB, GC), axis = 0)
    
    # ----------------------------------------------------------------
    # Finally, I will define the point of analysis for the facets.
    # For this I will use the bisetriz of  diagonal (could use any of
    # the two: makes no difference since these two points
    # are the same).
    HA = (A - H)/2 + H
    EB = (B - E)/2 + E
    FC = (C - F)/2 + F
    GD = (D - G)/2 + G
    BD = (D - B)/2 + B
    FH = (H - F)/2 + F
    
    grp_facets = np.concatenate((HA, EB, FC, GD, BD, FH), axis = 0)
    
    # ----------------------------------------------------------------
    # Finally, let's compute histograms! I think that things will go a little
    # bit faster if I can identify the superatoms that actually matter.
    CG_ZIF = np.zeros((1,5))
    CG_PVDF = np.zeros((1,5))
    
    all_points = np.concatenate((group, grp_edges, grp_facets), axis = 0)
    
    temporary_array = np.zeros((1,5))
    for it_1 in range (0, NA):
        position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
        for it_2 in range (0, len(all_points)):
            dist = ((position[0,0] - all_points[it_2,0])**2 + (position[0,1] - all_points[it_2,1])**2 + (position[0,2] - all_points[it_2,2])**2)**(1/2)
            if dist < 12:
                if output[it_1,1] <= 2:
                    temporary_array[0,:] = output[it_1,:]
                    CG_ZIF = np.concatenate((CG_ZIF, temporary_array), axis = 0)
                if output[it_1,1] > 2:
                    temporary_array[0,:] = output[it_1,:]
                    CG_PVDF = np.concatenate((CG_PVDF, temporary_array), axis = 0)
                break
    
    CG_ZIF = np.delete(CG_ZIF, 0, 0)
    CG_PVDF = np.delete(CG_PVDF, 0, 0)
    
    # I will define this CG_all below which will be useful for a neat code of density
    # profile.
    CG_all = np.concatenate((CG_ZIF, CG_PVDF), axis = 0)
    
    # Let's order the array CG_ZIF with ascending order of superatom type.
    for it_1 in range (0, len(CG_ZIF)):
        for it_2 in range (it_1 + 1, len(CG_ZIF)):
            if CG_ZIF[it_1,1] > CG_ZIF[it_2,1]:
                temporary_array[0,:] = CG_ZIF[it_1,:]
                CG_ZIF[it_1,:] = CG_ZIF[it_2,:]
                CG_ZIF[it_2,:] = temporary_array[0,:]
    
    N1_CG = 0
    N2_CG = 0
    # Let's count the superatoms of eah type
    for it_1 in range (0, len(CG_ZIF)):
        if CG_ZIF[it_1,1] == 1:
            N1_CG = N1_CG + 1
        if CG_ZIF[it_1,1] == 2:
            N2_CG = N2_CG + 1
    
    # Building the histogram for the verteces 
    for it_1 in range (0, len(group)):
        reference = np.array([[group[it_1,0], group[it_1,1], group[it_1,2]]])
        for it_2 in range (0, len(CG_PVDF)):
            position = np.array([[CG_PVDF[it_2,2], CG_PVDF[it_2,3], CG_PVDF[it_2,4]]])
            dist = ((position[0,0] - reference[0,0])**2 + (position[0,1] - reference[0,1])**2 + (position[0,2] - reference[0,2])**2)**(1/2)
            if dist < 12:
                vector_u = np.array([[com_mof[0,0] - reference[0,0], com_mof[0,1] - reference[0,1], com_mof[0,2] - reference[0,2]]])
                vector_v = np.array([[position[0,0] - reference[0,0], position[0,1] - reference[0,1], position[0,2] - reference[0,2]]])
                mod_u = (vector_u[0,0]**2 + vector_u[0,1]**2 + vector_u[0,2]**2)**(1/2)
                mod_v = (vector_v[0,0]**2 + vector_v[0,1]**2 + vector_v[0,2]**2)**(1/2)
                cos_value = (vector_u[0,0]*vector_v[0,0] + vector_u[0,1]*vector_v[0,1] + vector_u[0,2]*vector_v[0,2])/(mod_u*mod_v)
                # Angle value in degrees
                angle_value = math.acos(cos_value)*(180/math.pi)
                
                x_value = int(dist/bin_size)
                y_value = int(angle_value/bin_size_angle)
                # ------------------------------------------------
                if y_value == number_of_bins_angle:
                    y_value = number_of_bins_angle - 1
                # ------------------------------------------------
                rho_grp[x_value,y_value] = rho_grp[x_value,y_value] + 1
                
    # Building the histogram for the edges
    for it_1 in range (0, len(grp_edges)):
        reference = np.array([[grp_edges[it_1,0], grp_edges[it_1,1], grp_edges[it_1,2]]])
        for it_2 in range (0, len(CG_PVDF)):
            position = np.array([[CG_PVDF[it_2,2], CG_PVDF[it_2,3], CG_PVDF[it_2,4]]])
            dist = ((position[0,0] - reference[0,0])**2 + (position[0,1] - reference[0,1])**2 + (position[0,2] - reference[0,2])**2)**(1/2)
            if dist < 12:
                vector_u = np.array([[com_mof[0,0] - reference[0,0], com_mof[0,1] - reference[0,1], com_mof[0,2] - reference[0,2]]])
                vector_v = np.array([[position[0,0] - reference[0,0], position[0,1] - reference[0,1], position[0,2] - reference[0,2]]])
                mod_u = (vector_u[0,0]**2 + vector_u[0,1]**2 + vector_u[0,2]**2)**(1/2)
                mod_v = (vector_v[0,0]**2 + vector_v[0,1]**2 + vector_v[0,2]**2)**(1/2)
                cos_value = (vector_u[0,0]*vector_v[0,0] + vector_u[0,1]*vector_v[0,1] + vector_u[0,2]*vector_v[0,2])/(mod_u*mod_v)
                # Angle value in degrees
                angle_value = math.acos(cos_value)*(180/math.pi)
                
                x_value = int(dist/bin_size)
                y_value = int(angle_value/bin_size_angle)
                # ------------------------------------------------
                if y_value == number_of_bins_angle:
                    y_value = number_of_bins_angle - 1
                # ------------------------------------------------
                rho_edges[x_value,y_value] = rho_edges[x_value,y_value] + 1
                
    # Building the histogram for the facets
    for it_1 in range (0, len(grp_facets)):
        reference = np.array([[grp_facets[it_1,0], grp_facets[it_1,1], grp_facets[it_1,2]]])
        for it_2 in range (0, len(CG_PVDF)):
            position = np.array([[CG_PVDF[it_2,2], CG_PVDF[it_2,3], CG_PVDF[it_2,4]]])
            dist = ((position[0,0] - reference[0,0])**2 + (position[0,1] - reference[0,1])**2 + (position[0,2] - reference[0,2])**2)**(1/2)
            if dist < 12:
                vector_u = np.array([[com_mof[0,0] - reference[0,0], com_mof[0,1] - reference[0,1], com_mof[0,2] - reference[0,2]]])
                vector_v = np.array([[position[0,0] - reference[0,0], position[0,1] - reference[0,1], position[0,2] - reference[0,2]]])
                mod_u = (vector_u[0,0]**2 + vector_u[0,1]**2 + vector_u[0,2]**2)**(1/2)
                mod_v = (vector_v[0,0]**2 + vector_v[0,1]**2 + vector_v[0,2]**2)**(1/2)
                cos_value = (vector_u[0,0]*vector_v[0,0] + vector_u[0,1]*vector_v[0,1] + vector_u[0,2]*vector_v[0,2])/(mod_u*mod_v)
                # Angle value in degrees
                angle_value = math.acos(cos_value)*(180/math.pi)
                
                x_value = int(dist/bin_size)
                y_value = int(angle_value/bin_size_angle)
                # ------------------------------------------------
                if y_value == number_of_bins_angle:
                    y_value = number_of_bins_angle - 1
                # ------------------------------------------------
                rho_facets[x_value,y_value] = rho_facets[x_value,y_value] + 1
                
rho_grp = rho_grp/number_of_configurations_considered
rho_edges = rho_edges/number_of_configurations_considered
rho_facets = rho_facets/number_of_configurations_considered

ofi6 = open("rho_grp.dat", 'w')   
for it_1 in range(len(rho_grp)):
    for it_2 in range(0, number_of_bins_angle):
            ofi6.write(str(rho_grp[it_1,it_2]))
            ofi6.write('\t')
    ofi6.write('\n')
ofi6.close()

ofi6 = open("rho_edges.dat", 'w')   
for it_1 in range(len(rho_edges)):
    for it_2 in range(0, number_of_bins_angle):
            ofi6.write(str(rho_edges[it_1,it_2]))
            ofi6.write('\t')
    ofi6.write('\n')
ofi6.close()

ofi6 = open("rho_faces.dat", 'w')   
for it_1 in range(len(rho_facets)):
    for it_2 in range(0, number_of_bins_angle):
            ofi6.write(str(rho_facets[it_1,it_2]))
            ofi6.write('\t')
    ofi6.write('\n')
ofi6.close()
