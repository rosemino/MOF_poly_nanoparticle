import numpy as np

NA = 68592

number_of_configurations = 200

C = np.zeros((1,5))
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
    C = np.concatenate((C, add), axis = 0)
 
C = np.delete(C,0, axis=0)
box_size = np.delete(box_size, 0, 0)

# -------------------------------------------------------------------------------------
# I will need to read this for later. It was built from simply copying pasting the bond
# section of the initial file used to start the CG-MD simulation of the nanoparticle.
number_of_bonds = 5952

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
# I will set up a value of more or less 0.5 angs to the radial bin size. For this, I will cacu-
# late the average_lz and find out more or less the required number of bins. Evidently, the 
# density_profile will have a fixed number of bins, but since the dynamics is carried out
# in the npt ensemble and therefore each configuration reached has a different volume, a
# given bin will always have a different size. Hopefully, assuming the volume doesnt fluctuate
# wildly, this difference should be slight. The bins will always be identified with a same
# ID.
# Note that suffices to calculate the average only of one of the box dimensions since the
# other two are exactly teh asme in each and all configurations (and therefore the average
# is also the same).

numerator = 0
for it_1 in range (0, number_of_configurations):
    numerator = numerator + box_size[int((it_1*3) + 2),1] - box_size[int((it_1*3) + 2),0]

average_lz = numerator/number_of_configurations

number_of_bins = int((average_lz/2 - 3)/0.5)
# Actual bin size (should not be so different than 0.5)
bin_size = (average_lz/2 - 3)/number_of_bins

# ------------------------------------------------------------------------------------
# I will build the r_axis to be centered in the middle of the unit length the given
# bin concerns. Data will be acumulated in each bin as long as it is smaller than than the
# upper value that marks its range and higher than the lower value that marks its range.
density_profile_r_axis = np.zeros((1, number_of_bins))
for it_1 in range (0,number_of_bins):
    density_profile_r_axis[0,it_1] = it_1*bin_size + bin_size/2
    
density_profile_mof = np.zeros((1, number_of_bins))
density_profile_polymer = np.zeros((1, number_of_bins))
# This array is because I want to find out the average density in g/cm3 later,
# and this requires me actually counting the mass instead of atomic density (with
# the latter I will have no control afterwards in the amount of mass (no way to determine
# it)).
mass_polymer = np.zeros((1, number_of_bins))

M1 = 65.38
M2 = 81.0

counter_atom_mof = 0
counter_atom_polymer = 0

number_of_configurations_considered = 200

for loop_0 in range (0, number_of_configurations_considered):
    
    output = np.zeros((NA,5))
    output[:,:] = C[int((loop_0)*NA):int((loop_0 + 1)*NA),:]
    
    i = np.zeros((1,3))
    i[0,0] = box_size[int(loop_0*3),1] - box_size[int(loop_0*3),0]
    mod_i = ((i[0,0])**2 + (i[0,1])**2 + (i[0,2])**2)**(1/2)
    j = np.zeros((1,3))
    j[0,1] = box_size[int((loop_0*3) + 1),1] - box_size[int((loop_0*3) + 1),0]
    mod_j = ((j[0,0])**2 + (j[0,1])**2 + (j[0,2])**2)**(1/2)
    k = np.zeros((1,3))
    k[0,2] = box_size[int((loop_0*3) + 2),1] - box_size[int((loop_0*3) + 2),0]
    mod_k = ((k[0,0])**2 + (k[0,1])**2 + (k[0,2])**2)**(1/2)
    
    # Here I am using k but could be i or j since the box is cubic.
    length_r = (k[0,2]/2 - 3)
    microstate_bin_size = length_r/number_of_bins
    
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
    
    # # I will now translate the entire simulation box so that com_mof is in the very
    # # center of the simulation box.
    # center = np.array([[i[0,0]/2, j[0,1]/2, k[0,2]/2]])
    # translator = center - com_mof
    
    # for it_1 in range (0, len(output)):
    #     output[it_1,2] = output[it_1,2] + translator[0,0]
    #     output[it_1,3] = output[it_1,3] + translator[0,1]
    #     output[it_1,4] = output[it_1,4] + translator[0,2]
    # # The new com is now the center of the simulation box also necessarily.
    
    # --------------------------------------------------------------------------
    # Finally, now it suffices to do the desired calculation of density.
    
    for it_1 in range (0, NA):
        atom_ID = output[it_1,0]
        # ----------------------------------------------------------------------
        # Checking quickly the mass (only necessary for PVDF)
        atom_type = output[it_1,1]
        if (atom_type == 3):
            atom_mass = 128.0
        if (atom_type == 4):
            atom_mass = 65.0
        if (atom_type == 5):
            atom_mass = 83.0
        if (atom_type == 6):
            atom_mass = 64.0
        # ------------------------------------------------------------------------
        atom_position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
        dist_x = atom_position[0,0] - com_mof[0,0]
        dist_y = atom_position[0,1] - com_mof[0,1]
        dist_z = atom_position[0,2] - com_mof[0,2]
        
        # Shifting to get the minimal distance (i.e., to put the polymer in the simulation
        # box domain that is the closest to the com of the MOF)
        if (abs(dist_x) > mod_i/2) & (dist_x > 0):
            atom_position[0,0] = atom_position[0,0] - i[0,0]
        if (abs(dist_x) > mod_i/2) & (dist_x < 0):
            atom_position[0,0] = atom_position[0,0] + i[0,0]
        if (abs(dist_y) > mod_j/2) & (dist_y > 0):
            atom_position[0,1] = atom_position[0,1] - j[0,1]
        if (abs(dist_y) > mod_j/2) & (dist_y < 0):
            atom_position[0,1] = atom_position[0,1] + j[0,1]
        if (abs(dist_z) > mod_k/2) & (dist_z > 0):
            atom_position[0,2] = atom_position[0,2] - k[0,2]
        if (abs(dist_z) > mod_k/2) & (dist_z < 0):
            atom_position[0,2] = atom_position[0,2] + k[0,2]
        
        r_coordinate = ((atom_position[0,0] - com_mof[0,0])**2 + (atom_position[0,1] - com_mof[0,1])**2 + (atom_position[0,2] - com_mof[0,2])**2)**(1/2)
        corresponding_bin = int(r_coordinate/microstate_bin_size)
        # -----------------------------------------------------------
        # The condition below is quite possible and means that this atom is not part of the
        # domain I am investigating.
        # I am dismissing the possibility of being right in the edge of the bin also (i.e.,
        # corresponding_bin == number_of_bins as this would be a different treatment compa-
        # red to the other bins, which encapsulate only the bin lower edge and I noticed a
        # very high count if I treat the last bin differently (probably here it makes a lot
        # of difference because the box is big and thus I hve high prob. of having a larger
        # count due to this mere detail)).
        if corresponding_bin >= number_of_bins:
            continue
        # -----------------------------------------------------------
        if output[it_1,1] <= 2:
            density_profile_mof[0, corresponding_bin] = density_profile_mof[0, corresponding_bin] + 1
            counter_atom_mof = counter_atom_mof + 1
        if output[it_1,1] > 2:
            density_profile_polymer[0, corresponding_bin] = density_profile_polymer[0, corresponding_bin] + 1
            mass_polymer[0, corresponding_bin] = mass_polymer[0, corresponding_bin] + atom_mass
            counter_atom_polymer = counter_atom_polymer + 1
    
density_profile_r_axis = density_profile_r_axis.transpose()
density_profile_mof = density_profile_mof.transpose()
density_profile_polymer = density_profile_polymer.transpose()
mass_polymer = mass_polymer.transpose()

for it_1 in range (0, number_of_bins):
     volume = (4/3)*3.1415*((bin_size*(it_1+1))**3 - (bin_size*it_1)**3)
     density_profile_mof[it_1,0] = density_profile_mof[it_1,0]/volume
     density_profile_polymer[it_1,0] = density_profile_polymer[it_1,0]/volume
     mass_polymer[it_1,0] = mass_polymer[it_1,0]/volume

density_profile_mof = density_profile_mof/number_of_configurations_considered
density_profile_polymer = density_profile_polymer/number_of_configurations_considered
mass_polymer = mass_polymer/number_of_configurations_considered

output_data = np.concatenate((density_profile_r_axis, density_profile_mof, density_profile_polymer, mass_polymer), axis = 1)

ofi = open("density_profile_r.dat", 'w')   
for it_1 in range(len(output_data)):
    for it_2 in range(4):
            ofi.write(str(output_data[it_1,it_2]))
            ofi.write('\t')
    ofi.write('\n')

information = np.array([[bin_size, counter_atom_mof, counter_atom_polymer]])
ending = np.savetxt('information.txt', (information))
