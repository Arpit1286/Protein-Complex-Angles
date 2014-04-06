"""
Author: Arpit Tandon
Contributer:  Raul Mendez-Giraldez
Calculate the angle between the planes of the principal axes vector
"""
import __main__
__main__.pymol_argv = ['pymol'] #quiet and no GUI __main__ is now pymol

import numpy
#import os, sys
import pymol
pymol.finish_launching()


# user input for structures
#dock_path = os.path.abspath(sys.argv[1])
#dock_name = dock_path.split('/')[-1].split('.')[0]
#actin_vt_path = os.path.abspath(sys.argv[2])
#actin_vt_name = actin_vt_path.split('/')[-1].split('.')[0]

#dock_file = open('model.000.00.pdb','r')
#actin_file = open('00095_Pradeep.pdb','r')

pymol.cmd.load('model.000.00.pdb')
pymol.cmd.load('00095_Pradeep.pdb')

pymol.cmd.copy('Pradeep_2','00095_Pradeep')

pymol.cmd.align('00095_Pradeep','rec' )
pymol.cmd.align('Pradeep_2','lig.000.00')

pymol.cmd.save('p1.pdb','00095_Pradeep and (not chain A)')
pymol.cmd.save('p2.pdb','Pradeep_2 and (not chain A)')
# quit pymol
pymol.cmd.quit()


#using the principal_axes script from Pierre Poulain http://bit.ly/Px695Q 

center_list =[]
axis1_list =[]

def read_pdb_xyz(pdb_name):
    """
reads atomic coordinates of C-alpha atoms in a .pdb file
returns:
[[x1 y1 z1]
 [x2 y2 z2]
 [.. .. ..] 
 [xn yn zn]]
    """
    xyz = []
    pdb_file = open(pdb_name, 'r')
    for line in pdb_file:
        if line.startswith("ATOM"):
            # extract x, y, z coordinates for carbon alpha atoms
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            if line[12:16].strip() == "CA":
                xyz.append([x, y, z])
    pdb_file.close()
    return xyz 

pdb_name = ['protien1.pdb','protien2.pdb']

for i in pdb_name:
    xyz = read_pdb_xyz(i)
    print "%d CA atomes found if %s" %(len(xyz), i)

#create coordinates array
    coord = numpy.array(xyz, float)

# compute geometric center
    center = numpy.mean(coord, 0)
    print "Coordinates of the geometric center:\n", center
    center_list.append(center)
# center with geometric center
    coord = coord - center

# compute principal axis matrix
    inertia = numpy.dot(coord.transpose(), coord)
    e_values, e_vectors = numpy.linalg.eig(inertia)
# warning eigen values are not necessary ordered!
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html
    print "(Unordered) eigen values:"
    print e_values
    print "(Unordered) eigen vectors:"
    print e_vectors

#--------------------------------------------------------------------------
# order eigen values (and eigen vectors)
#
# axis1 is the principal axis with the biggest eigen value (eval1)
# axis2 is the principal axis with the second biggest eigen value (eval2)
# axis3 is the principal axis with the smallest eigen value (eval3)
#--------------------------------------------------------------------------
    for i in xrange(len(e_values)):
	# find biggest eigen value : This is the axis we need
         if e_values[i] == max(e_values):
             eval1 = e_values[i]
             axis1 = e_vectors[:,i]
             axis1_list.append(axis1)
	# find smallest eigen value
         elif e_values[i] == min(e_values):
             eval3 = e_values[i]
             axis3 = e_vectors[:,i]
	# middle eigen value
         else:
            eval2 = e_values[i]
            axis2 = e_vectors[:,i]


# Vector between the origin, take care of directions
c_vector1 = center_list[0] - center_list[1] #betwenn p1_axis1 and center
c_vector2 = center_list[1] - center_list[0] #between p2_axis1 and center, points in reverse direction

d1 = numpy.cross(axis1_list[0],c_vector1) 
d2 = numpy.cross(axis1_list[1],c_vector2)

# calculate the unit vector

def unit_vector(vector):
    unit_vector = vector / numpy.linalg.norm(vector)
    return unit_vector
    
# convert radian to degrees    
def angle_deg(angle_rad):
    angle_deg = numpy.degrees(angle_rad)
    return angle_deg

# calculate the angle
vector_1 = unit_vector(d1)
vector_2 = unit_vector(d2)
vector_d1 = unit_vector(axis1_list[0])
vector_d2 = unit_vector(axis1_list[1])

angle_axis1 = numpy.arccos(numpy.dot(vector_d1,vector_d2))
angle_planes = numpy.arccos(numpy.dot(vector_1,vector_2))
angle_planes_deg = angle_deg(angle_planes)
angle_axis1_deg = angle_deg(angle_axis1)

print "the angle between the two planes is %8.3f" % (180-angle_planes_deg)
print "the angle between the principal vectors of two actins is %8.3f" % (angle_axis1_deg)
