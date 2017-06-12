import math
import random
import re
import numpy as np
#epsilon = 0.039 #kcal/mol
#sigma   = 2.934 #angstrom
epsilon = 4.32
rmin   = 2.88358
r = 2.88358
step = 0.01
f1 = open("test.xyz","w+")
array_sites = []
initial_coord = []
with open("Au.pdb", "r") as file:
	array_pdb = []
	count = 0
	for line in file:
		coord = []
		line = re.sub( '\s+', ' ',  line).strip()
		temp_array = line.split(" ")
		if(len(temp_array) >= 6 and count!=0): 


			coord.append(float(temp_array[6]))
			coord.append(float(temp_array[7]))
			coord.append(float(temp_array[8]))

			array_pdb.append(coord)
        
           #print coord
		count=count+1


def find_sites(array_pdb):
	global array_sites
	global initial_coord
	global r
	array_sites = []
	x = []
	y = []
	z = []

	for i,atom in enumerate(array_pdb):
		x.append(atom[0])
		y.append(atom[1])
		z.append(atom[2])


	count = len(x)


	x1 = np.zeros((count,350))
	y1 = np.zeros((count,350))
	z1 = np.zeros((count,350))
	count = 0
	for i,j in np.ndindex(x1.shape):
		
  
 		u = random.uniform(0,1)
 		v = random.uniform(0,1)
		theta = 2*math.pi*u
		phi = math.acos(2*v - 1)
		x_coor = math.cos(theta)*math.sin(phi)
		y_coor = math.sin(theta)*math.sin(phi)
		z_coor = math.cos(phi)
		x1[i][j] = x[i] + r*x_coor
		y1[i][j] = y[i] + r*y_coor
		z1[i][j] = z[i] + r*z_coor

	initial_coord = [] #inittial coordinates for steepest descent

	for i,j in np.ndindex(x1.shape):
		count = 0	
		for m in range(len(x)):
			val = math.sqrt((math.pow((x1[i][j]-x[m]),2) + math.pow((y1[i][j]-y[m]),2) + math.pow((z1[i][j]-z[m]),2))) - r

			if(val >= 0):
				count = count + 1
			else:
				break

		if(count == len(x)):
			temp = [x1[i][j],y1[i][j],z1[i][j]]
			initial_coord.append(temp)


	#print initial_coord
			
	for i , val in enumerate(initial_coord):
    #print val
		x_site = val[0]
		y_site = val[1]
		z_site = val[2]
		for z in range(1000):
			fx_sum = 0.0
			fy_sum = 0.0
			fz_sum = 0.0
			for j , val1 in enumerate(array_pdb):
				x_pdb = val1[0]
				y_pdb = val1[1]
				z_pdb = val1[2]
				r = math.sqrt((x_site - x_pdb)*(x_site - x_pdb) + (y_site - y_pdb)*(y_site - y_pdb) + (z_site - z_pdb)*(z_site - z_pdb))
				#force = -4*epsilon*(-12.00*(math.pow((sigma/dis),12))*(1.00/dis) + 6.00*(math.pow((sigma/dis),6))*(1.00/dis))
				force = 12.00*epsilon*(-1.00*(math.pow((rmin/r),12))*1.00/r + (math.pow((rmin/r),6))*1.00/r)
				fx_sum = fx_sum + force*((x_site - x_pdb)/r)
				fy_sum = fy_sum + force*((y_site - y_pdb)/r)
				fz_sum = fz_sum + force*((z_site - z_pdb)/r)

			force = math.sqrt(fx_sum*fx_sum + fy_sum*fy_sum + fz_sum*fz_sum)
			
			if(abs(force) < 0.001):
				#print force
				coord = []
				coord = [x_site,y_site,z_site]
				#print coord
				array_sites.append(coord)
				break
	 		x_site = x_site - step*fx_sum
			y_site = y_site - step*fy_sum
			z_site = z_site - step*fz_sum

	count = 0
	remove_sites = {} # dictionary which contains index to remove sites
	for i in range(len(array_sites)):
		remove_sites[i] = 0
	for i in range(len(array_sites)): 
		if remove_sites[i]:
			continue   
    	for j in range(len(array_sites)):
			if(j > i and remove_sites[j] == 0):
				x = array_sites[i][0] - array_sites[j][0]
				y = array_sites[i][1] - array_sites[j][1]
				z = array_sites[i][2] - array_sites[j][2] 
				r = math.sqrt(x*x + y*y + z*z)
				if(r < 0.1):
					count = count + 1
					remove_sites[j] = 1

	print remove_sites				
		
	count = 0
	for key,val in remove_sites.iteritems():
		if(val == 1):
			del array_sites[key-count]
		count = count + 1


		#return array_sites
		

find_sites(array_pdb)


s = str(len(array_sites)-1)
f1.write(s)
f1.write("\n")
f1.write("\n")


for i,val in enumerate(array_sites):
	s = ""
	s = "Au" + "\t" + str(val[0]) + "\t" + str(val[1]) + "\t" + str(val[2])
	f1.write(s)
	f1.write("\n")