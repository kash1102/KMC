import re
import math
#from sortedcontainers import SortedDict
import random
import numpy as np
r = 2.88358

no_of_atoms = 108
radius_of_atom = 1.4417907268393706
# defining constants
Avagadro = 6.022*math.pow(10,23)
nearest_neighbour_radius = 2.8500419592663904
epsilon = 0.039 #kcal/mol
sigma   = 2.934 #angstrom
k_ads = 0.272
k_red  = 2.822 * math.pow(10,-3)
Au_plus = 12.5*math.pow(10,-3)
Reducing_agent = 10*12.5*math.pow(10,-3)
volume_of_solution = 8.329*math.pow(10,5)*math.exp(-24)
N_plus = (Au_plus*volume_of_solution)*Avagadro
Au_zero = 0
R =  0.00198588  #universal gas constant
A = 40*math.exp(12) # pre-exponential factor

T = 296.15
k_b = 0.002
#k_diff = A*math.exp(-1*E_adiff/(R*T))


#r_diff = k_diff*math.exp((-E_adiff - (E_i - E_f))/(2*k_b*T))
#r_des = k_red*math.exp((-E_adiff - E_i)/(2*k_b*T))
#r_ads = k_ads*Au_zero

ads = {} # rate mapping to site at which adsorption will take place
des = {} # rate mapping to coordinate at which desorption will take place
diff_site = {} # rate  mapping where after surface diffusion atom will move
diff_atom = {} # rate  mapping on which atom surface diffusion would occur
event_list = {} # on this rate wheter adsorption , desorption or surface diffusion would occur 
rate_list = {} # index of rate 0,1,2,3......
sites_hash = {}
atoms_hash = {}

t = 0.0 # global variable time
r_ads = 0.0 # initial rate adsorption
r_des = 0.0 # initial rate desorption
r_diff = 0.0 # initial rate diffusion
r_total = 0.0
rate_count = 0 # possible rates 0,1,2,3,..... for all the possible transitions
rate_list[rate_count] = 0
rate_count = rate_count + 1

file1 = open("time.txt","a+")
file2 = open("diameter.txt","a+") 

with open("Au.pdb", "r") as file: # all the coordinates of atom stored in array_pdb
    array_pdb = []
    count = 0
    for line in file:
    	coord = []
    	line = re.sub( '\s+', ' ',  line).strip()
    	temp_array = line.split(" ")
        if(len(temp_array) >= 6 and count!=0): 
    	   #print temp_array[6]
    	   
    	   coord.append(float(temp_array[6]))
           coord.append(float(temp_array[7]))
           coord.append(float(temp_array[8]))
           #atoms_hash[coord]=0

           array_pdb.append(coord)
        
           #print coord
        count=count+1

with open("step_descent.xyz","r") as file: # all the coordinates of sites are stored in array_sites
    array_sites = []
    for line in file:
        coord = []
        line = re.sub( '\s+', ' ',  line).strip()
        temp_array = line.split(" ")
        if(len(temp_array) >= 4):
        	coord.append(float(temp_array[1]))
        	coord.append(float(temp_array[2]))
        	coord.append(float(temp_array[3]))
        
        #sites_hash[temp_array]=0
       		array_sites.append(coord)

# lennard jones potentiol

def Concentration(step):
	global Reducing_agent
	global Au_plus
	global Au_zero
	Au_plus = ((N_plus/Avagadro)/volume_of_solution)
	Reducing_agent = Reducing_agent + (-1*k_red*Au_plus*Reducing_agent)*step
	Au_plus = Au_plus + (-1*k_red*Au_plus*Reducing_agent)*step
	Au_zero = Au_zero + k_red*Au_plus*Reducing_agent*step + (-r_ads*step + r_des*step)

def Calculate_Diameter(number_of_atoms):
	x_c = 0.0 
	y_c = 0.0
	z_c = 0.0
	Radius = 0.0
	for i,val in enumerate(array_pdb):
		x_c = x_c + val[0]
		y_c = y_c + val[1]
		z_c = z_c + val[2]
	x_c = (x_c)/len(array_pdb)
	y_c = (y_c)/len(array_pdb)
	z_c = (z_c)/len(array_pdb)
	c_o_m = [x_c,y_c,z_c]
	for i,val in enumerate(array_pdb):
		d = distance(val,c_o_m)
		d = d*d
		Radius = Radius + d
	Radius = Radius/(len(array_pdb))
	Radius = math.sqrt(Radius)	
	return 2*math.sqrt(5/3)*Radius

def surface_atoms():
	x = []
	y = []
	z = []
	list_of_surface_atoms = []
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

	#print x1,y1,z1
	
	for i,j in np.ndindex(x1.shape):
		if(j==0):
			flag_surface = 0
		
		count = 0
		
		for m in range(len(x)):
			val = math.sqrt((math.pow((x1[i][j]-x[m]),2) + math.pow((y1[i][j]-y[m]),2) + math.pow((z1[i][j]-z[m]),2))) - r

			if(val >= 0):
				count = count + 1
      			#print count

			else:
				count = 0
      			#print count
				break
      		
			#print count
		if(count == len(x)):
			flag_surface = flag_surface + 1
			#print flag_surface
  		
		if(j==349 and flag_surface != 0):
  			
  			temp = []
			temp.append(x[i])
			temp.append(y[i])
			temp.append(z[i])
			list_of_surface_atoms.append(temp)
			"""
			s = ""
			s = "Au" + "\t" + str(x[i]) + "\t" + str(y[i]) + "\t" + str(z[i])
			f1.write(s)
			f1.write("\n")

			"""
    
	return list_of_surface_atoms  			
	
		
def potential(r):
	v = 4*epsilon*(math.pow((sigma/r),12) - math.pow((sigma/r),6))
	#if(v > 60.0):
	#	print r
	return v

# distance between 2 coordinates point
def distance(r1,r2):
	x = r1[0] - r2[0]
	y = r1[1] - r2[1]
	z = r1[2] - r2[2]
	return math.sqrt(x*x + y*y + z*z)

# potential energy at particular position
def energy(pos):
	E_initial = 0.0
	for i, val in enumerate(array_pdb):
		r = distance(pos,val)
		if(r > 1.4):
			
			E_initial = E_initial + potential(r)

	return (E_initial)


def energy_neighbour(curr_pos):
	nbours = []
	neighbours_count = 0
	dist = []
	for i, atom in enumerate(array_pdb):
		r = distance(atom,curr_pos)
		if(r > 0.1):
			dist.append(r)

	dist.sort()		
	start = dist[0]
	
	for i,val in enumerate(dist):	
		if(val/start > 1.1):
			break
		else:
			neighbours_count = neighbours_count + 1
	
	print neighbours_count		
	return neighbours_count	


def nearest_neighbour(curr_pos):
	nbours = []
	dist = {}
	for i,site in enumerate(array_sites):
		r = distance(site,curr_pos)
		dist[r] = i
	
	keylist = dist.keys()
	keylist.sort()
	
	#dist = SortedDict(dist)
	#start = list(dist)[0]
	if(len(keylist) > 0):
		start = keylist[0]
  #print dist,
  	"""
	for i,val in dist.iteritems():	
		if(i/start > 1.1):
			break
		else:
			nbours.append(val)
	"""
	for i,val in enumerate(keylist):
		if (val/start > 1.1):
			break
		else:
			nbours.append(dist[val])
	#print len(nbours),"neighbours"
	return nbours

def adsorption():
	global r_ads
	global ads
	global event_list
	global rate_list
	global rate_count
	global Au_zero
	r_ads = 0
	ads.clear()
	event_list.clear()
	rate_count = 1
	rate_list.clear()
	rate_list[0] = 0
	#Au_zero = ((N_zero/Avagadro)/volume_of_solution)
	#Concentration()
	for i , site in enumerate(array_sites):
		r_ads = r_ads + k_ads*Au_zero
		#print r_ads
		ads[r_ads] = site  # basically on this rate, adsorption on this particular site 
		event_list[r_ads] = 'adsorption'
		rate_list[rate_count] = r_ads
		rate_count=rate_count+1
	
	print r_ads
	print "adsorption done"
	#event_list = SortedDict(event_list)
	#print event_list

def desorption():
	global event_list
	global rate_list
	global rate_count
	global r_des
	global des
	E_adiff = 42
	#E_adiff = 0
	des.clear()
	r_des = r_ads
	list_of_surface_atoms = surface_atoms()
	#k_diff = A*math.exp(-1*(E_adiff)/(R*T))
	#print k_diff
	for j, atom in enumerate(list_of_surface_atoms):
		E_initial = energy(atom)
		#print E_initial
		#E_initial  = 2.303*energy_neighbour(atom)
		r_des = r_des + math.pow(10,13)*math.exp((-1*E_adiff - E_initial)/(2*k_b*T))
		#print r_des
		des[r_des] = atom # basically on this rate, desorption of this particular atom
		event_list[r_des]='desorption'
		rate_list[rate_count] = r_des
		rate_count = rate_count + 1 

	print r_des	
	print "desorption done"

def surface_diffusion():
	global event_list
	global rate_listth
	global rate_count
	global r_diff
	global diff_site
	global diff_atom
	#E_adiff = -2*epsilon
	E_adiff = 42
	diff_site.clear()
	diff_atom.clear()
	list_of_surface_atoms = surface_atoms()
	r_diff = r_des
	for i, atom in enumerate(list_of_surface_atoms):
		nbors = [] # finding nearest neighbours where surface diffusion could occur
		nbors = nearest_neighbour(atom)
		
		E_initial = energy(atom)
		#E_initial = 2.306*energy_neighbour(atom)

		for i,site in enumerate(nbors):

			E_final = energy(array_sites[site])
			"""
			E_final = 4*energy_neighbour(array_sites[site])
			if(E_final > E_initial):
				continue
			temp = math.pow(10,13)*math.exp((-1*E_adiff - (E_initial - E_final))/(2*k_b*T))
			if(temp > 60):
				print E_initial,E_final
			"""
			#E_final = 0
			r_diff = r_diff + math.pow(10,13)*math.exp((-E_adiff -(E_initial - E_final))/(2*k_b*T))
			#print r_diff
			diff_site[r_diff] = array_sites[site]
			diff_atom[r_diff] = atom
			event_list[r_diff] = 'difusion'
			rate_list[rate_count] = r_diff
			rate_count = rate_count + 1 

	print r_diff		
	print "surface diffusion done"		


def find_event():
	global r_total
	r_total = r_diff
	bot = 0
	top = len(rate_list)
	random_num = random.uniform(0,r_total) # now finding range where it lies
	#print random_num
	#mid = 1
	while(bot <= top):
		mid = (bot + top)/2
		if(rate_list[mid] <= random_num and rate_list[mid+1] >= random_num ):
			break
		elif (rate_list[mid] > random_num):
			top = mid - 1
		else:
			bot = mid + 1
	#print mid		
	return rate_list[mid+1] 

def main():
	"""
	ads = {} # rate mapping to site at which adsorption will take place
	des = {} # rate mapping to coordinate at which desorption will take place
	diff_site = {} # rate  mapping where after surface diffusion atom will move
	diff_atom = {} # rate  mapping on which atom surface diffusion would occur
	event_list = {} # on this rate wheter adsorption , desorption or surface diffusion would occur 
	rate_list = {} # index of rate 0,1,2,3......
	sites_hash = {}
	atoms_hash = {}

	r_ads = 0 # initial rate adsorption
	r_des = 0 # initial rate desorption
	r_diff = 0 # initial rate diffusion
	rate_count = 0 # possible rates 0,1,2,3,..... for all the possible transitions
	rate_list[rate_count] = 0
	rate_count = rate_count + 1
	"""
	global t
	global array_pdb
	global array_sites
	global N_zero
	global N_plus
	global no_of_atoms
	
	while t < 3600:
		file1 = open("time.txt","a+")
		file2 = open("diameter.txt","a+") 
		adsorption()
		desorption()
		surface_diffusion()
		rate_event = find_event()
		if(event_list[rate_event] == 'adsorption'):
			N_plus = N_plus - 1
			no_of_atoms = no_of_atoms + 1
			array_pdb.append(ads[rate_event])
			array_sites.remove(ads[rate_event])
		elif(event_list[rate_event] == 'desorption'):
			N_plus = N_plus + 1
			no_of_atoms = no_of_atoms - 1
			array_pdb.remove(des[rate_event])
			array_sites.append(des[rate_event])
		else:

			array_pdb.append(diff_site[rate_event])
			array_pdb.remove(diff_atom[rate_event])
			array_sites.append(diff_atom[rate_event])
			array_sites.remove(diff_site[rate_event])

		r = random.uniform(0,1)
		step = (-1*(math.log(r)))/r_total
		#print r_total#, step 
		t = t + step
		file1.write(str(t))
		file1.write(str("\n"))
		diameter = Calculate_Diameter(no_of_atoms)
		file2.write(str(diameter))
		file2.write(str("\n"))
		print event_list[rate_event],t,diameter
		Concentration(step)

		file1.close()
		file2.close()	
	#adsorption()
	#desorption()
	#surface_diffusion()
	#rate_event = find_event()
	#print event_list[rate_event],diff_site[rate_event],diff_atom[rate_event]
	
	
	

main()	