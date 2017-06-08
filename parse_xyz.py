import re
import math
f=open("xyzpdb.txt", "w+")
f1 = open("step_descent.xyz","w+")

epsilon = 4.32
rmin   = 2.88358
step  = 0.001

with open("Au.pdb", "r") as file:
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
           s = "";
           s = temp_array[6] + ' ' + temp_array[7] +' ' +temp_array[8]
           f.write(s)
           f.write("\n")
           array_pdb.append(coord)
        
           #print coord
        count=count+1

with open("surfacepoints2.txt","r") as file:
    array_sites = []
    for line in file:
        coord = []
        line = re.sub( '\s+', ' ',  line).strip()
        temp_array = line.split(" ")
        temp_array[0] = float(temp_array[0])
        temp_array[1] = float(temp_array[1])
        temp_array[2] = float(temp_array[2])
        array_sites.append(temp_array)
    #print len(array_sites)



final_sites = []

for i , val in enumerate(array_sites):
    #print val
    x_site = val[0]
    y_site = val[1]
    z_site = val[2]
    for z in range(10000):
        fx_sum = 0.0
        fy_sum = 0.0
        fz_sum = 0.0
        for j , val1 in enumerate(array_pdb):
            x_pdb = val1[0]
            y_pdb = val1[1]
            z_pdb = val1[2]
            r = math.sqrt((x_site - x_pdb)*(x_site - x_pdb) + (y_site - y_pdb)*(y_site - y_pdb) + (z_site - z_pdb)*(z_site - z_pdb))
            force = 12.00*epsilon*(-1.00*(math.pow((rmin/r),12))*1.00/r + (math.pow((rmin/r),6))*1.00/r)
            fx_sum = fx_sum + force*((x_site - x_pdb)/r)
            fy_sum = fy_sum + force*((y_site - y_pdb)/r)
            fz_sum = fz_sum + force*((z_site - z_pdb)/r)

        force = math.sqrt(fx_sum*fx_sum + fy_sum*fy_sum + fz_sum*fz_sum)
        if(abs(force) < 0.0001):
            coord = []
            coord = [x_site,y_site,z_site]
            final_sites.append(coord)
            break
        x_site = x_site - step*fx_sum
        y_site = y_site - step*fy_sum
        z_site = z_site - step*fz_sum


count = 0
remove_sites = {} # dictionary which contains index to remove sites
for i in range(len(final_sites)):
    remove_sites[i] = 0
for i in range(len(final_sites)): 
    if remove_sites[i]:
        continue   
    for j in range(len(final_sites)):
        if(j > i and remove_sites[j] == 0):
            x = final_sites[i][0] - final_sites[j][0]
            y = final_sites[i][1] - final_sites[j][1]
            z = final_sites[i][2] - final_sites[j][2] 
            r = math.sqrt(x*x + y*y + z*z)
            if(r < 0.1):
                count = count + 1
                remove_sites[j] = 1
               

#print count

count = 0
for key,val in remove_sites.iteritems():
    if(val == 1):
        del final_sites[key-count]
        count = count + 1


print len(final_sites)

s = str(len(final_sites)-1)
f1.write(s)
f1.write("\n")



for i,val in enumerate(final_sites):
    s = ""
    s = "Au" + "\t" + str(val[0]) + "\t" + str(val[1]) + "\t" + str(val[2])
    f1.write(s)
    f1.write("\n")