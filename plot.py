import numpy as np
import matplotlib.pyplot as plt

with open('diameter.txt') as f:
    for line in f:
    	line = re.sub( '\s+', ' ',  line).strip()
    	y.append(flaot(line))

count = 0;
with open('time.txt') as f:
    for line in f:
    	line = re.sub( '\s+', ' ',  line).strip()
    	x.append(flaot(line))

plt.plot(x,y)
plt.show()
