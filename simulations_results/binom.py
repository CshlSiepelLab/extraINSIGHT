import sys
import numpy as np
import scipy
import random
f=sys.argv[1]
with open(f,'r') as fopen:
	for t in fopen:
		x=random.uniform(0,1)
		y=float(t.split()[6])
		if x<y:
			print('1',t.split()[4],"{:.8f}".format(y))
		else:
			print('0',t.split()[4],"{:.8f}".format(y))
