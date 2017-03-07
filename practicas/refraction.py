import matplotlib.pyplot as plt
import numpy as np
from math import *

def calc_teta_t(teta_i, n1, n2):
	teta_t = []
	for i in teta_i:
		teta_t.append(degrees(asin(n1/n2 * sin (radians(i)))))
	return teta_t

teta_i = range(0,90)
n_index = {'diamond':2.417, 'air':1.00, 'glass':1.490}
fig=plt.figure()

ax1 = fig.add_subplot(311)
n1 = n_index['air']
n2 = n_index['diamond']
teta_t = []
for i in teta_i:
	teta_t.append(degrees(asin(n1/n2 * sin (radians(i)))))

ax1.set_title('Air -> Diamond')
ax1.scatter(teta_i, teta_t)

ax2 = fig.add_subplot(312)
n1 = n_index['diamond']
n2 = n_index['air']
teta_t = []
for i in teta_i:
	try:
		teta_t.append(degrees(asin(n1/n2 * sin (radians(i)))))
	except:
		teta_t.append(0)
		print i

ax2.set_title('Diamond -> Air')
ax2.scatter(teta_i, teta_t)

ax3 = fig.add_subplot(313)
n1 = n_index['air']
n2 = n_index['glass']
teta_t = []
for i in teta_i:
	teta_t.append(degrees(asin(n1/n2 * sin (radians(i)))))

ax3.set_title('Air -> Glass')
ax3.scatter(teta_i, teta_t)
plt.show()

