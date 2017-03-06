import numpy as np
from math import *
import astropy.constants as apc
import astropy.units as u
import scipy.integrate as integrate
from scipy.optimize import fsolve
import multiprocessing
import warnings
warnings.filterwarnings("ignore")



def theorical_sagita(r,s,k,units_ref=u.um):
	if(type(r) == type(s)):
		r = r.to(units_ref).value
	c = 1./r
	s = s.to(units_ref).value
	z = ( (c*s**2) / (1+(1-(k+1)*c**2*s**2)**(0.5)) )
	return z * units_ref

def sagita(r, s, k, z, units_ref = u.um):
	if(type(r) == type(s)):
		r = r.to(units_ref).value
	sol = z.to(units_ref) - ( theorical_sagita(r,s,k,units_ref))
	return sol.to(units_ref).value


#Articulo Manuesl => diameters = [5.,5.5] rings= [45, 64, 72]
diameters = [4.]
for diameter in diameters:
	diameter = diameter*u.cm
	unit_ref = u.um
	rings= [6]
	r = 500.*u.m
	s = diameter * 0.5
	wavelength = 546.1 * u.nm 
	k=0

	print "Radius: %s\nWavelength: %s" % (s, wavelength)
	print "%s |%s |%s |%s |%s" % ('Rings Number'.ljust(20), 'Sagita'.ljust(20), 'Curvature Radius'.ljust(30), 'Theorical Sagita'.ljust(20), 'Dif')
	for rings_number in rings:
		z = (rings_number-1) * (0.5*wavelength)
		radius = fsolve(sagita, r.to(unit_ref).value, args=(s,k,z))
		#z = theorical_sagita(radius[0], s, k, unit_ref)
		print "%s |%s |%s |%s |%s" % (str(rings_number).ljust(20), str(z.to(u.um)).ljust(20), str((radius[0] * unit_ref).to(u.m)).ljust(30), str(theorical_sagita(radius[0]*unit_ref, s, k, unit_ref)).ljust(20),str(sagita(radius[0]*unit_ref, s, k, z, unit_ref)))


print theorical_sagita(98.8*u.cm, 134.05*u.mm*0.5,0).to(u.mm)

