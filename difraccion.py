import matplotlib.pyplot as plt
import numpy as np
from math import *
import astropy.units as u

lines = [150., 300., 400., 600., 1200.]
wave = 633 * u.nm
for line in lines:
	a =  1. / (line / (1*u.mm))
	teta = asin(wave.to(u.m).value/a.to(u.m).value) * u.radian
	print "Long.de onda= %s\t\ta = %s\t\tteta = %s (%s)" % (wave.to(u.nm),a.to(u.nm), teta.to(u.degree), teta.to(u.arcsec))

print "=========="

waves = [400 * u.nm, 700 * u.nm, 1 * u.nm]
line = 300
for wave in waves:
	a =  1. / (line / (1*u.mm))
	teta = asin(wave.to(u.m).value/a.to(u.m).value) * u.radian
	print "Long.de onda= %s\t\ta = %s\t\tteta = %s (%s)" % (wave.to(u.nm),a.to(u.nm), teta.to(u.degree), teta.to(u.arcsec))
