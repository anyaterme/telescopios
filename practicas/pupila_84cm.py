import astropy.units as u
import numpy as np
from math import *

f_efec = {"84cm":10323.06 * u.mm, "1.5m":19575 * u.mm, "2.1 (conf a)": 15824*u.mm , "2.1 (conf b)": 28816.8*u.mm, "2.1 (conf c)": 63460 * u.mm}
diam = {"84cm":84 * u.cm, "1.5m":1.5 * u.m, "2.1 (conf a)": 2.1*u.m , "2.1 (conf b)": 2.1*u.m, "2.1 (conf c)": 2.1* u.m}

S_0 = 2029.7 * u.mm 
h=840 * u.mm
Rs = -1555 * u.mm
Rp = 5287 * u.mm
f2 = Rs / 2.
f1 = Rp / 2.

deltaArray = [0, 1, 2.54]

for delta in deltaArray:
	S_0 = S_0 + delta * u.cm

	print "Desplazamiento de %s" % (delta*u.cm)

	S_i = f2 * S_0 / (S_0 -f2)
	m= -S_i/S_0
	h_2 = m*h

	print "S_i = %s, m = %s, pupila_out = %s" % (S_i, m, h_2)

	f_efe = (f1*f2)  / (f1 +f2 - S_0)

	Ep = (1 * u.rad).to(u.arcsec) / f_efe

	print "F_efe = %s, Escala de placa = %s" % (f_efe, Ep)

	print "Num F = %f" % (f_efe / h)

	pixels = 512

	SSD = 512 * 15 * u.micron
	print "Campo=%s x %s; Diag=%s; Por pix=%s" % ((SSD.to(u.mm) * Ep).to(u.arcmin),(SSD.to(u.mm) * Ep).to(u.arcmin),(2*((SSD.to(u.mm) * Ep).to(u.arcmin) )**2)**(0.5), (SSD.to(u.mm) * Ep).to(u.arcsec)/pixels)


	lambdas = [300,500,800]
	for long_onda in lambdas:
		print "Resolucion [%s] = %s" % (long_onda * u.nm, ((1.22 * long_onda * u.nm / h).cgs.value * u.rad).to(u.arcsec))

	t_exp = 81
	print "Magnitud maxima con %lf minutos = %lf " % (t_exp, 4 + 5 * np.log10 (h.cgs.value) + 2.5 * log10 (t_exp))

