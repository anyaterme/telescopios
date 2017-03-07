import astropy.units as u
import numpy as np
from math import *

f_efec = {"84cm":10323.06 * u.mm, "1.5m":19575 * u.mm, "2.1 (conf a)": 15824*u.mm , "2.1 (conf b)": 28816.8*u.mm, "2.1 (conf c)": 63460 * u.mm}
diam = {"84cm":84 * u.cm, "1.5m":1.5 * u.m, "2.1 (conf a)": 2.1*u.m , "2.1 (conf b)": 2.1*u.m, "2.1 (conf c)": 2.1* u.m}
lambdas = [500,300,800]

for key in f_efec:
	print "======= %s ========" % key

	h = diam[key]
	f_efe = f_efec[key]
	Ep = (1 * u.rad).to(u.arcsec) / f_efe
	print "F_efe = %s, Escala de placa = %s" % (f_efe, Ep)
	print "Num F = %f" % (f_efe / h)

	pixels = 512
	SSD = pixels * 13.5 * u.micron
	print "Campo=%s x %s; Diag=%s; Por pix=%s" % ((SSD.to(u.mm) * Ep).to(u.arcmin),(SSD.to(u.mm) * Ep).to(u.arcmin),(2*((SSD.to(u.mm) * Ep).to(u.arcmin) )**2)**(0.5), (SSD.to(u.mm) * Ep).to(u.arcsec)/pixels)


	for long_onda in lambdas:
		print "Resolucion [%s] = %s" % (long_onda * u.nm, ((1.22 * long_onda * u.nm / h).cgs.value * u.rad).to(u.arcsec))

	t_exp = 81
	print "Magnitud maxima con %lf minutos = %lf " % (t_exp, 4 + 5 * np.log10 (h.cgs.value) + 2.5 * log10 (t_exp))

