#!/usr/bin/env python2.7
#vim:set ts=3:
#vim:set sw=3:
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import glob
import math
import re
from astropy.modeling import models, fitting
from datetime import datetime

mask = np.zeros((2048,2048))
mask[245,71] = 1
mask[219,52] = 1
mask[146,80] = 1
mask[119,104] = 1
mask[338,21] = 1
mask[311,17] = 1
mask[192,49] = 1
mask[146,80] = 1
mask[053,88] = 1

class ApogeeFile:
	def __init__(self, filename, medium_wavelength, profile, fwhm, fwhm_line, medium_wl_gauss, profile_gauss, fwhm_gauss, fwhm_gauss_line, error=0.01, show="y"):
		self.filename = filename
		self.medium_wl = medium_wavelength
		self.medium_wl_gauss = medium_wl_gauss
		self.profile = profile
		self.fwhm = fwhm
		self.fwhm_line = fwhm_line
		self.profile_gauss = profile_gauss
		self.fwhm_gauss = fwhm_gauss
		self.fwhm_gauss_line = fwhm_gauss_line
		self.error = error
		if (show == "y"):
			self.show = True
		else:
			self.show = False
	


def apply_mask(data, mask, row):
	right_pixels = mask[row,0:len(data)] == 0
	wrong_pixels = mask[row,0:len(data)] == 1
	data[wrong_pixels] = np.max(data[right_pixels])
	return data

def get_FWHM(data, delta_px, to_unit):
    difference = np.max(data) - np.min(data)
    HM = difference * limit
    pos_extremum = data.argmax()  
    nearest_above = (np.abs(data[pos_extremum:-1] - HM)).argmin() + pos_extremum
    nearest_below = (np.abs(data[0:pos_extremum] - HM)).argmin()
    FWHM = nearest_above - nearest_below
    FWHM_units = (FWHM*delta_px).to(to_unit)
    return (FWHM_units, data[nearest_above], nearest_above - FWHM * 0.5)

files = []
listDirs = ["./", "./fits/"]
for mydir in listDirs:
    filestemp = glob.glob("%s*.fits" % mydir)
    filestemp = sorted(filestemp)
    files = files + filestemp
i = 0
lineFilename = ""
for filename in files:
	i+=1
	strFilename = "[%2d] %s" % (i, filename)
	lineFilename = lineFilename + strFilename.ljust(50)
	if i % 2 == 0:
		print lineFilename
		lineFilename = ""

files_selected = raw_input("\nPlease, select files (index separated by colon, ex: 1,4,5): ")
files_index = files_selected.split(',')
initial_wavelength= raw_input("\nPlease, introduce the center wavelength (AA): ")
result = []
resultLines = []
for selection in files_index:
	for index in range(int(selection.split('-')[0]), int(selection.split('-')[-1])+1):
		try:
			filename = files[int(index)-1]
			print "\n\n================ [%s] ====================" % filename
			hdulist = pyfits.open(filename)
			perfil = []

			for i in range(len(hdulist[0].data)):
				data = hdulist[0].data[i]
				median = np.median(data)
				data = apply_mask(data, mask, i)
				data = data - median
				suma = np.sum(data)
				perfil.append(suma)

			x = np.arange(len(perfil))
			perfil = np.asarray(perfil)
			if (min(perfil) < 0.):
				perfil += abs(min(perfil))

			FWHM_perfil = np.zeros(len(perfil))

			limit = 0.5
			valid_answers = ["y","n",""]
			error = 0.01
			answer = "__"
			while answer not in valid_answers:
				answer = raw_input("Please, you must specify the noise level (0-1.) (default value = 0.01): ") 
				try:
					if (answer != ""):
						error = float(answer)
						answer = ""
				except:
					pass
			low_values = perfil < max(perfil) * error
			perfil[low_values] = np.median(perfil[low_values])
			FWHM, value, medium_value = get_FWHM(perfil, 0.87*u.AA, u.nm)

			perfil_gauss = None
			try:
				g_init = models.Gaussian1D(amplitude=max(perfil), mean=np.mean(perfil), stddev=np.std(perfil))
				fit_g = fitting.LevMarLSQFitter()
				g = fit_g(g_init, x, perfil)
				perfil_gauss = g(x)

				FWHM_gauss, value_gauss, medium_value_gauss = get_FWHM(perfil_gauss, 0.87*u.AA, u.nm)
			except Exception as e:
				FWHM_gauss, value_gauss, medium_value_gauss = 0,0,0
				perfil_gauss = None
				print "\tNo es posible ajustar (%s)" % str(e)

			print "\tDimension: (%d,%d)" % (len(hdulist[0].data[0]), len (hdulist[0].data))
			print "\tCentral Wavelength: %s" % (float(initial_wavelength)*u.AA)
			initial_window = (float(initial_wavelength)*u.AA - (len(perfil)-medium_value) * 0.87 * u.AA)
			final_window = initial_window + len(perfil) * 0.87 * u.AA
			print "\tWindow Values : %s - %s" % (initial_window, final_window)
			print "\tFWHM: %s" % FWHM
			initial_window_gauss= (float(initial_wavelength)*u.AA - (len(perfil)-medium_value_gauss) * 0.87 * u.AA)
			print "\tWindow Values (Gaussian): %s - %s" % (initial_window_gauss, initial_window_gauss + (len(perfil) * 0.87 * u.AA))
			print "\tFWHM (gaussian): %s" % (FWHM_gauss)
			print "\tExposure time: %s" % hdulist[0].header["EXPTIME"]
			strLine = "%s;(%d-%d);%E;%E;%E;%E;%d" % (filename.split('/')[-1],len(hdulist[0].data[0]), len (hdulist[0].data), float(initial_wavelength)*u.AA.cgs.value, initial_window.cgs.value, final_window.cgs.value, FWHM.cgs.value, hdulist[0].header["EXPTIME"])
			resultLines.append(strLine)

			valid_answers = ["y","n",""]
			answer = "__"
			while answer not in valid_answers:
				answer = raw_input("Do you want show graphic? y/N: ") 
			FWHM_perfil = np.zeros(len(perfil))
			FWHM_perfil += value
			FWHM_perfil_gauss = np.zeros(len(perfil))
			if (perfil_gauss is not None):
				FWHM_perfil_gauss = np.zeros(len(perfil_gauss))
				FWHM_perfil_gauss += value_gauss
				result.append(ApogeeFile(filename, medium_value, perfil, FWHM, FWHM_perfil, medium_value_gauss, perfil_gauss, FWHM_gauss, FWHM_perfil_gauss, error, answer))
		except Exception as e:
			print e


num_graphs = len(result)
if (num_graphs > 0):
	for i in range(num_graphs):
		fig = plt.figure()
		graph = result[i]
		initial_window= (float(initial_wavelength)*u.AA - (len(perfil)-graph.medium_wl) * 0.87 * u.AA)
		final_window = initial_window + len(perfil) * 0.87 * u.AA
		x = np.linspace(initial_window.value, final_window.value, len(perfil))

		ax = fig.add_subplot(1, 1,1)
		if (graph.profile is not None):
			graph.profile = list(graph.profile)
			graph.profile.reverse()
			ax.plot(x,graph.profile, 'r', label='Data')
			ax.plot(x,graph.fwhm_line, 'b', label='FWHM data')
		if (graph.profile_gauss is not None):
			graph.profile_gauss = list(graph.profile_gauss)
			graph.profile_gauss.reverse()
			ax.plot(x,graph.profile_gauss, 'g', label='Gaussian')
			ax.plot(x,graph.fwhm_gauss_line, 'black', label='FWHM Gauss')
		ax.set_title(r'%s, FWHM = %s, FWHM_gauss = %s' % (graph.filename, graph.fwhm, graph.fwhm_gauss))
		handles, labels = ax.get_legend_handles_labels()
		ax.legend(handles, labels, fontsize=10)
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()): 
			item.set_fontsize(10)
		if (graph.show):
			plt.show()
		fig.savefig("./images/%s.png" % graph.filename.split('/')[-1])
		plt.close(fig)

filename = "%s_calibration.csv" % datetime.now().strftime('%Y%m%d%H%M%S')
f = open(filename, "w")
f.write("filename;dimension;central_wavelength;initial_wavelength;final_wavelength;fwhm;exptime(s)\n")
for line in resultLines:
	f.write("%s\n" % line)
f.close()
