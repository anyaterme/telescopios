#!python
#Author: Daniel Jacobo Diaz Gonzalez
#Description
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import glob
from astropy.modeling import models, fitting

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
	def __init__(self, filename, profile, fwhm, fwhm_line, profile_gauss, fwhm_gauss, fwhm_gauss_line, error=0.01, show='n'):
		self.filename = filename
		self.profile = profile
		self.fwhm = fwhm
		self.fwhm_line = fwhm_line
		self.profile_gauss = profile_gauss
		self.fwhm_gauss = fwhm_gauss
		self.fwhm_gauss_line = fwhm_gauss_line
		self.error = error
		self.show = False
		if (show == "y"):
			self.show = True
	


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
	return (FWHM_units, data[nearest_above])

listDirs = ["./", "./fits/"]
totalfiles = []
for mydir in listDirs:
	files = glob.glob("%s*.fits" % mydir)
	files=sorted(files)
	totalfiles += files
i = 0
for filename in files:
	i+=1
	print "[%2d] %s" % (i, filename)

files_selected = raw_input("\nPlease, select files (index separated by colon, ex: 1,4,5): ")

files_index = files_selected.split(',')
result = []
for index in files_index:
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

		#perfil = (perfil * 1.) / max(perfil)

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
		FWHM, value = get_FWHM(perfil, 0.87*u.AA, u.nm)

		perfil_gauss = None
		try:
			g_init = models.Gaussian1D(amplitude=max(perfil), mean=np.mean(perfil), stddev=np.std(perfil))
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, x, perfil)
			perfil_gauss = g(x)
			FWHM_gauss, value_gauss = get_FWHM(perfil_gauss, 0.87*u.AA, u.nm)

		except:
			FWHM_gauss, value_gauss = 0,0
			perfil_gauss = None
			print "\tNo es posible ajustar"

		print "\tDimension: (%d,%d)" % (len(hdulist[0].data[0]), len (hdulist[0].data))
		print "\tFWHM: %s" % FWHM
		print "\tFWHM (gaussian): %s" % FWHM_gauss
		print "\tExposure time: %s" % hdulist[0].header["EXPTIME"]

		valid_answers = ["y","n",""]
		answer = "__"
		while answer not in valid_answers:
			answer = raw_input("Do you want show graphic? y/N: ") 
		FWHM_perfil = np.zeros(len(perfil))
		FWHM_perfil += value
		FWHM_perfil_gauss = np.zeros(len(perfil))
		FWHM_perfil_gauss += value_gauss
		result.append(ApogeeFile(filename, perfil, FWHM, FWHM_perfil, perfil_gauss, FWHM_gauss, FWHM_perfil_gauss, error, answer))
	except Exception as e:
		print e

num_graphs = len(result)
for i in range(num_graphs):
	fig = plt.figure()
	graph = result[i]
	x = np.arange(len(graph.profile))
	ax = fig.add_subplot(num_graphs, 1,1)
	if (graph.profile is not None):
		ax.plot(graph.profile, 'r', label='Data')
		ax.plot(graph.fwhm_line, 'b', label='FWHM data')
	if (graph.profile_gauss is not None):
		ax.plot(graph.profile_gauss, 'g', label='Gaussian')
		ax.plot(graph.fwhm_gauss_line, 'black', label='FWHM Gauss')
	ax.set_title(r'%s, FWHM = %s, FWHM_gauss = %s' % (graph.filename, graph.fwhm, graph.fwhm_gauss))
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels, fontsize=10)
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()): 
		item.set_fontsize(10)
	if graph.show :
		plt.show()
	filename = graph.filename.split('/')[-1]
	plt.savefig('./images/%s.png' % filename)
exit()
