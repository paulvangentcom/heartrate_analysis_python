from datetime import datetime
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import butter, lfilter

__author__ = "Paul van Gent"
__version__ = "Version 0.9"

measures = {}
working_data = {}

#Data handling
def get_data(filename, delim = ',', column_name = 'None'):
	file_ext = filename.split('.')[-1]
	if file_ext == 'csv' or file_ext == 'txt':
		if (column_name != 'None'):
			hrdata = np.genfromtxt(filename, delimiter=delim, names=True, dtype=None)
			try:
				hrdata = hrdata[column_name]
			except:
				print('\nError loading column "%s" from file "%s". Is column name specified correctly?' %(column_name, filename))
		elif (column_name == 'None'):
			hrdata = np.genfromtxt(filename, delimiter=delim, dtype = np.float64)
		else:
			print('\nError: column name "%s" not found in header of "%s".\n' %(column_name, filename))
	elif file_ext == 'mat':
		print('getting matlab file')
		import scipy.io
		data = scipy.io.loadmat(filename)
		if (column_name != "None"):
			hrdata = np.array(data[column_name][:,0], dtype=np.float64)
		else:
			print("\nError: column name required for Matlab .mat files\n\n")
	else:
		print('unknown file format')
		hrdata = np.nan
	return hrdata

#Preprocessing
def get_samplerate_mstimer(timerdata):
	fs = ((len(timerdata) / (timerdata[-1]-timerdata[0]))*1000)
	working_data['fs'] = fs
	return fs

def get_samplerate_datetime(datetimedata, timeformat='%H:%M:%S.%f'):
	elapsed = ((datetime.strptime(datetimedata[-1], timeformat) - datetime.strptime(datetimedata[0], timeformat)).total_seconds())
	fs = (len(datetimedata) / elapsed)
	working_data['fs'] = fs
	return fs

def rollwindow(x, window):
	shape = x.shape[:-1] + (x.shape[-1] - window + 1, window)
	strides = x.strides + (x.strides[-1],)
	return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides)

def rolmean(hrdata, hrw, fs):
	avg_hr = (np.mean(hrdata)) 
	hrarr = np.array(hrdata)
	rol_mean = np.mean(rollwindow(hrarr, int(hrw*fs)), axis=1)
	ln = np.array([avg_hr for i in range(0,abs(len(hrarr)-len(rol_mean))/2)])
	rol_mean = np.insert(rol_mean, 0, ln)
	rol_mean = np.append(rol_mean, ln)
	rol_mean = rol_mean * 1.1
	return rol_mean

def butter_lowpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	normal_cutoff = cutoff / nyq
	b, a = butter(order, normal_cutoff, btype='low', analog=False)
	return b, a

def butter_lowpass_filter(hrdata, cutoff, fs, order):
	b, a = butter_lowpass(cutoff, fs, order=order)
	y = lfilter(b, a, hrdata)
	return y    

def filtersignal(hrdata, cutoff, fs, order):
	hr = np.power(np.array(hrdata), 3)
	hrfiltered = butter_lowpass_filter(hr, cutoff, fs, order)
	return hrfiltered

#Peak detection
def detect_peaks(hrdata, rol_mean, ma_perc, fs):
	rm = np.array(rol_mean)
	rolmean = rm+((rm/100)*ma_perc)
	peaksx = np.where((hrdata > rolmean))[0]
	peaksy = hrdata[np.where((hrdata > rolmean))[0]]
	peakedges = np.concatenate((np.array([0]), (np.where(np.diff(peaksx) > 1)[0]), np.array([len(peaksx)])))
	peaklist = []
	ybeat = []

	for i in range(0, len(peakedges)-1):
		try:
			y = peaksy[peakedges[i]:peakedges[i+1]].tolist()
			peaklist.append(peaksx[peakedges[i] + y.index(max(y))])
		except:
			pass
	
	working_data['peaklist'] = peaklist
	working_data['ybeat'] = [hrdata[x] for x in peaklist]
	working_data['rolmean'] = rolmean
	calc_RR(fs)
	if len(working_data['RR_list']):
		working_data['rrsd'] = np.std(working_data['RR_list'])
	else:
		working_data['rrsd'] = np.inf

def fit_peaks(hrdata, rol_mean, fs):
	ma_perc_list = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300]
	rrsd = []
	valid_ma = []
	for x in ma_perc_list:
		detect_peaks(hrdata, rol_mean, x, fs)
		bpm = ((len(working_data['peaklist'])/(len(working_data['hr'])/fs))*60)
		rrsd.append([working_data['rrsd'], bpm, x])

	for x,y,z in rrsd:
		if ((x > 1) and ((y > 40) and (y < 150))):
			valid_ma.append([x, z])
	
	working_data['best'] = min(valid_ma, key = lambda t: t[0])[1]
	detect_peaks(hrdata, rol_mean, min(valid_ma, key = lambda t: t[0])[1], fs)

def check_peaks():
	RR_arr = np.array(working_data['RR_list'])
	peaklist = np.array(working_data['peaklist'])
	peaklist2 = working_data['peaklist']
	ybeat = np.array(working_data['ybeat'])
	upper_threshold = np.mean(RR_arr) + 300
	lower_threshold = np.mean(RR_arr) - 300
	working_data['RR_list_cor'] = RR_arr[np.where((RR_arr > lower_threshold) & (RR_arr < upper_threshold))]
	peaklist_cor = peaklist[np.where((RR_arr > lower_threshold) & (RR_arr < upper_threshold))[0]+1]
	working_data['peaklist_cor'] = np.insert(peaklist_cor, 0, peaklist[0])
	working_data['removed_beats'] = peaklist[np.where((RR_arr <= lower_threshold) | (RR_arr >= upper_threshold))[0]+1]
	working_data['removed_beats_y'] = removed_beats_y = ybeat[np.where((RR_arr <= lower_threshold) | (RR_arr >= upper_threshold))[0]+1]

#Calculating all measures
def calc_RR(fs):
	peaklist = np.array(working_data['peaklist'])
	RR_list = (np.diff(peaklist) / fs) * 1000.0
	RR_diff = np.abs(np.diff(RR_list))
	RR_sqdiff = np.power(RR_diff, 2)
	working_data['RR_list'] = RR_list
	working_data['RR_diff'] = RR_diff
	working_data['RR_sqdiff'] = RR_sqdiff
	
def calc_ts_measures():
	RR_list = working_data['RR_list_cor']
	RR_diff = working_data['RR_diff']
	RR_sqdiff = working_data['RR_sqdiff']
	measures['bpm'] = 60000 / np.mean(RR_list)
	measures['ibi'] = np.mean(RR_list)
	measures['sdnn'] = np.std(RR_list)
	measures['sdsd'] = np.std(RR_diff)
	measures['rmssd'] = np.sqrt(np.mean(RR_sqdiff))
	NN20 = [x for x in RR_diff if (x>20)]
	NN50 = [x for x in RR_diff if (x>50)]
	measures['nn20'] = NN20
	measures['nn50'] = NN50
	measures['pnn20'] = float(len(NN20)) / float(len(RR_diff))
	measures['pnn50'] = float(len(NN50)) / float(len(RR_diff))
	measures['hr_mad'] = np.median(np.abs(RR_list-np.median(RR_list)))

def calc_fd_measures(hrdata, fs):
	peaklist = working_data['peaklist_cor']
	RR_list = working_data['RR_list_cor']
	RR_x = peaklist[1:]
	RR_y = RR_list
	RR_x_new = np.linspace(RR_x[0],RR_x[-1],RR_x[-1])
	f = interp1d(RR_x, RR_y, kind='cubic')
	n = len(hrdata)
	frq = np.fft.fftfreq(len(hrdata), d=((1/fs)))
	frq = frq[range(n/2)]
	Y = np.fft.fft(f(RR_x_new))/n
	Y = Y[range(n/2)]
	measures['lf'] = np.trapz(abs(Y[(frq>=0.04) & (frq<=0.15)]))
	measures['hf'] = np.trapz(abs(Y[(frq>=0.16) & (frq<=0.5)]))
	measures['lf/hf'] = measures['lf'] / measures['hf']

#Plotting it
def plotter(show=True, title='Heart Rate Signal Peak Detection'):
	peaklist = working_data['peaklist']
	ybeat = working_data['ybeat']
	rejectedpeaks = working_data['removed_beats']
	rejectedpeaks_y = working_data['removed_beats_y']
	plt.title(title)
	plt.plot(working_data['hr'], alpha=0.5, color='blue', label='heart rate signal')
	plt.scatter(peaklist, ybeat, color='green', label='BPM:%.2f' %(measures['bpm']))
	plt.scatter(rejectedpeaks, rejectedpeaks_y, color='red', label='rejected peaks')
	plt.legend(loc=4, framealpha=0.6) 
	if show == True:
		plt.show() 
	else:
		return plt

#Wrapper function
def process(hrdata, fs, hrw = 0.75):
	t1 = time.clock()
	hrdata = filtersignal(hrdata, 4, fs, 5)
	working_data['hr'] = hrdata
	rol_mean = rolmean(hrdata, hrw, fs)
	fit_peaks(hrdata, rol_mean, fs)
	calc_RR(fs)
	check_peaks()
	calc_ts_measures()
	calc_fd_measures(hrdata, fs)
	print('\nFinished in %.4fsec' %(time.clock()-t1))
	return measures
