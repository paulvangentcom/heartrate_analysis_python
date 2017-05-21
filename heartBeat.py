from datetime import datetime
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.signal import butter, lfilter

measures = {}

#Preprocessing
def get_data(filename):
	print "getting file"
	dataset = pd.read_csv(filename)
	print "file gotten"
	return dataset

def get_samplerate_mstimer(dataset):
	sampletimer = [x for x in dataset.timer]
	fs = ((len(sampletimer) / (sampletimer[-1]-sampletimer[0]))*1000)
	measures['fs'] = fs
	return fs

def get_samplerate_datetime(dataset):
	times = [x for x in dataset.datetime]
	elapsed = ((datetime.strptime(times[-1].split(" ")[-1], "%H:%M:%S.%f") - datetime.strptime(times[0].split(" ")[-1], "%H:%M:%S.%f")).total_seconds())
	fs = (len(times) / elapsed)
	measures['fs'] = fs
	return fs

def rolmean(dataset, hrw, fs):
	avg_hr = (np.mean(dataset.hr)) 
	mov_avg = dataset.rolling(window=int(hrw*fs), center=False)['hr'].mean().fillna(avg_hr)
	mov_avg = [x*1.1 for x in mov_avg]
	dataset['hr_rollingmean'] = mov_avg

def butter_lowpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	normal_cutoff = cutoff / nyq
	b, a = butter(order, normal_cutoff, btype='low', analog=False)
	return b, a

def butter_lowpass_filter(data, cutoff, fs, order):
	b, a = butter_lowpass(cutoff, fs, order=order)
	y = lfilter(b, a, data)
	return y    

def filtersignal(hrdata, cutoff, fs, order):
	hr = np.power(np.array(hrdata), 3)
	hrfiltered = butter_lowpass_filter(hr, cutoff, fs, order)
	return hrfiltered

#Peak detection
def detect_peaks(dataset, ma_perc, fs):
	rm = np.array(dataset.hr_rollingmean)
	rolmean = rm+((rm/100)*ma_perc)
	hr = np.array(dataset.hr)
	peaksx = np.where((hr > rolmean))[0]
	peaksy = hr[np.where((hr > rolmean))[0]]
	peakedges = np.concatenate((np.array([0]), (np.where(np.diff(peaksx) > 1)[0]), np.array([len(peaksx)])))
	peaklist = []
	ybeat = []
	
	for i in range(0, len(peakedges)-1):
		try:
			y = peaksy[peakedges[i]:peakedges[i+1]].tolist()
			peaklist.append(peaksx[peakedges[i] + y.index(max(y))])
		except:
			pass
	
	measures['peaklist'] = peaklist
	measures['ybeat'] = [dataset.hr[x] for x in peaklist]
	measures['rolmean'] = rolmean
	calc_RR(dataset, fs)
	measures['rrsd'] = np.std(measures['RR_list'])

def fit_peaks(dataset, fs):
	ma_perc_list = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300]
	rrsd = []
	valid_ma = []
	for x in ma_perc_list:
		detect_peaks(dataset, x, fs)
		bpm = ((len(measures['peaklist'])/(len(dataset.hr)/fs))*60)
		rrsd.append([measures['rrsd'], bpm, x])

	for x,y,z in rrsd:
		if ((x > 1) and ((y > 40) and (y < 150))):
			valid_ma.append([x, z])
	
	measures['best'] = min(valid_ma, key = lambda t: t[0])[1]
	detect_peaks(dataset, min(valid_ma, key = lambda t: t[0])[1], fs)

def check_peaks(dataset):
	RR_arr = np.array(measures['RR_list'])
	peaklist = np.array(measures['peaklist'])
	peaklist2 = measures['peaklist']
	ybeat = np.array(measures['ybeat'])
	upper_threshold = np.mean(RR_arr) + 300
	lower_threshold = np.mean(RR_arr) - 300
	measures['RR_list_cor'] = RR_arr[np.where((RR_arr > lower_threshold) & (RR_arr < upper_threshold))]
	peaklist_cor = peaklist[np.where((RR_arr > lower_threshold) & (RR_arr < upper_threshold))[0]+1]
	measures['peaklist_cor'] = np.insert(peaklist_cor, 0, peaklist[0])
	measures['removed_beats'] = peaklist[np.where((RR_arr <= lower_threshold) | (RR_arr >= upper_threshold))[0]+1]
	measures['removed_beats_y'] = removed_beats_y = ybeat[np.where((RR_arr <= lower_threshold) | (RR_arr >= upper_threshold))[0]+1]

#Calculating all measures
def calc_RR(dataset, fs):
	peaklist = np.array(measures['peaklist'])
	RR_list = (np.diff(peaklist) / fs) * 1000.0
	RR_diff = np.abs(np.diff(RR_list))
	RR_sqdiff = np.power(RR_diff, 2)
	measures['RR_list'] = RR_list
	measures['RR_diff'] = RR_diff
	measures['RR_sqdiff'] = RR_sqdiff
	
def calc_ts_measures(dataset):
	RR_list = measures['RR_list_cor']
	RR_diff = measures['RR_diff']
	RR_sqdiff = measures['RR_sqdiff']
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

def calc_fd_measures(dataset, fs):
	peaklist = measures['peaklist_cor']
	RR_list = measures['RR_list_cor']
	RR_x = peaklist[1:]
	RR_y = RR_list
	RR_x_new = np.linspace(RR_x[0],RR_x[-1],RR_x[-1])
	f = interp1d(RR_x, RR_y, kind='cubic')
	n = len(dataset.hr)
	frq = np.fft.fftfreq(len(dataset.hr), d=((1/fs)))
	frq = frq[range(n/2)]
	Y = np.fft.fft(f(RR_x_new))/n
	Y = Y[range(n/2)]
	measures['lf'] = np.trapz(abs(Y[(frq>=0.04) & (frq<=0.15)]))
	measures['hf'] = np.trapz(abs(Y[(frq>=0.16) & (frq<=0.5)]))
	measures['lf/hf'] = measures['lf'] / measures['hf']

#Plotting it
def plotter(dataset):
	peaklist = measures['peaklist']
	ybeat = measures['ybeat']
	rejectedpeaks = measures['removed_beats']
	rejectedpeaks_y = measures['removed_beats_y']
	plt.title("Heart Rate Signal Peak Detection")
	plt.plot(dataset.hr, alpha=0.5, color='blue', label="heart rate signal")
	plt.plot(measures['rolmean'], color ='grey', label="moving average")
	plt.scatter(peaklist, ybeat, color='green', label="BPM:%.2f" %(measures['bpm']))
	plt.scatter(rejectedpeaks, rejectedpeaks_y, color='red', label="rejected peaks")
	plt.legend(loc=4, framealpha=0.6) 
	plt.show() 

#Wrapper function
def process(dataset, hrw, fs):
	dataset['hr'] = filtersignal(dataset.hr, 2.5, fs, 5)
	rolmean(dataset, hrw, fs)
	fit_peaks(dataset, fs)
	calc_RR(dataset, fs)
	check_peaks(dataset)
	calc_ts_measures(dataset)
	calc_fd_measures(dataset, fs)
	return measures