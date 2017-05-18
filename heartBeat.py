import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.signal import butter, lfilter, detrend
from datetime import datetime

measures = {}

#Preprocessing
def get_data(filename):
	print "getting file"
	dataset = pd.read_csv(filename)
	print "file gotten"
	return dataset

def get_samplerate_mstimer(dataset):
	sampletimer = [x for x in dataset.datetime]
	fs = ((len(sampletimer) / sampletimer[-1])*1000)
	measures['fs'] = fs
	return fs

def get_samplerate_datetime(dataset):
	times = [x for x in dataset.datetime]
	elapsed = ((datetime.strptime(times[-1].split(" ")[-1], "%H:%M:%S.%f") - datetime.strptime(times[0].split(" ")[-1], "%H:%M:%S.%f")).total_seconds())
	fs = (len(times) / elapsed)
	measures['fs'] = fs
	return fs

def rolmean(dataset, hrw, fs):
	mov_avg = dataset.rolling(window=int(hrw*fs), center=False)['hr'].mean()
	avg_hr = (np.mean(dataset.hr)) 
	mov_avg = [avg_hr if math.isnan(x) else x for x in mov_avg]
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

def filtersignal(data, cutoff, fs, order):
	hr = [math.pow(x, 3) for x in dataset.hr]
	hrfiltered = butter_lowpass_filter(hr, cutoff, fs, order)
	dataset['hr'] = hrfiltered

#Peak detection
def detect_peaks(dataset, ma_perc, fs):
	rolmean = [(x+((x/100)*ma_perc)) for x in dataset.hr_rollingmean]
	window = []
	peaklist = []
	listpos = 0 
	for datapoint in dataset.hr:
		rollingmean = rolmean[listpos]
		if (datapoint <= rollingmean) and (len(window) <= 1):
			listpos += 1
		elif (datapoint > rollingmean):
			window.append(datapoint)
			listpos += 1
		else:
			maximum = max(window)
			beatposition = listpos - len(window) + (window.index(max(window)))
			peaklist.append(beatposition)
			window = []
			listpos += 1
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
		print x, y, z
		if ((x > 1) and ((y > 30) and (y < 130))):
			valid_ma.append([x, z])
	
	measures['best'] = min(valid_ma, key = lambda t: t[0])[1]
	detect_peaks(dataset, min(valid_ma, key = lambda t: t[0])[1], fs)

def check_peaks(dataset):
	RR_list = measures['RR_list']
	peaklist = measures['peaklist']
	ybeat = measures['ybeat']
	
	upper_threshold = np.mean(RR_list) + 300
	lower_threshold = np.mean(RR_list) - 300
	
	removed_beats = []
	removed_beats_y = []
	RR_list_cor = []
	peaklist_cor = [peaklist[0]]
	
	cnt = 0
	while cnt < len(RR_list):
		if (RR_list[cnt] < upper_threshold) and (RR_list[cnt] > lower_threshold):
			RR_list_cor.append(RR_list[cnt])
			peaklist_cor.append(peaklist[cnt+1]) 
			cnt += 1
		else:	
			removed_beats.append(peaklist[cnt+1])
			removed_beats_y.append(ybeat[cnt+1])
			cnt += 1

	measures['RR_list_cor'] = RR_list_cor
	measures['peaklist_cor'] = peaklist_cor
	measures['removed_beats'] = removed_beats
	measures['removed_beats_y'] = removed_beats_y

#Calculating all measures
def calc_RR(dataset, fs):
	peaklist = measures['peaklist']
	RR_list = []
	cnt = 0
	while (cnt < (len(peaklist)-1)):
		RR_interval = (peaklist[cnt+1] - peaklist[cnt])
		ms_dist = ((RR_interval / fs) * 1000.0)
		RR_list.append(ms_dist)
		cnt += 1

	RR_diff = []
	RR_sqdiff = []
	cnt = 0 
	while (cnt < (len(RR_list)-1)): 
		RR_diff.append(abs(RR_list[cnt] - RR_list[cnt+1])) 
		RR_sqdiff.append(math.pow(RR_list[cnt] - RR_list[cnt+1], 2))
		cnt += 1
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
	filtersignal(dataset, 2.5, fs, 5)
	rolmean(dataset, hrw, fs)
	fit_peaks(dataset, fs)
	calc_RR(dataset, fs)
	check_peaks(dataset)
	calc_ts_measures(dataset)
	calc_fd_measures(dataset, fs)


if __name__ == "__main__":
	dataset = get_data("data2.csv")
	#dataset = dataset[4500:].reset_index(drop=True)
	process(dataset, 0.75, get_samplerate_mstimer(dataset))
	plotter(dataset)