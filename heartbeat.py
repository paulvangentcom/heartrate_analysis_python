'''Noise-resistant heart rate analysis module for Python

#Reference:
<Article submitted, awaiting publication>

See also:
http://www.paulvangent.com/2016/03/15/analyzing-a-discrete-heart-rate-signal-using-python-part-1/
http://www.paulvangent.com/2016/03/21/analyzing-a-discrete-heart-rate-signal-using-python-part-2/
http://www.paulvangent.com/2016/03/30/analyzing-a-discrete-heart-rate-signal-using-python-part-3/
<part 4 to follow after publication>
'''

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
def get_data(filename, delim=',', column_name='None'):
    '''Loads data from a .CSV or .MAT file into numpy array.

    Keyword Arguments:
    filename -- absolute or relative path to the file object to read
    delim -- the delimiter used if CSV file passed, (default ',')
    column_name -- for CSV files with header: specify column that contains the data
                   for matlab files it specifies the table name that contains the data
                   (default 'None')
    '''
    file_ext = filename.split('.')[-1]
    if file_ext == 'csv' or file_ext == 'txt':
        if column_name != 'None':
            hrdata = np.genfromtxt(filename, delimiter=delim, names=True, dtype=None)
            try:
                hrdata = hrdata[column_name]
            except Exception as error:
                print('\nError loading column "%s" from file "%s". \
                Is column name specified correctly?\n' %(column_name, filename))
                print('------\nError message: ' + str(error) + '\n------')
        elif column_name == 'None':
            hrdata = np.genfromtxt(filename, delimiter=delim, dtype=np.float64)
        else:
            print('\nError: column name "%s" not found in header of "%s".\n'
                  %(column_name, filename))
    elif file_ext == 'mat':
        print('getting matlab file')
        import scipy.io
        data = scipy.io.loadmat(filename)
        if column_name != "None":
            hrdata = np.array(data[column_name][:, 0], dtype=np.float64)
        else:
            print("\nError: column name required for Matlab .mat files\n\n")
    else:
        print('unknown file format')
        return None
    return hrdata

#Preprocessing
def get_samplerate_mstimer(timerdata):
    '''Determines sample rate of data from ms-based timer.

    Keyword arguments:
    timerdata -- array containing values of a timer, in ms
    '''
    sample_rate = ((len(timerdata) / (timerdata[-1]-timerdata[0]))*1000)
    working_data['sample_rate'] = sample_rate
    return sample_rate

def get_samplerate_datetime(datetimedata, timeformat='%H:%M:%S.%f'):
    '''Determines sample rate of data from datetime-based timer.

    Keyword arguments:
    timerdata -- array containing values of a timer, datetime strings
    timeformat -- the format of the datetime-strings in datetimedata
    default('%H:%M:%S.f', 24-hour based time including ms: 21:43:12.569)
    '''
    elapsed = ((datetime.strptime(datetimedata[-1], timeformat) -
                datetime.strptime(datetimedata[0], timeformat)).total_seconds())
    sample_rate = (len(datetimedata) / elapsed)
    working_data['sample_rate'] = sample_rate
    return sample_rate

def rollwindow(data, windowsize):
    '''Returns rolling window of size 'window' over dataset 'data'.

    Keyword arguments:
    data -- 1-dimensional numpy array
    window -- window size
    '''
    shape = data.shape[:-1] + (data.shape[-1] - windowsize + 1, windowsize)
    strides = data.strides + (data.strides[-1],)
    return np.lib.stride_tricks.as_strided(data, shape=shape, strides=strides)

def rolmean(data, windowsize, sample_rate):
    '''Calculates the rolling mean over passed data.

    Keyword arguments:
    data -- 1-dimensional numpy array or list
    windowsize -- the window size to use, in seconds (calculated as windowsize * sample_rate)
    sample_rate -- the sample rate of the data set
    '''
    avg_hr = (np.mean(data))
    data_arr = np.array(data)
    rol_mean = np.mean(rollwindow(data_arr, int(windowsize*sample_rate)), axis=1)
    missing_vals = np.array([avg_hr for i in range(0, int(abs(len(data_arr) - len(rol_mean))/2))])
    rol_mean = np.insert(rol_mean, 0, missing_vals)
    rol_mean = np.append(rol_mean, missing_vals)
    rol_mean = rol_mean * 1.1
    return rol_mean

def butter_lowpass(cutoff, sample_rate, order=2):
    '''Defines standard Butterworth lowpass filter.

    use 'butter_lowpass_filter' to call the filter.
    '''
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, sample_rate, order):
    '''Applies the Butterworth lowpass filter

    Keyword arguments:
    data -- 1-dimensional numpy array or list containing the to be filtered data
    cutoff -- the cutoff frequency of the filter
    sample_rate -- the sample rate of the data set
    order -- the filter order (default 2)
    '''
    b, a = butter_lowpass(cutoff, sample_rate, order=order)
    filtered_data = lfilter(b, a, data)
    return filtered_data

def filtersignal(data, cutoff, sample_rate, order):
    '''Filters the given signal using a Butterworth lowpass filter.

    Keyword arguments:
    data -- 1-dimensional numpy array or list containing the to be filtered data
    cutoff -- the cutoff frequency of the filter
    sample_rate -- the sample rate of the data set
    order -- the filter order (default 2)
    '''
    data = np.power(np.array(data), 3)
    filtered_data = butter_lowpass_filter(data, cutoff, sample_rate, order)
    return filtered_data

#Peak detection
def detect_peaks(hrdata, rol_mean, ma_perc, sample_rate):
    '''Detects heartrate peaks in the given dataset.

    Keyword arguments:
   hr data -- 1-dimensional numpy array or list containing the heart rate data
    rol_mean -- 1-dimensional numpy array containing the rolling mean of the heart rate signal
    ma_perc -- the percentage with which to raise the rolling mean,
    used for fitting detection solutions to data
    sample_rate -- the sample rate of the data set
    '''
    rmean = np.array(rol_mean)
    rol_mean = rmean+((rmean/100)*ma_perc)
    peaksx = np.where((hrdata > rol_mean))[0]
    peaksy = hrdata[np.where((hrdata > rol_mean))[0]]
    peakedges = np.concatenate((np.array([0]),
                                (np.where(np.diff(peaksx) > 1)[0]),
                                np.array([len(peaksx)])))
    peaklist = []

    for i in range(0, len(peakedges)-1):
        try:
            y_values = peaksy[peakedges[i]:peakedges[i+1]].tolist()
            peaklist.append(peaksx[peakedges[i] + y_values.index(max(y_values))])
        except:
            pass

    working_data['peaklist'] = peaklist
    working_data['ybeat'] = [hrdata[x] for x in peaklist]
    working_data['rolmean'] = rol_mean
    calc_rr(sample_rate)
    if len(working_data['RR_list']):
        working_data['rrsd'] = np.std(working_data['RR_list'])
    else:
        working_data['rrsd'] = np.inf

def fit_peaks(hrdata, rol_mean, sample_rate):
    '''Runs variations in peak detection given a noisy heart rate signal

    Keyword arguments:
    hrdata - 1-dimensional numpy array or list containing the heart rate data
    rol_mean -- 1-dimensional numpy array containing the rolling mean of the heart rate signal
    sample_rate -- the sample rate of the data set
    '''
    ma_perc_list = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300]
    rrsd = []
    valid_ma = []
    for ma_perc in ma_perc_list:
        detect_peaks(hrdata, rol_mean, ma_perc, sample_rate)
        bpm = ((len(working_data['peaklist'])/(len(working_data['hr'])/sample_rate))*60)
        rrsd.append([working_data['rrsd'], bpm, ma_perc])

    for _rrsd, _bpm, _ma_perc in rrsd:
        if (_rrsd > 1) and ((_bpm > 40) and (_ma_perc < 150)):
            valid_ma.append([_rrsd, _ma_perc])

    working_data['best'] = min(valid_ma, key=lambda t: t[0])[1]
    detect_peaks(hrdata, rol_mean, min(valid_ma, key=lambda t: t[0])[1], sample_rate)

def check_peaks():
    '''Determines the best fit for peak detection variations run by fit_peaks().'''
    rr_arr = np.array(working_data['RR_list'])
    peaklist = np.array(working_data['peaklist'])
    ybeat = np.array(working_data['ybeat'])
    upper_threshold = np.mean(rr_arr) + 300
    lower_threshold = np.mean(rr_arr) - 300
    working_data['RR_list_cor'] = rr_arr[np.where((rr_arr > lower_threshold) &
                                                  (rr_arr < upper_threshold))]
    peaklist_cor = peaklist[np.where((rr_arr > lower_threshold) &
                                     (rr_arr < upper_threshold))[0]+1]
    working_data['peaklist_cor'] = np.insert(peaklist_cor, 0, peaklist[0])
    working_data['removed_beats'] = peaklist[np.where((rr_arr <= lower_threshold) |
                                                      (rr_arr >= upper_threshold))[0]+1]
    working_data['removed_beats_y'] = ybeat[np.where((rr_arr <= lower_threshold) |
                                                     (rr_arr >= upper_threshold))[0]+1]

#Calculating all measures
def calc_rr(sample_rate):
    '''Calculates the R-R (peak-peak) data required for further analysis.

    Uses calculated measures stored in the working_data{} dict to calculate
    all required peak-peak datasets. Stores results in the working_data{} dict.

    Keyword arguments:
    sample_rate -- the sample rate of the data set
    '''
    peaklist = np.array(working_data['peaklist'])
    rr_list = (np.diff(peaklist) / sample_rate) * 1000.0
    rr_diff = np.abs(np.diff(rr_list))
    rr_sqdiff = np.power(rr_diff, 2)
    working_data['RR_list'] = rr_list
    working_data['RR_diff'] = rr_diff
    working_data['RR_sqdiff'] = rr_sqdiff

def calc_ts_measures():
    '''Calculates the time-series measurements.

    Uses calculated measures stored in the working_data{} dict to calculate
    the time-series measurements of the heart rate signal.
    Stores results in the measures{} dict object.
    '''
    rr_list = working_data['RR_list_cor']
    rr_diff = working_data['RR_diff']
    rr_sqdiff = working_data['RR_sqdiff']
    measures['bpm'] = 60000 / np.mean(rr_list)
    measures['ibi'] = np.mean(rr_list)
    measures['sdnn'] = np.std(rr_list)
    measures['sdsd'] = np.std(rr_diff)
    measures['rmssd'] = np.sqrt(np.mean(rr_sqdiff))
    nn20 = [x for x in rr_diff if x > 20]
    nn50 = [x for x in rr_diff if x > 50]
    measures['nn20'] = nn20
    measures['nn50'] = nn50
    measures['pnn20'] = float(len(nn20)) / float(len(rr_diff))
    measures['pnn50'] = float(len(nn50)) / float(len(rr_diff))
    measures['hr_mad'] = np.median(np.abs(rr_list-np.median(rr_list)))

def calc_fd_measures(hrdata, sample_rate):
    '''Calculates the frequency-domain measurements.

    Uses calculated measures stored in the working_data{} dict to calculate
    the frequency-domain measurements of the heart rate signal.
    Stores results in the measures{} dict object.
    '''
    peaklist = working_data['peaklist_cor']
    rr_list = working_data['RR_list_cor']
    rr_x = peaklist[1:]
    rr_y = rr_list
    rr_x_new = np.linspace(rr_x[0], rr_x[-1], rr_x[-1])
    interpolated_func = interp1d(rr_x, rr_y, kind='cubic')
    datalen = len(hrdata)
    frq = np.fft.fftfreq(len(hrdata), d=((1/sample_rate)))
    frq = frq[range(int(datalen/2))]
    Y = np.fft.fft(interpolated_func(rr_x_new))/datalen
    Y = Y[range(int(datalen/2))]
    measures['lf'] = np.trapz(abs(Y[(frq >= 0.04) & (frq <= 0.15)]))
    measures['hf'] = np.trapz(abs(Y[(frq >= 0.16) & (frq <= 0.5)]))
    measures['lf/hf'] = measures['lf'] / measures['hf']

#Plotting it
def plotter(show=True, title='Heart Rate Signal Peak Detection'):
    '''Plots the analysis results.

    Uses calculated measures and data stored in the working_data{} and measures{}
    dict objects to visualise the fitted peak detection solution.

    Keyword arguments:
    show -- whether to display the plot (True) or return a plot object (False) (default True)
    title -- the title used in the plot
    '''
    peaklist = working_data['peaklist']
    ybeat = working_data['ybeat']
    rejectedpeaks = working_data['removed_beats']
    rejectedpeaks_y = working_data['removed_beats_y']
    plt.title(title)
    plt.plot(working_data['hr'], alpha=0.5, color='blue', label='heart rate signal')
    plt.scatter(peaklist, ybeat, color='green', label='BPM:%.2f' %(measures['bpm']))
    plt.scatter(rejectedpeaks, rejectedpeaks_y, color='red', label='rejected peaks')
    plt.legend(loc=4, framealpha=0.6)
    if show:
        plt.show()
    else:
        return plt

#Wrapper function
def process(hrdata, sample_rate, windowsize=0.75):
    '''Processed the passed heart rate data. Returns measures{} dict containing results.

    Keyword arguments:
    hrdata -- 1-dimensional numpy array or list containing heart rate data
    sample_rate -- the sample rate of the heart rate data
    windowsize -- the window size to use, in seconds (calculated as windowsize * sample_rate)
    '''
    t1 = time.clock()
    working_data['hr'] = hrdata
    rol_mean = rolmean(hrdata, windowsize, sample_rate)
    fit_peaks(hrdata, rol_mean, sample_rate)
    calc_rr(sample_rate)
    check_peaks()
    calc_ts_measures()
    calc_fd_measures(hrdata, sample_rate)
    print('\nFinished in %.8s sec' %(time.clock()-t1))
    return measures
