'''
Experimental module for detecting first and second heart sounds on audio, 
and deriving measures from this.
'''

print('\n-------------------Warning-------------------\nThis is a highly \
experimental module. \n\nPlease be aware of this and do not use for any \
critical scientific work in its current form.\n \
---------------------------------------------')

from scipy.io import wavfile
from scipy.signal import resample
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../')
import heartbeat as hb

__author__ = "Paul van Gent"
__version__ = "Version 0.9"

def read_wav(file):
    '''Function reads file and outputs sample rate and data array'''
    samplerate, wv = wavfile.read(file)
    wv_data = np.array(wv, dtype=np.float32)
    return samplerate, wv_data

def mark_audio_peaks(peaklist, samplerate):
    '''
    Marks peaks in audio signal, implements crude type detection

    Keyword arguments:
    peaklist -- 1-dimensional list or array containing detected peak locations
    samplerate -- sample rate of the original signal
    '''

    hb.calc_rr(samplerate)
    first_heartsounds = []
    second_heartsounds = []
    rejected = []

    cnt = 0
    for rr in hb.working_data['RR_list']:
        if rr <= 0.6*np.mean(hb.working_data['RR_list']):
            rejected.append(peaklist[cnt])
        elif rr >= 0.9*np.mean(hb.working_data['RR_list']):
            first_heartsounds.append(peaklist[cnt])
        else:
            second_heartsounds.append(peaklist[cnt])
        cnt += 1

    hb.working_data['peaklist'] = first_heartsounds
    hb.working_data['second_heartsounds'] = second_heartsounds
    hb.working_data['removed_beats'] = rejected

def plotter(show=True, title='Heart Rate Audio Signal Peak Detection'):
    '''Alternative plotter function to plot signals'''

    plt.title(title)
    plt.plot(abs(hb.working_data['hr']))
    plt.plot(hb.working_data['rolmean'])
    plt.scatter(hb.working_data['peaklist'], [hb.working_data['hr'][x] for x in hb.working_data['peaklist']], color='green', label='first heartsounds')
    plt.scatter(hb.working_data['second_heartsounds'], [hb.working_data['hr'][x] for x in hb.working_data['second_heartsounds']], color='black', label='second heartsounds')
    plt.scatter(hb.working_data['removed_beats'], [hb.working_data['hr'][x] for x in hb.working_data['removed_beats']], color='red', label='rejected peaks')
    plt.legend()
    plt.show()

def process(filename):
    '''Processes the wav file passed to it

    Keyword arguments:
    filename -- absolute or relative path to WAV file
    '''

    print('resampling signal to 1000Hz')
    samplerate, wv_data = read_wav(filename)
    wv_data = abs(wv_data)
    wavlength = len(wv_data) / samplerate

    wv_data = resample(wv_data, int(1000*wavlength))
    new_samplerate = len(wv_data) / wavlength

    print('detecting peaks')
    hb.working_data['hr'] = wv_data
    hb.working_data['hr'] = hb.butter_lowpass_filter(hb.working_data['hr'], 5, new_samplerate, 2)
    rol_mean = hb.rolmean(hb.working_data['hr'], 0.75, new_samplerate)
    hb.working_data['hr'] = np.resize(hb.working_data['hr'], rol_mean.shape)
    hb.fit_peaks(hb.working_data['hr'], rol_mean, new_samplerate)
    mark_audio_peaks(hb.working_data['peaklist'], new_samplerate)
    hb.working_data['peaklist_cor'] = hb.working_data['peaklist']
    hb.calc_rr(new_samplerate)
    hb.working_data['RR_list_cor'] = hb.working_data['RR_list']

    print('calculating measures')
    hb.calc_ts_measures()
    hb.calc_fd_measures(hb.working_data['hr'], samplerate)
    return hb.measures
    #plotter()

if __name__ == '__main__':
    filename = 'heartbeat.wav'
    measures = process(filename)
    print(measures)

