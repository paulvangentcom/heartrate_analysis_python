'''
Functions that handle computation of heart rate (HR) and
heart rate variability (HRV) measures.
'''

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.signal import welch, periodogram

from .datautils import MAD, rolling_mean, outliers_iqr_method, outliers_modified_z
import heartpy as hp


__all__ = ['calc_rr',
           'update_rr',
           'calc_rr_segment',
           'clean_rr_intervals',
           'calc_ts_measures',
           'calc_fd_measures',
           'calc_breathing']


def calc_rr(peaklist, sample_rate, working_data={}):
    '''calculates peak-peak intervals

    Function that calculates the peak-peak data required for 
    further analysis. Stores results in the working_data{} dict.

    Parameters
    ----------
    peaklist : 1d list or array
        list or array containing detected peak positions

    sample_rate : int or float
        the sample rate with which the heart rate signal is collected

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    working_data : dict
        working_data dictionary object containing all of heartpy's temp objects

    Examples
    --------
    Let's assume we detected peaks at these positions in the signal:

    >>> peaklist = [200, 280, 405, 501, 615]
    
    It is then easy to call calc_rr to compute what we need:

    >>> wd = calc_rr(peaklist, sample_rate = 100.0)
    >>> wd['RR_list']
    array([ 800., 1250.,  960., 1140.])
    >>> wd['RR_diff']
    array([450., 290., 180.])
    >>> wd['RR_sqdiff']
    array([202500.,  84100.,  32400.])

    Note that the list of peak-peak intervals is of length len(peaks) - 1
    the length of the differences is of length len(peaks) - 2
    '''
    peaklist = np.array(peaklist) #cast numpy array to be sure or correct array type

    #delete first peak if within first 150ms (signal might start mid-beat after peak)
    if len(peaklist) > 0:
        if peaklist[0] <= ((sample_rate / 1000.0) * 150):
            peaklist = np.delete(peaklist, 0)
            working_data['peaklist'] = peaklist
            working_data['ybeat'] = np.delete(working_data['ybeat'], 0)

    rr_list = (np.diff(peaklist) / sample_rate) * 1000.0
    rr_diff = np.abs(np.diff(rr_list))
    rr_sqdiff = np.power(rr_diff, 2)
    working_data['RR_list'] = rr_list
    working_data['RR_diff'] = rr_diff
    working_data['RR_sqdiff'] = rr_sqdiff
    return working_data


def update_rr(working_data={}):
    '''updates differences between adjacent peak-peak distances
    
    Function that updates RR differences and RR squared differences 
    based on corrected RR list

    Parameters
    ----------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects
        will be created if not passed to function

    Returns
    -------
    out : dict
        working_data dictionary object containing all of heartpy's temp objects

    Examples
    --------
    Let's assume we detected peaks at these positions in the signal:
    
    >>> peaklist = [200, 280, 405, 410, 501, 615]
    
    And we subsequently ran further analysis:
    
    >>> wd = calc_rr(peaklist, sample_rate = 100.0)

    The peak at position 410 is likely an incorrect detection and 
    will be marked as such by other heartpy functions. This is indicated
    by an array 'binary_peaklist' in working_data. Binary peaklist is of
    the same length as peaklist, and is formatted as a mask:

    For now let's set it manually, normally this is done by the check_peaks()
    function from HeartPy's peakdetection module.
    
    >>> wd['binary_peaklist'] = [1, 1, 1, 0, 1, 1]

    Rejected peaks are marked with a zero and accepted with a 1.

    By now running update_rr(), heartpy will update all associated measures
    and will only compute peak-peak intervals between two accepted peaks.

    >>> wd = update_rr(wd)

    This will have generated a corrected RR_list object in the dictionary:
    
    >>> wd['RR_list_cor']
    [800.0, 1250.0, 1140.0]

    As well as updated the lists RR_diff (differences between adjacent peak-peak intervals) and 
    RR_sqdiff (squared differences between adjacent peak-peak intervals).
    '''
    rr_source = working_data['RR_list']
    b_peaklist = working_data['binary_peaklist']
    rr_list = [rr_source[i] for i in range(len(rr_source)) if b_peaklist[i] + b_peaklist[i+1] == 2]
    rr_mask = [0 if (b_peaklist[i] + b_peaklist[i+1] == 2) else 1 for i in range(len(rr_source))]
    rr_masked = np.ma.array(rr_source, mask=rr_mask)
    #TODO: don't do diff between rr-intervals where masked entry/entries lie in between
    rr_diff = np.abs(np.diff(rr_masked))
    rr_diff = rr_diff[~rr_diff.mask]
    rr_sqdiff = np.power(rr_diff, 2)
    
    working_data['RR_masklist'] = rr_mask
    working_data['RR_list_cor'] = rr_list
    working_data['RR_diff'] = rr_diff
    working_data['RR_sqdiff'] = rr_sqdiff

    return working_data


def calc_rr_segment(rr_source, b_peaklist):
    '''calculates peak-peak differences for segmentwise processing

    Function that calculates rr-measures when analysing segmentwise 
    in the 'fast' mode.

    Parameters
    ----------
    rr_source : 1d list or array
        list or array containing peak-peak intervals.

    b_peaklist : 1d list or array
        list or array containing mask for peaklist.

    Returns
    -------
    rr_list : array
        array containing peak-peak intervals.

    rr_diff : array
        array containing differences between adjacent peak-peak intervals

    rr_sqdiff : array
        array containing squared differences between adjacent peak-peak intervals

    Examples
    --------
    The function works in the same way as update_rr, except it returns
    three separate objects. It's used by process_segmentwise.
    Revert to doc on update_rr for more information.

    >>> rr, rrd, rrsd = calc_rr_segment(rr_source = [ 800., 1250.,   50.,  910., 1140., 1002., 1142.], 
    ... b_peaklist = [1, 1, 1, 0, 1, 1, 1, 1])
    >>> print(rr)
    [800.0, 1250.0, 1140.0, 1002.0, 1142.0]
    >>> print(rrd)
    [450.0 138.0 140.0]
    >>> print(rrsd)
    [202500.  19044.  19600.]
    '''
    rr_list = [rr_source[i] for i in range(len(rr_source)) if b_peaklist[i] + b_peaklist[i+1] == 2]
    rr_mask = [0 if (b_peaklist[i] + b_peaklist[i+1] == 2) else 1 for i in range(len(rr_source))]
    rr_masked = np.ma.array(rr_source, mask=rr_mask)
    rr_diff = np.abs(np.diff(rr_masked))
    rr_diff = rr_diff[~rr_diff.mask]
    rr_sqdiff = np.power(rr_diff, 2)
    
    return rr_list, rr_diff, rr_sqdiff


def clean_rr_intervals(sample_rate, working_data, method='iqr', calc_freq=False, freq_method='welch'):
    '''detects and rejects outliers in peak-peak intervals

    Function that detects and rejects outliers in the peak-peak intervals. It updates
    the RR_list_cor in the working data dict

    Parameters
    ----------
    sample_rate : int or float
        the sample rate with which the original signal was recorded in Hz.

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        Needs to contain RR_list_cor, meaning one analysis cycle has already completed.

    method : str
        which method to use for outlier rejection, included are 'iqr', based on
        the inter-quartile range, and 'z-score', which uses the modified z-score method.
        default : iqr

    calc_freq : bool
        whether to compute time-series measurements 
        default : False

    freq_method : str
        method used to extract the frequency spectrum. Available: 'fft' (Fourier Analysis), 
        'periodogram', and 'welch' (Welch's method). 
        default : 'welch'

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function.

    Returns
    -------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Examples
    --------
    Let's load some data

    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(1)
    >>> sample_rate = hp.get_samplerate_mstimer(timer)
    
    Run at least one analysis cycle first so that the dicts are populated

    >>> wd, m = hp.process(data, sample_rate)
    >>> wd = clean_rr_intervals(sample_rate, working_data = wd)
    >>> ['%.3f' %x for x in wd['RR_list_cor'][0:5]]
    ['897.470', '811.997', '829.091', '965.849', '803.449']

    You can also specify the outlier rejection method to be used:

    >>> wd = clean_rr_intervals(sample_rate, working_data = wd, method = 'z-score')
    >>> ['%.3f' %x for x in wd['RR_list_cor'][0:5]]
    ['897.470', '811.997', '829.091', '965.849', '803.449']
    '''

    #clean rr-list
    if method.lower() == 'iqr':
        rr_cleaned = outliers_iqr_method(working_data['RR_list_cor'])
    elif method.lower() == 'z-score':
        rr_cleaned = outliers_modified_z(working_data[ 'RR_list_cor'])

    rr_diff = np.diff(rr_cleaned)
    rr_sqdiff = np.power(rr_diff, 2)
    working_data['RR_list_cor'] = rr_cleaned
    working_data['RR_diff'] = rr_diff
    working_data['RR_sqdiff'] = rr_sqdiff

    #TODO: update rejected peaks for plotting!

    return working_data


def calc_ts_measures(rr_list, rr_diff, rr_sqdiff, measures={}, working_data={}):
    '''calculates standard time-series measurements.
    
    Function that calculates the time-series measurements for HeartPy.

    Parameters
    ----------
    rr_list : 1d list or array
        list or array containing peak-peak intervals

    rr_diff : 1d list or array
        list or array containing differences between adjacent peak-peak intervals

    rr_sqdiff : 1d list or array
        squared rr_diff

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function.

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.

    measures : dict
        dictionary object used by heartpy to store computed measures.

    Examples
    --------
    Normally this function is called during the process pipeline of HeartPy. It can
    of course also be used separately.

    Assuming we have the following peak-peak distances:

    >>> import numpy as np
    >>> rr_list = [1020.0, 990.0, 960.0, 1000.0, 1050.0, 1090.0, 990.0, 900.0, 900.0, 950.0, 1080.0]
    
    we can then compute the other two required lists by hand for now:

    >>> rr_diff = np.diff(rr_list)
    >>> rr_sqdiff = np.power(rr_diff, 2)
    >>> wd, m = calc_ts_measures(rr_list, rr_diff, rr_sqdiff)

    All output measures are then accessible from the measures object through
    their respective keys:

    >>> print('%.3f' %m['bpm'])
    60.384
    >>> print('%.3f' %m['rmssd'])
    67.082
    '''
    
    measures['bpm'] = 60000 / np.mean(rr_list)
    measures['ibi'] = np.mean(rr_list)

    ##TODO: 
    measures['sdnn'] = np.std(rr_list)
    measures['sdsd'] = np.std(rr_diff)
    measures['rmssd'] = np.sqrt(np.mean(rr_sqdiff))
    nn20 = rr_diff[np.where(rr_diff > 20.0)]
    nn50 = rr_diff[np.where(rr_diff > 50.0)]
    working_data['nn20'] = nn20
    working_data['nn50'] = nn50
    try:
        measures['pnn20'] = float(len(nn20)) / float(len(rr_diff))
    except:
        measures['pnn20'] = np.nan
    try:
        measures['pnn50'] = float(len(nn50)) / float(len(rr_diff))
    except:
        measures['pnn50'] = np.nan
    measures['hr_mad'] = MAD(rr_list)

    return working_data, measures


def calc_fd_measures(method='welch', square_spectrum=True, measures={}, working_data={}):
    '''calculates the frequency-domain measurements.

    Function that calculates the frequency-domain measurements for HeartPy.

    Parameters
    ----------
    method : str
        method used to compute the spectrogram of the heart rate.
        available methods: fft, periodogram, and welch
        default : welch

    square_spectrum : bool
        whether to square the power spectrum returned.
        default : true

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function.

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.

    measures : dict
        dictionary object used by heartpy to store computed measures.

    Examples
    --------
    Normally this function is called during the process pipeline of HeartPy. It can
    of course also be used separately.

    Let's load an example and get a list of peak-peak intervals

    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(0)
    >>> wd, m = hp.process(data, 100.0)
    
    wd now contains a list of peak-peak intervals that has been cleaned of
    outliers ('RR_list_cor'). Calling the function then is easy

    >>> wd, m = calc_fd_measures(method = 'periodogram', measures = m, working_data = wd)
    >>> print('%.3f' %m['lf/hf'])
    1.368

    Available methods are 'fft', 'welch' and 'periodogram'. To set another method, do:

    >>> wd, m = calc_fd_measures(method = 'fft', measures = m, working_data = wd)
    >>> print('%.3f' %m['lf/hf'])
    1.368
    '''
    rr_list = working_data['RR_list_cor']
    rr_x = []
    pointer = 0
    for x in rr_list:
        pointer += x
        rr_x.append(pointer)
    rr_x_new = np.linspace(rr_x[0], rr_x[-1], rr_x[-1])
    interpolated_func = UnivariateSpline(rr_x, rr_list, k=3)
    
    if method=='fft':
        datalen = len(rr_x_new)
        frq = np.fft.fftfreq(datalen, d=((1/1000.0)))
        frq = frq[range(int(datalen/2))]
        Y = np.fft.fft(interpolated_func(rr_x_new))/datalen
        Y = Y[range(int(datalen/2))]
        psd = np.power(Y, 2)
    elif method=='periodogram':
        frq, psd = periodogram(interpolated_func(rr_x_new), fs=1000.0)
    elif method=='welch':
        frq, psd = welch(interpolated_func(rr_x_new), fs=1000.0, nperseg=len(rr_x_new) - 1)
    else:
        raise ValueError("specified method incorrect, use 'fft', 'periodogram' or 'welch'")
    
    working_data['frq'] = frq
    working_data['psd'] = psd
    measures['lf'] = np.trapz(abs(psd[(frq >= 0.04) & (frq <= 0.15)]))
    measures['hf'] = np.trapz(abs(psd[(frq >= 0.16) & (frq <= 0.5)]))
    measures['lf/hf'] = measures['lf'] / measures['hf']
    working_data['interp_rr_function'] = interpolated_func
    working_data['interp_rr_linspace'] = (rr_x[0], rr_x[-1], rr_x[-1])
    return working_data, measures


def calc_breathing(rrlist, hrdata, sample_rate, measures={}, working_data={}):
    '''estimates breathing rate

    Function that estimates breathing rate from heart rate signal. 
    Upsamples the list of detected rr_intervals by interpolation then 
    tries to extract breathing peaks in the signal.

    Parameters
    ----------
    rr_list : 1d list or array
        list or array containing peak-peak intervals

    hrdata : 1d array or list
        sequence containing raw heart rate data

    sample_rate : int or float
        sample rate with which the heart rate data was measured.

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function.

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    measures : dict
        dictionary object used by heartpy to store computed measures.

    Examples
    --------
    Normally this function is called during the process pipeline of HeartPy. It can
    of course also be used separately.

    Let's load an example and get a list of peak-peak intervals

    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(0)
    >>> wd, m = hp.process(data, 100.0)

    Breathing is then computed with the function

    >>> m = calc_breathing(wd['RR_list_cor'], data, sample_rate = 100.0, measures = m, working_data = wd)
    >>> round(m['breathingrate'], 3)
    0.161

    There we have it, .16Hz, or about one breathing cycle in 6.25 seconds.
    '''
    x = np.linspace(0, len(rrlist), len(rrlist))
    x_new = np.linspace(0, len(rrlist), len(rrlist)*10)
    interp = UnivariateSpline(x, rrlist, k=3)
    breathing = interp(x_new)
    breathing_rolling_mean = rolling_mean(breathing, 0.75, sample_rate)
    peaks, working_data = hp.peakdetection.detect_peaks(breathing, breathing_rolling_mean, 1, sample_rate, 
                                                        update_dict=False)
    
    if len(peaks) > 1:
        signaltime = len(hrdata) / sample_rate
        measures['breathingrate'] = len(peaks) / signaltime
    else:
        measures['breathingrate'] = np.nan # pragma: no cover

    return measures