'''
Functions that handle computation of heart rate (HR) and
heart rate variability (HRV) measures.
'''

import warnings

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.signal import welch, periodogram

from .datautils import MAD, rolling_mean, outliers_iqr_method, outliers_modified_z
from .filtering import quotient_filter
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
    rr_indices = [(peaklist[i], peaklist[i+1]) for i in range(len(peaklist) - 1)]
    rr_diff = np.abs(np.diff(rr_list))
    rr_sqdiff = np.power(rr_diff, 2)
    working_data['RR_list'] = rr_list
    working_data['RR_indices'] = rr_indices
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


def clean_rr_intervals(working_data, method='quotient-filter'):
    '''detects and rejects outliers in peak-peak intervals

    Function that detects and rejects outliers in the peak-peak intervals. It updates
    the RR_list_cor in the working data dict

    Parameters
    ----------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        Needs to contain RR_list_cor, meaning one analysis cycle has already completed.

    method : str
        which method to use for outlier rejection, included are:
        - 'quotient-filter', based on the work in "Piskorki, J., Guzik, P. (2005), Filtering Poincare plots",
        - 'iqr', which uses the inter-quartile range, 
        - 'z-score', which uses the modified z-score method.
        default : quotient-filter

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
    >>> wd = clean_rr_intervals(working_data = wd)
    >>> ['%.3f' %x for x in wd['RR_list_cor'][0:5]]
    ['897.470', '811.997', '829.091', '777.807', '803.449']

    You can also specify the outlier rejection method to be used, for example using 
    the z-score method:

    >>> wd = clean_rr_intervals(working_data = wd, method = 'z-score')
    >>> ['%.3f' %x for x in wd['RR_list_cor'][0:5]]
    ['897.470', '811.997', '829.091', '777.807', '803.449']

    Or the inter-quartile range (iqr) based method:
    
    >>> wd = clean_rr_intervals(working_data = wd, method = 'iqr')
    >>> ['%.3f' %x for x in wd['RR_list_cor'][0:5]]
    ['897.470', '811.997', '829.091', '965.849', '803.449']
    '''

    #generate RR_list_cor indices relative to RR_list
    RR_cor_indices = [i for i in range(len(working_data['RR_masklist']))
                       if working_data['RR_masklist'][i] == 0]

    #clean rr-list
    if method.lower() == 'iqr':
        rr_cleaned, replaced_indices = outliers_iqr_method(working_data['RR_list_cor'])
        rr_mask = working_data['RR_masklist']
        for i in replaced_indices:
            rr_mask[RR_cor_indices[i]] = 1

    elif method.lower() == 'z-score':
        rr_cleaned, replaced_indices = outliers_modified_z(working_data['RR_list_cor'])
        rr_mask = working_data['RR_masklist']
        for i in replaced_indices:
            rr_mask[RR_cor_indices[i]] = 1

    elif method.lower() == 'quotient-filter':
        rr_mask = quotient_filter(working_data['RR_list'], working_data['RR_masklist'])
        rr_cleaned = [x for x,y in zip(working_data['RR_list'], rr_mask) if y == 0]


    else:
        raise ValueError('Incorrect method specified, use either "iqr", "z-score" or "quotient-filtering". \
Nothing to do!')

    rr_masked = np.ma.array(working_data['RR_list'], mask=rr_mask)
    rr_diff = np.abs(np.diff(rr_masked))
    rr_diff = rr_diff[~rr_diff.mask]
    rr_sqdiff = np.power(rr_diff, 2)
    working_data['RR_masked'] = rr_masked
    working_data['RR_list_cor'] = np.asarray(rr_cleaned)
    working_data['RR_diff'] = rr_diff
    working_data['RR_sqdiff'] = rr_sqdiff

    try:
        removed_beats = [x for x in working_data['removed_beats']]
        removed_beats_y = [x for x in working_data['removed_beats_y']]
        peaklist = working_data['peaklist']
        ybeat = working_data['ybeat']

        for i in range(len(rr_mask)):
            if rr_mask[i] == 1 and peaklist[i] not in removed_beats:
                removed_beats.append(peaklist[i])
                removed_beats_y.append(ybeat[i])

        working_data['removed_beats'] = np.asarray(removed_beats)
        working_data['removed_beats_y'] = np.asarray(removed_beats_y)
    except:
        pass

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
    >>> data, timer = hp.load_exampledata(2)
    >>> sample_rate = hp.get_samplerate_datetime(timer, timeformat='%Y-%m-%d %H:%M:%S.%f')
    >>> wd, m = hp.process(data, sample_rate)
    
    wd now contains a list of peak-peak intervals that has been cleaned of
    outliers ('RR_list_cor'). Calling the function then is easy

    >>> wd, m = calc_fd_measures(method = 'periodogram', measures = m, working_data = wd)
    >>> print('%.3f' %m['lf/hf'])
    4.964

    Available methods are 'fft', 'welch' and 'periodogram'. To set another method, do:

    >>> wd, m = calc_fd_measures(method = 'fft', measures = m, working_data = wd)
    >>> print('%.3f' %m['lf/hf'])
    4.964

    If there are no valid peak-peak intervals specified, returned measures are NaN:
    >>> wd['RR_list_cor'] = []
    >>> wd, m = calc_fd_measures(working_data = wd)
    >>> np.isnan(m['lf/hf'])
    True

    If there are rr-intervals but not enough to reliably compute frequency measures, a
    warning is raised:

    --------------
    RuntimeWarning: Short signal.
    ---------Warning:---------
    too few peak-peak intervals for (reliable) frequency domain measure computation,
    frequency output measures are still computed but treat them with caution!

    HF is usually computed over a minimum of 1 minute of good signal.
    LF is usually computed over a minimum of 2 minutes of good signal.
    The LF/HF ratio is usually computed over minimum 24 hours, although an
    absolute minimum of 5 min has also been suggested.

    For more info see: \nShaffer, F., Ginsberg, J.P. (2017).
    An Overview of Heart Rate Variability Metrics and Norms.

    Task Force of Pacing and Electrophysiology (1996), Heart Rate Variability
    in: European Heart Journal, vol.17, issue 3, pp354-381
    
    
    This warning will not repeat'
    --------------
    '''
    rr_list = working_data['RR_list_cor']

    if len(rr_list) <= 1:
        working_data['frq'] = np.nan
        working_data['psd'] = np.nan
        measures['lf'] = np.nan
        measures['hf'] = np.nan
        measures['lf/hf'] = np.nan
        return working_data, measures
    elif np.sum(rr_list) <= 300000: # pragma: no cover
        #warn if signal is short
        msg = ''.join(('Short signal.\n',
                       '\n---------Warning:---------\n',
                       'too few peak-peak intervals for (reliable) frequency domain measure computation, ',
                       'frequency output measures are still computed but treat them with caution!\n\n',
                       'HF is usually computed over a minimum of 1 minute of good signal. ',
                       'LF is usually computed over a minimum of 2 minutes of good signal.',
                       'The LF/HF ratio is usually computed over minimum 24 hours, although an ',
                       'absolute minimum of 5 min has also been suggested.\n\n',
                       'For more info see: \nShaffer, F., Ginsberg, J.P. (2017), ',
                       'An Overview of Heart Rate Variability Metrics and Norms.\n\n',
                       'Task Force of Pacing and Electrophysiology (1996), Heart Rate Variability, ',
                       'in: European Heart Journal, vol.17, issue 3, pp354-381'
                       '\n\nThis warning will not repeat'))
        warnings.warn(msg, UserWarning)

    rr_x = []
    pointer = 0
    for x in rr_list:
        pointer += x
        rr_x.append(pointer)
    rr_x_new = np.linspace(int(rr_x[0]), int(rr_x[-1]), int(rr_x[-1]))
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


def calc_breathing(rrlist, hrdata, sample_rate, method='fft',
                   measures={}, working_data={}):
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

    >>> m, wd = calc_breathing(wd['RR_list_cor'], data, sample_rate = 100.0, measures = m, working_data = wd)
    >>> round(m['breathingrate'], 3)
    0.128

    There we have it, .16Hz, or about one breathing cycle in 6.25 seconds.
    '''

    #resample RR-list to 1000Hz
    x = np.linspace(0, len(rrlist), len(rrlist))
    x_new = np.linspace(0, len(rrlist), np.sum(rrlist, dtype=np.int32))
    interp = UnivariateSpline(x, rrlist, k=3)
    breathing = interp(x_new)

    if method.lower() == 'fft':
        datalen = len(breathing)
        frq_ = np.fft.fftfreq(datalen, d=((1/1000.0)))
        frq_ = frq_[range(int(datalen/2))]
        Y = np.fft.fft(breathing)/datalen
        Y = Y[range(int(datalen/2))]
        psd_ = np.power(np.abs(Y), 2)
    elif method.lower() == 'welch':
        frq_, psd_ = welch(breathing, fs=1000, nperseg=len(breathing))
    else:
        raise ValueError('Breathing rate extraction method not understood! Must be \'welch\' or \'fft\'!')

    #take out lowest peak
    frq = frq_[frq_ >= 0.04]
    psd = psd_[frq_ >= 0.04]
    
    #find max
    measures['breathingrate'] = frq[np.argmax(psd)]
    working_data['breathing_signal'] = breathing
    working_data['breathing_psd'] = psd
    working_data['breathing_frq'] = frq

    return measures, working_data


def calc_poincare(rr_list, rr_mask=[], measures={}, working_data={}):
    '''computes poincare parameters

    Function that takes peak-peak intervals and computes poincare parameters:
    [0] standard deviation perpendicular to identity line (SD1)
    [1] standard deviation along identity line (SD2)
    [2] area of ellipse described by SD1 and SD2
    [3] SD1/SD2 ratio

    Based on:
    "Shaffer, F., Ginsberg, J.P. (2017), An Overview of Heart Rate
    Variability Metrics and Norms"

    Parameters
    ----------
    rr_list : 1d array or list
        list or array containing peak-peak intervals

    rr_mask : 1d array or list
        list or array containing mask for rejected peak-peak intervals

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
        poincare values are appended to measures['poincare']
    '''

    #generate vectors of adjacent peak-peak intervals
    x_plus = []
    x_minus = []

    for i in range(len(working_data['RR_masklist']) - 1):
        if working_data['RR_masklist'][i] + working_data['RR_masklist'][i + 1] == 0:
            #only add adjacent RR-intervals that are not rejected
            x_plus.append(working_data['RR_list'][i])
            x_minus.append(working_data['RR_list'][i + 1])
        else:
            pass
    
    #cast to arrays so we can do numerical work easily
    x_plus = np.asarray(x_plus)
    x_minus = np.asarray(x_minus)

    #compute parameters and append to dict
    x_one = (x_plus - x_minus) / np.sqrt(2)
    x_two = (x_plus + x_minus) / np.sqrt(2)
    sd1 = np.sqrt(np.var(x_one)) #compute stdev perpendicular to identity line
    sd2 = np.sqrt(np.var(x_two)) #compute stdev parallel to identity line
    s = np.pi * sd1 * sd2 #compute area of ellipse

    #write computed measures to dicts
    measures['sd1'] = sd1
    measures['sd2'] = sd2 
    measures['s'] = s
    measures['sd1/sd2'] = sd1 / sd2

    working_data['poincare'] = {}
    working_data['poincare']['x_plus'] = x_plus
    working_data['poincare']['x_minus'] = x_minus
    working_data['poincare']['x_one'] = x_one
    working_data['poincare']['x_two'] = x_two

    return measures
