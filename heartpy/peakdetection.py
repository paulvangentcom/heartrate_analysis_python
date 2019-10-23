'''
functions for peak detection and related tasks
'''

import numpy as np
from scipy.signal import resample

from .analysis import calc_rr, update_rr
from .exceptions import BadSignalWarning
from .filtering import quotient_filter


__all__ = ['make_windows',
           'append_dict',
           'detect_peaks',
           'fit_peaks',
           'check_peaks',
           'check_binary_quality',
           'interpolate_peaks']


def make_windows(data, sample_rate, windowsize=120, overlap=0, min_size=20):
    '''slices data into windows

    Funcion that slices data into windows for concurrent analysis. 
    Used by process_segmentwise wrapper function.
    
    Parameters
    ----------
    data : 1-d array 
        array containing heart rate sensor data

    sample_rate : int or float
        sample rate of the data stream in 'data'

    windowsize : int 
        size of the window that is sliced in seconds

    overlap : float
        fraction of overlap between two adjacent windows: 0 <= float < 1.0

    min_size : int 
        the minimum size for the last (partial) window to be included. Very short windows
        might not stable for peak fitting, especially when significant noise is present. 
        Slightly longer windows are likely stable but don't make much sense from a 
        signal analysis perspective.
    
    Returns
    -------
    out : array
        tuples of window indices

    Examples
    --------
    Assuming a given example data file:

    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(1)
    
    We can split the data into windows:

    >>> indices = make_windows(data, 100.0, windowsize = 30, overlap = 0.5, min_size = 20)
    >>> indices.shape
    (9, 2)

    Specifying min_size = -1 will include the last window no matter what:

    >>> indices = make_windows(data, 100.0, windowsize = 30, overlap = 0.5, min_size = -1)
    '''
    ln = len(data)
    window = windowsize * sample_rate
    stepsize = (1 - overlap) * window
    start = 0
    end = window
    
    slices = []
    while end < len(data):
        slices.append((start, end))
        start += stepsize
        end += stepsize
    
    if min_size == -1: 
        slices[-1] = (slices[-1][0], len(data))
    elif (ln - start) / sample_rate >= min_size:
        slices.append((start, ln))
        
    return np.array(slices, dtype=np.int32)


def append_dict(dict_obj, measure_key, measure_value):
    '''appends data to keyed dict.
    
    Function that appends key to continuous dict, creates if doesn't exist. EAFP

    Parameters
    ----------
    dict_obj : dict
        dictionary object that contains continuous output measures

    measure_key : str 
        key for the measure to be stored in continuous_dict
    
    measure_value : any data container
        value to be appended to dictionary

    Returns
    -------
    dict_obj : dict
        dictionary object passed to function, with specified data container appended

    Examples
    --------
    Given a dict object 'example' with some data in it:

    >>> example = {}
    >>> example['call'] = ['hello']

    We can use the function to append it:

    >>> example = append_dict(example, 'call', 'world')
    >>> example['call']
    ['hello', 'world']

    A new key will be created if it doesn't exist:

    >>> example = append_dict(example, 'different_key', 'hello there!')
    >>> sorted(example.keys())
    ['call', 'different_key']
    ''' 
    try:
        dict_obj[measure_key].append(measure_value)
    except KeyError:
        dict_obj[measure_key] = [measure_value]
    return dict_obj


def detect_peaks(hrdata, rol_mean, ma_perc, sample_rate, update_dict=True, working_data={}):
    '''detect peaks in signal

    Function that detects heartrate peaks in the given dataset.

    Parameters
    ----------
    hr data : 1-d numpy array or list 
        array or list containing the heart rate data

    rol_mean : 1-d numpy array 
        array containing the rolling mean of the heart rate signal

    ma_perc : int or float
        the percentage with which to raise the rolling mean,
        used for fitting detection solutions to data

    sample_rate : int or float
        the sample rate of the provided data set

    update_dict : bool
        whether to update the peak information in the module's data structure
        Settable to False to allow this function to be re-used for example by 
        the breath analysis module.
        default : True

    Examples
    --------
    Normally part of the peak detection pipeline. Given the first example data
    it would work like this:

    >>> import heartpy as hp
    >>> from heartpy.datautils import rolling_mean, _sliding_window
    >>> data, _ = hp.load_exampledata(0)
    >>> rol_mean = rolling_mean(data, windowsize = 0.75, sample_rate = 100.0)
    >>> wd = detect_peaks(data, rol_mean, ma_perc = 20, sample_rate = 100.0)

    Now the peaklist has been appended to the working data dict. Let's look
    at the first five peak positions:

    >>> wd['peaklist'][0:5]
    [63, 165, 264, 360, 460]
    '''
    rmean = np.array(rol_mean)

    #rol_mean = rmean + ((rmean / 100) * ma_perc)
    mn = np.mean(rmean / 100) * ma_perc
    rol_mean = rmean + mn

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

    if update_dict:
        working_data['peaklist'] = peaklist
        working_data['ybeat'] = [hrdata[x] for x in peaklist]
        working_data['rolling_mean'] = rol_mean
        working_data = calc_rr(working_data['peaklist'], sample_rate,
                               working_data=working_data)
        if len(working_data['RR_list']) > 0:
            working_data['rrsd'] = np.std(working_data['RR_list'])
        else:
            working_data['rrsd'] = np.inf
        return working_data
    else:
        return peaklist, working_data


def fit_peaks(hrdata, rol_mean, sample_rate, bpmmin=40, bpmmax=180, working_data={}):
    '''optimize for best peak detection

    Function that runs fitting with varying peak detection thresholds given a 
    heart rate signal.
    
    Parameters
    ----------
    hrdata : 1d array or list 
        array or list containing the heart rate data

    rol_mean : 1-d array 
        array containing the rolling mean of the heart rate signal

    sample_rate : int or float
        the sample rate of the data set

    bpmmin : int
        minimum value of bpm to see as likely 
        default : 40

    bpmmax : int 
        maximum value of bpm to see as likely 
        default : 180

    Returns
    -------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Examples
    --------
    Part of peak detection pipeline. Uses moving average as a peak detection
    threshold and rises it stepwise. Determines best fit by minimising 
    standard deviation of peak-peak distances as well as getting a bpm that
    lies within the expected range.

    Given included example data let's show how this works

    >>> import heartpy as hp
    >>> from heartpy.datautils import rolling_mean, _sliding_window
    >>> data, _ = hp.load_exampledata(0)
    >>> rol_mean = rolling_mean(data, windowsize = 0.75, sample_rate = 100.0)

    We can then call this function and let the optimizer do its work:

    >>> wd = fit_peaks(data, rol_mean, sample_rate = 100.0)

    Now the wd dict contains the best fit paramater(s):

    >>> wd['best']
    20

    This indicates the best fit can be obtained by raising the moving average
    with 20%.
    
    The results of the peak detection using these parameters are included too.
    To illustrate, these are the first five detected peaks:

    >>> wd['peaklist'][0:5]
    [63, 165, 264, 360, 460]

    and the corresponding peak-peak intervals:

    >>> wd['RR_list'][0:4]
    array([1020.,  990.,  960., 1000.])
    '''
    ma_perc_list = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300]
    rrsd = []
    valid_ma = []
    for ma_perc in ma_perc_list:
        working_data = detect_peaks(hrdata, rol_mean, ma_perc, sample_rate, 
                                    update_dict=True, working_data=working_data)
        bpm = ((len(working_data['peaklist'])/(len(hrdata)/sample_rate))*60)
        rrsd.append([working_data['rrsd'], bpm, ma_perc])

    for _rrsd, _bpm, _ma_perc in rrsd:
        if (_rrsd > 0.1) and ((bpmmin <= _bpm <= bpmmax)):
            valid_ma.append([_rrsd, _ma_perc])

    if len(valid_ma) > 0:
        working_data['best'] = min(valid_ma, key=lambda t: t[0])[1]
        working_data = detect_peaks(hrdata, rol_mean, min(valid_ma, key=lambda t: t[0])[1], 
                                    sample_rate, update_dict=True, working_data=working_data)
        return working_data
    else:
        raise BadSignalWarning('\n----------------\nCould not determine best fit for \
given signal. Please check the source signal.\n Probable causes:\n- detected heart rate falls \
outside of bpmmin<->bpmmax constraints\n- no detectable heart rate present in signal\n\
- very noisy signal (consider filtering and scaling)\nIf you\'re sure the signal contains heart\
rate data, consider filtering and/or scaling first.\n----------------\n')


def check_peaks(rr_arr, peaklist, ybeat, quotient_filter=False, reject_segmentwise=False, 
                working_data={}):
    '''find anomalous peaks.

    Funcion that checks peaks for outliers based on anomalous peak-peak distances and corrects
    by excluding them from further analysis.

    Parameters
    ----------
    rr_arr : 1d array or list
        list or array containing peak-peak intervals

    peaklist : 1d array or list
        list or array containing detected peak positions

    ybeat : 1d array or list
        list or array containing corresponding signal values at 
        detected peak positions. Used for plotting functionality
        later on.

    reject_segmentwise : bool
        if set, checks segments per 10 detected peaks. Marks segment
        as rejected if 30% of peaks are rejected.
        default : False

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    working_data : dict
        working_data dictionary object containing all of heartpy's temp objects

    Examples
    --------
    Part of peak detection pipeline. No standalone examples exist. See docstring
    for hp.process() function for more info
    '''
    rr_arr = np.array(rr_arr)
    peaklist = np.array(peaklist)
    ybeat = np.array(ybeat)

    mean_rr = np.mean(rr_arr)
    upper_threshold = mean_rr + 300 if (0.3 * mean_rr) <= 300 else mean_rr + (0.3 * mean_rr)
    lower_threshold = mean_rr - 300 if (0.3 * mean_rr) <= 300 else mean_rr - (0.3 * mean_rr)
    
    working_data['removed_beats'] = peaklist[np.where((rr_arr <= lower_threshold) |
                                                      (rr_arr >= upper_threshold))[0]+1]
    working_data['removed_beats_y'] = ybeat[np.where((rr_arr <= lower_threshold) |
                                                     (rr_arr >= upper_threshold))[0]+1]
    working_data['binary_peaklist'] = np.asarray([0 if x in working_data['removed_beats'] 
                                                  else 1 for x in peaklist])

    if reject_segmentwise: 
        working_data = check_binary_quality(peaklist, working_data['binary_peaklist'],
                                            working_data=working_data)

    working_data = update_rr(working_data=working_data)
    return working_data


def check_binary_quality(peaklist, binary_peaklist, maxrejects=3, working_data={}):
    '''checks signal in chunks of 10 beats. 
    
    Function that checks signal in chunks of 10 beats. It zeros out chunk if 
    number of rejected peaks > maxrejects. Also marks rejected segment coordinates 
    in tuples (x[0], x[1] in working_data['rejected_segments']
    
    Parameters
    ----------
    peaklist : 1d array or list
        list or array containing detected peak positions

    binary_peaklist : 1d array or list
        list or array containing mask for peaklist, coding which peaks are rejected

    maxjerects : int
        maximum number of rejected peaks per 10-beat window 
        default : 3

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    working_data : dict
        working_data dictionary object containing all of heartpy's temp objects

    Examples
    --------
    Part of peak detection pipeline. No standalone examples exist. See docstring
    for hp.process() function for more info

    Given some peaklist and binary mask:
    >>> peaklist = [30, 60, 90, 110, 130, 140, 160, 170, 200, 220]
    >>> binary_peaklist = [0, 1, 1, 0, 0, 1, 0, 1, 0, 0]
    >>> wd = check_binary_quality(peaklist, binary_peaklist)
    >>> wd['rejected_segments']
    [(30, 220)]

    The whole segment is rejected as it contains more than the specified 3 rejections 
    per 10 beats.
    '''
    idx = 0
    working_data['rejected_segments'] = []
    for i in range(int(len(binary_peaklist) / 10)):
        if np.bincount(binary_peaklist[idx:idx + 10])[0] > maxrejects:
            binary_peaklist[idx:idx + 10] = [0 for i in range(len(binary_peaklist[idx:idx+10]))]
            if idx + 10 < len(peaklist): 
                working_data['rejected_segments'].append((peaklist[idx], peaklist[idx + 10]))
            else:
                working_data['rejected_segments'].append((peaklist[idx], peaklist[-1]))
        idx += 10
    return working_data


def interpolate_peaks(data, peaks, sample_rate, desired_sample_rate=1000.0, working_data={}):
    '''interpolate detected peak positions and surrounding data points

    Function that enables high-precision mode by taking the estimated peak position,
    then upsampling the peak position +/- 100ms to the specified sampling rate, subsequently
    estimating the peak position with higher accuracy.

    Parameters
    ----------
    data : 1d list or array 
        list or array containing heart rate data

    peaks : 1d list or array 
        list or array containing x-positions of peaks in signal
    
    sample_rate : int or float
        the sample rate of the signal (in Hz)

    desired_sampled-rate : int or float
        the sample rate to which to upsample. 
        Must be sample_rate < desired_sample_rate

    Returns
    -------
    working_data : dict
        working_data dictionary object containing all of heartpy's temp objects

    Examples
    --------
    Given the output of a normal analysis and the first five peak-peak intervals:

    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(0)
    >>> wd, m = hp.process(data, 100.0)
    >>> wd['peaklist'][0:5]
    [63, 165, 264, 360, 460]

    Now, the resolution is at max 10ms as that's the distance between data points. 
    We can use the high precision mode for example to approximate a more precise
    position, for example if we had recorded at 1000Hz:
    
    >>> wd = interpolate_peaks(data = data, peaks = wd['peaklist'], 
    ... sample_rate = 100.0, desired_sample_rate = 1000.0, working_data = wd)
    >>> wd['peaklist'][0:5]
    [63.5, 165.4, 263.6, 360.4, 460.2]

    As you can see the accuracy of peak positions has increased.
    Note that you cannot magically upsample nothing into something. Be reasonable.
    '''
    assert desired_sample_rate > sample_rate, "desired sample rate is lower than actual sample rate \
this would result in downsampling which will hurt accuracy."
    
    num_samples = int(0.1 * sample_rate)
    ratio = sample_rate / desired_sample_rate
    interpolation_slices = [(x - num_samples, x + num_samples) for x in peaks]
    peaks = []

    for i in interpolation_slices:
        slice = data[i[0]:i[1]]
        resampled = resample(slice, int(len(slice) * (desired_sample_rate / sample_rate)))
        peakpos = np.argmax(resampled)
        peaks.append((i[0] + (peakpos * ratio)))

    working_data['peaklist'] = peaks

    return working_data
