'''
Functions to handle common preprocessing tasks

Data scaling
------------
- 'scale_data' -- scales given array to be between specified values
                  default range(0,1024)
- 'scale_section' -- applies scale_data within a moving window function
- 'enhance_peaks' -- function to apply some linear scaling in order
                     to bring out peaks.
- 'flip_signal' -- inverts the signal. Useful when using certain types of ECG


Clipping tools
--------------
- 'interpolate_clipping' -- function able to find and correct clipping sections of
                            the passed heart rate signal.

Hidden helper functions
-----------------------
- 'mark_clipping' -- marks beginning and end of clipping segments. Used by
                     interpolate_clipping() to find clipping segments.
'''

import numpy as np
from scipy.interpolate import UnivariateSpline

from .filtering import filter_signal


__all__ = ['scale_data',
           'scale_sections',
           'enhance_peaks'
           'interpolate_clipping',
           'flip_signal']


def scale_data(data, lower=0, upper=1024):
    '''scales passed sequence between thresholds

    Function that scales passed data so that it has specified lower 
    and upper bounds.
    
    Parameters
    ----------
    data : 1-d array or list
        Sequence to be scaled

    lower : int or float
        lower threshold for scaling
        default : 0

    upper : int or float
        upper threshold for scaling
        default : 1024

    Returns
    -------
    out : 1-d array
        contains scaled data

    Examples
    --------
    >>> x = [2, 3, 4, 5]

    Only passing data to the function means it scales 0-1024
    >>> scale_data(x)
    array([   0.        ,  341.33333333,  682.66666667, 1024.        ])

    Or we can specify our own range.
    >>> scale_data(x, lower = 50, upper = 124)
    array([ 50.        ,  74.66666667,  99.33333333, 124.        ])
    '''

    rng = np.max(data) - np.min(data)
    minimum = np.min(data)
    data = (upper - lower) * ((data - minimum) / rng) + lower
    return data


def scale_sections(data, sample_rate, windowsize=2.5, lower=0, upper=1024):
    '''scales data using sliding window approach

    Function that scales the data within the defined sliding window between 
    the defined lower and upper bounds.

    Parameters
    ----------
    data : 1-d array or list
        Sequence to be scaled

    sample_rate : int or float
        Sample rate of the passed signal

    windowsize : int or float
        size of the window within which signal is scaled, in seconds
        default : 2.5

    lower : int or float
        lower threshold for scaling. Passed to scale_data.
        default : 0

    upper : int or float
        upper threshold for scaling. Passed to scale_data.
        default : 1024

    Returns
    -------
    out : 1-d array
        contains scaled data

    Examples
    --------
    >>> x = [20, 30, 20, 30, 70, 80, 20, 30, 20, 30]
    >>> scale_sections(x, sample_rate=1, windowsize=2, lower=20, upper=30)
    array([20., 30., 20., 30., 20., 30., 20., 30., 20., 30.])
    '''

    total_length = len(data) / sample_rate
    window_dimension = int(windowsize * sample_rate)
    
    data_start = 0
    data_end = window_dimension
    
    output = np.empty(len(data))
    
    while data_end <= len(data):
        sliced = data[data_start:data_end]
        sliced = np.power(sliced, 2)
        scaled = scale_data(sliced, lower, upper)
        
        output[data_start:data_end] = scaled
        data_start += window_dimension
        data_end += window_dimension
        
    return np.array(output)


def mark_clipping(data, threshold=1020):
    '''marks clipping sections
    
    Function that marks start and end of clipping part
    it detects the start and end of clipping segments and returns them

    Parameters
    ----------
    data : 1-d numpy array
        Sequence to be scaled

    threshold: int or float
        the threshold for clipping, recommended to
        be a few data points below ADC or sensor max value, 
        to compensate for signal noise 
        default : 1020

    Returns
    -------
    out : list of tuples
        the output is a list of tuples. Each tuple marks the start
        and endpoint of the detected clipping segment

    Examples
    --------
    >>> from heartpy import datautils
    >>> data, _ = datautils.load_exampledata(example=2)
    >>> x = data[2000:3000]
    >>> mark_clipping(x, threshold=970)
    [(369, 375), (426, 437), (486, 493), (544, 552), (604, 610), (663, 665), \
(721, 722), (776, 781), (831, 836), (883, 891), (995, 999)]
    '''

    clip_binary = np.where(data > threshold)
    clipping_edges = np.where(np.diff(clip_binary) > 1)[1]

    clipping_segments = []

    for i in range(0, len(clipping_edges)):
        if i == 0: #if first clipping segment
            clipping_segments.append((clip_binary[0][0], 
                                      clip_binary[0][clipping_edges[0]]))
        elif i == len(clipping_edges) - 1:
            #append last entry
            clipping_segments.append((clip_binary[0][clipping_edges[i]+1],
                                      clip_binary[0][-1]))    
        else:
            clipping_segments.append((clip_binary[0][clipping_edges[i-1] + 1],
                                      clip_binary[0][clipping_edges[i]]))

    return clipping_segments


def interpolate_clipping(data, sample_rate, threshold=1020):
    '''interpolate peak waveform

    Function that interpolates peaks between the clipping segments using 
    cubic spline interpolation. It takes the clipping start +/- 100ms to 
    calculate the spline.
    
    Parameters
    ----------
    data : 1d list or numpy array
        data section to be evaluated 

    sample_rate : int or float
        sample rate with which the data array is sampled

    threshold : int or float
        the threshold for clipping, recommended to
        be a few data points below ADC or sensor max value, 
        to compensate for signal noise 
        default : 1020

    Returns
    -------
    out : array 
        the output is an array with clipping segments replaced
        by interpolated segments

    Examples
    --------
    >>> from heartpy import datautils
    >>> data, _ = datautils.load_exampledata(example=2)
    >>> x = data[2000:3000]
    >>> x[425:445]
    array([948, 977, 977, 977, 977, 978, 978, 977, 978, 977, 977, 977, 977,
           914, 820, 722, 627, 536, 460, 394])
    >>> intp = interpolate_clipping(x, sample_rate=117, threshold=970)
    >>> intp[425:445]
    array([ 972, 1043, 1098, 1138, 1163, 1174, 1173, 1159, 1134, 1098, 1053,
            998,  934,  848,  747,  646,  552,  470,  402,  348])
    '''

    clipping_segments = mark_clipping(data, threshold)
    num_datapoints = int(0.1 * sample_rate)
    newx = []
    newy = []
    
    i = 0
    
    for segment in clipping_segments:
        if segment[0] < num_datapoints: 
            #if clipping is present at start of signal, skip.
            #We cannot interpolate accurately when there is insufficient data prior to clipping segment.
            pass
        else: 
            antecedent = data[segment[0] - num_datapoints : segment[0]]
            consequent = data[segment[1] : segment[1] + num_datapoints]
            segment_data = np.concatenate((antecedent, consequent))
        
            interpdata_x = np.concatenate(([x for x in range(segment[0] - num_datapoints, segment[0])],
                                            [x for x in range(segment[1], segment[1] + num_datapoints)]))
            x_new = np.linspace(segment[0] - num_datapoints,
                                segment[1] + num_datapoints,
                                ((segment[1] - segment[0]) + (2 * num_datapoints)))
        
            try:
                interp_func = UnivariateSpline(interpdata_x, segment_data, k=3)
                interp_data = interp_func(x_new)
        
                data[segment[0] - num_datapoints :
                     segment[1] + num_datapoints] = interp_data
            except:
                #pass over failed interpolation: leave original data alone
                pass
       
    return data


def flip_signal(data, enhancepeaks=False, keep_range=True):
    '''invert signal waveforms.

    Function that flips raw signal with negative mV peaks to normal ECG.
    Required for proper peak finding in case peaks are expressed as
    negative dips.

    Parameters
    ----------
    data : 1d list or numpy array
        data section to be evaluated 
    
    enhance_peaks : bool
        whether to apply peak accentuation 
        default : False

    keep_range : bool
        whether to scale the inverted data so that the original
        range is maintained

    Returns
    -------
    out : 1d array

    Examples
    --------
    Given an array of data
    >>> x = [200, 300, 500, 900, 500, 300, 200]

    We can call the function. If keep_range is False, the signal
    will be inverted relative to its mean.
    >>> flip_signal(x, keep_range=False)
    array([628.57142857, 528.57142857, 328.57142857, -71.42857143,
           328.57142857, 528.57142857, 628.57142857])

    However, by specifying keep_range, the inverted signal will be
    put 'back in place' in its original range.
    >>> flip_signal(x, keep_range=True)
    array([900., 800., 600., 200., 600., 800., 900.])

    It's also possible to use the enhance_peaks function:
    >>> flip_signal(x, enhancepeaks=True)
    array([1024.        ,  621.75746332,  176.85545623,    0.        ,
            176.85545623,  621.75746332, 1024.        ])
    '''
    data_mean = np.mean(data)
    data_min = np.min(data)
    data_max = np.max(data)

    #invert signal
    data = (data_mean - data) + data_mean
    
    if keep_range:
        #scale data so original range is maintained
        data = scale_data(data, lower = data_min, upper = data_max)
    if enhancepeaks:
        data = enhance_peaks(data)
    return data


def enhance_peaks(hrdata, iterations=2):
    '''enhances peak amplitude relative to rest of signal
    
    Function thta attempts to enhance the signal-noise ratio by accentuating 
    the highest peaks. Note: denoise first
    
    Parameters
    ----------
    data : 1-d numpy array or list 
        sequence containing heart rate data
    iterations : int
        the number of scaling steps to perform 
        default : 2

    Returns
    -------
    out : 1-d numpy array
        array containing enhanced peaks

    Examples
    --------
    Given an array of data, the peaks can be enhanced using the function
    >>> x = [200, 300, 500, 900, 500, 300, 200]
    >>> enhance_peaks(x)
    array([   0.        ,    4.31776016,   76.16528926, 1024.        ,
             76.16528926,    4.31776016,    0.        ])
    '''
    scale_data(hrdata)
    for i in range(iterations):
        hrdata = np.power(hrdata, 2)
        hrdata = scale_data(hrdata)
    return hrdata  