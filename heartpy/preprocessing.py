import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d

np.seterr(divide='ignore') #disable div by zero warnings
np.seterr(invalid='ignore')

from .filtering import filter_signal


__all__ = ['scale_data',
           'scale_sections',
           'enhance_peaks',
           'enhance_ecg_peaks',
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
    When passing data without further arguments to the function means it scales 0-1024
    
    >>> x = [2, 3, 4, 5]
    >>> scale_data(x)
    array([   0.        ,  341.33333333,  682.66666667, 1024.        ])

    Or you can specify a range:

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
    Import heartpy and load example data

    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(example=2)

    Let's slice a part of the data that I know contains clipping

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
    First let's load some example data:

    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(example=2)
    >>> x = data[2000:3000]
    >>> x[425:445]
    array([948, 977, 977, 977, 977, 978, 978, 977, 978, 977, 977, 977, 977,
           914, 820, 722, 627, 536, 460, 394])

    And interpolate any clipping segments as such:

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
    hrdata : 1-d numpy array or list 
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


def enhance_ecg_peaks(hrdata, sample_rate, iterations=4, aggregation='mean',
                      notch_filter=True):
    '''enhances ECG peaks

    Function that convolves synthetic QRS templates with the signal, leading
    to a strong increase signal-to-noise ratio. Function ends with an optional
    Notch filterstep (default : true) to reduce noise from the iterating
    convolution steps.

    Parameters
    ----------
    hrdata : 1-d numpy array or list 
        sequence containing heart rate data

    sample_rate : int or float
        sample rate with which the data is sampled 

    iterations : int
        how many convolutional iterations should be run. More will result in
        stronger peak enhancement, but over a certain point (usually between 12-16)
        overtones start appearing in the signal. Only increase this if the peaks
        aren't amplified enough.
        default : 4

    aggregation : str
        how the data from the different convolutions should be aggregated.
        Can be either 'mean' or 'median'.
        default : mean

    notch_filter : bool
        whether to apply a notch filter after the last convolution to get rid of
        remaining low frequency noise.
        default : true

    Returns
    -------
    output : 1d array
        The array containing the filtered data with enhanced peaks

    Examples
    --------
    First let's import the module and load the data

    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(1)
    >>> sample_rate = hp.get_samplerate_mstimer(timer)

    After loading the data we call the function like so:

    >>> filtered_data = enhance_ecg_peaks(data, sample_rate, iterations = 3)

    By default the module uses the mean to aggregate convolutional outputs. It
    is also possible to use the median.

    >>> filtered_data = enhance_ecg_peaks(data, sample_rate, iterations = 3,
    ... aggregation = 'median', notch_filter = False)

    In the last example we also disabled the notch filter.
    '''

    #assign output
    output = np.copy(hrdata)

    #generate synthetic QRS complexes
    templates = generate_ecg_templates(sample_rate)

    for i in range(int(iterations)):
        convolved = denoise_convolutions(output, sample_rate, templates)
        if aggregation == 'mean':
            output = np.nanmean(convolved, axis=0)
        elif aggregation == 'median':
            output = np.nanmedian(convolved, axis=0)

    #offset signal shift (shifts 1 datapoint for every iteration after the first)
    output = output[int(iterations) - 1:-int(iterations)]

    if notch_filter:
        output = filter_signal(output, 0.05, sample_rate, filtertype='notch')

    return output


def generate_ecg_templates(sample_rate, widths=[50, 60, 70, 80, 100],
                           presets=[[0, 2, 2.5, 3, 3.5, 5],
                                    [0, 1, 1.5, 2, 2.5, 3],
                                    [0, 3, 3.5, 4, 4.5, 6]],
                           amplitude = [0, -0.1, 1, -0.5, 0, 0]):
    '''helper function for enhance_ecg_peaks
    
    Helper function that generates synthetic QRS complexes of varying sizes
    to convolve with the signal.
    '''

    templates = []

    for i in presets:
        for j in widths:
            #duration < 120ms
            duration = (j / 1000) * sample_rate
            step = duration / len(i)
            new_t = [int(step * x) for x in i]
            new_x = np.linspace(new_t[0], new_t[-1], new_t[-1])
            #interpolate peak to fit in correct sampling rate
            interp_func = interp1d(new_t, amplitude, kind='linear')
            templates.append(interp_func(new_x))
    
    return templates


def denoise_convolutions(data, sample_rate, templates):
    '''helper function for enhance_ecg_peaks
    
    Helper function that convolves the generated synthetic QRS templates
    with the provided signal.
    '''

    convolutions = []

    for i in range(len(templates)):
            convolved = np.convolve(data, templates[i], mode='same')
            convolutions.append(convolved)
                         
    return np.asarray(convolutions)