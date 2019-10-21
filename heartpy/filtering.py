'''
Functions for data filtering tasks.
'''

from scipy.signal import butter, filtfilt, iirnotch, savgol_filter
import numpy as np

from .datautils import MAD

__all__ = ['filter_signal',
           'hampel_filter',
           'hampel_correcter',
           'smooth_signal']

def butter_lowpass(cutoff, sample_rate, order=2):
    '''standard lowpass filter.

    Function that defines standard Butterworth lowpass filter

    Parameters
    ----------
    cutoff : int or float
        frequency in Hz that acts as cutoff for filter.
        All frequencies above cutoff are filtered out.

    sample_rate : int or float
        sample rate of the supplied signal

    order : int
        filter order, defines the strength of the roll-off
        around the cutoff frequency. Typically orders above 6
        are not used frequently.
        default: 2
    
    Returns
    -------
    out : tuple
        numerator and denominator (b, a) polynomials
        of the defined Butterworth IIR filter.

    Examples
    --------
    >>> b, a = butter_lowpass(cutoff = 2, sample_rate = 100, order = 2)
    >>> b, a = butter_lowpass(cutoff = 4.5, sample_rate = 12.5, order = 5)
    '''
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_highpass(cutoff, sample_rate, order=2):
    '''standard highpass filter.

    Function that defines standard Butterworth highpass filter

    Parameters
    ----------
    cutoff : int or float
        frequency in Hz that acts as cutoff for filter.
        All frequencies below cutoff are filtered out.

    sample_rate : int or float
        sample rate of the supplied signal

    order : int
        filter order, defines the strength of the roll-off
        around the cutoff frequency. Typically orders above 6
        are not used frequently.
        default : 2
    
    Returns
    -------
    out : tuple
        numerator and denominator (b, a) polynomials
        of the defined Butterworth IIR filter.

    Examples
    --------
    we can specify the cutoff and sample_rate as ints or floats.

    >>> b, a = butter_highpass(cutoff = 2, sample_rate = 100, order = 2)
    >>> b, a = butter_highpass(cutoff = 4.5, sample_rate = 12.5, order = 5)
    '''
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a


def butter_bandpass(lowcut, highcut, sample_rate, order=2):
    '''standard bandpass filter.
    Function that defines standard Butterworth bandpass filter.
    Filters out frequencies outside the frequency range
    defined by [lowcut, highcut].

    Parameters
    ----------
    lowcut : int or float
        Lower frequency bound of the filter in Hz

    highcut : int or float
        Upper frequency bound of the filter in Hz

    sample_rate : int or float
        sample rate of the supplied signal

    order : int
        filter order, defines the strength of the roll-off
        around the cutoff frequency. Typically orders above 6
        are not used frequently.
        default : 2
    
    Returns
    -------
    out : tuple
        numerator and denominator (b, a) polynomials
        of the defined Butterworth IIR filter.

    Examples
    --------
    we can specify lowcut, highcut and sample_rate as ints or floats.

    >>> b, a = butter_bandpass(lowcut = 1, highcut = 6, sample_rate = 100, order = 2)
    >>> b, a = butter_bandpass(lowcut = 0.4, highcut = 3.7, sample_rate = 72.6, order = 2)
    '''
    nyq = 0.5 * sample_rate
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def filter_signal(data, cutoff, sample_rate, order=2, filtertype='lowpass',
                  return_top = False):
    '''Apply the specified filter

    Function that applies the specified lowpass, highpass or bandpass filter to
    the provided dataset.

    Parameters
    ----------
    data : 1-dimensional numpy array or list 
        Sequence containing the to be filtered data

    cutoff : int, float or tuple
        the cutoff frequency of the filter. Expects float for low and high types
        and for bandpass filter expects list or array of format [lower_bound, higher_bound]

    sample_rate : int or float
        the sample rate with which the passed data sequence was sampled

    order : int
        the filter order 
        default : 2

    filtertype : str
        The type of filter to use. Available:
        - lowpass : a lowpass butterworth filter
        - highpass : a highpass butterworth filter
        - bandpass : a bandpass butterworth filter
        - notch : a notch filter around specified frequency range
        both the highpass and notch filter are useful for removing baseline wander. The notch
        filter is especially useful for removing baseling wander in ECG signals.


    Returns
    -------
    out : 1d array
        1d array containing the filtered data

    Examples
    --------
    >>> import numpy as np
    >>> import heartpy as hp

    Using standard data provided

    >>> data, _ = hp.load_exampledata(0)

    We can filter the signal, for example with a lowpass cutting out all frequencies
    of 5Hz and greater (with a sloping frequency cutoff)

    >>> filtered = filter_signal(data, cutoff = 5, sample_rate = 100.0, order = 3, filtertype='lowpass')
    >>> print(np.around(filtered[0:6], 3))
    [530.175 517.893 505.768 494.002 482.789 472.315]

    Or we can cut out all frequencies below 0.75Hz with a highpass filter:

    >>> filtered = filter_signal(data, cutoff = 0.75, sample_rate = 100.0, order = 3, filtertype='highpass')
    >>> print(np.around(filtered[0:6], 3))
    [-17.975 -28.271 -38.609 -48.992 -58.422 -67.902]

    Or specify a range (here: 0.75 - 3.5Hz), outside of which all frequencies
    are cut out.

    >>> filtered = filter_signal(data, cutoff = [0.75, 3.5], sample_rate = 100.0, 
    ... order = 3, filtertype='bandpass')
    >>> print(np.around(filtered[0:6], 3))
    [-12.012 -23.159 -34.261 -45.12  -55.541 -65.336]

    A 'Notch' filtertype is also available (see remove_baseline_wander).
    
    >>> filtered = filter_signal(data, cutoff = 0.05, sample_rate = 100.0, filtertype='notch')

    Finally we can use the return_top flag to only return the filter response that
    has amplitute above zero. We're only interested in the peaks, and sometimes
    this can improve peak prediction:

    >>> filtered = filter_signal(data, cutoff = [0.75, 3.5], sample_rate = 100.0, 
    ... order = 3, filtertype='bandpass', return_top = True)
    >>> print(np.around(filtered[48:53], 3))
    [ 0.     0.     0.409 17.088 35.673]
    '''
    if filtertype.lower() == 'lowpass':
        b, a = butter_lowpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'highpass':
        b, a = butter_highpass(cutoff, sample_rate, order=order)
    elif filtertype.lower() == 'bandpass':
        assert type(cutoff) == tuple or list or np.array, 'if bandpass filter is specified, \
cutoff needs to be array or tuple specifying lower and upper bound: [lower, upper].'
        b, a = butter_bandpass(cutoff[0], cutoff[1], sample_rate, order=order)
    elif filtertype.lower() == 'notch':
        b, a = iirnotch(cutoff, Q = 0.005, fs = sample_rate)
    else:
        raise ValueError('filtertype: %s is unknown, available are: \
lowpass, highpass, bandpass, and notch' %filtertype)

    filtered_data = filtfilt(b, a, data)
    
    if return_top:
        return np.clip(filtered_data, a_min = 0, a_max = None)
    else:
        return filtered_data


def remove_baseline_wander(data, sample_rate, cutoff=0.05):
    '''removes baseline wander

    Function that uses a Notch filter to remove baseline
    wander from (especially) ECG signals

    Parameters
    ----------
    data : 1-dimensional numpy array or list 
        Sequence containing the to be filtered data

    sample_rate : int or float
        the sample rate with which the passed data sequence was sampled

    cutoff : int, float 
        the cutoff frequency of the Notch filter. We recommend 0.05Hz.
        default : 0.05

    Returns
    -------
    out : 1d array
        1d array containing the filtered data

    Examples
    --------
    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(0)

    baseline wander is removed by calling the function and specifying
    the data and sample rate.

    >>> filtered = remove_baseline_wander(data, 100.0)
    '''

    return filter_signal(data = data, cutoff = cutoff, sample_rate = sample_rate,
                         filtertype='notch')


def hampel_filter(data, filtsize=6):
    '''Detect outliers based on hampel filter
    
    Funcion that detects outliers based on a hampel filter.
    The filter takes datapoint and six surrounding samples.
    Detect outliers based on being more than 3std from window mean.
    See:
    https://www.mathworks.com/help/signal/ref/hampel.html
    
    Parameters
    ----------
    data : 1d list or array
        list or array containing the data to be filtered

    filtsize : int
        the filter size expressed the number of datapoints
        taken surrounding the analysed datapoint. a filtsize
        of 6 means three datapoints on each side are taken.
        total filtersize is thus filtsize + 1 (datapoint evaluated)

    Returns
    -------
    out :  array containing filtered data

    Examples
    --------
    >>> from .datautils import get_data, load_exampledata
    >>> data, _ = load_exampledata(0)
    >>> filtered = hampel_filter(data, filtsize = 6)
    >>> print('%i, %i' %(data[1232], filtered[1232]))
    497, 496
    '''

    #generate second list to prevent overwriting first
    #cast as array to be sure, in case list is passed
    output = np.copy(np.asarray(data)) 
    onesided_filt = filtsize // 2
    for i in range(onesided_filt, len(data) - onesided_filt - 1):
        dataslice = output[i - onesided_filt : i + onesided_filt]
        mad = MAD(dataslice)
        median = np.median(dataslice)
        if output[i] > median + (3 * mad):
            output[i] = median
    return output


def hampel_correcter(data, sample_rate):
    '''apply altered version of hampel filter to suppress noise.

    Function that returns te difference between data and 1-second 
    windowed hampel median filter. Results in strong noise suppression 
    characteristics, but relatively expensive to compute.

    Result on output measures is present but generally not large. However,
    use sparingly, and only when other means have been exhausted.

    Parameters
    ----------
    data : 1d numpy array
        array containing the data to be filtered

    sample_rate : int or float
        sample rate with which data was recorded
       
    Returns
    -------
    out : 1d numpy array
        array containing filtered data

    Examples
    --------
    >>> from .datautils import get_data, load_exampledata
    >>> data, _ = load_exampledata(1)
    >>> filtered = hampel_correcter(data, sample_rate = 116.995)

    '''

    return data - hampel_filter(data, filtsize=int(sample_rate))


def quotient_filter(RR_list, RR_list_mask = [], iterations=2):
    '''applies a quotient filter

    Function that applies a quotient filter as described in
    "Piskorki, J., Guzik, P. (2005), Filtering Poincare plots"

    Parameters
    ----------
    RR_list - 1d array or list
        array or list of peak-peak intervals to be filtered

    RR_list_mask - 1d array or list
        array or list containing the mask for which intervals are 
        rejected. If not supplied, it will be generated. Mask is 
        zero for accepted intervals, one for rejected intervals.

    iterations - int
        how many times to apply the quotient filter. Multipled
        iterations have a stronger filtering effect
        default : 2

    Returns
    -------
    RR_list_mask : 1d array
        mask for RR_list, 1 where intervals are rejected, 0 where
        intervals are accepted.

    Examples
    --------
    Given some example data let's generate an RR-list first
    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(1)
    >>> sample_rate = hp.get_samplerate_mstimer(timer)
    >>> wd, m = hp.process(data, sample_rate)
    >>> rr = wd['RR_list']
    >>> rr_mask = wd['RR_masklist']

    Given this data we can use this function to further clean the data:
    >>> new_mask = quotient_filter(rr, rr_mask)

    Although specifying the mask is optional, as you may not always have a
    pre-computed mask available:
    >>> new_mask = quotient_filter(rr)
    
    '''

    if len(RR_list_mask) == 0:
        RR_list_mask = np.zeros((len(RR_list)))
    else:
        assert len(RR_list) == len(RR_list_mask), \
        'error: RR_list and RR_list_mask should be same length if RR_list_mask is specified'

    for iteration in range(iterations):
        for i in range(len(RR_list) - 1):
            if RR_list_mask[i] + RR_list_mask[i + 1] != 0:
                pass #skip if one of both intervals is already rejected
            elif 0.8 <= RR_list[i] / RR_list[i + 1] <= 1.2:
                pass #if R-R pair seems ok, do noting
            else: #update mask
                RR_list_mask[i] = 1
                #RR_list_mask[i + 1] = 1

    return np.asarray(RR_list_mask)


def smooth_signal(data, sample_rate, window_length=None, polyorder=3):
    '''smooths given signal using savitzky-golay filter

    Function that smooths data using savitzky-golay filter using default settings.

    Functionality requested by Eirik Svendsen. Added since 1.2.4

    Parameters
    ----------
    data : 1d array or list
        array or list containing the data to be filtered

    sample_rate : int or float
        the sample rate with which data is sampled

    window_length : int or None
        window length parameter for savitzky-golay filter, see Scipy.signal.savgol_filter docs.
        Must be odd, if an even int is given, one will be added to make it uneven.
        default : 0.1  * sample_rate

    polyorder : int
        the order of the polynomial fitted to the signal. See scipy.signal.savgol_filter docs.
        default : 3

    Returns
    -------
    smoothed : 1d array
        array containing the smoothed data

    Examples
    --------
    Given a fictional signal, a smoothed signal can be obtained by smooth_signal():

    >>> x = [1, 3, 4, 5, 6, 7, 5, 3, 1, 1]
    >>> smoothed = smooth_signal(x, sample_rate = 2, window_length=4, polyorder=2)
    >>> np.around(smoothed[0:4], 3)
    array([1.114, 2.743, 4.086, 5.   ])

    If you don't specify the window_length, it is computed to be 10% of the 
    sample rate (+1 if needed to make odd)
    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(0)
    >>> smoothed = smooth_signal(data, sample_rate = 100)

    '''

    if window_length == None:
        window_length = sample_rate // 10
        
    if window_length % 2 == 0 or window_length == 0: window_length += 1

    smoothed = savgol_filter(data, window_length = window_length,
                             polyorder = polyorder)

    return smoothed