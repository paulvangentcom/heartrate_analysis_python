'''
Functions for data filtering tasks.

Frequency filtering
-------------------
- 'filter_singnal' -- filter specified frequency range from signal
- 'hampel_filter' -- a type of median filter used to remove outliers
- 'hampel_correcter' -- filter based on variant of hampel_filter that 
                        has good noise suppression characteristics.

Hidden helper functions
-----------------------
- 'butter_lowpass' -- returns polynomials for lowpass IIR butterworth filter
- 'butter_highpass' -- returns polynomials for highpass IIR butterworth filter
- 'butter_bandpass' -- returns polynomials for bandpass IIR butterworth filter
'''

from scipy.signal import butter, filtfilt
import numpy as np

from .datautils import MAD

__all__ = ['filter_signal',
           'hampel_filter',
           'hampel_correcter']

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
    # we can specify the cutoff and sample_rate as ints or floats.
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
    # we can specify lowcut, highcut and sample_rate as ints or floats.
    >>> b, a = butter_bandpass(lowcut = 1, highcut = 6, sample_rate = 100, order = 2)
    >>> b, a = butter_bandpass(lowcut = 0.4, highcut = 3.7, sample_rate = 72.6, order = 2)
    '''
    nyq = 0.5 * sample_rate
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def filter_signal(data, cutoff, sample_rate, order, filtertype='lowpass',
                  return_top = False):
    '''Apply the specified filed

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

    Returns
    -------
    out : tuple
        numerator and denominator (b, a) polynomials
        of the defined Butterworth IIR filter.

    Examples
    --------
    >>> import numpy as np
    >>> import heartpy as hp

    Assuming standard data provided
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

    Finally we can use the return_top flag to only return the filter response that
    has amplitute above zero. We're only interested in the peaks, and sometimes
    this can improve peak prediction:
    >>> filtered = filter_signal(data, cutoff = [0.75, 3.5], sample_rate = 100.0, 
    ... order = 3, filtertype='bandpass', return_top = True)
    >>> filtered[48:53]
    array([ 0.        ,  0.        ,  0.40944503, 17.08776349, 35.67256091])
    '''
    if filtertype == 'lowpass':
        b, a = butter_lowpass(cutoff, sample_rate, order=order)
    elif filtertype == 'highpass':
        b, a = butter_highpass(cutoff, sample_rate, order=order)
    elif filtertype == 'bandpass':
        assert type(cutoff) == tuple or list or np.array, 'if bandpass filter is specified, \
cutoff needs to be array or tuple specifying lower and upper bound: [lower, upper].'
        b, a = butter_bandpass(cutoff[0], cutoff[1], sample_rate, order=order)
    filtered_data = filtfilt(b, a, data)
    if return_top:
        return np.clip(filtered_data, a_min = 0, a_max = None)
    else:
        return filtered_data


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