'''
Functions for loading and slicing data
'''

from datetime import datetime
from pkg_resources import resource_filename

import numpy as np
from scipy.io import loadmat

__all__ = ['get_data',
           'get_samplerate_mstimer',
           'get_samplerate_datetime',
           '_sliding_window',
           'rolling_mean',
           'outliers_iqr_method',
           'outliers_modified_z',
           'MAD',
           'load_exampledata']


def get_data(filename, delim=',', column_name='None', encoding=None, 
             ignore_extension=False):
    '''load data from file

    Function to load data from a .CSV or .MAT file into numpy array.
    File can be accessed from local disk or url.

    Parameters
    ----------
    filename : string
        absolute or relative path to the file object to read

    delim : string
        the delimiter used if CSV file passed
        default : ','

    column_name : string
        for CSV files with header: specify column that contains the data
        for matlab files it specifies the table name that contains the data
        default : 'None'

    ignore_extension : bool
        if True, extension is not tested, use for example for files where
        the extention is not .csv or .txt but the data is formatted as if
        it is.
        default : False

    Returns
    -------
    out : 1-d numpy array
        array containing the data from the requested column of the specified file

    Examples
    --------
    As an example, let's load two example data files included in the package
    For this we use pkg_resources for automated testing purposes, you don't need
    this when using the function.

    >>> from pkg_resources import resource_filename
    >>> filepath = resource_filename(__name__, 'data/data.csv')

    So, assuming your file lives at 'filepath', you open it as such:

    >>> get_data(filepath)
    array([530., 518., 506., ..., 492., 493., 494.])

    Files with multiple columns can be opened by specifying the 'column_name' where
    the data resides:

    >>> filepath = resource_filename(__name__, 'data/data2.csv')

    Again you don't need the above. It is there for automated testing.

    >>> get_data(filepath, column_name='timer')
    array([0.00000000e+00, 8.54790319e+00, 1.70958064e+01, ...,
           1.28192904e+05, 1.28201452e+05, 1.28210000e+05])

    You can open matlab files in much the same way by specifying the column
    where the data lives:

    >>> filepath = resource_filename(__name__, 'data/data2.mat')

    Again you don't need the above. It is there for automated testing.
    Open matlab file by specifying the column name as well:

    >>> get_data(filepath, column_name='hr')
    array([515., 514., 514., ..., 492., 494., 496.])

    You can any csv formatted text file no matter the extension if you
    set ignore_extension to True:

    >>> filepath = resource_filename(__name__, 'data/data.log')
    >>> get_data(filepath, ignore_extension = True)
    array([530., 518., 506., ..., 492., 493., 494.])

    You can specify column names in the same way when using ignore_extension

    >>> filepath = resource_filename(__name__, 'data/data2.log')
    >>> data = get_data(filepath, column_name = 'hr', ignore_extension = True)
    '''
    file_ext = filename.split('.')[-1]
    if file_ext == 'csv' or file_ext == 'txt':
        if column_name != 'None':
            hrdata = np.genfromtxt(filename, delimiter=delim, names=True, dtype=None, encoding=None)
            try:
                hrdata = hrdata[column_name]
            except Exception as error:
                raise LookupError('\nError loading column "%s" from file "%s". \
Is column name specified correctly?\n The following error was provided: %s' 
                                 %(column_name, filename, error))
        elif column_name == 'None':
            hrdata = np.genfromtxt(filename, delimiter=delim, dtype=np.float64)
        else: # pragma: no cover
            raise LookupError('\nError: column name "%s" not found in header of "%s".\n'
                              %(column_name, filename))
    elif file_ext == 'mat':
        data = loadmat(filename)
        if column_name != "None":
            hrdata = np.array(data[column_name][:, 0], dtype=np.float64)
        else: # pragma: no cover
            raise LookupError('\nError: column name required for Matlab .mat files\n\n')
    else:
        if ignore_extension:
            if column_name != 'None':
                hrdata = np.genfromtxt(filename, delimiter=delim, names=True, dtype=None, encoding=None)
                try:
                    hrdata = hrdata[column_name]
                except Exception as error:
                    raise LookupError('\nError loading column "%s" from file "%s". \
Is column name specified correctly?\n' 
                                      %(column_name, filename))
            elif column_name == 'None': # pragma: no cover
                hrdata = np.genfromtxt(filename, delimiter=delim, dtype=np.float64)
            else: # pragma: no cover
                raise LookupError('\nError: column name "%s" not found in header of "%s".\n'
                                  %(column_name, filename))
        else:
            raise IncorrectFileType('unknown file format')
            return None 
    return hrdata


def get_samplerate_mstimer(timerdata):
    '''detemine sample rate based on ms timer

    Function to determine sample rate of data from ms-based timer list or array.

    Parameters
    ----------
    timerdata : 1d numpy array or list
        sequence containing values of a timer, in ms

    Returns
    -------
    out : float
        the sample rate as determined from the timer sequence provided
        
    Examples
    --------
    first we load a provided example dataset

    >>> data, timer = load_exampledata(example = 1)
    
    since it's a timer that counts miliseconds, we use this function.
    Let's also round to three decimals

    >>> round(get_samplerate_mstimer(timer), 3)
    116.996

    of course if another time unit is used, converting it to ms-based
    should be trivial.
    '''
    sample_rate = ((len(timerdata) / (timerdata[-1]-timerdata[0]))*1000)
    return sample_rate


def get_samplerate_datetime(datetimedata, timeformat='%H:%M:%S.%f'):
    '''determine sample rate based on datetime

    Function to determine sample rate of data from datetime-based timer
    list or array.

    Parameters
    ----------
    timerdata : 1-d numpy array or list
        sequence containing datetime strings

    timeformat : string
        the format of the datetime-strings in datetimedata
        default : '%H:%M:%S.f' (24-hour based time including ms: e.g. 21:43:12.569)

    Returns
    -------
    out : float
        the sample rate as determined from the timer sequence provided

    Examples
    --------
    We load the data like before

    >>> data, timer = load_exampledata(example = 2)
    >>> timer[0]
    '2016-11-24 13:58:58.081000'

    Note that we need to specify the timeformat used so that datetime understands
    what it's working with:

    >>> round(get_samplerate_datetime(timer, timeformat = '%Y-%m-%d %H:%M:%S.%f'), 3)
    100.42
    '''
    datetimedata = np.asarray(datetimedata, dtype='str') #cast as str in case of np.bytes type
    elapsed = ((datetime.strptime(datetimedata[-1], timeformat) -
                datetime.strptime(datetimedata[0], timeformat)).total_seconds())
    sample_rate = (len(datetimedata) / elapsed)
    return sample_rate


def _sliding_window(data, windowsize):
    '''segments data into windows

    Function to segment data into windows for rolling mean function.
    Function returns the data segemented into sections.

    Parameters
    ----------
    data : 1d array or list
        array or list containing data over which sliding windows are computed

    windowsize : int
        size of the windows to be created by the function

    Returns
    -------
    out : array of arrays
        data segmented into separate windows.

    Examples
    --------
    >>> import numpy as np
    >>> data = np.array([1, 2, 3, 4, 5])
    >>> windows = _sliding_window(data, windowsize = 3)
    >>> windows.shape
    (3, 3)
    '''
    shape = data.shape[:-1] + (data.shape[-1] - windowsize + 1, windowsize)
    strides = data.strides + (data.strides[-1],)
    return np.lib.stride_tricks.as_strided(data, shape=shape, strides=strides)


def rolling_mean(data, windowsize, sample_rate):
    '''calculates rolling mean

    Function to calculate the rolling mean (also: moving average) over the passed data.

    Parameters
    ----------
    data : 1-dimensional numpy array or list
        sequence containing data over which rolling mean is to be computed

    windowsize : int or float 
        the window size to use, in seconds 
        calculated as windowsize * sample_rate

    sample_rate : int or float
        the sample rate of the data set

    Returns
    -------
    out : 1-d numpy array
        sequence containing computed rolling mean

    Examples
    --------
    >>> data, _ = load_exampledata(example = 1)
    >>> rmean = rolling_mean(data, windowsize=0.75, sample_rate=100)
    >>> rmean[100:110]
    array([514.49333333, 514.49333333, 514.49333333, 514.46666667,
           514.45333333, 514.45333333, 514.45333333, 514.45333333,
           514.48      , 514.52      ])
    '''
    avg_hr = (np.mean(data))
    data_arr = np.array(data)
    rol_mean = np.mean(_sliding_window(data_arr, int(windowsize*sample_rate)), axis=1)
    missing_vals = np.array([avg_hr for i in range(0, int(abs(len(data_arr) - len(rol_mean))/2))])
    rol_mean = np.insert(rol_mean, 0, missing_vals)
    rol_mean = np.append(rol_mean, missing_vals)

    #only to catch length errors that sometimes unexplicably occur. 
    ##Generally not executed, excluded from testing and coverage
    if len(rol_mean) != len(data): # pragma: no cover
        lendiff = len(rol_mean) - len(data)
        if lendiff < 0:
            rol_mean = np.append(rol_mean, 0)
        else:
            rol_mean = rol_mean[:-1]            
    return rol_mean


def outliers_iqr_method(hrvalues):
    '''removes outliers

    Function that removes outliers based on the interquartile range method and
    substitutes them for the median
    see: https://en.wikipedia.org/wiki/Interquartile_range

    Parameters
    ----------
    hrvalues : 1-d numpy array or list 
        sequence of values, from which outliers need to be identified

    Returns
    -------
    out : tuple
        [0] cleaned sequence with identified outliers substituted for the median
        [1] list of indices that have been replaced in the original array or list

    Examples
    --------
    >>> x = [2, 4, 3, 4, 6, 7, 35, 2, 3, 4]
    >>> outliers_iqr_method(x)
    ([2, 4, 3, 4, 6, 7, 4.0, 2, 3, 4], [6])
    '''
    med = np.median(hrvalues)
    q1, q3 = np.percentile(hrvalues, [25, 75])
    iqr = q3 - q1
    lower = q1 - (1.5 * iqr)
    upper = q3 + (1.5 * iqr)
    output = []
    replaced_indices = []
    for i in range(0,len(hrvalues)):
        if hrvalues[i] < lower or hrvalues[i] > upper:
            output.append(med)
            replaced_indices.append(i)
        else:
            output.append(hrvalues[i])
    return output, replaced_indices


def outliers_modified_z(hrvalues):
    '''removes outliers

    Function that removes outliers based on the modified Z-score metric and
    substitutes them for the median

    Parameters
    ----------
    hrvalues : 1-d numpy array or list 
        sequence of values, from which outliers need to be identified

    Returns
    -------
    out : tuple
        [0] cleaned sequence with identified outliers substituted for the median
        [1] list of indices that have been replaced in the original array or list

    Examples
    --------
    >>> x = [2, 4, 3, 4, 6, 7, 35, 2, 3, 4]
    >>> outliers_modified_z(x)
    ([2, 4, 3, 4, 6, 7, 4.0, 2, 3, 4], [6])
    '''
    hrvalues = np.array(hrvalues)
    threshold = 3.5
    med = np.median(hrvalues)
    mean_abs_dev = MAD(hrvalues)
    modified_z_result = 0.6745 * (hrvalues - med) / mean_abs_dev
    output = []
    replaced_indices = []
    for i in range(0, len(hrvalues)):
        if np.abs(modified_z_result[i]) <= threshold:
            output.append(hrvalues[i])
        else:
            output.append(med)
            replaced_indices.append(i)
    return output, replaced_indices


def MAD(data):
    '''computes median absolute deviation

    Function that compute median absolute deviation of data slice
    See: https://en.wikipedia.org/wiki/Median_absolute_deviation
    
    Parameters
    ----------
    data : 1-dimensional numpy array or list
        sequence containing data over which to compute the MAD

    Returns
    -------
    out : float
        the Median Absolute Deviation as computed

    Examples
    --------
    >>> x = [2, 4, 3, 4, 6, 7, 35, 2, 3, 4]
    >>> MAD(x)
    1.5
    '''
    med = np.median(data)
    return np.median(np.abs(data - med))


def load_exampledata(example=0):
    '''loads example data

    Function to load one of the example datasets included in HeartPy
    and used in the documentation.

    Parameters
    ----------
    example : int (0, 1, 2)
        selects example data used in docs of three datafiles.
        Available (see github repo for source of files):
        0 : data.csv
        1 : data2.csv
        2 : data3.csv
        default : 0

    Returns
    -------
    out : tuple of two arrays
        Contains the data and timer column. If no timer data is
        available, such as in example 0, an empty second
        array is returned.

    Examples
    --------
    This function can load one of the three example data files provided
    with HeartPy. It returns both the data and a timer if that is present

    For example:

    >>> data, _ = load_exampledata(0)
    >>> data[0:5]
    array([530., 518., 506., 494., 483.])

    And another example:

    >>> data, timer = load_exampledata(1)
    >>> [round(x, 2) for x in timer[0:5]]
    [0.0, 8.55, 17.1, 25.64, 34.19]
    '''

    timer = []
    
    if example == 0:
        path = path = 'data/data.csv'
        filepath = resource_filename(__name__, path)
        data = get_data(filepath)
    elif example == 1:
        path = path = 'data/data2.csv'
        filepath = resource_filename(__name__, path)
        data = get_data(filepath, column_name = 'hr')
        timer = get_data(filepath, column_name = 'timer')
    elif example == 2:
        path = path = 'data/data3.csv'
        filepath = resource_filename(__name__, path)
        data = get_data(filepath, column_name = 'hr')
        timer = get_data(filepath, column_name = 'datetime')
    else:
        raise ValueError('Incorrect data file specified.\
available datafiles are data.csv (0), data2.csv(1), data3.csv(2).')

    return data, timer
