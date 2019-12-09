'''
main module for HeartPy.
'''

from datetime import datetime
import time
import os

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.signal import butter, filtfilt, welch, periodogram, resample_poly, resample

from . import exceptions
from .datautils import get_data, get_samplerate_mstimer, get_samplerate_datetime,\
                       rolling_mean, outliers_iqr_method, outliers_modified_z, \
                       load_exampledata
from .preprocessing import scale_data, scale_sections, interpolate_clipping, \
                           flip_signal, enhance_peaks, enhance_ecg_peaks
from .filtering import filter_signal, hampel_filter, hampel_correcter, \
                       remove_baseline_wander, smooth_signal
from .peakdetection import make_windows, append_dict, fit_peaks, check_peaks, \
                           check_binary_quality, interpolate_peaks
from .visualizeutils import plotter, segment_plotter, plot_poincare, plot_breathing
from .analysis import calc_rr, calc_rr_segment, clean_rr_intervals, calc_ts_measures, \
                      calc_fd_measures, calc_breathing, calc_poincare

from . import config
config.init() #initialize global conf vars

__all__ = ['enhance_peaks',
           'enhance_ecg_peaks',
           'get_data',
           'get_samplerate_mstimer',
           'get_samplerate_datetime',
           'hampel_correcter',
           'hampel_filter',
           'load_exampledata',
           'plotter',
           'plot_breathing',
           'plot_poincare',
           'process',
           'process_rr',
           'process_segmentwise',
           'flip_signal',
           'remove_baseline_wander',
           'scale_data',
           'scale_sections',
           'segment_plotter',
           'smooth_signal',
           'filter_signal',
           'run_tests']


def process(hrdata, sample_rate, windowsize=0.75, report_time=False, 
            calc_freq=False, freq_method='welch', freq_square=True,
            interp_clipping=False, clipping_scale=False, interp_threshold=1020, 
            hampel_correct=False, bpmmin=40, bpmmax=180, reject_segmentwise=False, 
            high_precision=False, high_precision_fs=1000.0, breathing_method='fft',
            clean_rr=False, clean_rr_method='quotient-filter', measures={}, working_data={}):
    '''processes passed heart rate data.
    
    Processes the passed heart rate data. Returns measures{} dict containing results.

    Parameters
    ----------
    hrdata : 1d array or list 
        array or list containing heart rate data to be analysed

    sample_rate : int or float
        the sample rate with which the heart rate data is sampled

    windowsize : int or float
        the window size in seconds to use in the calculation of the moving average.
        Calculated as windowsize * sample_rate
        default : 0.75

    report_time : bool
        whether to report total processing time of algorithm 
        default : True

    calc_freq : bool
        whether to compute time-series measurements 
        default : False

    freq_method : str
        method used to extract the frequency spectrum. Available: 'fft' (Fourier Analysis), 
        'periodogram', and 'welch' (Welch's method). 
        default : 'welch'

    freq_square : bool
        whether to square the power spectrum returned when computing frequency measures
        default : true

    interp_clipping : bool 
        whether to detect and interpolate clipping segments of the signal 
        default : True

    clipping_scale : bool
        whether to scale the data prior to clipping detection. Can correct errors 
        if signal amplitude has been affected after digitization (for example through 
        filtering). Not recommended by default. 
        default : False

    interp_threshold : int or float
        threshold to use to detect clipping segments. Recommended to be a few
        datapoints below the sensor or ADC's maximum value (to account for
        slight data line noise). 
        default : 1020, 4 below max of 1024 for 10-bit ADC

    hampel_correct : bool 
        whether to reduce noisy segments using large median filter. Disabled by
        default due to computational complexity and (small) distortions induced
        into output measures. Generally it is not necessary.
        default : False

    bpmmin : int or float 
        minimum value to see as likely for BPM when fitting peaks
        default : 40

    bpmmax : int or float
        maximum value to see as likely for BPM when fitting peaks
        default : 180

    reject_segmentwise : bool
        whether to reject segments with more than 30% rejected beats. 
        By default looks at segments of 10 beats at a time. 
        default : False

    high_precision : bool 
        whether to estimate peak positions by upsampling signal to sample rate
        as specified in high_precision_fs
        default : false

    high_precision_fs : int or float 
        the sample rate to which to upsample for more accurate peak position estimation 
        default : 1000 Hz

    breathing_method : str
        method to use for estimating breathing rate, should be 'welch' or 'fft'
        default : fft

    clean_rr : bool
        if true, the RR_list is further cleaned with an outlier rejection pass
        default : false

    clean_rr_method: str
        how to find and reject outliers. Available methods are ' quotient-filter', 
        'iqr' (interquartile range), and 'z-score'.
        default : 'quotient-filter'

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function.

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    working_data : dict
        dictionary object used to store temporary values.
    
    measures : dict
        dictionary object used by heartpy to store computed measures.

    Examples
    --------
    There's example data included in HeartPy to help you get up to speed. Here are
    provided two examples of how to approach heart rate analysis.

    The first example contains noisy sections and comes with a timer column that
    counts miliseconds since start of recording. 

    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(1)
    >>> sample_rate = hp.get_samplerate_mstimer(timer)
    >>> '%.3f' %sample_rate
    '116.996'

    The sample rate is one of the most important characteristics during the
    heart rate analysis, as all measures are relative to this.
    
    With all data loaded and the sample rate determined, nalysis is now easy:

    >>> wd, m = hp.process(data, sample_rate = sample_rate)

    The measures ('m') dictionary returned contains all determined measures

    >>> '%.3f' %m['bpm']
    '62.376'
    >>> '%.3f' %m['rmssd']
    '57.070'

    Using a slightly longer example:

    >>> data, timer = hp.load_exampledata(2)
    >>> print(timer[0])
    2016-11-24 13:58:58.081000

    As you can see something is going on here: we have a datetime-based timer.
    HeartPy can accomodate this and determine sample rate nontheless:

    >>> sample_rate = hp.get_samplerate_datetime(timer, timeformat = '%Y-%m-%d %H:%M:%S.%f')
    >>> '%.3f' %sample_rate
    '100.420'

    Now analysis can proceed. Let's also compute frequency domain data and interpolate clipping.
    In this segment the clipping is visible around amplitude 980 so let's set that as well:

    >>> wd, m = hp.process(data, sample_rate = sample_rate, calc_freq = True, 
    ... interp_clipping = True, interp_threshold = 975)
    >>> '%.3f' %m['bpm']
    '97.270'
    >>> '%.3f' %m['rmssd']
    '34.743'
    >>> '%.3f' %m['lf/hf']
    '4.960'

    High precision mode will upsample 200ms of data surrounding detected peak
    and attempt to estimate the peak's real position with higher accuracy.
    Use high_precision_fs to set the virtual sample rate to which the peak
    will be upsampled (e.g. 1000Hz gives an estimated 1ms accuracy)

    >>> wd, m = hp.process(data, sample_rate = sample_rate, calc_freq = True, 
    ... high_precision = True, high_precision_fs = 1000.0)

    Finally setting reject_segmentwise will reject segments with more than 30% rejected beats
    See check_binary_quality in the peakdetection.py module.

    >>> wd, m = hp.process(data, sample_rate = sample_rate, calc_freq = True, 
    ... reject_segmentwise = True)

    Final test for code coverage, let's turn all bells and whistles on that haven't been
    tested yet

    >>> wd, m = hp.process(data, sample_rate = 100.0, calc_freq = True, 
    ... interp_clipping = True, clipping_scale = True, hampel_correct = True,
    ... reject_segmentwise = True, clean_rr = True)
    '''
    t1 = time.perf_counter()

    assert np.asarray(hrdata).ndim == 1, 'error: multi-dimensional data passed to process() \
function. Please supply a 1d array or list containing heart rate signal data. \n\nDid you perhaps \
include an index column?'

    if interp_clipping:
        if clipping_scale:
            hrdata = scale_data(hrdata)
        hrdata = interpolate_clipping(hrdata, sample_rate, threshold=interp_threshold)

    if hampel_correct:
        hrdata = enhance_peaks(hrdata)
        hrdata = hampel_correcter(hrdata, sample_rate)

    working_data['hr'] = hrdata
    rol_mean = rolling_mean(hrdata, windowsize, sample_rate)

    working_data = fit_peaks(hrdata, rol_mean, sample_rate, bpmmin=bpmmin,
                             bpmmax=bpmmax, working_data=working_data)
    
    if high_precision:
        working_data = interpolate_peaks(hrdata, working_data['peaklist'], sample_rate=sample_rate, 
                                         desired_sample_rate=high_precision_fs, working_data=working_data)

    working_data = calc_rr(working_data['peaklist'], sample_rate, working_data=working_data)
    working_data = check_peaks(working_data['RR_list'], working_data['peaklist'], working_data['ybeat'],
                               reject_segmentwise, working_data=working_data)

    if clean_rr:
        working_data = clean_rr_intervals(working_data, method = clean_rr_method)

    working_data, measures = calc_ts_measures(working_data['RR_list_cor'], working_data['RR_diff'],
                                              working_data['RR_sqdiff'], measures=measures, 
                                              working_data=working_data)
    
    measures = calc_poincare(working_data['RR_list'], working_data['RR_masklist'], measures = measures,
                             working_data = working_data)

    try:
        measures, working_data = calc_breathing(working_data['RR_list_cor'], hrdata, sample_rate, 
                                                method = breathing_method, measures=measures, 
                                                working_data=working_data)
    except:
        measures['breathingrate'] = np.nan

    if calc_freq:
        working_data, measures = calc_fd_measures(method=freq_method, measures=measures,
                                                  working_data = working_data)
    
    #report time if requested. Exclude from tests, output is untestable.
    if report_time: # pragma: no cover
        print('\nFinished in %.8s sec' %(time.perf_counter()-t1))

    return working_data, measures


def process_segmentwise(hrdata, sample_rate, segment_width=120, segment_overlap=0,
                        segment_min_size=20, replace_outliers=False, outlier_method='iqr',
                        mode='full', **kwargs):
    '''processes passed heart rate data with a windowed function

    Analyses a long heart rate data array by running a moving window 
    over the data, computing measures in each iteration. Both the window width
    and the overlap with the previous window location are settable.

    Parameters
    ----------
    hrdata : 1d array or list 
        array or list containing heart rate data to be analysed

    sample_rate : int or float
        the sample rate with which the heart rate data is sampled

    segment_width : int or float
        width of segments in seconds
        default : 120

    segment_overlap: float
        overlap fraction of adjacent segments.
        Needs to be 0 <= segment_overlap < 1.
        default : 0 (no overlap)

    segment_min_size : int
        often a tail end of the data remains after segmenting into segments.
        segment_min_size indicates the minimum length (in seconds) the tail 
        end needs  to be in order to be included in analysis. It is discarded 
        if it's shorter.
        default : 20

    replace_outliers : bool
        whether to detct and replace outliers in the segments. Will iterate over
        all computed measures and evaluate each.

    outlier_method : str
        what method to use to detect outlers. Available are 'iqr', which uses the
        inter-quartile range, and 'z-score', which uses the modified z-score approach.

    mode : str
        'full' or 'fast'

    Keyword arguments:
    ------------------
    hrdata -- 1-dimensional numpy array or list containing heart rate data
    sample_rate -- the sample rate of the heart rate data
    segment_width -- the width of the segment, in seconds, within which all measures 
                     will be computed.
    segment_overlap -- the fraction of overlap of adjacent segments, 
                       needs to be 0 <= segment_overlap < 1
    segment_min_size -- After segmenting the data, a tail end will likely remain that is shorter than the specified
                        segment_size. segment_min_size sets the minimum size for the last segment of the 
                        generated series of segments to still be included. Default = 20.
    replace_outliers -- bool, whether to replace outliers (likely caused by peak fitting
                        errors on one or more segments) with the median.
    outlier_method -- which  method to use to detect outliers. Available are the
                      'interquartile-range' ('iqr') and the 'modified z-score' ('z-score') methods.

    Returns
    -------
    working_data : dict
        dictionary object used to store temporary values.
    
    measures : dict
        dictionary object used by heartpy to store computed measures.
        
    Examples
    --------
    Given one of the included example datasets we can demonstrate this function:

    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(2)
    >>> sample_rate = hp.get_samplerate_datetime(timer, timeformat = '%Y-%m-%d %H:%M:%S.%f')
    >>> wd, m = hp.process_segmentwise(data, sample_rate, segment_width=120, segment_overlap=0.5)
    >>> len(m['bpm'])
    11

    The function has split the data into 11 segments and analysed each one. Every key in the
    measures (m) dict now contains a list of that measure for each segment.

    >>> [round(x, 1) for x in m['bpm']]
    [100.0, 96.8, 97.2, 97.9, 96.7, 96.8, 96.8, 95.0, 92.9, 96.7, 99.2]

    Specifying mode = 'fast' will run peak detection once and use detections
    to compute measures over each segment. Useful for speed ups, but typically
    the full mode has better results.

    >>> wd, m = hp.process_segmentwise(data, sample_rate, segment_width=120, segment_overlap=0.5, 
    ... mode = 'fast', replace_outliers = True)

    You can specify the outlier detection method ('iqr' - interquartile range, or 'z-score' for 
    modified z-score approach).
    
    >>> wd, m = hp.process_segmentwise(data, sample_rate, segment_width=120, segment_overlap=0.5, 
    ... mode = 'fast', replace_outliers = True, outlier_method = 'z-score')

    '''

    assert 0 <= segment_overlap < 1.0, 'value error: segment_overlap needs to be \
0 <= segment_overlap < 1.0!'

    assert outlier_method in ['iqr', 'z-score'], 'Unknown outlier detection method specified, \
use either \'iqr\' or \'z-score\''

    s_measures={}
    s_working_data={}

    slice_indices = make_windows(hrdata, sample_rate, segment_width, segment_overlap, segment_min_size)

    if mode == 'full':
        for i, ii in slice_indices:
            try:
                working_data, measures = process(hrdata[i:ii], sample_rate, **kwargs)
                for k in measures.keys():
                    s_measures = append_dict(s_measures, k, measures[k])
                for k in working_data.keys():
                    s_working_data = append_dict(s_working_data, k, working_data[k])
                s_measures = append_dict(s_measures, 'segment_indices', (i, ii))
                s_working_data = append_dict(s_working_data, 'segment_indices', (i, ii))
            except exceptions.BadSignalWarning:
                pass

    elif mode == 'fast':
        working_data, measures = process(hrdata, sample_rate, **kwargs)
        peaklist = np.asarray(working_data['peaklist'])
        for i, ii in slice_indices:
            #pks = [x for x in peaklist if i <= x < ii]
            pks = peaklist[np.where((peaklist >= i) & (peaklist < ii))]
            pks_b = working_data['binary_peaklist'][np.int(np.where(peaklist == pks[0])[0]):
                                                    np.int(np.where(peaklist == pks[-1])[-1]) + 1]
            rr_list = (np.diff(pks) / sample_rate) * 1000.0
            rr_list, rr_diff, rr_sqdiff = calc_rr_segment(rr_list, pks_b)
            _, tmp = calc_ts_measures(rr_list, rr_diff, rr_sqdiff)
            for k in tmp.keys():
                s_measures = append_dict(s_measures, k, tmp[k])
            s_measures = append_dict(s_measures, 'segment_indices', (i, ii))
            s_working_data = append_dict(s_working_data, 'segment_indices', (i, ii))
            s_working_data = append_dict(s_working_data, 'rr_list', rr_list)
            s_working_data = append_dict(s_working_data, 'rr_diff', rr_diff)
            s_working_data = append_dict(s_working_data, 'rr_sqdiff', rr_sqdiff)
            s_working_data = append_dict(s_working_data, 'peaklist', peaklist)

    else:
        raise ValueError('mode not understood! Needs to be either \'fast\' or \'full\', passed: %s' %mode)

    if replace_outliers:
        if outlier_method.lower() == 'iqr':
            for k in s_measures.keys():
                if k not in ['nn20', 'nn50', 'interp_rr_function', 
                             'interp_rr_linspace', 'segment_indices']: #skip these measures
                    s_measures[k], _ = outliers_iqr_method(s_measures[k])
        elif outlier_method.lower() == 'z-score':
            for k in s_measures.keys():
                if k not in ['nn20', 'nn50', 'interp_rr_function', 
                             'interp_rr_linspace', 'segment_indices']: #skip these measures
                    s_measures[k], _ = outliers_modified_z(s_measures[k])

    return s_working_data, s_measures


def process_rr(rr_list, threshold_rr=False, clean_rr=False, 
               clean_rr_method='quotient-filter', calc_freq=False, 
               freq_method='welch', square_spectrum=True, 
               measures={}, working_data={}):
    '''process rr-list

    Function that takes and processes a list of peak-peak intervals.
    Computes all measures as computed by the regular process() function, and
    sets up all dicts required for plotting poincare plots.
    
    Several filtering methods are available as well.

    Parameters
    ----------
    rr_list : 1d array or list
        list or array containing peak-peak intervals (in ms).

    threshold_rr : bool
        if true, the peak-peak intervals are cleaned using a threshold filter, which
        rejects all intervals that differ 30% from the mean peak-peak interval, with
        a minimum of 300ms. 
        default : false

    clean_rr : bool
        if true, the RR_list is further cleaned with an outlier rejection pass. This pass
        is performed after threshold_rr, if that is specified.
        default : false

    clean_rr_method: str
        how to find and reject outliers. Available methods are ' quotient-filter', 
        'iqr' (interquartile range), and 'z-score'.
        default : 'quotient-filter'

    calc_freq : bool
        whether to compute time-series measurements 
        default : False

    freq_method : str
        method used to extract the frequency spectrum. Available: 'fft' (Fourier Analysis), 
        'periodogram', and 'welch' (Welch's method). 
        default : 'welch'

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
        dictionary object used to store temporary values.
    
    measures : dict
        dictionary object used by heartpy to store computed measures.

    Examples
    --------
    Let's generate an RR-list first.

    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(2)
    >>> sample_rate = hp.get_samplerate_datetime(timer, timeformat = '%Y-%m-%d %H:%M:%S.%f')
    >>> wd, m = hp.process(data, sample_rate)
    >>> rr_list = wd['RR_list']

    Using only the RR-list (in ms!) we can now call this function, and let's put the results
    into a differently named container so we're sure all measures are unique:
    >>> wd2, m2 = process_rr(rr_list, threshold_rr = True, clean_rr = True, calc_freq = True)
    >>> '%.3f' %m2['rmssd']
    '45.641'

    If you want to, you can turn off all filters and rejection features:
    >>> wd2, m2 = process_rr(rr_list, threshold_rr = False, clean_rr = False)
    >>> '%.3f' %m2['rmssd']
    '162.645'

    In this case it seems the filtering was necessary: without the RMSSD lies outside the
    range expected in healthy humans.
    '''

    working_data['RR_list'] = rr_list

    if threshold_rr:
        #do thresholding pass
        mean_rr = np.mean(rr_list)
        upper_threshold = mean_rr + 300 if (0.3 * mean_rr) <= 300 else mean_rr + (0.3 * mean_rr)
        lower_threshold = mean_rr - 300 if (0.3 * mean_rr) <= 300 else mean_rr - (0.3 * mean_rr)
        rr_list_cor = [x for x in rr_list if x > lower_threshold and x < upper_threshold]
        rr_mask = [1 if x <= lower_threshold or x >= upper_threshold else 0 for x in rr_list]
        working_data['RR_list_cor'] = rr_list_cor
        working_data['RR_masklist'] = rr_mask

    if clean_rr:
        #do clean_rr pass
        working_data = clean_rr_intervals(working_data = working_data, method = clean_rr_method)

    if not threshold_rr and not clean_rr:
        working_data['RR_list_cor'] = rr_list
        working_data['RR_masklist'] = [0 for i in range(len(rr_list))]
        rr_diff = np.abs(np.diff(rr_list))
        rr_sqdiff = np.power(rr_diff, 2)
    else:
        rr_diff = np.abs(np.diff(working_data['RR_list_cor']))
        rr_sqdiff = np.power(rr_diff, 2)


    #compute ts measures
    working_data, measures = calc_ts_measures(rr_list = working_data['RR_list_cor'], rr_diff = rr_diff, 
                                              rr_sqdiff = rr_sqdiff, measures = measures, 
                                              working_data = working_data)

    measures = calc_poincare(rr_list = working_data['RR_list'], rr_mask = working_data['RR_masklist'], 
                             measures = measures, working_data = working_data)
    if calc_freq:
        #compute freq measures
        working_data, measures = calc_fd_measures(method = freq_method, square_spectrum = square_spectrum,
                                                  measures = measures, working_data = working_data)
        
    return working_data, measures


def run_tests():
    '''
    function to run doctest on all of HeartPy
    '''

    from . import analysis, datautils, filtering, peakdetection, preprocessing, visualizeutils, config
    import doctest
    
    succeeded = 0

    print('testing config')
    results = doctest.testmod(config)
    if results.failed == 0: # pragma: no cover
        print('success!')
        succeeded += 1

    print('testing analysis')
    results = doctest.testmod(analysis)
    if results.failed == 0: # pragma: no cover
        print('success!')
        succeeded += 1
        
    print('testing datautils')
    results = doctest.testmod(datautils)
    if results.failed == 0: # pragma: no cover
        print('success!')
        succeeded += 1

    print('testing filtering')
    results = doctest.testmod(filtering)
    if results.failed == 0: # pragma: no cover
        print('success!') 
        succeeded += 1

    print('testing peakdetection')
    results = doctest.testmod(peakdetection)
    if results.failed == 0: # pragma: no cover
        print('success!')
        succeeded += 1

    print('testing preprocessing')
    results = doctest.testmod(preprocessing)
    if results.failed == 0: # pragma: no cover
        print('success!')
        succeeded += 1

    print('testing visualization utils')
    results = doctest.testmod(visualizeutils)
    if results.failed == 0: # pragma: no cover
        print('success!')
        succeeded += 1

    print('testing main processing pipeline')
    from . import heartpy as hptester
    results = doctest.testmod(hptester)
    if results.failed == 0: # pragma: no cover
        print('success!')
        succeeded += 1

    if succeeded == 8: # pragma: no cover
        print('all tests passed, ready to go!')
    else:
        print('some tests failed...')