'''
Functions that help visualize results
'''

import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

from . import config

__all__ = ['plotter',
           'segment_plotter',
           'plot_poincare',
           'plot_breathing']

def plotter(working_data, measures, show=True, figsize=None,
            title='Heart Rate Signal Peak Detection', moving_average=False): # pragma: no cover
    '''plots the analysis results.

    Function that uses calculated measures and data stored in the working_data{} and measures{}
    dict objects to visualise the fitted peak detection solution.

    Parameters
    ----------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function

    show : bool
        when False, function will return a plot object rather than display the results.
        default : True

    figsize: tuple
        Set dimensions of image in inches like in matplotlib. figsize=(x, y)
        default: None => (6.4, 4.8)

    title : string
        title for the plot.
        default : "Heart Rate Signal Peak Detection"

    moving_average : bool
        whether to display the moving average on the plot.
        The moving average is used for peak fitting.
        default: False

    Returns
    -------
    out : matplotlib plot object
        only returned if show == False.

    Examples
    --------
    First let's load and analyse some data to visualise

    >>> import heartpy as hp
    >>> data, _ = hp.load_exampledata(0)
    >>> wd, m = hp.process(data, 100.0)

    Then we can visualise

    >>> plot_object = plotter(wd, m, show=False, title='some awesome title')

    This returns a plot object which can be visualized or saved or appended.
    See matplotlib API for more information on how to do this.

    A matplotlib plotting object is returned. This can be further processed and saved
    to a file.

    '''
    #get color palette
    colorpalette = config.get_colorpalette_plotter()

    # create plot x-var
    fs = working_data['sample_rate']
    plotx = np.arange(0, len(working_data['hr'])/fs, 1/fs)
    #check if there's a rounding error causing differing lengths of plotx and signal
    diff = len(plotx) - len(working_data['hr'])
    if diff < 0:
        #add to linspace
        plotx = np.append(plotx, plotx[-1] + (plotx[-2] - plotx[-1]))
    elif diff > 0:
        #trim linspace
        plotx = plotx[0:-diff]
        
    peaklist = working_data['peaklist']
    ybeat = working_data['ybeat']
    rejectedpeaks = working_data['removed_beats']
    rejectedpeaks_y = working_data['removed_beats_y']

    fig, ax = plt.subplots(figsize=figsize)

    ax.set_title(title)
    ax.plot(plotx, working_data['hr'], color=colorpalette[0], label='heart rate signal', zorder=-10)
    ax.set_xlabel('Time (s)')

    if moving_average:
        ax.plot(plotx, working_data['rolling_mean'], color='gray', alpha=0.5)

    ax.scatter(np.asarray(peaklist)/fs, ybeat, color=colorpalette[1], label='BPM:%.2f' %(measures['bpm']))
    ax.scatter(rejectedpeaks/fs, rejectedpeaks_y, color=colorpalette[2], label='rejected peaks')

    #check if rejected segment detection is on and has rejected segments
    try:
        if len(working_data['rejected_segments']) >= 1:
            for segment in working_data['rejected_segments']:
                ax.axvspan(segment[0], segment[1], facecolor='red', alpha=0.5)
    except:
        pass

    ax.legend(loc=4, framealpha=0.6)

    if show:
        fig.show()
    else:
        return fig

def segment_plotter(working_data, measures, title='Heart Rate Signal Peak Detection',
                    figsize=(6, 6), path='', start=0, end=None, step=1): # pragma: no cover
    '''plots analysis results

    Function that plots the results of segmentwise processing of heart rate signal
    and writes all results to separate files at the path provided.

    Parameters
    ----------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function

    title : str
        the title used in the plot

    figsize : tuple
        figsize tuple to be passed to matplotlib

    path : str
        the path where the files will be stored, folder must exist.

    start : int
        what segment to start plotting with
        default : 0

    end : int
        last segment to plot. Must be smaller than total number of segments
        default : None, will plot until end

    step : int
        stepsize used when iterating over plots every step'th segment will be plotted
        default : 1

    Returns
    -------
        None

    Examples
    --------
    This function has no examples. See documentation of heartpy for more info.
    '''
    #sanity check
    assert 0 < step < len(working_data['hr']), 'step must be larger than zero and smaller than total number of segments'

    #set endpoint if not explicitly defined
    if end == None:
        end = len(working_data['hr'])
    else:
        #make sure it is defined within boundary conditions
        assert end <= len(working_data['hr']), 'defined "end" endpoint is larger than number of segments'

    #add trailing path slash if user omitted it
    if not (path.endswith('/') or path.endswith('\\')) and len(path) > 0:
        path += '/'
        #create path if it doesn't exist
        if not os.path.isdir(path):
            os.makedirs(path)

    #make plots
    filenum = 0
    for i in range(start, end, step):
        wd_segment = {}
        m_segment = {}
        #assign values to sub-object for plotting purposes
        wd_segment['peaklist'] = working_data['peaklist'][i]
        wd_segment['ybeat'] = working_data['ybeat'][i]
        wd_segment['removed_beats'] = working_data['removed_beats'][i]
        wd_segment['removed_beats_y'] = working_data['removed_beats_y'][i]
        wd_segment['hr'] = working_data['hr'][i]
        wd_segment['rolling_mean'] = working_data['rolling_mean'][i]
        wd_segment['sample_rate'] = working_data['sample_rate'][i]
        m_segment['bpm'] = measures['bpm'][i]
        try:
            wd_segment['rejected_segments'] = working_data['rejected_segments'][i]
        except:
            pass

        #plot it using built-in plotter
        plt.figure(figsize = figsize)
        p = plotter(wd_segment, m_segment, show=False)
        p.savefig('%s%i.png' %(path, filenum))
        plt.close('all')
        filenum += 1


def plot_poincare(working_data, measures, show = True, figsize=None,
                  title='Poincare plot'): # pragma: no cover
    '''visualize poincare plot

    function that visualises poincare plot.

    Parameters
    ----------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function

    show : bool
        whether to show the plot right away, or return a matplotlib object for
        further manipulation

    figsize: tuple
        Set dimensions of image in inches like in matplotlib. figsize=(x, y)
        default: None => (6.4, 4.8)

    title : str
        the title used in the plot

    Returns
    -------
    out : matplotlib plot object
        only returned if show == False.

    Examples
    --------
    This function has no examples. See documentation of heartpy for more info.
    '''

    #get color palette
    colorpalette = config.get_colorpalette_poincare()

    #get values from dict
    x_plus = working_data['poincare']['x_plus']
    x_minus = working_data['poincare']['x_minus']
    sd1 = measures['sd1']
    sd2 = measures['sd2']

    #define figure
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'}, figsize=figsize)

    #plot scatter
    ax.scatter(x_plus, x_minus, color = colorpalette[0],
                alpha = 0.75, label = 'peak-peak intervals')

    #plot identity line
    mins = np.min([x_plus, x_minus])
    maxs = np.max([x_plus, x_minus])
    identity_line = np.linspace(np.min(mins), np.max(maxs))
    ax.plot(identity_line, identity_line, color='black', alpha=0.5,
             label = 'identity line')

    #rotate SD1, SD2 vectors 45 degrees counterclockwise
    sd1_xrot, sd1_yrot = rotate_vec(0, sd1, 45)
    sd2_xrot, sd2_yrot = rotate_vec(0, sd2, 45)

    #plot rotated SD1, SD2 lines
    ax.plot([np.mean(x_plus), np.mean(x_plus) + sd1_xrot],
             [np.mean(x_minus), np.mean(x_minus) + sd1_yrot],
             color = colorpalette[1], label = 'SD1')
    ax.plot([np.mean(x_plus), np.mean(x_plus) - sd2_xrot],
             [np.mean(x_minus), np.mean(x_minus) + sd2_yrot],
             color = colorpalette[2], label = 'SD2')

    #plot ellipse
    xmn = np.mean(x_plus)
    ymn = np.mean(x_minus)
    el = Ellipse((xmn, ymn), width = sd2 * 2, height = sd1 * 2, angle = 45.0)
    ax.add_artist(el)
    el.set_edgecolor((0,0,0))
    el.fill = False

    ax.set_xlabel(r'RRi$_n$ (ms)')
    ax.set_ylabel(r'RRi$_{n+1}$ (ms)')
    ax.legend(loc=4, framealpha=0.6)
    ax.set_title(title)

    if show:
        fig.show()
    else:
        return fig


def rotate_vec(x, y, angle):
    '''rotates vector around origin point

    Function that takes vector and angle, and rotates around origin point
    with given amount of degrees.

    Helper function for poincare plotting

    Parameters
    ----------
    x : int or float
        vector x coordinate

    y : int or float
        vector y coordinate

    angle: int or float
        the angle of rotation applied to the vecftor

    Returns
    -------
    x_rot : float
        new x coordinate with rotation applied

    y_rot : float
        new x coordinate with rotation applied

    Examples
    --------
    Given a vector (0,1), if we apply a rotation of 90 degrees clockwise
    we expect to get (1,0). Let's test

    >>> x_new, y_new = rotate_vec(0, 1, -90)
    >>> print('%.3f, %.3f' %(x_new, y_new))
    1.000, 0.000
    '''
    theta = np.radians(angle)

    cs = np.cos(theta)
    sn = np.sin(theta)

    x_rot = (x * cs) - (y * sn)
    y_rot = (x * sn) + (y * cs)

    return x_rot, y_rot


def plot_breathing(working_data, measures, show=True, figsize=None): # pragma: no cover
    '''plots extracted breathing signal and spectrogram

    Function that plots the breathing signal extracted from RR-intervals alongside
    its computed spectrogram representation.

    Parameters
    ----------
    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function

    show : bool
        whether to show the plot right away, or return a matplotlib object for
        further manipulation

    figsize: tuple
        Set dimensions of image in inches like in matplotlib. figsize=(x, y)
        default: None => (6.4, 4.8)

    Returns
    -------
    out : matplotlib plot object
        only returned if show == False.

    Examples
    --------
    This function has no examples. See documentation of heartpy for more info.
    '''

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

    ax1.plot(working_data['breathing_signal'], label='breathing signal')
    ax1.set_xlabel('ms')
    ax1.set_title('breathing signal extracted from RR-intervals')

    ax2.plot(working_data['breathing_frq'], working_data['breathing_psd'], label='spectrogram')
    ax2.set_xlim(0, 1)
    ax2.set_xlabel('Hz')
    ax2.set_title('spectrogram extracted from breathing rate signal')

    ax2.legend()
    plt.tight_layout()

    if show:
        fig.show()
    else:
        return fig

