'''
Functions that help visualize results
'''

import os

import matplotlib.pyplot as plt


__all__ = ['plotter',
           'segment_plotter']

def plotter(working_data, measures, show=True, title='Heart Rate Signal Peak Detection'): # pragma: no cover
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

    title : string
        title for the plot.
        default : "Heart Rate Signal Peak Detection"

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
    peaklist = working_data['peaklist']
    ybeat = working_data['ybeat']
    rejectedpeaks = working_data['removed_beats']
    rejectedpeaks_y = working_data['removed_beats_y']
    plt.title(title)
    plt.plot(working_data['hr'], alpha=0.5, color='blue', label='heart rate signal')
    plt.plot(working_data['rolling_mean'])
    plt.scatter(peaklist, ybeat, color='green', label='BPM:%.2f' %(measures['bpm']))
    plt.scatter(rejectedpeaks, rejectedpeaks_y, color='red', label='rejected peaks')
    
    #check if rejected segment detection is on and has rejected segments
    try:
        if len(working_data['rejected_segments']) >= 1:
            for segment in working_data['rejected_segments']:
                plt.axvspan(segment[0], segment[1], facecolor='red', alpha=0.5)
    except:
        pass

    plt.legend(loc=4, framealpha=0.6)
    if show:
        plt.show()
    else:
        return plt

def segment_plotter(working_data, measures, title='Heart Rate Signal Peak Detection',
                    path = '', start=0, end=None, step=1): # pragma: no cover
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
    if not path.endswith('/') or path.endswith('\\'):
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
        m_segment['bpm'] = measures['bpm'][i]
        try:
            wd_segment['rejected_segments'] = working_data['rejected_segments'][i]
        except:
            pass

        #plot it using built-in plotter
        p = plotter(wd_segment, m_segment, show=False)
        p.savefig('%s%i.png' %(path, filenum))
        p.close()
        filenum += 1