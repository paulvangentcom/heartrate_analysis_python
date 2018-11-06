.. _quickstart:

****************
Quickstart Guide
****************

Basic Example
=============
Import the `HeartPy` module and load a file


.. code-block:: python

    import heartpy as hp

    hrdata = hp.get_data('data.csv')


This returns a :code:`numpy.ndarray`.

Analysis requires the sampling rate for your data. If you know this *a priori*, supply it when calling the `process()` function, which returns a `dict{}` object containing all measures:

.. code-block:: python

    import heartpy as hp

    data = hp.get_data('data.csv') #data.csv is sampled at 100Hz
    working_data, measures = hp.process(data, 100.0)


**process(dataset, fs, windowsize=0.75, report_time=False,
calc_freq=False, freq_method='welch', interp_clipping=True, 
clipping_scale=True, interp_threshold=1020, hampel_correct=False, 
bpmmin=40, bpmmax=180, reject_segmentwise=False, measures = {},
working_data = {})**
               
requires two arguments:

* **dataset:** An 1-dimensional list, numpy array or array-like object containing the heart rate data;
* **fs**: The samplerate of the signal in Hz;

Several optional arguments are available:

* **windowsize:** _optional_ `windowsize` is the window size used for the calculation of the moving average. The windowsize is defined as `windowsize * samplerate`. Default windowsize=0.75.
* **report_time:** _optional_ whether to report total processing time of process() loop.
* **calc_fft:** _optional_ whether to calculate frequency domain measures. Default = false Note: can cause slowdowns in some cases.
* **calc_freq:** _optional_ whether to calculate frequency domain measures. Default = false Note: can cause slowdowns in some cases.
* **freq_method:** _optional_ method used to extract the frequency spectrum. Available: 'fft' (Fourier Analysis), 'periodogram', and 'welch' (Welch's method), Default = 'welch'
* **interp_clipping:** if True, clipping parts of the signal are identified and the implied peak shape is interpolated. Default=True
* **clipping_scale:** whether to scale the data priod to clipping detection. Can correct errors if signal amplitude has been affected after digitization (for example through filtering). Default = false
* **interp_threshold**: the amplitude threshold beyond which will be checked for clipping. Recommended is to take this as the maximum value of the ADC with some margin for signal noise (default 1020, default ADC max 1024) 
* **hampel_correct:** whether to reduce noisy segments using large median filter. Disabled by default due to computational complexity, and generally it is not necessary. Default = false.
* **bpmmin:** minimum value to see as likely for BPM when fitting peaks. Default = 40
* **bpmmax:** maximum value to see as likely for BPM when fitting peaks. Default = 180
* **reject_segmentwise:** whether to reject segments with more than 30% rejected beats. By default looks at segments of 10 beats at a time. Default = false.
* **measures:** measures dict in which results are stored. Custom dictionary can be passed, otherwise one is created and returned.
* **working_data:** working_data dict in which results are stored. Custom dictionary can be passed, otherwise one is created and returned.

Two :code:`dict{}` objects are returned: one working data dict, and one containing all measures. Access as such:

.. code-block:: python

    import heartpy as hp

    data = hp.get_data('data.csv') 
    fs = 100.0 #example file 'data.csv' is sampled at 100.0 Hz

    working_data, measures = hp.process(data, fs, report_time=True)

    print(measures['bpm']) #returns BPM value
    print(measures['rmssd']) # returns RMSSD HRV measure

    #You can also use Pandas if you so desire
    import pandas as pd
    df = pd.read_csv("data.csv", names=['hr'])
    #note we need calc_freq if we want frequency-domain measures
    working_data, measures = hp.process(df['hr'].values, fs, calc_freq=True)
    print(measures['bpm'])
    print(measures['lf/hf'])

    
Getting Data From Files
=======================
The toolkit has functionality to open and parse delimited .csv and .txt files, as well as matlab .mat files. Opening a file is done by the :code:`get_data()` function:

.. code-block:: python

    import heartpy as hp

    data = hp.get_data('data.csv')

This returns a 1-dimensional :code:`numpy.ndarray` containing the heart rate data.

:code:`get_data(filename, delim = ',', column_name = 'None')` requires one argument:

* **filename:** absolute or relative path to a valid (delimited .csv/.txt or matlab .mat) file;

Several optional arguments are available:

* **delim** _optional_: when loading a delimited .csv or .txt file, this specifies the delimiter used. Default delim = ',';
* **column_name** _optional_: In delimited files with header: specifying column_name will return data from that column. Not specifying column_name for delimited files will assume the file contains only numerical data, returning np.nan values where data is not numerical. For matlab files: column_name specifies the table name in the matlab file.


Examples:

.. code-block:: python

    import heartpy as hp

    #load data from a delimited file without header info
    headerless_data = hp.get_data('data.csv')

    #load data from column labeles 'hr' in a delimited file with header info
    headered_data = hp.get_data('data2.csv', column_name = 'hr')

    #load matlab file
    matlabdata = hp.get_data('data2.mat', column_name = 'hr')
    #note that the column_name here represents the table name in the matlab file
        

Estimating Sample Rate
======================
The toolkit has a simple built-in sample-rate detection. It can handle ms-based timers and datetime-based timers.

.. code-block:: python

    import heartpy as hp

    #if you have a ms-based timer:
	mstimer_data = hp.get_data('data2.csv', column_name='timer')
    fs = hp.get_samplerate_mstimer(mstimer_data)
	print(fs)

    #if you have a datetime-based timer:
	datetime_data = hp.get_data('data3.csv', column_name='datetime')
    fs = hp.get_samplerate_datetime(datetime_data, timeformat='%Y-%m-%d %H:%M:%S.%f')
	print(fs)


:code:`get_samplerate_mstimer(timerdata)` requires one argument:

* **timerdata:** a list, numpy array or array-like object containing ms-based timestamps (float or int).


:code:`get_samplerate_datetime(datetimedata, timeformat = '%H:%M:%S.f')` requires one argument:

* **datetimedata:** a list, numpy array or array-like object containing datetime-based timestamps (string);

One optional argument is available:

* **timeformat** _optional_: the format of the datetime-strings in your dataset. Default timeformat='%H:%M:%S.f', 24-hour based time including ms: 21:43:12.569.


Plotting Results
================
A plotting function is included. It plots the original signal and overlays the detected peaks and the rejected peaks (if any were rejected). 

Example with the included `data.csv` example file (recorded at 100.0Hz):

.. code-block:: python

    import heartpy as hp

    data = hp.get_data('data.csv')
    working_data, measures = hp.process(data, 100.0)
    hp.plotter(working_data, measures)

This returns:

.. image:: images/output1.jpeg

:code:`plotter(working_data, measures, show = True, title = 'Heart Rate Signal Peak Detection')` has two required arguments:

* **working_data** The working data :code:`dict{}` container returned by the :code:`process()` function.
* **measures** The measures :code:`dict{}` container returned by the :code:`process()` function.

Several optional arguments are available:

* **show** _optional_: if set to True a plot is visualised, if set to False a matplotlib.pyplot object is returned. Default show = True;
* **title** _optional_: Sets the title of the plot. If not specified, default title is used.

**Examples:**

.. code-block:: python

    import heartpy as hp
    hrdata = hp.get_data('data2.csv', column_name='hr')
    timerdata = hp.get_data('data2.csv', column_name='timer')

    working_data, measures = hp.process(hrdata, hp.get_samplerate_mstimer(timerdata))

    #plot with different title
    hp.plotter(working_data, measures, title='Heart Beat Detection on Noisy Signal')


.. image:: images/output2.jpeg

Measures are only calculated for non-rejected peaks and intervals between two non-rejected peaks. Rejected detections do not influence the calculated measures.

By default a plot is visualised when plotter() is called. The function returns a matplotlib.pyplot object if the argument show=False is passed:

.. code-block:: python

    working_data, measures = hp.process(hrdata, hp.get_samplerate_mstimer(timerdata))
    plot_object = hp.plotter(working_data, measures, show=False)

This returns:

.. code-block:: python

    <module 'matplotlib.pyplot' [...]>

Object can then be saved, appended to, or visualised:

.. code-block:: python

    working_data, measures = hp.process(hrdata, hp.get_samplerate_mstimer(timerdata))
    plot_object = hp.plotter(working_data, measures, show=False)

    plot_object.savefig('plot_1.jpg') #saves the plot as JPEG image.

    plot_object.show() #displays plot 
      