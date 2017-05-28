# Python Heart Rate Analysis Toolkit

The **Python Heart Rate Analysis Toolkit** is a module for heart rate analysis in Python. It started as pure-python implementation to analyse physiological data taken in naturalistisch driving and cycling experiments.

The module is described in my tutorial series:

* [Analyzing a Discrete Heart Rate Signal Using Python - Part 1](http://www.paulvangent.com/2016/03/15/analyzing-a-discrete-heart-rate-signal-using-python-part-1/)
* [Analyzing a Discrete Heart Rate Signal Using Python - Part 2](http://www.paulvangent.com/2016/03/21/analyzing-a-discrete-heart-rate-signal-using-python-part-2/)
* [Analyzing a Discrete Heart Rate Signal Using Python - Part 3](http://www.paulvangent.com/2016/03/30/analyzing-a-discrete-heart-rate-signal-using-python-part-3/)
* Analyzing a Discrete Heart Rate Signal Using Python - Part 4: in development


The module takes a discrete heart rate signal and outputs time-domain and frequency-domain measures often found in scientific literature:


Time domain:
* beats per minute, BPM
* interbeat interval, IBI
* standard deviation  if intervals between adjacent beats, SDNN
* standard deviation of successive differences between adjacent R-R intervals, SDSD
* root mean square of successive differences between adjacend R-R intervals, RMSSD
* proportion of differences between R-R intervals greater than 20ms, 50ms, pNN20, pNN50
* median absolute deviation, MAD


Frequency domain
* low frequency component (0.04-0.15Hz), LF
* high frequency component (0.16-0.5Hz), HF
* lf/hf ratio, Lf/HF


## Basic analysis example

Import the `heartBeat` module and load a file


```python
import heartBeat as hb

hrdata = hb.get_data('yourdata.csv', column_name = 'hr')
```

This returns a numpy.ndarray.

Analysis requires the sampling rate for your data. If you know this _a priori_, supply it when calling the `process()` function, which returns a `dict{}` object containing all measures:

```python
import heartBeat as hb

data = hb.get_data('yourdata.csv')
measures = hb.process(data, 100.0)
```

`process(dataset, fs, hrw = 0.75)` requires two arguments:
* **dataset:** An 1-dimensional list, numpy array or array-like object containing the heart rate data;
* **fs**: The samplerate of the signal in Hz;
* **hrw:** _optional_ `hrw` is the window size used for the calculation of the moving average. The windowsize is defined as `hrw * samplerate`. Default hrw = 0.75.

A `dict{}` object is returned containing all measures. The object is also stored in the module. Access as such:

```python
import heartBeat as hb

data = hb.get_data('data.csv') 
fs = 100.0 #example file 'data.csv' is sampled at 100.0 Hz

measures = hb.process(data, fs)

print(measures['bpm']) #returns BPM value
print(measures['lf/hf'] # returns LF:HF ratio

#Alternatively, use dictionary stored in module:
print(hb.measures['bpm']) #returns BPM value
print(hb.measures['lf/hf'] # returns LF:HF ratio

#You can also use Pandas if you so desire
import pandas as pd
df = pd.read_csv("data.csv")
measures = hb.process(df['hr'].values, fs)
print("measures['bpm'])
print("measures['lf/hf'])
```



## Getting data from files

The toolkit has functionality to open and parse delimited .csv and .txt files, as well as matlab .mat files. Opening a file is done by the `get_data()` function:

```python
import heartBeat as hb

data = hb.get_data('data.csv')
```

This returns a 1-dimensional `numpy.ndarray` containing the heart rate data.

`get_data(filename, delim = ',', column_name = 'None')` requires one argument:
* **filename:** absolute or relative path to a valid (delimited .csv/.txt or matlab .mat) file;
* **delim** _optional_: when loading a delimited .csv or .txt file, this specifies the delimiter used. Default delim = ',';
* **column_name** _optional_: In delimited files with header: specifying column_name will return data from that column. Not specifying column_name for delimited files will assume the file contains only numerical data, returning np.nan values where data is not numerical. For matlab files: column_name specifies the table name in the matlab file.


Examples:
```python
import heartBeat as hb

#load data from a delimited file without header info
headerless_data = hb.get_data('data.csv')

#load data from column labeles 'hr' in a delimited file with header info
headered_data = hb.get_data('data.csv', column_name = 'hr')

#load matlab file
matlabdata = hb.get_data('data2.mat', column_name = 'hr')
```


## Estimating Sample Rate

The toolkit has a simple built-in sample-rate detection. It can handle ms-based timers and datetime-based timers.

```python
import heartBeat as hb

#if you have a ms-based timer:
fs = hb.get_samplerate_mstimer(mstimer_data)

#if you have a datetime-based timer:
fs = hb.get_samplerate_datetime(datetime_data, timeformat='%Y-%m-%d %H:%M:%s.%f')
```


`get_samplerate_mstimer(timerdata)` requires one argument:
* **timerdata:** a list, numpy array or array-like object containing ms-based timestamps (float or int).


`get_samplerate_datetime(datetimedata, timeformat = '%H:%M:%S.f')` requires one argument:
* **datetimedata:** a list, numpy array or array-like object containing datetime-based timestamps (string);
* **timeformat** _optional_: the format of the datetime-strings in datetimedata. Default timeformat = '%H:%M:%S.f', 24-hour based time including ms: 21:43:12.569.

## Plotting your signal

A plotting function is included. It plots the original signal and overlays the detected peaks and the rejected peaks (if any were rejected). 

Example with the included `data.csv` example file (recorded at 100.0Hz):

```python
import heartBeat as hb

data = hb.get_data('data.csv')
measures = hb.process(data, 100.0)
hb.plotter()
```
This returns:

![output 1 of HR analysis](https://github.com/paulvangentcom/heartrate_analysis_python/blob/master/images/output1.jpeg)

`plotter(show = True, title = 'Heart Rate Signal Peak Detection')` has two optional arguments:
* **show** _optional_: if set to True a plot is visualised, if set to False a matplotlib.pyplot object is returned. Default show = True;
* **title** _optional_: Sets the title of the plot. If not specified, default title is used.

Examples:

```python
import heartBeat as hb
hrdata = hb.get_data('data2.csv', column_name = 'hr')
timerdata = hb.get_data('data2.csv., column_name = 'timer')

hb.process(dataset, hb.get_samplerate_mstimer(timerdata))

#plot with different title
hb.plotter(title='Heart Beat Detection on Noisy Signal')
```

![output 1 of HR analysis](https://github.com/paulvangentcom/heartrate_analysis_python/blob/master/images/output2.jpeg)

Measures are only calculated for non-rejected peaks and intervals between two non-rejected peaks. Rejected detections do not influence the calculated measures.

By default a plot is visualised when plotter() is called. The function returns a matplotlib.pyplot object if the argument show=False is passed:

```python
hb.process(dataset, 0.75, get_samplerate_mstimer(dataset))
hb.plotter(show=False)
```
This returns:
```
<module 'matplotlib.pyplot' from '<pythonpath>\lib\site-packages\matplotlib\pyplot.pyc'>
```

Object can then be saved or visualised:
```python
hb.process(dataset, 0.75, get_samplerate_mstimer(dataset))
plot_object = hb.plotter(show=False)

plot_object.savefig('plot_1.jpg') #saves the plot as JPEG image.

plt.object.show() #displays plot
```



## License
The module is licensed under the [GNU General Public License Version3, GPL-v3](https://opensource.org/licenses/GPL-3.0)

## To-do

The module is still in active development. The to-do for the coming months is:

- [X] Replace numerical work with numpy functions, to increase speed of processing
- [X] Drop dependency on pandas
- [X] Implement data handler function, recognising most used formats and parsing correctly
- [ ] Make values for inaccurate bpm rejection settable
- [X] Increase versatility of sampling rate detection
- [ ] Improve accuracy of peak detection/rejection with an FFT-based implementation.
- [X] Add MAD time-domain measure
