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


Frequency domain
* low frequency component (0.04-0.15Hz), LF
* high frequency component (0.16-0.5Hz), HF
* lf/hf ratio, Lf/HF


## Basic analysis example

Import the `heartBeat` module and load a file


```python
import heartBeat as hb

data = hb.get_data("yourdata.csv")
```

This creates a pandas dataframe containing the data from your file.

Analysis requires the sampling rate for your data. If you know this _a priori_, supply it when calling the `process()` function, which returns a `dict{}` object containing all measures:

```python
import heartBeat as hb

data = hb.get_data("yourdata.csv")
measures = hb.process(data, 0.75, fs)
```

`process(dataset, hrw, fs)` requires three arguments:
* **dataset:** the dataset you imported using `get_data`. Alternatively, you can supply a standard pandas dataframe. The heart rate data column should be labeled `hr`. If you wish to use the built-in sample rate detection, the time column should be labeled `datetime`;
* **hrw:** the algorithm uses a moving average during peak detection. `hrw` is the window size used for the calculation. The windowsize is `hrw * samplerate`;
* **fs**: The samplerate of the signal in Hz.

A `dict{}` object is returned containing all measures, and stored in the module. Access as such:

```python
import heartBeat as hb

data = gb.get_data("yourdata.csv")
measures = hb.process(data, 0.75, fs)

print(measures['bpm']) #returns BPM value
print(measures['lf/hf'] # returns LF:HF ratio

#Alternatively, use dictionary stored in module:
print(hb.measures['bpm']) #returns BPM value
print(hb.measures['lf/hf'] # returns LF:HF ratio
```


## Estimating Sample Rate

The toolkit has a simple built-in sample-rate detection. It can handle ms-based timers ([0.0, 10.0, 20.0, N] for a 100Hz signal), and datetime-based timers (in the format: `2016-03-06 09:14:40.650000`, or `yyyy-mm-dd HH:MM:SS.f).
The ms-based timer expects the timer column to be labeled "timer", the datetime timer expects the column named "datetime"

```python
import heartBeat as hb

data = hb.get_data("yourdata.csv")

#if you have a ms-based timer:
fs = hb.get_samplerate_mstimer(data)

#if you have a datetime-based timer:
fs = hb.get_samplerate_datetime(data)
```
In addition to being returned, the samplerate is also stored in the module measures `dict{}`: `print(hb.measures['fs'])`

**Please note:** When using a ms-based timer, 

## Plotting your signal
A basic plotting function is included. It plots the original signal, the moving average, the detected peaks and the rejected peaks (if any were rejected). Usage example with the included `data.csv` example file (recorded at 100Hz):

```python
import heartBeat as hb

data = hb.get_data("data.csv")
measures = hb.process(data, 0.75, 100.0)
hb.plotter(data)
```
This returns:

![output 1 of HR analysis](http://www.paulvangent.com/github/output1.jpeg)

Using the much more noisy data2.csv (containing  only noise for the first 40 sec, then a few noise segments between the signal):

```python
hb.process(dataset, 0.75, get_samplerate_mstimer(dataset))
hb.plotter(dataset)
```

![output 1 of HR analysis](http://www.paulvangent.com/github/output2.jpeg)

Measures are only calculated for non-rejected peaks and intervals between two non-rejected peaks. Rejected detections do not influence the calculated measures.

## License
The module is licensed under the [GNU General Public License Version3, GPL-v3](https://opensource.org/licenses/GPL-3.0)

## To-do

The module is still in active development. The to-do for the coming months is:

1. ~~Replace pandas data handling with numpy data handling, to increase speed of processing~~
2. ~~Drop dependency on pandas~~
3. Implement data handler function, recognising most used formats and parsing correctly
4. Increase versatility of sampling rate detection
5. Improve accuracy of peak detection/rejection with an FFT-based implementation."
6. Add MAD time-domain measure
