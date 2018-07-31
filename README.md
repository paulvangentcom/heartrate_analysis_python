# A note to our friends from the Humanist 2018 conference

This is the package presented at the Humanist 2018 conference in the Hague, Netherlands. Note that it is always in active development. Whenever including it in your research, pull the latest version from this GitHub. Support is always available at: P.vanGent@tudelft.nl.

**Thank you for your interest**

# Official Documentation

The official documentation is online! [You can find the official documentation here](https://python-heart-rate-analysis-toolkit.readthedocs.io)

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

## Documentation

[You can find the official documentation here](https://python-heart-rate-analysis-toolkit.readthedocs.io)

## License
The module is licensed under the [GNU General Public License Version3, GPL-v3](https://opensource.org/licenses/GPL-3.0)

## To-do

The module is still in active development. See the changelog for past changes. The to-do for the coming months is:

to do before V1.1
- [ ] Validate performance on Physionet datasets
- [ ] Add R-position interpolation ('high accuracy' mode)
- [ ] Handle cases where strong T-peak in signal is present