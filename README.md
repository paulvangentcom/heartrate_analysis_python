# HeartPy - Python Heart Rate Analysis Toolkit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1324311.svg)](https://doi.org/10.5281/zenodo.1324311) [![Build Status](https://travis-ci.org/paulvangentcom/heartrate_analysis_python.svg?branch=master)](https://travis-ci.org/paulvangentcom/heartrate_analysis_python) [![codecov](https://codecov.io/gh/paulvangentcom/heartrate_analysis_python/branch/master/graph/badge.svg)](https://codecov.io/gh/paulvangentcom/heartrate_analysis_python) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/heartpy)



# Structural update

HeartPy V1.2 has landed! The structure of the package has been reworked to be in separate modules now in preparation of the next big update, which will feature many analysis expansions and the first steps towards a GUI for HeartPy. HeartPy has been growing steadily and had reached the point where it became cluttered and unwieldy to keep in a single file. The API remains unchanged.

An 'Examples' folder has been added to the repo which will be expanded soon. Now there's two notebooks explaining how to analyse ppg signals from smartwatches and smart rings.

# Installation
```
python -m setup.py install
```

Alternatively, we're also on PIP:
```
python -m pip install heartpy
```

That's it!

# Official Documentation

The official documentation is online! [You can find the official documentation here](https://python-heart-rate-analysis-toolkit.readthedocs.io)

# Python 2.7
The module compiles and and runs fine on Python 2.7, **but** as of now the unit tests fail due to different formatting of outputs between Python 2 and Python 3. You can still install and use HeartPy on Python 2.7 if you want.

# More information
**HeartPy**, the **Python Heart Rate Analysis Toolkit** is a module for heart rate analysis in Python. It started as pure-python implementation to analyse physiological data taken in naturalistic driving and cycling experiments.

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

**When using the package in your research, please cite**:

van Gent, P., Farah, H., van Nes, N., & van Arem, B. (2018). Heart Rate Analysis for Human Factors: Development and Validation of an Open Source Toolkit for Noisy Naturalistic Heart Rate Data. In Proceedings of the 6th HUMANIST Conference (pp. 173–178).

van Gent, P. Van, Farah, H., Nes, N. Van, & Arem, B. Van. (manuscript submitted). Analysing Noisy Driver Physiology Real-Time Using Off-the-Shelf Sensors: Heart Rate Analysis Software from the Taking the Fast Lane Project.. Doi: doi.org/10.13140/RG.2.2.24895.56485


## Documentation

[You can find the official documentation here](https://python-heart-rate-analysis-toolkit.readthedocs.io)

The module is also to some extent described in my tutorial series:

* [Analyzing a Discrete Heart Rate Signal Using Python - Part 1](http://www.paulvangent.com/2016/03/15/analyzing-a-discrete-heart-rate-signal-using-python-part-1/)
* [Analyzing a Discrete Heart Rate Signal Using Python - Part 2](http://www.paulvangent.com/2016/03/21/analyzing-a-discrete-heart-rate-signal-using-python-part-2/)
* [Analyzing a Discrete Heart Rate Signal Using Python - Part 3](http://www.paulvangent.com/2016/03/30/analyzing-a-discrete-heart-rate-signal-using-python-part-3/)
* Analyzing a Discrete Heart Rate Signal Using Python - Part 4: in development


## License
The module is licensed under the [GNU General Public License Version3, GPL-v3](https://opensource.org/licenses/GPL-3.0)

## Validation
Initial results of the validation have been reported in [1, 2]. Updates here are soon to follow once the papers are published.


[1]van Gent, P., Farah, H., van Nes, N., & van Arem, B. (2018). Heart Rate Analysis for Human Factors: Development and Validation of an Open Source Toolkit for Noisy Naturalistic Heart Rate Data. In Proceedings of the 6th HUMANIST Conference (pp. 173–178).

[2] van Gent, P. Van, Farah, H., Nes, N. Van, & Arem, B. Van. (manuscript submitted for publication). A Novel Heart Rate Algorithm for the Analysis of Noisy Signals.



## To-do

The module is still in active development. See the changelog for past changes. The to-do for the coming months is:

to do before V1.2
- [X] Validate performance on Physionet datasets
- [ ] Add several extra filtering options (Savitzky-Golay, Notch)
- [ ] Add convolutional pre-processing pipeline for both PPG and ECG
- [ ] Mark which RR-intervals are adjacent (no rejection gaps in between) to improve things like breathing rate detection
- [ ] Update threshold parameter optimization to handle cases of double peaks and large portion of signal containing little to no hr
- [ ] Report validation performance on repo (published paper + key-points document once published)
- [X] Add R-position interpolation ('high accuracy' mode)
- [X] Handle cases where strong T-peak in signal is present
- [X] Include example data on PIP release
- [ ] Change peak fitting to optimize locally rather than globally, to handle long signals with noisy intermediate sections better
- [X] Change segmented analysis method to run peak detection once on whole signal, then do segmentwise computation from working_data object (drastic speed up)
