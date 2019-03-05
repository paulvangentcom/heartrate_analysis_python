# A note to our friends from the Humanist 2018 conference

This is the package presented at the Humanist 2018 conference in the Hague, Netherlands. Note that it is always in active development. Whenever including it in your research, pull the latest version from this GitHub. Support is always available at: P.vanGent@tudelft.nl.

**Thank you for your interest**

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

# HeartPy - Python Heart Rate Analysis Toolkit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1324311.svg)](https://doi.org/10.5281/zenodo.1324311)

**HeartPy**, the **Python Heart Rate Analysis Toolkit** is a module for heart rate analysis in Python. It started as pure-python implementation to analyse physiological data taken in naturalistisch driving and cycling experiments.

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
- [ ] Update threshold parameter optimization to handle cases of double peaks and large portion of signal containing little to no hr
- [ ] Report validation performance on repo (published paper + key-points document)
- [ ] Add R-position interpolation ('high accuracy' mode)
- [ ] Handle cases where strong T-peak in signal is present
- [ ] Include example data on PIP release
