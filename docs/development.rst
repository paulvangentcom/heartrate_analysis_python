***********
Development
***********

Release Notes
=============

V0.8.1
~~~~~~

- Added changelog to repository
- Implemented clipping detection and interpolation functionality
- Changed FFT calculation flag to default False, as in some cases the FFT takes very long to compute. Possible causes and fixes to be investigated
- Pushed readthedocs.io documentation source structure to repository
- Added encoding argument to get_data function, per the NumPy deprecation of not using encoding. For more info: https://docs.scipy.org/doc/numpy-1.14.0/release.html#encoding-argument-for-text-io-functions

V0.8.2
~~~~~~

- RR_difference interval no longer taken into account when RR-intervals are not technically adjacent due to rejected peak presence in between
- Moved matplotlib import statement so that it is no longer necessary unless calling the plot functionality, reduces need to install irrelevant dependencies when plotting functionality not needed
- Added Hampel Filter with settable filtersize
- Added method to suppress noisy segments called 'Hampel Corrector', called such as it's simply a Hampel Filter with large window size. Computationally on the expensive side so disabled by default, but very good at suppressing noisy segments without influencing peak positions in the rest of the signal.
- Added breathing rate extraction method. Stores estimated breathing rate in measures['breathingrate']
- Made BPM threshold values settable
- Added Periodogram- and Welch-based PSD estimation
- Added support for edge case where clipping segment starts early in signal, meaning there is insufficient data to interpolate accurately.

V1.0
~~~~
- Released Version 1.0
- Added flag to disable scaling when interpolating clipping segments. Useful for data with large amplitude variations.
- Added marking of rejected segments in plotter
- Added automatic peak rejection when first peak occurs within 150ms, since the signal might start just after a peak, which creates slight inaccuracy.
- Added segment rejection based on percentage of incorrect peaks.

V1.0.1
~~~~~~
- Changed segmentwise rejection API to simplify plotting

V1.1
~~~~
- We are now officially called HeartPy
- Changed overall structure to get rid of global dicts, allows for modular or multithreaded use easier.
- Changed docs to reflect changes

V1.1.2
~~~~~~
- Added high-pass and band-pass Butterworth filters
- Fixed case where no peak-peak differences over 20ms in a signal caused an exception
- Fixed case where intermittent noisy signals resulted in exception when calculating breathing rate
- Added scale_sections() function that uses local scaling rather than global
- Added preprocess_ecg(data, sample_rate) function that attempts to preprocess ecg data. Note: doubles sampling rate of returned data.
- Added highpass and bandpass filtering options. Called through filtersignal function with argument filtertype= lowpass/highpass/bandpass.
- Changed way peak fitting works by adding extra peak validation step between fitting phases

V1.1.3
~~~~~~
- Added functions to allow for continous measure output
- Added make_windows() function to divide input data into evenly sized segments with settable windowsize and settable overlap
- Added two functions to remove outliers from continous set: outliers_modified_z(), and outliers_iqr_method(). Both take a list of numpy array of one continous measure and remove outliers due to incorrectly analysed sections, if any outliers are persent.

V1.1.4
~~~~~~
- Added wrapper function 'process_segmentwise()' that splits hrdata in sections (overlap between sections is settable), and analyses each section separately. Returns two dict objects with all outputs.
- Changed rolling mean function to no longer flatten off the more it is raised, turned out to be more robust.
- Removed peak validation step implemented in V1.1.2 -> after further testing it turned out to be detrimental to many cases.
- Updated docs to reflect the changes to the codebase.

V1.1.5
~~~~~~
- Adapted `make_windows()` to accept tail end of data. Updated `process_segmentise()` to make use of this.
- Updated docs explaining the new functionality
- Fixed error where 'fast' segmentwise method returned single instance in dict rather than sequence
- Fixed error where 'fast' segmentwise method returned empty working_data object
- Started properly structuring module.

V1.1.6
~~~~~~
- moved nn20/nn50 temp containers to 'working_data' dict in stead of output measures (see issue #15).
- fixed wrongly unpacked kwargs in process_segmentwise
- deprecated doctests.txt for now - no longer functional and need updating.
- process_segmentwise() now returns the indices of the slices made in the original data, these are appended to both the returned measures{} and working_data{}. Closes #14
- updated process_segmentwise to allow control over last (likely incomplete) segment. Setting min_size to -1 now puts the tail end of the data into the last bin, making it larger than the others. Closes #16
- fixed sample_rate not being passed to rolmean() when esitmating the breathing rate

V1.1.7
~~~~~~
- added peak interpolation (high precision mode) method 'interpolate_peaks' that allows more accurate estimation of peak positions in signal of low sampling frequency
- in segmentwise processing, fixed bug where the interquartile-range was also used when modified z-score approach was requested.
- fixed mistake in argument order in process_segmentwise function docstring
- implemented 'segment_plotter()' function. This will plot segments and save plots to folder after running 'process_segmentwise()'.
- updated docs to include new functionality.

V1.1.7a
~~~~~~~
- hotfix for process_segmentwise issue where multiple copies of the same index range were placed in the output.

V1.2
~~~~
- Changed organisation HeartPy, it is now split into multiple modules to keep the growing library ordered. This opens the way to the  planned addition of a GUI.
- Added examples that also function as doctests to all functions
- Added extensive documentation docstrings to all functions
- Added function load_exampledata() that loads the available example data files directly from github.
- Added several jupyter notebooks in Examples folder, illustrating how to work with different types of data.
- Added function to reject outliers in RR-list and compute measures based on cleaned list. See: clean_rr_intervals()


Questions
=========
contact me at P.vanGent@tudelft.nl