Welcome to HeartPy - Python Heart Rate Analysis Toolkit's documentation!
========================================================================

.. image:: images/output_4.jpg

Welcome to the documentation of the HeartPy, Python Heart Rate Analysis Toolkit. The toolkit is designed to handle (noisy) PPG data collected with either PPG or camera sensors.

* The toolkit was presented at the Humanist 2018 conference in The Hague (`see paper here <https://www.researchgate.net/publication/325967542_Heart_Rate_Analysis_for_Human_Factors_Development_and_Validation_of_an_Open_Source_Toolkit_for_Noisy_Naturalistic_Heart_Rate_Data>`_ ). 

* A technical paper about the functionality `is available here <http://doi.org/10.13140/RG.2.2.24895.56485>`_

**Please cite one or both of these papers when using the toolkit in your research!**

The documentation will help you get up to speed quickly. Follow the :ref:`quickstart` guide for a general overview of how to use the toolkit in only a few lines of code. For a more in-depth review of the module's functionality you can refer to the papers mentioned above, or the :ref:`heart rate analysis` overview.

Example Notebooks are available!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you're looking for a few hands-on examples on how to get started with HeartPy, have a look at the links below! These notebooks show how to handle various analysis tasks with HeartPy, from smartwatch data, smart ring data, regular PPG, and regular (and very noisy) ECG. The notebooks sometimes don't render through the github engine, so either open them locally, or use an online viewer like [nbviewer](https://nbviewer.jupyter.org/).

We recommend you follow the notebooks in order:
- [1. Analysing a PPG signal](https://github.com/paulvangentcom/heartrate_analysis_python/blob/master/examples/1_regular_PPG/Analysing_a_PPG_signal.ipynb), a notebook for starting out with HeartPy using built-in examples.
- [2. Analysing an ECG signal](https://github.com/paulvangentcom/heartrate_analysis_python/blob/master/examples/1_regular_ECG/Analysing_a_regular_ECG_signal.ipynb), a notebook for working with HeartPy and typical ECG data.
- [3. Analysing smartwatch data](https://github.com/paulvangentcom/heartrate_analysis_python/blob/master/examples/smartwatch_data/Analysing_Smartwatch_Data.ipynb), a notebook on analysing low resolution PPG data from a smartwatch.
- [4. Analysing smart ring data](https://github.com/paulvangentcom/heartrate_analysis_python/blob/master/examples/smartring_data/Analysing_Smart_Ring_Data.ipynb), a notebook on analysing smart ring PPG data.
- [5. Analysing noisy ECG data](https://github.com/paulvangentcom/heartrate_analysis_python/blob/master/examples/noisy_ECG/Analysing_Noisy_ECG.ipynb), an advanced notebook on working with very noisy ECG data, using data from the MIT-BIH noise stress test dataset.


Note on using it in scientific research
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Support is available at P.vanGent@tudelft.nl. When using the toolkit in your scientific work: please include me in the process. I can help you implement the toolkit, and the collaboration will also help improve the toolkit so that it can handle more types of data in the future.

Index
======

.. toctree::
   :maxdepth: 3
   :caption: .
   

   quickstart
   apiref
   heartrateanalysis
   algorithmfunctioning
   development
   

