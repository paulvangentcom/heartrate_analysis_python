"""Load example data for experimenting and training."""

import os

import numpy as np

from heartpy.heartpy import get_data


__all__ = [
    'load_example_basic',
    'load_example_milliseconds',
    'load_example_datetime'
]


def _get_data_dir():
    """Return the location of the data folder."""
    return os.path.dirname(os.path.realpath(__file__))


def load_example_basic():
    """Load signal without time reference.

    Returns
    -------
    numpy.ndarray
        Return amplitudes in array.

    Example
    -------

    ```
    > from heartpy.datasets import load_example_basic
    amplitude = load_example_basic()
    ```
    """
    return get_data(os.path.join(_get_data_dir(), "data", "data_basic.csv"))


def load_example_milliseconds():
    """Load signal and milliseconds sequence.

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Return milliseconds and magnitude in arrays.

    Example
    -------

    ```
    > from heartpy.datasets import load_example_milliseconds
    amplitude, time_milliseconds = load_example_milliseconds()
    ```
    """
    fp = os.path.join(_get_data_dir(), "data", "data_milliseconds.csv")

    data = np.loadtxt(fp, delimiter=",", skiprows=1, dtype=np.float64)

    return data[:, 1], data[:, 0]


def load_example_datetime():
    """Load signal and datetime.

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Return milliseconds and magnitude in arrays.

    Example
    -------

    ```
    > from heartpy.datasets import load_example_datetime
    amplitude, time_datetime = load_example_datetime()
    ```
    """
    fp = os.path.join(_get_data_dir(), "data", "data_datetime.csv")

    dt_values = np.loadtxt(
        fp, delimiter=",", skiprows=1, usecols=0, dtype=np.datetime64)
    amplitudes = np.loadtxt(
        fp, delimiter=",", skiprows=1, usecols=1, dtype=np.float64)

    return amplitudes, dt_values
