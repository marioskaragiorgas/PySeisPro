"""
Trace_analysis.py

This module provides a set of functions to analyze seismic traces using various signal processing techniques.
The primary analyses included are:
  - Periodogram: To estimate the power spectral density of a trace.
  - Welch Periodogram: An improved periodogram method to reduce noise by averaging.
  - Wavelet Transform: Continuous wavelet transform using the Ricker wavelet.
  - Spectrogram: Time-frequency analysis showing how spectral content changes over time.

Functions:
    trace_periodogram(trace, fs): Computes the periodogram of a seismic trace.
    trace_welch_periodogram(trace, fs): Computes the Welch periodogram of a seismic trace.
    trace_wavelet_transform(trace, widths): Computes the continuous wavelet transform (CWT) of a seismic trace.
    trace_spectrogram(trace, fs): Computes the spectrogram (time-frequency representation) of a seismic trace.
"""

import numpy as np
from scipy.signal import periodogram, welch, cwt, ricker, spectrogram
from scipy.fft import fft, fftfreq

def trace_periodogram(trace, fs):

    """
    Compute the periodogram of a seismic trace.

    The periodogram estimates the power spectral density (PSD) of the signal, 
    showing the distribution of power into frequency components composing the signal.

    Parameters:
        trace (ndarray): The seismic trace to analyze.
        fs (float): The sampling frequency of the trace.

    Returns:
        tuple: Contains:
            - f (ndarray): Array of sample frequencies.
            - Pxx (ndarray): Power spectral density of the trace.
    """

    f, Pxx = periodogram(trace, fs)
    return f, Pxx

def trace_welch_periodogram(trace, fs):
    
    """
    Compute the Welch periodogram of a seismic trace.

    The Welch method is an improvement over the standard periodogram by splitting 
    the signal into overlapping segments, computing the periodogram of each segment, 
    and then averaging them. This helps reduce noise in the power spectral density estimate.

    Parameters:
        trace (ndarray): The seismic trace to analyze.
        fs (float): The sampling frequency of the trace.

    Returns:
        tuple: Contains:
            - f (ndarray): Array of sample frequencies.
            - Pxx (ndarray): Power spectral density of the trace using Welch's method.
    """

    f, Pxx = welch(trace, fs)
    return f, Pxx

def trace_wavelet_transform(trace, widths):
    
    """
    Compute the continuous wavelet transform (CWT) of a seismic trace using the Ricker wavelet.

    Wavelet transform provides a time-frequency representation of the signal. 
    The Ricker wavelet, also known as the "Mexican hat" wavelet, is commonly used 
    for seismic analysis because of its similarity to seismic wavelets.

    Parameters:
        trace (ndarray): The seismic trace to analyze.
        widths (ndarray): Widths of the wavelet. This determines the frequency scale.

    Returns:
        ndarray: CWT matrix where each row corresponds to a wavelet transform at a different width.
    """
    
    return cwt(trace, ricker, widths)


def trace_spectrogram(trace, fs):
    
    """
    Compute the spectrogram of a seismic trace.

    The spectrogram is a time-frequency representation that shows how the 
    frequency content of the trace changes over time. It is computed using 
    short-time Fourier transform (STFT).

    Parameters:
        trace (ndarray): The seismic trace to analyze.
        fs (float): The sampling frequency of the trace.

    Returns:
        tuple: Contains:
            - f (ndarray): Array of sample frequencies.
            - t (ndarray): Array of time segments.
            - Sxx (ndarray): Spectrogram of the trace.
    """

    f, t, Sxx = spectrogram(trace, fs)
    return f, t, Sxx