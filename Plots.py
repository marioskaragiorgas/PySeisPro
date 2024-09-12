"""
Plots.py

This module contains functions for visualizing seismic data. It provides a variety of plotting methods, 
including plotting seismic traces, periodograms, Welch periodograms, wavelet transforms, spectrograms, 
and seismic images. The functions utilize Matplotlib for rendering the visualizations.

Functions:
    plot_trace(ax, trace, trace_number, delta): Plots a seismic trace.
    plot_periodogram(ax, f, Pxx, trace_number): Plots the periodogram of a seismic trace.
    plot_welch_periodogram(ax, f, Pxx, trace_number): Plots the Welch periodogram of a seismic trace.
    plot_wavelet_transform(ax, cwt_matrix, widths, trace_number): Plots the wavelet transform of a seismic trace.
    plot_spectrogram(ax, f, t, Sxx, trace_number): Plots the spectrogram of a seismic trace.
    plot_seismic_image(ax, seismic_data): Displays a seismic section as an image.
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_trace(ax, trace, trace_number, delta):
    
    """
    Plot the seismic trace on the given axis.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the trace.
        trace (ndarray): The seismic trace data (amplitude values).
        trace_number (int): The index of the seismic trace to be plotted.
        delta (float): The time interval between samples (sampling period).

    Returns:
        None: The trace is plotted on the provided Matplotlib axis.
    """

    time = np.arange(len(trace)) * delta
    ax.plot(time, trace, color='black')
    ax.set_title(f"Seismic Trace {trace_number}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")

def plot_periodogram(ax, f, Pxx, trace_number):
    
    """
    Plot the periodogram of a seismic trace.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the periodogram.
        f (ndarray): Array of sample frequencies.
        Pxx (ndarray): Power spectral density (PSD) values for each frequency.
        trace_number (int): The index of the seismic trace being analyzed.

    Returns:
        None: The periodogram is plotted on the provided Matplotlib axis.
    """

    ax.semilogy(f, Pxx, color = 'black')
    ax.set_title(f'Seismic Trace {trace_number} Periodogram')
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power/Frequency (dB/Hz)")

def plot_welch_periodogram(ax, f, Pxx, trace_number):
    
    """
    Plot the Welch periodogram of a seismic trace.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the Welch periodogram.
        f (ndarray): Array of sample frequencies.
        Pxx (ndarray): Power spectral density (PSD) values for each frequency using Welch's method.
        trace_number (int): The index of the seismic trace being analyzed.

    Returns:
        None: The Welch periodogram is plotted on the provided Matplotlib axis.
    """

    ax.semilogy(f, Pxx, color = 'black')
    ax.set_title(f"Seismic Trace {trace_number} Welch Periodogram")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power/Frequency (dB/Hz)")

def plot_wavelet_transform(ax, cwt_matrix, widths, trace_number):
    
    """
    Plot the wavelet transform of a seismic trace.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the wavelet transform.
        cwt_matrix (ndarray): The continuous wavelet transform (CWT) matrix. Each row corresponds to a different scale.
        widths (ndarray): Array of wavelet widths (scales) used for the transform.
        trace_number (int): The index of the seismic trace being analyzed.

    Returns:
        None: The wavelet transform is displayed as an image on the provided Matplotlib axis.
    """

    ax.imshow(np.abs(cwt_matrix), aspect='auto', extent=[0, len(cwt_matrix[0]), min(widths), max(widths)])
    ax.set_title(f"Seismic Trace {trace_number} Wavelet Transform")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Scale")

def plot_spectrogram(ax, f, t, Sxx, trace_number):
    
    """
    Plot the spectrogram of a seismic trace.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot the spectrogram.
        f (ndarray): Array of sample frequencies.
        t (ndarray): Array of time segments.
        Sxx (ndarray): Spectrogram matrix, representing the power at each time-frequency point.
        trace_number (int): The index of the seismic trace being analyzed.

    Returns:
        None: The spectrogram is displayed as a color map on the provided Matplotlib axis.
    """

    ax.pcolormesh(t, f, 10 * np.log10(Sxx))
    ax.set_title(f"Seismic Trace {trace_number} Spectrogram")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")

def plot_seismic_image(ax, seismic_data, twt, n_traces):

    """
    Display a seismic section as an image.

    This function takes seismic data in the form of a 2D array (traces x samples) and displays it 
    as an image where each pixel represents a seismic amplitude.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to display the seismic image.
        seismic_data (ndarray): 2D array representing seismic traces (rows) and samples (columns).

    Returns:
        None: The seismic image is displayed on the provided Matplotlib axis.
    """

    ax.imshow(seismic_data, cmap='seismic', aspect='auto', extent=[0, n_traces, twt[-1], twt[0]])
    ax.set_xlabel('Trace Number')
    ax.set_ylabel('Two Way Time (ms)')