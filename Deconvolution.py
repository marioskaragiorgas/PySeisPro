"""
deconvolution.py

Deconvolution techniques are essential in seismic data processing for enhancing the temporal resolution of seismic signals. 
The primary goal of deconvolution is to compress the seismic wavelet into a spike, thereby retrieving the earth's reflectivity series. 
This script implements several key deconvolution methods, including spiking deconvolution, predictive deconvolution, 
and Wiener deconvolution. These methods help remove the effects of the seismic wavelet and attenuate multiples and reverberations, 
yielding a clearer image of subsurface structures.

Key concepts:
- **Convolutional Model:** The recorded seismic trace is a convolution of the earth’s impulse response with the seismic wavelet.
- **Wavelet Compression:** Deconvolution aims to compress the wavelet, ideally leaving only the earth's reflectivity.
- **Minimum Phase Assumption:** Deconvolution typically assumes that the seismic wavelet is minimum phase.
"""

import numpy as np
from scipy.signal import lfilter, correlate, wiener

class Deconvolution:

    @staticmethod
    def spiking_deconvolution(trace, wavelet, noise_level=0.001):
        
        """
        Spiking Deconvolution (also known as least-squares inverse filtering) aims to compress the seismic wavelet 
        into a spike, effectively enhancing the resolution of the seismic data. This method assumes that the seismic 
        wavelet is minimum phase and that the trace can be represented as a convolution of this wavelet with the 
        earth's reflectivity series.

        Parameters:
        - trace: 1D numpy array, the seismic trace to be deconvolved.
        - wavelet: 1D numpy array, the estimated seismic wavelet.
        - noise_level: A small constant added to stabilize the inverse filter (default=0.001).

        Returns:
        - deconvolved_trace: 1D numpy array, the trace after spiking deconvolution.
        """

        autocorr = correlate(wavelet, wavelet, mode='full')
        mid_point = len(autocorr) // 2
        autocorr = autocorr[mid_point:]

        # Adding a small noise level to stabilize the inverse filter
        autocorr[0] += noise_level
        
        inverse_filter = np.linalg.inv(np.toeplitz(autocorr)).dot(np.eye(len(autocorr))[:, 0])
        deconvolved_trace = lfilter(inverse_filter, [1.0], trace)
        
        return deconvolved_trace
    
    @staticmethod
    def predictive_deconvolution(trace, prediction_distance, filter_length, noise_level=0.001):
        
        """
        Predictive Deconvolution uses a prediction error filter to remove periodic components (such as multiples) 
        from the seismic trace. The technique aims to predict the primary reflections by filtering out predictable 
        (repeated) components, thus enhancing the primary signal.

        Parameters:
        - trace: 1D numpy array, the seismic trace to be deconvolved.
        - prediction_distance: The lag distance for prediction.
        - filter_length: Length of the prediction error filter.
        - noise_level: A small constant added to stabilize the inverse filter (default=0.001).

        Returns:
        - deconvolved_trace: 1D numpy array, the trace after predictive deconvolution.
        """

        autocorr = correlate(trace, trace, mode='full')
        mid_point = len(autocorr) // 2
        autocorr = autocorr[mid_point:]

        # Creating the prediction error filter
        R = np.toeplitz(autocorr[:filter_length])
        R[:, 0] += noise_level
        p = autocorr[prediction_distance:prediction_distance + filter_length]
        prediction_error_filter = np.linalg.solve(R, p)
        
        # Deconvolution using the prediction error filter
        deconvolved_trace = lfilter(prediction_error_filter, [1.0], trace)
        
        return deconvolved_trace
    
    @staticmethod
    def wiener_deconvolution(seismic_data, window_size, noise_power):
        
        """
        Wiener Deconvolution aims to minimize the mean square error between the desired output and the actual output.
        This filter is optimal in the least-squares sense and can be designed to convert the seismic wavelet into any desired shape,
        typically a spike. Unlike spiking deconvolution, Wiener deconvolution can balance between wavelet compression and noise attenuation.

        Parameters:
        - seismic data: N-dimensional numpy array, the seismic data to be deconvolved.
        - window_size: int or array_like, optional. A scalar or an N-length
            list giving the size of the Wiener filter window in each dimension. 
            Elements of mysize should be odd. If mysize is a scalar, then this 
            scalar is used as the size in each dimension.
        - noise_power: float, optional. The noise-power to use. If None, 
            then noise is estimated as the average of the local variance of the input.

        Returns:
        - deconvolved_data: Wiener filtered result with the same shape as the imput data.
        """
        
        deconvolved_data = np.zeros_like(seismic_data)
        for i, trace in enumerate(seismic_data):
            deconvolved_data[i] = wiener(trace, mysize=window_size, noise = noise_power)
        return deconvolved_data