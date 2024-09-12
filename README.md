# PySeisPro
A small python GUI application for processing and interpreting seismic data from SEG-Y files.

## Overview
This application is a GUI tool for processing and interpreting seismic data from SEG-Y files. Built with PyQt5 and leveraging segyio for seismic data handling, the tool provides several features such as applying filters, plotting seismic images, and performing trace analysis.

## Features
- **SEG-Y File Handling**: Load and visualize SEG-Y seismic data.
- **Data Filtering**: Apply FIR and IIR filters to the seismic data.
- **Seismic Interpretation**: Interactive window for seismic data interpretation.
- **Seismic Trace Analysis**: Periodogram, Welch periodogram, wavelet transform, and spectrogram analysis on seismic traces.
- **Data vizualization**: Display seismic images, trace plots, spectrograms, and more using Matplotlib.

## Installation

1. Clone the repository:

```bash
git clone https://github.com/marioskaragiorgas/PySeisPro.git
cd PySeisPro
```

2. Install dependencies:

```bash
pip install -r requirements.txt
```

PySeisPro depends on the following Python libraries:

1.PyQt5

2.Segyio

3.Matplotlib

4.Numpy

5.Scipy

6.Scikit-image

7.Pandas

8.Obspy

9.PyWavelets

**Custom modules**: Data_handling, Filters, Mute, Gains, Plots, Trace_analysis, Deconvolution, Seismic_Interpretation

## Usage

Run the application:

In the PySeisPro main directory open a terminal and type:

``` bash
python PySeisPro.py
```
or 

``` bash
python3 PySeisPro.py
```

## In the PySeisPro GUI:

1. Use the File Menu to import a SEG-Y file. 
2. Navigate through various processing options, including applying filters and plotting seismic data.
3. Open the Seismic Interpretation window for more detailed analysis and interpretation of seismic images.

## Key Functionalities

**SEG-Y File Loading**

The application uses segyio to load SEG-Y files. Generally segyio can load 2D/3D seismic data and visualize it using Matplotlib. At this phase the user can load only 2D data.

**Filters**

IIR and FIR Filters: The tool includes FIR and IIR filters, such as the FK filter, that can be applied to seismic data.

**Plotting Options**

- Seismic Image: Plot seismic sections in various styles.
- Spectrogram, Wavelet Transform: Visualize seismic traces with time-frequency analysis methods.
- Periodograms: Generate frequency-domain representations of seismic traces.

**Data Muting Techniques**

This tool Implements muting functions for specific sections of seismic traces.

**Gains**

This tool includes functions for automatic gain control (AGC), time-varying gain (TVG), and constant gain, that can be applied to seismic data.

**Seismic Deconvolution**

This tool implements deconvolution techniques for seismic processing.

**Seismic Interpretation**

The Seismic Interpretation Window allows users to interactively interpret seismic lines, providing a high level of flexibility in data analysis.

## Custom Modules
The following custom modules extend the functionality of the application:

- **Data_handling.py**: Manages SEG-Y file loading and basic processing.

- **Filters.py**: Contains filter functions like FIR, IIR, and FK filtering.

- **Mute.py**: Implements muting functions for specific sections of seismic traces.

- **Gains.py**: Provides functions for automatic gain control (AGC), time-varying gain (TVG), and constant gain application.

- **Plots.py**: Manages all plotting functions, including seismic images, periodograms, and spectrograms.

- **Trace_analysis.py**: Performs analysis like periodograms and wavelet transforms on individual seismic traces.

- **Deconvolution.py**: Implements deconvolution techniques for seismic processing.

- **Seismic_Interpretation**: Manages the seismic interpretation window where users can analyze the seismic data interactively.

## Future Enhancements
- **File Format extension**: Extent the file formats for seismic data import (Seismic Unix SU, SAC) and export (as images not only SEG-Y/SGY).
  
- **Basemap Integration**: Add functionality to plot seismic lines on a geographic map using Cartopy or other geospatial tools.

- **Improved Filter Set and Deconvolution techniques**: Enhance the available filters with additional seismic-specific options.

- **Deep Learning techniques for automatic horizon and fault exctraction**: Use of pretrained Convolutional Neural Networks (CNNs) specialized on horizon and fault exctraction.
