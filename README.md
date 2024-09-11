# PySeisPro
A small python GUI application for processing and interpreting seismic data from SEG-Y files

Overview
This application is a GUI tool for processing and interpreting seismic data from SEG-Y files. Built with PyQt5 and leveraging segyio for seismic data handling, the tool provides several features such as applying filters, plotting seismic images, and performing trace analysis.

Features
SEG-Y File Handling: Load and visualize SEG-Y seismic data.
Filtering: Apply FIR and IIR filters to the seismic data.
Seismic Interpretation: Interactive window for seismic data interpretation.
Trace Analysis: Periodogram, Welch periodogram, wavelet transform, and spectrogram analysis on seismic traces.
Plotting: Display seismic images, trace plots, spectrograms, and more using Matplotlib.

Installation
Clone the repository:

bash
Copy code
git clone https://github.com/your-repo-url
cd seismic-gui
Install dependencies:

bash
Copy code
pip install -r requirements.txt
Make sure the following Python libraries are installed:

PyQt5
segyio
Matplotlib
Numpy
Custom modules: Data_handling, Filters, Mute, Gains, Plots, Trace_analysis, Deconvolution, Interpretation_new
Usage
Run the application:

bash
Copy code
python main_gui_app_segyio.py
In the GUI:

Use the File Menu to import a SEG-Y file.
Navigate through various processing options, including applying filters and plotting seismic data.
Open the Seismic Interpretation window for more detailed analysis and interpretation of seismic lines.
Key Functionalities
SEG-Y File Loading
The application uses segyio to load SEG-Y files. You can load 2D/3D seismic data and visualize it using Matplotlib.
Filters
IIR and FIR Filters: The tool includes FIR and IIR filters, such as the FK filter, that can be applied to seismic data.
Plotting Options
Seismic Image: Plot seismic sections in various styles.
Spectrogram, Wavelet Transform: Visualize seismic traces with time-frequency analysis methods.
Periodograms: Generate frequency-domain representations of seismic traces.
Seismic Interpretation
The Seismic Interpretation Window allows users to interactively interpret seismic lines, providing a high level of flexibility in data analysis.
Custom Modules
The following custom modules extend the functionality of the application:

Data_handling: Manages SEG-Y file loading and basic processing.
Filters: Contains filter functions like FIR, IIR, and FK filtering.
Mute: Implements muting functions for specific sections of seismic traces.
Gains: Provides functions for automatic gain control (AGC), time-varying gain (TVG), and constant gain application.
Plots: Manages all plotting functions, including seismic images, periodograms, and spectrograms.
Trace_analysis: Performs analysis like periodograms and wavelet transforms on individual seismic traces.
Deconvolution: Implements deconvolution techniques for seismic processing.
Interpretation_new: Manages the seismic interpretation window where users can analyze the seismic data interactively.
Future Enhancements
Basemap Integration: Add functionality to plot seismic lines on a geographic map using Cartopy or other geospatial tools.
Improved Filter Set: Enhance the available filters with additional seismic-specific options.
