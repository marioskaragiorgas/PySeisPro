# main_gui_app_segyio.py

import sys
import os
import numpy as np
import segyio
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QLabel, QVBoxLayout, QWidget, 
    QAction, QSizePolicy, QInputDialog, QMessageBox, QToolBar
)
from PyQt5.QtCore import pyqtSlot, QSize, QThread, pyqtSignal
from PyQt5.QtGui import QIcon
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib import pyplot as plt
from Data_handling import SEGYHandler
from Filters import IIR_Filters, FIR_Filters
from Mute import Mute, PredefinedMute
from Gains import agc_gain, tvg_gain, constant_gain
from Plots import plot_trace, plot_periodogram, plot_welch_periodogram, plot_wavelet_transform, plot_spectrogram, plot_seismic_image
from Trace_analysis import trace_periodogram, trace_welch_periodogram, trace_wavelet_transform, trace_spectrogram
from Deconvolution import Deconvolution
from Seismic_Interpretation import SeismicInterpretationWindow

import logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s %(message)s'
)

class SeismicFilterWorker(QThread):
    finished = pyqtSignal(object)
    error = pyqtSignal(str)

    def __init__(self, filter_func, data, *args):
        super().__init__()
        self.filter_func = filter_func
        self.data = data
        self.args = args

    def run(self):
        try:
            if self.filter_func == FIR_Filters.fk_filter:
                # Apply the F-K filter to the entire data (2D array)
                result = self.filter_func(self.data, *self.args)
            else:
                # Apply other filters trace by trace    
                result = np.array([self.filter_func(trace, *self.args) for trace in self.data])
            self.finished.emit(result)
        except Exception as e:
            logging.error(f"Filter application failed: {e}")
            self.error.emit(str(e))

class SeismicDataApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.init_ui()
        self.create_menus()
        self.create_toolbar()
        self.segy_handler = SEGYHandler()
        self.segy_file = None
        self.spec = None
        self.data = None  # Raw seismic data
        self.n_traces = None
        self.twt = None
        self.processed_data = None  # Processed seismic data
        self.sample_interval = None 
        self.sample_rate = None
        self.interpretation_window = None


    def init_ui(self):
        """Initialize the main UI components."""
        self.setWindowTitle("PySeisPro")

        # Create central widget and layout
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)

        # Create GUI components
        self.data_info_label = QLabel("")
        self.layout.addWidget(self.data_info_label)

        # Data name label with text wrapping
        self.data_name_label = QLabel("No Data Loaded")
        self.layout.addWidget(self.data_name_label)

        # Matplotlib figure for displaying plots
        self.figure, (self.ax_raw, self.ax_top) = plt.subplots(2, 1, figsize=(15, 20))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.layout.addWidget(self.canvas)

        # Add Matplotlib Navigation Toolbar for zooming and panning
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout.addWidget(self.toolbar)

    def create_menus(self):
        """Create application menus."""
        self.menu_bar = self.menuBar()

        # File menu
        file_menu = self.menu_bar.addMenu("File")
        file_submenu = file_menu.addMenu("Import...")
        file_submenu.addAction(self.create_action("Import Data", self.import_data))
        file_menu.addAction(self.create_action("Export as...", self.save_processed_segy))
        file_menu.addAction(self.create_action("Close", self.close_segy))

        # Trace Analysis menu
        self.trace_menu = self.menu_bar.addMenu("Trace Analysis")
        self.add_trace_menu_action("Trace Plot", self.plot_trace)
        self.add_trace_menu_action("Trace Periodogram", self.plot_trace_periodogram)
        self.add_trace_menu_action("Trace Welch Periodogram", self.plot_trace_periodogram_welch)
        self.add_trace_menu_action("Trace Wavelet Transform", self.plot_trace_wavelet)
        self.add_trace_menu_action("Trace Spectrogram", self.plot_trace_spectrogram)

        # Filter menu
        filter_menu = self.menu_bar.addMenu("Filters")
        
        IIR_filter_submenu = filter_menu.addMenu("IIR")
        
        butterworth_menu = IIR_filter_submenu.addMenu("Butterworth")
        butterworth_menu.addAction(self.create_action("Bandpass Filter", self.apply_butter_bandpass_filter))
        butterworth_menu.addAction(self.create_action("Highpass Filter", self.apply_butter_highpass_filter))
        butterworth_menu.addAction(self.create_action("Lowpass Filter", self.apply_butter_lowpass_filter))
        
        cheby_menu = IIR_filter_submenu.addMenu("Chebyshev")
        cheby_menu.addAction(self.create_action("Bandpass Filter", self.apply_cheby_bandpass_filter))
        cheby_menu.addAction(self.create_action("Highpass Filter", self.apply_cheby_highpass_filter))
        cheby_menu.addAction(self.create_action("Lowpass Filter", self.apply_cheby_lowpass_filter))

        FIR_filter_submenu = filter_menu.addMenu("FIR")

        FIR_filter_submenu.addAction(self.create_action("Bandpass Filter", self.apply_fir_bandpass_filter))
        FIR_filter_submenu.addAction(self.create_action("Highpass Filter", self.apply_fir_highpass_filter))
        FIR_filter_submenu.addAction(self.create_action("Lowpass Filter", self.apply_fir_lowpass_filter))
        FIR_filter_submenu.addAction(self.create_action("F-K Filter", self.apply_fk_filter))
        FIR_filter_submenu.addAction(self.create_action("Zero Phase Filter", self.apply_zero_phase_filter))
        FIR_filter_submenu.addAction(self.create_action("Wavelet Filter", self.apply_wavelet_filter))

        # Gain menu
        gain_menu = self.menu_bar.addMenu("Gains")
        gain_menu.addAction(self.create_action("AGC Gain", self.agc_gain))
        gain_menu.addAction(self.create_action("TVG Gain", self.tvg_gain))
        gain_menu.addAction(self.create_action("Constant Gain", self.const_gain))

        # Mute menu
        mute_menu = self.menu_bar.addMenu("Mute")
        mute_functions_submenu = mute_menu.addMenu("Mutting Functions")
        mute_functions_submenu.addAction(self.create_action("Top Mute", self.apply_top_mute))
        mute_functions_submenu.addAction(self.create_action("Bottom Mute", self.apply_bottom_mute))
        #mute_functions_submenu.addAction(self.create_action("Offset Mute", self.apply_offset_mute))
        mute_functions_submenu.addAction(self.create_action("Time Variant Mute", self.apply_time_variant_mute))
        mute_functions_submenu.addAction(self.create_action("Interactive Mute", self.apply_interactive_mute))
        mute_predifined_functions_submenu = mute_menu.addMenu("Predifined Mutting Functions")
        mute_predifined_functions_submenu.addAction(self.create_action("Shallow Zone Mute ", self.apply_SZ_mute))
        mute_predifined_functions_submenu.addAction(self.create_action("Deep Zone Mute ", self.apply_DZ_mute))
        mute_predifined_functions_submenu.addAction(self.create_action("Direct Wave Mute ", self.apply_DW_mute))

        # Deconvolution menu 
        deconvolution_menu = self.menu_bar.addMenu("Deconvolution")
        #deconvolution_menu.addAction(self.create_action("Spiking", self.apply_spiking_dec))
        #deconvolution_menu.addAction(self.create_action("Predictive", self.apply_predictive_dec))
        deconvolution_menu.addAction(self.create_action("Wiener", self.apply_wiener_dec))

        # Interpretation menu
        interpretation_menu = self.menu_bar.addMenu("Interpretation")
        interpretation_menu.addAction(self.create_action("Horizon Picking", self.apply_Horizon_pick))

    def create_action(self, name, method):
        """Helper function to create a QAction."""
        action = QAction(name, self)
        action.triggered.connect(method)
        return action

    def add_trace_menu_action(self, name, method):
        """Add an action to the Trace Analysis menu."""
        action = self.create_action(name, method)
        self.trace_menu.addAction(action)

    def create_toolbar(self):
        """Create application toolbar."""
        self.toolbar = QToolBar()
        self.toolbar.setIconSize(QSize(16, 16))
        self.toolbar.setFixedHeight(36)

        # Specify the path to your icon file
        open_file_icon_path = os.path.join(os.path.dirname(__file__), "Images\mIconFolderOpen.png")
        
        # Create the QIcon object
        open_file_icon = QIcon(open_file_icon_path)

        # Create the QAction with the icon and connect it to the function
        import_file_action = QAction(open_file_icon, "Import file", self)
        import_file_action.triggered.connect(self.import_data)

        # Add the action to the toolbar
        self.toolbar.addAction(import_file_action)
        self.addToolBar(self.toolbar)

    def show_error(self, title, message):
        QMessageBox.critical(self, title, message)    

    @pyqtSlot()
    def import_data(self):
        logging.info("Starting load_data_segyio function")
        """Load seismic data from file using Segyio."""
        file_path, _ = QFileDialog.getOpenFileName(self, "Open Data File", "", "SEG-Y Files (*.sgy *.segy)")
        if file_path:
            try:
                self.segy_file, self.spec, self.data, self.n_traces, n_samples, self.twt, data_format, self.sample_interval, self.sample_rate = self.segy_handler.load_data_segyio(file_path)
                self.data_info_label.setText(f"Loaded data from {file_path} with Segyio")
                self.plot_raw_seismic_image()
            except Exception as e:
                self.data_info_label.setText(f"Error loading data: {e}")
        logging.info("Ending load_data_segyio function")

    @pyqtSlot()
    def save_processed_segy(self):
        """Save processed data as SEG-Y file."""
        if self.processed_data is None:
            QMessageBox.warning(self, "Warning", "No processed data available. Saving raw data instead.")
            data_to_save = self.data
        else:
            data_to_save = self.processed_data
        
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Data File", "", "SEG-Y Files (*.sgy *.segy)")
        
        if file_path:
            try:
                self.segy_handler.save_segy_file(file_path, self.spec, data_to_save)
                self.data_info_label.setText(f"Saved data to {file_path}")
            except Exception as e:
                self.data_info_label.setText(f"Error saving data: {e}")

    @pyqtSlot()
    def close_segy(self):
        """Close the currently open SEG-Y file and clear any displayed data."""
        try:
            self.segy_handler.close_file()
            self.data_info_label.setText("Files closed successfully.")

            # Clear any processed or raw data
            self.data = None
            self.processed_data = None

            # Clear the plots
            self.ax_top.clear()
            self.ax_raw.clear()
            self.canvas.draw()

        except Exception as e:
            self.data_info_label.setText(f"Error closing file: {e}")


    def plot_raw_seismic_image(self):
        """Plot the seismic image from the raw data."""
        if self.data is None:
            QMessageBox.critical(self, "Error", "No data loaded.")
            return
        self.ax_top.clear()
        plot_seismic_image(self.ax_top, self.data.T, self.twt, self.n_traces)
        self.canvas.draw()

    def plot_processed_seismic_image(self):
        """Plot the seismic image from the processed data."""
        if self.processed_data is None:
            QMessageBox.critical(self, "Error", "No processed data available.")
            return
        self.ax_top.clear()
        plot_seismic_image(self.ax_top, self.processed_data.T, self.twt, self.n_traces)
        self.canvas.draw()

    def plot_trace(self):
        """Plot the seismic trace."""
        if self.data is None:
            QMessageBox.critical(self, "Error", "No data loaded.")
            return
        
        trace_number, ok = QInputDialog.getInt(self, "Enter Trace Number", "Trace Number:", 0, 0, self.data.shape[0] - 1, 1)
        
        if ok:
            trace = self.processed_data[trace_number] if self.processed_data is not None else self.data[trace_number]
            self.ax_raw.clear()
            plot_trace(self.ax_raw, trace, trace_number, delta=self.sample_interval)
            self.canvas.draw()

    def plot_trace_periodogram(self):
        """Plot the periodogram of the seismic trace."""
        if self.data is None:
            QMessageBox.critical(self, "Error", "No data loaded.")
            return

        trace_number, ok = QInputDialog.getInt(self, "Enter Trace Number", "Trace Number:", 0, 0, self.data.shape[0] - 1, 1)
        if ok:
            trace = self.processed_data[trace_number] if self.processed_data is not None else self.data[trace_number]
            f, Pxx = trace_periodogram(trace, fs=self.sample_rate)
            self.ax_raw.clear()
            plot_periodogram(self.ax_raw, f, Pxx, trace_number)
            self.canvas.draw()

    def plot_trace_periodogram_welch(self):
        """Plot the Welch periodogram of the seismic trace."""
        if self.data is None:
            QMessageBox.critical(self, "Error", "No data loaded.")
            return

        trace_number, ok = QInputDialog.getInt(self, "Enter Trace Number", "Trace Number:", 0, 0, self.data.shape[0] - 1, 1)
        if ok:
            trace = self.processed_data[trace_number] if self.processed_data is not None else self.data[trace_number]
            fs = 1.0  # Assuming sampling rate of 1.0 Hz, adjust as needed
            f, Pxx = trace_welch_periodogram(trace, fs=self.sample_rate)
            self.ax_raw.clear()
            plot_welch_periodogram(self.ax_raw, f, Pxx, trace_number)
            self.canvas.draw()

    def plot_trace_wavelet(self):
        """Plot the wavelet transform of the seismic trace."""
        if self.data is None:
            QMessageBox.critical(self, "Error", "No data loaded.")
            return

        trace_number, ok = QInputDialog.getInt(self, "Enter Trace Number", "Trace Number:", 0, 0, self.data.shape[0] - 1, 1)
        if ok:
            trace = self.processed_data[trace_number] if self.processed_data is not None else self.data[trace_number]
            widths = np.arange(1, 128)
            cwt_matrix = trace_wavelet_transform(trace, widths)
            self.ax_raw.clear()
            plot_wavelet_transform(self.ax_raw, cwt_matrix, widths, trace_number)
            self.canvas.draw()

    def plot_trace_spectrogram(self):
        """Plot the spectrogram of the seismic trace."""
        if self.data is None:
            QMessageBox.critical(self, "Error", "No data loaded.")
            return

        trace_number, ok = QInputDialog.getInt(self, "Enter Trace Number", "Trace Number:", 0, 0, self.data.shape[0] - 1, 1)
        if ok:
            trace = self.processed_data[trace_number] if self.processed_data is not None else self.data[trace_number]
            fs = 1.0  # Assuming sampling rate of 1.0 Hz, adjust as needed
            f, t, Sxx = trace_spectrogram(trace, fs=self.sample_rate)
            self.ax_raw.clear()
            plot_spectrogram(self.ax_raw, f, t, Sxx, trace_number)
            self.canvas.draw()

    def get_filter_params(self, bandpass=False, filter_type='IIR'):
        order, ok = QInputDialog.getInt(self, "Enter Filter Order", "Filter Order:", 4, 1, 100, 1)
        if not ok:
            raise ValueError("No filter order provided")
        freq, ok = QInputDialog.getDouble(self, "Enter Critical Frequency (Hz)", "Critical Frequency (Hz):", 10, 0, 10000, 1)
        if not ok:
            raise ValueError("No critical frequency provided")
        if bandpass:
            freqmax, ok = QInputDialog.getDouble(self, "Enter Maximum Critical Frequency (Hz)", "Maximum Critical Frequency (Hz):", 100, 0, 10000, 1)
            if not ok:
                raise ValueError("No maximum critical frequency provided")
            self.validate_filter_params(order, freq, self.sample_rate, filter_type, bandpass=True, freqmax=freqmax)
            return order, freq, freqmax
        self.validate_filter_params(order, freq, self.sample_rate, filter_type)
        return order, freq

    def validate_filter_params(self, order, freq, sample_rate, filter_type, bandpass=False, freqmax=None):
        nyquist = sample_rate / 2
        if order < 1 or order > 100:
            raise ValueError("Filter order must be between 1 and 100.")
        if freq <= 0 or freq >= nyquist:
            raise ValueError(f"{filter_type.capitalize()} filter frequency must be between 0 and Nyquist frequency {nyquist}.")
        if bandpass and (freqmax is None or freqmax <= freq or freqmax >= nyquist):
            raise ValueError(f"Bandpass filter max frequency must be between {freq} and Nyquist frequency {nyquist}.")

    def apply_filter(self, filter_func, *args):
        """Apply the filter asynchronously to keep the GUI responsive."""
        data_to_filter = self.processed_data if self.processed_data is not None else self.data
        if data_to_filter is None:
            self.show_error("Error", "No data loaded.")
            return
        else:
            logging.info(f"Data shape before filtering: {data_to_filter.shape}")
            self.worker = SeismicFilterWorker(filter_func, data_to_filter, *args)
            self.worker.finished.connect(self.on_filter_finished)
            self.worker.error.connect(self.show_error)
            self.worker.start()

    def on_filter_finished(self, result):
        self.processed_data = result
        self.plot_processed_seismic_image()
        self.data_info_label.setText("Filter applied and processed data updated.")

    def show_error(self, title, message):
        QMessageBox.critical(self, title, message)

    @pyqtSlot()
    def apply_butter_bandpass_filter(self):
        try:
            order, freqmin, freqmax = self.get_filter_params(bandpass=True, filter_type='Butterworth')
            self.apply_filter(IIR_Filters.bandpass_filter, self.sample_rate, order, freqmin, freqmax)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply bandpass filter: {e}")

    @pyqtSlot()
    def apply_butter_highpass_filter(self):
        try:
            order, freq = self.get_filter_params(filter_type='Butterworth')
            self.apply_filter(IIR_Filters.highpass_filter, self.sample_rate, order, freq)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply highpass filter: {e}")

    @pyqtSlot()
    def apply_butter_lowpass_filter(self):
        try:
            order, freq = self.get_filter_params(filter_type='Butterworth')
            self.apply_filter(IIR_Filters.lowpass_filter, self.sample_rate, order, freq)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply lowpass filter: {e}")

    @pyqtSlot()
    def apply_cheby_bandpass_filter(self):
        try:
            order, freqmin, freqmax = self.get_filter_params(bandpass=True, filter_type='Chebyshev')
            ripple = self.get_chebyshev_ripple()
            self.apply_filter(IIR_Filters.cheby2_bandpass_filter, self.sample_rate, order, freqmin, freqmax, ripple)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply Chebyshev bandpass filter: {e}")

    @pyqtSlot()
    def apply_cheby_highpass_filter(self):
        try:
            order, freq = self.get_filter_params(filter_type='Chebyshev')
            ripple = self.get_chebyshev_ripple()
            self.apply_filter(IIR_Filters.cheby2_highpass_filter, self.sample_rate, order, freq, ripple)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply Chebyshev highpass filter: {e}")

    @pyqtSlot()
    def apply_cheby_lowpass_filter(self):
        try:
            order, freq = self.get_filter_params(filter_type='Chebyshev')
            ripple = self.get_chebyshev_ripple()
            self.apply_filter(IIR_Filters.cheby2_lowpass_filter, self.sample_rate, order, freq, ripple)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply Chebyshev lowpass filter: {e}")

    def get_chebyshev_ripple(self):
        ripple, ok = QInputDialog.getDouble(self, "Enter Ripple (dB)", "Ripple (dB):", 0.5, 0.01, 10.0, 2)
        if not ok:
            raise ValueError("No ripple provided")
        return ripple

    @pyqtSlot()
    def apply_fir_bandpass_filter(self):
        try:
            order, freqmin, freqmax = self.get_filter_params(bandpass=True, filter_type='FIR')
            window = self.get_fir_window()
            self.apply_filter(FIR_Filters.bandpass_filter, freqmin, freqmax, self.sample_rate, order, window)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply FIR bandpass filter: {e}")

    @pyqtSlot()
    def apply_fir_highpass_filter(self):
        try:
            order, freq = self.get_filter_params(filter_type='FIR')
            window = self.get_fir_window()
            self.apply_filter(FIR_Filters.highpass_filter, freq, self.sample_rate, order, window)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply FIR highpass filter: {e}")

    @pyqtSlot()
    def apply_fir_lowpass_filter(self):
        try:
            order, freq = self.get_filter_params(filter_type='FIR')
            window = self.get_fir_window()
            self.apply_filter(FIR_Filters.lowpass_filter, freq, self.sample_rate, order, window)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply FIR lowpass filter: {e}")

    @pyqtSlot()
    def apply_fk_filter(self):
        try:
            order, ok = QInputDialog.getInt(self, "Enter Filter Order", "Filter Order:", 10, 1, 100, 1)
            if not ok:
                raise ValueError("No filter order provided")
            self.apply_filter(FIR_Filters.fk_filter, self.sample_rate, order)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply F-K filter: {e}")

    @pyqtSlot()
    def apply_zero_phase_filter(self):
        try:
            order, freqmin, freqmax = self.get_filter_params(bandpass=True, filter_type='FIR')
            self.apply_filter(FIR_Filters.zero_phase_bandpass_filter, freqmin, freqmax, self.sample_rate, order)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply zero-phase bandpass filter: {e}")

    @pyqtSlot()
    def apply_wavelet_filter(self):
        try:
            wavelet_type, ok = QInputDialog.getText(self, "Enter Wavelet Type", "Wavelet Type:", text="db4")
            if not ok:
                raise ValueError("No wavelet type provided")
            level, ok = QInputDialog.getInt(self, "Enter Decomposition Level", "Decomposition Level:", 5, 1, 10, 1)
            if not ok:
                raise ValueError("No decomposition level provided")
            self.apply_filter(FIR_Filters.wavelet_filter, wavelet_type, level)
        except ValueError as e:
            self.show_error("Error", f"Unable to apply wavelet filter: {e}")

    def get_fir_window(self):
        window, ok = QInputDialog.getText(self, "Enter Window Type", "Window Type:", text="hamming")
        if not ok:
            raise ValueError("No window type provided")
        return window

    def agc_gain(self):
        param1, ok = QInputDialog.getInt(self, "Enter AGC Window Size", "Window Size:", 51, 1, 1000, 1)
        if ok:
            self.apply_gain('agc', param1)

    def tvg_gain(self):
        param1, ok = QInputDialog.getDouble(self, "Enter TVG time gradient", "Time gradient:", 2.0, 0.1, 10.0, 1)
        if ok:
            self.apply_gain('tvg', param1)

    def const_gain(self):
        param1, ok = QInputDialog.getDouble(self, "Enter Constant Gain Factor", "Gain Factor:", 2.0, 0.1, 10.0, 1)
        if ok:
            self.apply_gain('const', param1)

    def apply_gain(self, gain_type, param1=None):
        if self.data is None:
            QMessageBox.critical(self, "Error", "No data loaded.")
            return
        try:
            data_to_gain = self.processed_data if self.processed_data is not None else self.data
            if gain_type == 'agc':
                self.processed_data = agc_gain(data_to_gain, param1)
            elif gain_type == 'tvg':
                self.processed_data = tvg_gain(data_to_gain, param1)
            elif gain_type == 'const':
                self.processed_data = constant_gain(data_to_gain, param1)
            
            self.plot_processed_seismic_image()
            self.data_info_label.setText(f"Applied {gain_type.upper()} gain and updated processed data.")
        except ValueError as e:
            self.show_error("Error", f"Unable to apply {gain_type} gain: {e}")

    def apply_top_mute(self):
        """Slot to apply top mute using a user-specified mute time."""
        mute_time, ok = QInputDialog.getDouble(self, "Enter Mute Time", "Mute Time (seconds):", 0.1, 0, 10, 2)
        if ok:
            self.processed_data = Mute.top_mute(self.processed_data if self.processed_data is not None else self.data, mute_time, self.sample_interval)
            self.plot_processed_seismic_image()
            self.data_info_label.setText("Applied top mute and updated processed data.")
    
    def apply_bottom_mute(self):
        """Slot to apply bottom mute using a user-specified mute time."""
        mute_time, ok = QInputDialog.getDouble(self, "Enter Mute Time", "Mute Time (seconds):", 0.1, 0, 10, 2)
        if ok:
            self.processed_data = Mute.top_mute(self.processed_data if self.processed_data is not None else self.data, mute_time, self.sample_interval)
            self.plot_processed_seismic_image()
            self.data_info_label.setText("Applied bottom mute and updated processed data.")
    
    """
    def apply_offset_mute(self):
        #Slot to apply offset mute.

         if self.processed_data is not None:
            self.processed_data = offset_mute(self.processed_data, offsets, mute_offset)
       
         else:
            self.processed_data = offset_mute(self.data, offsets, mute_offset)

        self.plot_processed_seismic_image()
        self.data_info_label.setText("Applied offset mute and updated processed data.")   
    """
    def apply_time_variant_mute(self):
        """Slot to apply time-variant mute."""
        initial_time, ok1 = QInputDialog.getDouble(self, "Enter Initial Mute Time", "Initial Mute Time (seconds):", 0.1, 0, 10, 2)
        if not ok1: return
        final_time, ok2 = QInputDialog.getDouble(self, "Enter Final Mute Time", "Final Mute Time (seconds):", 1.0, 0, 10, 2)
        if not ok2: return
        
        self.processed_data = Mute.time_variant_mute(self.processed_data if self.processed_data is not None else self.data, initial_time, final_time, self.sample_interval)
        self.plot_processed_seismic_image()
        self.data_info_label.setText("Applied time-variant mute and updated processed data.")
    
    def apply_SZ_mute(self):
        """Slot to apply shallow mute."""
        self.processed_data = PredefinedMute.shallow_zone_mute(self.processed_data if self.processed_data is not None else self.data, self.sample_interval)
        self.plot_processed_seismic_image()
        self.data_info_label.setText("Applied shallow zone mute and updated processed data.")

    def apply_DW_mute(self):
        """Slot to apply shallow mute."""
        self.processed_data = PredefinedMute.marine_direct_wave_mute(self.processed_data if self.processed_data is not None else self.data, self.sample_interval)
        self.plot_processed_seismic_image()
        self.data_info_label.setText("Applied shallow zone mute and updated processed data.") 
    
    def apply_DZ_mute(self):
        """Slot to apply shallow mute."""
        self.processed_data = PredefinedMute.deep_zone_mute(self.processed_data if self.processed_data is not None else self.data, self.sample_interval)
        self.plot_processed_seismic_image()
        self.data_info_label.setText("Applied shallow zone mute and updated processed data.")       

    pyqtSlot()
    def apply_interactive_mute(self):
        """Slot to enable interactive mute."""
        self.processed_data = Mute.interactive_mute(self.ax_raw, self.data)
        self.plot_processed_seismic_image()
        self.data_info_label.setText("Applied interactive mute based on the drawn polygon.")        
    
    pyqtSlot()
    def apply_wiener_dec(self):
        """ Slot to enable wiener deconvolution"""
        self.processed_data = Deconvolution.wiener_deconvolution(self.processed_data if self.processed_data is not None else self.data, window_size=15, noise_power=None)
        self.plot_processed_seismic_image()
        self.data_info_label.setText("Applied wiener Deconvolution.")

    @pyqtSlot()
    def apply_Horizon_pick(self):
        if not hasattr(self, 'interpretation_window') or self.interpretation_window is None:
            self.interpretation_window = SeismicInterpretationWindow(
                seismic_data=self.processed_data if self.processed_data is not None else self.data
            )
        self.interpretation_window.show()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = SeismicDataApp()
    main_window.show()
    sys.exit(app.exec_())