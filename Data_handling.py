"""
Data_handling.py

This module is designed to handle the loading, processing, and saving of seismic data from SEG-Y files. 
It uses two primary libraries: `segyio` for efficient handling of SEG-Y files and `obspy` for reading 
seismic data formats. The module provides functionality to:
  - Load seismic data using either `segyio` or `obspy`.
  - Extract basic seismic metadata such as sample rate, number of traces, and two-way travel time (TWT).
  - Save processed seismic data back to SEG-Y format.
  - Close any open seismic files gracefully.

Classes:
    SEGYHandler: Handles loading, saving, and processing seismic data from SEG-Y files.
"""

import segyio
import obspy

class SEGYHandler:

    """
    SEGYHandler class provides methods to load, save, and manage seismic data from SEG-Y files.

    Attributes:
        segyio_file (segyio.SegyFile): A SEG-Y file object used by segyio for data handling.
        stream (obspy.Stream): A stream object used by obspy to handle seismic traces.

    Methods:
        load_data_segyio(file_path): Load seismic data using Segyio.
        load_data_obspy(file_path): Load seismic data using Obspy.
        close_file(): Close the opened SEG-Y file.
        save_segy_file(file_path, file_spec, save_data): Save processed data back into a SEG-Y file.
    """

    def __init__(self):

        """Initialize SEGYHandler with default attributes."""

        self.segyio_file = None
        self.stream = None

    def load_data_segyio(self, file_path):
    
        """
        Load seismic data from a SEG-Y file using Segyio.

        Parameters:
            file_path (str): The path to the SEG-Y file to be loaded.

        Returns:
            tuple: Contains the loaded SEG-Y file object, metadata specification, 
                   seismic data, number of samples, two-way travel time (TWT),
                   data format, sample interval, and sample rate.

        Raises:
            ValueError: If the file cannot be opened or parsed properly.
        """
        print("Opening the seismic file with Segyio")
        try:
            with segyio.open(file_path, 'r', ignore_geometry=True) as segyfile:
                self.segyio_file = segyfile
                spec = segyio.tools.metadata(segyfile)
                seismic_data = segyfile.trace.raw[:]
                print(seismic_data)

                # Extract basic metadata such as format, number of traces, sample rate, etc.
                n_traces = segyfile.tracecount
                data_format = segyfile.format
                
                sample_interval = segyio.tools.dt(segyfile) / 1e6  # Convert sample interval to seconds because is in μseconds
                sample_rate = 1 / sample_interval
                n_samples = segyfile.samples.size
                twt = segyfile.samples
                
                # Displaying extracted data
                print(f"Data Format: {data_format}")
                print(f"Number of Traces: {n_traces}")
                print(f"Number of Samples: {n_samples}")
                print(f"Two Way Time (TWT): {twt}")
                print(f"Sample Rate: {sample_rate}")

                return self.segyio_file, spec, seismic_data, n_traces, n_samples, twt, data_format, sample_interval, sample_rate
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error: Can't open the raw segyfile: {e}")
        
    def load_data_obspy(self, file_path):
        
        """
        Load seismic data from a SEG-Y file using Obspy.

        Parameters:
            file_path (str): The path to the SEG-Y file to be loaded.

        Returns:
            tuple: Contains the obspy stream object and the data type of the loaded seismic data.

        Raises:
            ValueError: If the file cannot be opened or parsed properly.
        """

        print("Opening the seismic file with Obspy")
        try:
            self.stream = obspy.read(file_path, format='segy')
            trace = self.stream[0]  # Accessing the first trace for metadata
            sample_rate = trace.stats.sampling_rate  # Get sample rate from the trace metadata
            samples = trace.stats.npts  # Number of samples in the trace
            raw_seismic_data_type = trace.data[0].dtype  # Data type of the seismic samples

            return self.stream, raw_seismic_data_type
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error: Can't open the raw segyfile with ObsPy: {e}")

    def close_file(self):
        
        """
        Close the opened SEG-Y file and clear the data.

        This method closes any SEG-Y file opened using Segyio and clears the 
        loaded stream data from Obspy if applicable.
        """

        if self.segyio_file is not None:
            print("Closing the Segyio file")
            self.segyio_file.close()
            self.segyio_file = None
            
            """
        if self.stream is not None:
            print("Clearing the Obspy stream")
            self.stream.clear()  # clears data but doesn't actually close a file, since ObsPy reads into memory
            self.stream = None
            """
        print("File closed successfully.")

    def save_segy_file(self, file_path, file_spec, save_data):
        
        """
        Save processed seismic data back to a SEG-Y file.

        Parameters:
            file_path (str): The output path where the SEG-Y file will be saved.
            file_spec (dict): The SEG-Y file specification (metadata) to be used for the saved file.
            save_data (ndarray): The seismic data that will be written into the SEG-Y file.

        Returns:
            segyio.SegyFile: The saved SEG-Y file object.

        Raises:
            ValueError: If the file cannot be saved due to an error.
        """
        try:
            with segyio.create(file_path, file_spec) as dst_file:
                self.dst_segyfile = dst_file
                dst_file.trace = save_data
                
                return self.dst_segyfile
            
        except (ValueError, IndexError, Exception) as e:
            raise ValueError(f"Error saving SEG-Y file: {e}")    