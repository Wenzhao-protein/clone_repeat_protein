import pandas as pd
import numpy as np
import h5py
import datetime # Retained for general date/time operations
import scipy.io

# %pip install h5py # Ensure this is run in a separate cell if needed, or managed in your environment

class MatFileProcessor:
    """
    Processes SEC data from a .mat file, supporting both HDF5 v7.3 and older formats.
    """
    FINAL_DF_COLUMNS = ['datetime', 'abs_wl', 'sample_name', 'time', 'signal']

    def __init__(self, mat_filename):
        self.mat_filename = mat_filename
        self.h_file = None # To store the h5py file object

    def _get_str_from_hdf5_char_array(self, item_ref):
        """Extracts a string from an HDF5 char array reference."""
        char_codes = self.h_file[item_ref][:].ravel()
        return "".join(chr(c) for c in char_codes)

    def _get_array_from_hdf5_dataset(self, item_ref):
        """Extracts a 1D numpy array from an HDF5 dataset reference."""
        return self.h_file[item_ref][:].ravel()

    def _determine_hdf5_access_method(self, hdf5_dataset_or_group_field):
        """Determines the number of elements and access method for HDF5 datasets."""
        num_elements = 0
        access_method = None

        if hdf5_dataset_or_group_field.ndim == 0: # Scalar
            num_elements = 1
            access_method = lambda ds, i: ds[()]
        elif hdf5_dataset_or_group_field.ndim == 1: # 1D array, shape (N,)
            num_elements = hdf5_dataset_or_group_field.shape[0]
            access_method = lambda ds, i: ds[i]
        elif hdf5_dataset_or_group_field.ndim == 2: # 2D array
            if hdf5_dataset_or_group_field.shape[1] == 1: # Nx1
                num_elements = hdf5_dataset_or_group_field.shape[0]
                access_method = lambda ds, i: ds[i, 0]
            elif hdf5_dataset_or_group_field.shape[0] == 1: # 1xN
                num_elements = hdf5_dataset_or_group_field.shape[1]
                access_method = lambda ds, i: ds[0, i]
            elif hdf5_dataset_or_group_field.shape == (1,1): # Single element 2D array
                num_elements = 1
                access_method = lambda ds, i: ds[0,0]
            else:
                raise ValueError(f"Dataset has shape {hdf5_dataset_or_group_field.shape}, not a simple vector of references.")
        else:
            raise ValueError(f"Dataset has unsupported ndim {hdf5_dataset_or_group_field.ndim}")
        return num_elements, access_method

    def _process_hdf5_ans_dataset(self, ans_dataset_of_refs):
        """Processes 'ans' when it's an HDF5 Dataset of object references."""
        data_list = []
        num_elements, access_method_ref = self._determine_hdf5_access_method(ans_dataset_of_refs)

        for i in range(num_elements):
            struct_ref = access_method_ref(ans_dataset_of_refs, i)
            if not isinstance(self.h_file[struct_ref], h5py.Group):
                raise ValueError(f"Element {i} in 'ans' dataset does not reference a group (struct).")
            struct_group = self.h_file[struct_ref]
            
            data_list.append({
                'datetime_str': self._get_str_from_hdf5_char_array(struct_group['DateTime'][()]),
                'abs_wl': self._get_str_from_hdf5_char_array(struct_group['SigDesc'][()]),
                'sample_name': self._get_str_from_hdf5_char_array(struct_group['SampleName'][()]),
                'time': self._get_array_from_hdf5_dataset(struct_group['Time'][()]),
                'signal': self._get_array_from_hdf5_dataset(struct_group['Signal'][()])
            })
        return data_list

    def _process_hdf5_ans_group(self, ans_group):
        """Processes 'ans' when it's an HDF5 Group where fields are datasets."""
        data_list = []
        required_fields = ['DateTime', 'SigDesc', 'SampleName', 'Time', 'Signal']
        for field in required_fields:
            if field not in ans_group:
                raise ValueError(f"'ans' group missing field dataset '{field}'.")
            if not isinstance(ans_group[field], h5py.Dataset):
                raise ValueError(f"Field '{field}' in 'ans' group is not a dataset.")

        dt_field_ds = ans_group['DateTime'] # Dataset of references
        num_elements, access_method_field = self._determine_hdf5_access_method(dt_field_ds)

        for i in range(num_elements):
            data_list.append({
                'datetime_str': self._get_str_from_hdf5_char_array(access_method_field(ans_group['DateTime'], i)),
                'abs_wl': self._get_str_from_hdf5_char_array(access_method_field(ans_group['SigDesc'], i)),
                'sample_name': self._get_str_from_hdf5_char_array(access_method_field(ans_group['SampleName'], i)),
                'time': self._get_array_from_hdf5_dataset(access_method_field(ans_group['Time'], i)),
                'signal': self._get_array_from_hdf5_dataset(access_method_field(ans_group['Signal'], i))
            })
        return data_list

    def _read_with_h5py(self):
        """Reads data using h5py for HDF5 v7.3 .mat files."""
        with h5py.File(self.mat_filename, 'r') as h_file:
            self.h_file = h_file # Store for helper methods
            if 'ans' not in h_file:
                raise ValueError("MAT file (HDF5 v7.3) does not contain 'ans' dataset/group.")

            ans_entity = h_file['ans']
            data_list = []

            if isinstance(ans_entity, h5py.Dataset):
                data_list = self._process_hdf5_ans_dataset(ans_entity)
            elif isinstance(ans_entity, h5py.Group):
                data_list = self._process_hdf5_ans_group(ans_entity)
            else:
                raise ValueError(f"'ans' in HDF5 file is neither a Dataset nor a Group, but {type(ans_entity)}.")
            
            self.h_file = None # Clear after use
            return self._finalize_dataframe(data_list)

    def _read_with_scipy(self):
        """Reads data using scipy.io.loadmat for older .mat files."""
        mat = scipy.io.loadmat(self.mat_filename, squeeze_me=True, struct_as_record=False)
        
        if 'ans' not in mat:
             raise ValueError("MAT file (scipy) does not contain 'ans' variable.")

        ans_data = mat['ans']

        if not isinstance(ans_data, (list, np.ndarray)):
            ans_data = [ans_data] 
        elif isinstance(ans_data, np.ndarray) and ans_data.ndim == 0:
            ans_data = [ans_data.item()]
        
        if len(ans_data) == 0 or (hasattr(ans_data, 'size') and ans_data.size == 0) :
             return pd.DataFrame(columns=self.FINAL_DF_COLUMNS)
            
        processed_data_list = []
        for s_obj in ans_data:
            if s_obj is None or not all(hasattr(s_obj, attr) for attr in ['DateTime', 'SigDesc', 'SampleName', 'Time', 'Signal']):
                print(f"Skipping an invalid/empty struct object: {s_obj}")
                continue

            processed_data_list.append({
                'datetime_str': s_obj.DateTime,
                'abs_wl': s_obj.SigDesc,
                'sample_name': s_obj.SampleName,
                'time': np.asarray(s_obj.Time).ravel() if hasattr(s_obj.Time, 'ravel') else np.asarray(s_obj.Time), 
                'signal': np.asarray(s_obj.Signal).ravel() if hasattr(s_obj.Signal, 'ravel') else np.asarray(s_obj.Signal)
            })
        
        return self._finalize_dataframe(processed_data_list)

    def _finalize_dataframe(self, data_list):
        """Converts a list of dicts to a DataFrame and processes datetime."""
        if not data_list:
            return pd.DataFrame(columns=self.FINAL_DF_COLUMNS)
        
        df = pd.DataFrame(data_list)
        df['datetime'] = pd.to_datetime(df['datetime_str'], format='%Y-%m-%d %H:%M:%S')
        return df[self.FINAL_DF_COLUMNS]

    def read_mat_file(self):
        """
        Reads SEC data from the .mat file.
        Attempts h5py for v7.3 files, falls back to scipy.io.loadmat.
        """
        is_hdf5_file = False
        try:
            is_hdf5_file = h5py.is_hdf5(self.mat_filename)
        except Exception:
            pass # is_hdf5_file remains False

        if is_hdf5_file:
            try:
                return self._read_with_h5py()
            except Exception as e_h5py:
                print(f"Failed to read MAT file '{self.mat_filename}' with h5py (v7.3 attempt) ({type(e_h5py).__name__}: {e_h5py}). Attempting with scipy.io.loadmat.")
                # Fall through to scipy.io.loadmat
        
        # Fallback or direct attempt for non-HDF5 files
        try:
            return self._read_with_scipy()
        except Exception as e_scipy:
            if is_hdf5_file: 
                print(f"Failed to read MAT file '{self.mat_filename}' with scipy.io.loadmat as well ({type(e_scipy).__name__}: {e_scipy}). Both methods failed.")
            else:
                print(f"Failed to read MAT file '{self.mat_filename}' with scipy.io.loadmat ({type(e_scipy).__name__}: {e_scipy}).")
            raise

def read_mat(mat_filename):
    """
    Reads SEC data from a .mat file using the MatFileProcessor.
    """
    processor = MatFileProcessor(mat_filename)
    return processor.read_mat_file()


def filter_by_sample_name(input_df):
    """
    Filters the DataFrame based on datetime, abs_wl, and sample_name criteria.

    Args:
        input_df (pd.DataFrame): The DataFrame to filter.

    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    
    df_filtered_step1 = input_df.copy()
    
    # 2. Filter by sample_name criteria, applied to the result of step 1
    
    # Condition: sample_name should not be NaN.
    condition_not_na = df_filtered_step1['sample_name'].notna()
    condition_not_empty_string = df_filtered_step1['sample_name'] != ''
    condition_not_literal_blank = df_filtered_step1['sample_name'].str.lower() != 'blank'
    condition_not_contain_wash = ~df_filtered_step1['sample_name'].str.lower().str.contains('wash', na=False)
    condition_length_greater_than_2 = df_filtered_step1['sample_name'].str.len() > 2

    # Apply all sample_name conditions together
    df_final_filtered = df_filtered_step1[
        condition_not_na &
        condition_not_empty_string &
        condition_not_literal_blank &
        condition_not_contain_wash &
        condition_length_greater_than_2
    ]
    
    return df_final_filtered


def filter_by_abs_wl(df, abs_wl_values=False, remove_fld=True):
    """
    Filters the DataFrame based on the specified abs_wl values.
    
    Args:
        df (pd.DataFrame): The DataFrame to filter.
        abs_wl_values (list): List of abs_wl values to filter by.
        
    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    if remove_fld:
        df = df[~(df['abs_wl'].str.contains('FLD'))]

    if abs_wl_values != False:
        pattern = '|'.join(map(str, abs_wl_values))
        df = df[df['abs_wl'].str.contains(pattern)]

    df['abs_wl'] = "A"+df['abs_wl'].str[11:14]

    return df


