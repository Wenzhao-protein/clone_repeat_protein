#!/usr/bin/env python3
"""
SCN File Metadata Reader
读取Bio-Rad SCN文件的元数据，包括创建时间等信息
"""

import struct
import re
from datetime import datetime
from pathlib import Path


def read_scn_metadata(filepath):
    """
    Read metadata from Bio-Rad SCN file
    
    Parameters:
    -----------
    filepath : str
        Path to the SCN file
    
    Returns:
    --------
    dict : Dictionary containing metadata
        - creation_time: datetime object or None
        - file_size: file size in bytes
        - mime_headers: dict of MIME headers
        - raw_metadata: raw metadata string
    """
    
    metadata = {
        'filepath': str(filepath),
        'filename': Path(filepath).name,
        'file_size': Path(filepath).stat().st_size,
        'creation_time': None,
        'mime_headers': {},
        'raw_metadata': '',
        'error': None
    }
    
    try:
        with open(filepath, 'rb') as f:
            # Read first part of file to find MIME headers
            header_data = f.read(10000)  # Read first 10KB
            
            try:
                # Try to decode as text to find MIME headers
                header_text = header_data.decode('latin-1')
                metadata['raw_metadata'] = header_text[:2000]  # Store first 2KB
                
                # Extract MIME headers (format: key: value)
                mime_pattern = r'^([A-Za-z\-]+):\s*(.+?)$'
                for line in header_text.split('\n')[:100]:  # Check first 100 lines
                    match = re.match(mime_pattern, line.strip())
                    if match:
                        key, value = match.groups()
                        metadata['mime_headers'][key] = value.strip()
                
                # Look for creation date/time in various formats
                # Common keys: Date, Creation-Date, Timestamp, X-LastSaveDate, etc.
                date_keys = ['X-LastSaveDate', 'Date', 'Creation-Date', 'Timestamp', 'Created', 'DateTime']
                
                for key in date_keys:
                    if key in metadata['mime_headers']:
                        date_str = metadata['mime_headers'][key]
                        # Try to parse date
                        creation_time = parse_date_string(date_str)
                        if creation_time:
                            metadata['creation_time'] = creation_time
                            metadata['date_source'] = f'mime_header:{key}'
                            break
                
                # If no date found in headers, look for date patterns in raw text
                if not metadata['creation_time']:
                    # Look for common date patterns
                    date_patterns = [
                        r'(\d{4}[-/]\d{2}[-/]\d{2}\s+\d{2}:\d{2}:\d{2})',  # 2022-12-12 14:30:45
                        r'(\d{2}[-/]\d{2}[-/]\d{4}\s+\d{2}:\d{2}:\d{2})',  # 12/12/2022 14:30:45
                        r'(\d{4}\d{2}\d{2}\s+\d{2}\d{2}\d{2})',  # 20221212 143045
                    ]
                    
                    for pattern in date_patterns:
                        match = re.search(pattern, header_text)
                        if match:
                            date_str = match.group(1)
                            creation_time = parse_date_string(date_str)
                            if creation_time:
                                metadata['creation_time'] = creation_time
                                metadata['date_source'] = 'pattern_match'
                                break
                
            except Exception as e:
                metadata['error'] = f"Error parsing metadata: {e}"
        
        # If still no creation time, use file system timestamp as fallback
        if not metadata['creation_time']:
            file_mtime = Path(filepath).stat().st_mtime
            metadata['creation_time'] = datetime.fromtimestamp(file_mtime)
            metadata['date_source'] = 'filesystem'
        else:
            if 'date_source' not in metadata:
                metadata['date_source'] = 'mime_header'
    
    except Exception as e:
        metadata['error'] = f"Error reading file: {e}"
    
    return metadata


def parse_date_string(date_str):
    """
    Try to parse date string in various formats
    
    Parameters:
    -----------
    date_str : str
        Date string to parse
    
    Returns:
    --------
    datetime or None
    """
    
    # List of common date formats
    date_formats = [
        '%Y.%m.%d.%H.%M.%S.%f',  # Bio-Rad format: 2022.12.13.13.26.11.391
        '%Y-%m-%d %H:%M:%S',
        '%Y/%m/%d %H:%M:%S',
        '%d-%m-%Y %H:%M:%S',
        '%d/%m/%Y %H:%M:%S',
        '%m/%d/%Y %H:%M:%S',
        '%Y%m%d %H%M%S',
        '%Y-%m-%d',
        '%Y/%m/%d',
        '%d-%m-%Y',
        '%d/%m/%Y',
        '%m/%d/%Y',
        '%Y%m%d',
        # RFC 2822 format (common in MIME headers)
        '%a, %d %b %Y %H:%M:%S %z',
        '%a, %d %b %Y %H:%M:%S',
        '%d %b %Y %H:%M:%S',
    ]
    
    for fmt in date_formats:
        try:
            return datetime.strptime(date_str.strip(), fmt)
        except ValueError:
            continue
    
    return None


def read_scn_batch(directory, pattern="*.scn"):
    """
    Read metadata from all SCN files in a directory
    
    Parameters:
    -----------
    directory : str
        Directory path
    pattern : str
        File pattern (default: "*.scn")
    
    Returns:
    --------
    list : List of metadata dictionaries
    """
    
    dir_path = Path(directory)
    scn_files = sorted(dir_path.glob(pattern))
    
    results = []
    for scn_file in scn_files:
        metadata = read_scn_metadata(scn_file)
        results.append(metadata)
    
    return results


def print_metadata(metadata):
    """Pretty print metadata"""
    
    print("=" * 80)
    print(f"File: {metadata['filename']}")
    print("=" * 80)
    print(f"Path: {metadata['filepath']}")
    print(f"Size: {metadata['file_size']:,} bytes ({metadata['file_size']/1024:.1f} KB)")
    
    if metadata['creation_time']:
        print(f"Creation Time: {metadata['creation_time'].strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Date Source: {metadata.get('date_source', 'unknown')}")
    else:
        print("Creation Time: Not found")
    
    if metadata['mime_headers']:
        print("\nMIME Headers:")
        for key, value in metadata['mime_headers'].items():
            print(f"  {key}: {value}")
    
    if metadata['error']:
        print(f"\nError: {metadata['error']}")
    
    print()


# Test code
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        # Test with specific file
        filepath = sys.argv[1]
        metadata = read_scn_metadata(filepath)
        print_metadata(metadata)
    else:
        # Test with current directory
        print("Testing with files in current directory...\n")
        results = read_scn_batch(".", "*.scn")
        
        if results:
            for metadata in results[:5]:  # Show first 5
                print_metadata(metadata)
        else:
            print("No SCN files found in current directory")
            print("\nUsage: python scn_metadata_reader.py <path_to_scn_file>")
