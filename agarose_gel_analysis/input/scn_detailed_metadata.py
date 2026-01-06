#!/usr/bin/env python3
"""
详细的SCN文件元数据读取器
使用email.parser来完整解析MIME结构，获取所有可能的元数据
"""

import re
from datetime import datetime
from email import policy
from email.parser import BytesParser
from pathlib import Path


def read_scn_detailed_metadata(filepath):
    """
    详细读取SCN文件的所有元数据
    
    Parameters:
    -----------
    filepath : str
        Path to the SCN file
    
    Returns:
    --------
    dict : Dictionary containing comprehensive metadata
    """
    
    metadata = {
        'filepath': str(filepath),
        'filename': Path(filepath).name,
        'file_size': Path(filepath).stat().st_size,
        'mime_parts': [],
        'all_headers': {},
        'xml_metadata': {},
        'dates': {},
        'creation_time': None,
        'error': None
    }
    
    try:
        # Read binary content
        with open(filepath, 'rb') as f:
            binary_content = f.read()
        
        # Parse MIME structure
        msg = BytesParser(policy=policy.default).parsebytes(binary_content)
        
        # Extract all headers from main message
        for key in msg.keys():
            metadata['all_headers'][key] = msg[key]
        
        # Walk through all MIME parts
        part_idx = 0
        for part in msg.walk():
            part_info = {
                'index': part_idx,
                'content_type': part.get_content_type(),
                'headers': {}
            }
            
            # Extract all headers from this part
            for key in part.keys():
                value = part[key]
                part_info['headers'][key] = value
                
                # Store in all_headers with part prefix
                metadata['all_headers'][f"Part{part_idx}_{key}"] = value
            
            # Check for Content-Description
            if 'Content-Description' in part_info['headers']:
                part_info['description'] = part_info['headers']['Content-Description']
            
            metadata['mime_parts'].append(part_info)
            part_idx += 1
        
        # Read as text to extract XML metadata
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            text_content = f.read()
        
        # Extract XML attributes
        xml_patterns = {
            'org_size': r'<org_size_pix width="(\d+)" height="(\d+)"',
            'size': r'<size_pix width="(\d+)" height="(\d+)"',
            'capture_time': r'<capture_time>([^<]+)</capture_time>',
            'created': r'<created>([^<]+)</created>',
            'modified': r'<modified>([^<]+)</modified>',
            'timestamp': r'<timestamp>([^<]+)</timestamp>',
            'date': r'<date>([^<]+)</date>',
            'acquisition_date': r'<acquisition_date>([^<]+)</acquisition_date>',
            'acquisition_time': r'<acquisition_time>([^<]+)</acquisition_time>',
        }
        
        for key, pattern in xml_patterns.items():
            match = re.search(pattern, text_content)
            if match:
                if 'size' in key:
                    metadata['xml_metadata'][key] = {
                        'width': int(match.group(1)),
                        'height': int(match.group(2))
                    }
                else:
                    metadata['xml_metadata'][key] = match.group(1)
        
        # Extract all date-related fields
        date_patterns = {
            'X-LastSaveDate': r'X-LastSaveDate:\s*([^\r\n]+)',
            'X-CreationDate': r'X-CreationDate:\s*([^\r\n]+)',
            'X-AcquisitionDate': r'X-AcquisitionDate:\s*([^\r\n]+)',
            'Date': r'Date:\s*([^\r\n]+)',
            'Timestamp': r'Timestamp:\s*([^\r\n]+)',
        }
        
        for key, pattern in date_patterns.items():
            match = re.search(pattern, text_content)
            if match:
                date_str = match.group(1).strip()
                metadata['dates'][key] = date_str
        
        # Also check for dates in XML
        for key, value in metadata['xml_metadata'].items():
            if any(date_word in key.lower() for date_word in ['time', 'date', 'created', 'modified']):
                metadata['dates'][f'XML_{key}'] = value
        
        # Try to parse all dates
        parsed_dates = {}
        for source, date_str in metadata['dates'].items():
            parsed = parse_date_string(date_str)
            if parsed:
                parsed_dates[source] = parsed
        
        # Try to extract date from filename
        filename_date = extract_date_from_filename(Path(filepath).name)
        if filename_date:
            metadata['dates']['filename'] = filename_date.strftime('%Y-%m-%d')
            metadata['parsed_dates_temp'] = {**parsed_dates}  # temp copy
            parsed_dates['filename'] = filename_date
        
        # Store all parsed dates
        metadata['parsed_dates'] = parsed_dates
        
        # Determine creation time priority
        # Priority: filename > X-CreationDate > X-AcquisitionDate > XML_capture_time > Date/X-LastSaveDate > filesystem
        # Note: X-LastSaveDate and Date appear to be the last save/modification time, not initial creation
        priority_sources = [
            'filename',  # Best guess - often contains the actual gel run date
            'X-CreationDate',
            'X-AcquisitionDate', 
            'XML_capture_time',
            'XML_acquisition_time',
            'XML_created',
            'Date',  # Likely same as X-LastSaveDate (last save time)
            'Timestamp',
            'XML_timestamp',
            'X-LastSaveDate',  # Last save/modification, not creation
            'XML_modified'
        ]
        
        for source in priority_sources:
            if source in parsed_dates:
                metadata['creation_time'] = parsed_dates[source]
                metadata['date_source'] = source
                break
        
        # Store all parsed dates
        metadata['parsed_dates'] = parsed_dates
        
        # Fallback to filesystem
        if not metadata['creation_time']:
            file_mtime = Path(filepath).stat().st_mtime
            metadata['creation_time'] = datetime.fromtimestamp(file_mtime)
            metadata['date_source'] = 'filesystem'
    
    except Exception as e:
        metadata['error'] = f"Error reading file: {e}"
        import traceback
        metadata['error_trace'] = traceback.format_exc()
    
    return metadata


def extract_date_from_filename(filename):
    """
    Extract date from filename
    
    Parameters:
    -----------
    filename : str
        Filename to parse
    
    Returns:
    --------
    datetime or None
    """
    import re
    
    # Try to match YYYY-MM-DD format at start of filename
    match = re.match(r'(\d{4})-(\d{2})-(\d{2})', filename)
    if match:
        try:
            return datetime(int(match.group(1)), int(match.group(2)), int(match.group(3)))
        except ValueError:
            pass
    
    # Try to match YYYYMMDD format
    match = re.match(r'(\d{4})(\d{2})(\d{2})', filename)
    if match:
        try:
            return datetime(int(match.group(1)), int(match.group(2)), int(match.group(3)))
        except ValueError:
            pass
    
    return None


def parse_date_string(date_str):
    """
    Try to parse date string in various formats
    """
    
    date_formats = [
        # Bio-Rad formats
        '%Y.%m.%d.%H.%M.%S.%f',  # 2022.12.13.13.26.11.391
        '%Y.%m.%d.%H.%M.%S',     # 2022.12.13.13.26.11
        '%Y.%m.%d',              # 2022.12.13
        # ISO formats
        '%Y-%m-%dT%H:%M:%S.%fZ',
        '%Y-%m-%dT%H:%M:%SZ',
        '%Y-%m-%dT%H:%M:%S',
        '%Y-%m-%d %H:%M:%S.%f',
        '%Y-%m-%d %H:%M:%S',
        '%Y-%m-%d',
        # Common formats
        '%Y/%m/%d %H:%M:%S',
        '%Y/%m/%d',
        '%d-%m-%Y %H:%M:%S',
        '%d/%m/%Y %H:%M:%S',
        '%m/%d/%Y %H:%M:%S',
        '%Y%m%d%H%M%S',
        '%Y%m%d',
        # RFC formats
        '%a, %d %b %Y %H:%M:%S %z',
        '%a, %d %b %Y %H:%M:%S',
        '%d %b %Y %H:%M:%S',
    ]
    
    for fmt in date_formats:
        try:
            return datetime.strptime(date_str.strip(), fmt)
        except (ValueError, AttributeError):
            continue
    
    return None


def get_scn_creation_time(filepath):
    """
    简洁版：获取SCN文件的创建时间
    
    Parameters:
    -----------
    filepath : str
        Path to the SCN file
    
    Returns:
    --------
    datetime : 创建时间（datetime对象）
        优先从文件名提取，其次从MIME头，最后使用文件系统时间
    """
    from pathlib import Path
    
    # 1. 尝试从文件名提取日期（最可靠）
    filename_date = extract_date_from_filename(Path(filepath).name)
    if filename_date:
        return filename_date
    
    # 2. 尝试从MIME头读取
    try:
        with open(filepath, 'rb') as f:
            content = f.read(5000)  # 只读前5KB
        
        content_str = content.decode('latin-1', errors='ignore')
        
        # 查找X-LastSaveDate或Date字段
        for pattern in [r'X-LastSaveDate:\s*([^\r\n]+)', r'Date:\s*([^\r\n]+)']:
            match = re.search(pattern, content_str)
            if match:
                date_str = match.group(1).strip()
                parsed = parse_date_string(date_str)
                if parsed:
                    return parsed
    except:
        pass
    
    # 3. 回退到文件系统时间
    return datetime.fromtimestamp(Path(filepath).stat().st_mtime)


def get_scn_creation_timestamp(filepath):
    """
    最简版：获取SCN文件的创建时间戳
    
    Parameters:
    -----------
    filepath : str
        Path to the SCN file
    
    Returns:
    --------
    float : Unix时间戳
    """
    return get_scn_creation_time(filepath).timestamp()


def print_detailed_metadata(metadata):
    """Pretty print detailed metadata"""
    
    print("=" * 100)
    print(f"SCN File Detailed Metadata Analysis")
    print("=" * 100)
    print(f"File: {metadata['filename']}")
    print(f"Path: {metadata['filepath']}")
    print(f"Size: {metadata['file_size']:,} bytes ({metadata['file_size']/1024:.1f} KB)")
    print()
    
    # Date information
    print("-" * 100)
    print("DATE INFORMATION")
    print("-" * 100)
    
    if metadata['creation_time']:
        print(f"✓ Best Guess Creation Time: {metadata['creation_time'].strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Source: {metadata.get('date_source', 'unknown')}")
    else:
        print("✗ No creation time found")
    
    print("\nAll Date Fields Found:")
    if metadata['dates']:
        for source, date_str in metadata['dates'].items():
            parsed = metadata['parsed_dates'].get(source)
            if parsed:
                print(f"  {source:.<40} {date_str:.<30} ✓ {parsed.strftime('%Y-%m-%d %H:%M:%S')}")
            else:
                print(f"  {source:.<40} {date_str:.<30} ✗ (parse failed)")
    else:
        print("  (none found)")
    
    # MIME structure
    print("\n" + "-" * 100)
    print("MIME STRUCTURE")
    print("-" * 100)
    print(f"Number of MIME parts: {len(metadata['mime_parts'])}")
    
    for part in metadata['mime_parts']:
        print(f"\n  Part {part['index']}:")
        print(f"    Content-Type: {part['content_type']}")
        if 'description' in part:
            print(f"    Description: {part['description']}")
        print(f"    Headers ({len(part['headers'])}):")
        for key, value in part['headers'].items():
            value_str = str(value)[:80]  # Truncate long values
            print(f"      {key}: {value_str}")
    
    # XML metadata
    print("\n" + "-" * 100)
    print("XML METADATA")
    print("-" * 100)
    
    if metadata['xml_metadata']:
        for key, value in metadata['xml_metadata'].items():
            print(f"  {key}: {value}")
    else:
        print("  (none found)")
    
    # All headers
    print("\n" + "-" * 100)
    print("ALL HEADERS")
    print("-" * 100)
    
    if metadata['all_headers']:
        for key, value in sorted(metadata['all_headers'].items()):
            value_str = str(value)[:100]  # Truncate long values
            print(f"  {key}: {value_str}")
    else:
        print("  (none found)")
    
    # Error information
    if metadata['error']:
        print("\n" + "-" * 100)
        print("ERROR")
        print("-" * 100)
        print(metadata['error'])
        if 'error_trace' in metadata:
            print("\nTraceback:")
            print(metadata['error_trace'])
    
    print("\n" + "=" * 100)


# Test code
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        # Test with specific file
        filepath = sys.argv[1]
        metadata = read_scn_detailed_metadata(filepath)
        print_detailed_metadata(metadata)
    else:
        print("Usage: python scn_detailed_metadata.py <path_to_scn_file>")
