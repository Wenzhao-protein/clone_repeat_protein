"""
SCN文件读取模块
用于读取和解析Bio-Rad凝胶成像仪的SCN格式文件
"""

import numpy as np
import re
from email import policy
from email.parser import BytesParser


def read_scn_file(filepath, priority='crop'):
    """
    读取SCN文件并提取图像数据和矩阵
    
    参数:
        filepath (str): SCN文件的完整路径
        priority (str): 尺寸优先级，'crop'使用裁剪尺寸，'whole'使用原始尺寸
    
    返回:
        tuple: (image_array, dimensions, error)
            - image_array (numpy.ndarray): 图像数据矩阵，如果失败则为None
            - dimensions (tuple): (width, height) 图像尺寸，如果失败则为None
            - error (str): 错误信息，如果成功则为None
    
    示例:
        >>> img_array, dims, error = read_scn_file('/path/to/file.scn')
        >>> if error is None:
        >>>     print(f"成功读取图像，尺寸: {dims}")
        >>>     print(f"矩阵形状: {img_array.shape}")
    """
    try:
        # 读取二进制内容
        with open(filepath, 'rb') as f:
            binary_content = f.read()
        
        # 同时读取文本内容（用于查找XML属性）
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            text_content = f.read()
        
        # 从XML属性提取尺寸
        org_match = re.search(r'<org_size_pix width="(\d+)" height="(\d+)"', text_content)
        size_match = re.search(r'<size_pix width="(\d+)" height="(\d+)"', text_content)
        
        if not org_match and not size_match:
            return None, None, "未找到尺寸信息"
        
        # 解析MIME结构
        msg = BytesParser(policy=policy.default).parsebytes(binary_content)
        
        # 查找ImageData部分
        image_data = None
        for part in msg.walk():
            if part.get('Content-Description') == 'ImageData':
                payload = part.get_payload(decode=True)
                
                # 有些文件在ImageData后附加了XML,需要分离
                xml_start = payload.find(b'<?xml')
                if xml_start != -1:
                    image_data = payload[:xml_start]
                else:
                    image_data = payload
                break
        
        if image_data is None:
            return None, None, "未找到ImageData"
        
        # 根据priority决定尝试顺序
        sizes_to_try = []
        if priority == 'whole':
            if org_match:
                sizes_to_try.append((int(org_match.group(1)), int(org_match.group(2))))
            if size_match:
                sizes_to_try.append((int(size_match.group(1)), int(size_match.group(2))))
        else:  # priority == 'crop' or default
            if size_match:
                sizes_to_try.append((int(size_match.group(1)), int(size_match.group(2))))
            if org_match:
                sizes_to_try.append((int(org_match.group(1)), int(org_match.group(2))))
        
        # 尝试所有尺寸配置
        for width, height in sizes_to_try:
            # 尝试不同的字节/像素配置
            for bpp in [3, 1, 2, 4]:
                expected = width * height * bpp
                
                # 允许一定的误差（可能有header）
                if abs(len(image_data) - expected) < 2000:
                    try:
                        if bpp == 1:
                            img_array = np.frombuffer(image_data[:width*height], dtype=np.uint8).reshape((height, width))
                        elif bpp == 2:
                            img_array = np.frombuffer(image_data[:width*height*2], dtype=np.uint16).reshape((height, width))
                        elif bpp == 3:
                            img_array = np.frombuffer(image_data[:width*height*3], dtype=np.uint8).reshape((height, width, 3))
                        else:  # bpp == 4
                            img_array = np.frombuffer(image_data[:width*height*4], dtype=np.uint8).reshape((height, width, 4))
                        
                        # 验证数据有效性
                        if img_array.max() > 0:
                            return img_array, (width, height), None
                    except:
                        continue
            
            # 如果精确匹配失败，尝试offset扫描（针对16-bit图像）
            expected_bytes = width * height * 2
            if len(image_data) >= expected_bytes:
                best_result, best_std = None, 0
                for offset in range(0, min(5000, len(image_data) - expected_bytes), 2):
                    if offset + expected_bytes <= len(image_data):
                        try:
                            img_test = np.frombuffer(image_data[offset:offset+expected_bytes], 
                                                    dtype=np.uint16).reshape(height, width)
                            if img_test.std() > best_std and img_test.max() > 0:
                                best_std = img_test.std()
                                best_result = img_test
                        except:
                            continue
                if best_result is not None:
                    return best_result, (width, height), None
        
        return None, None, "无法解析图像数据"
        
    except Exception as e:
        return None, None, f"处理失败: {e}"


def get_image_info(img_array):
    """
    获取图像数组的基本信息
    
    参数:
        img_array (numpy.ndarray): 图像数据矩阵
    
    返回:
        dict: 包含图像信息的字典
    """
    if img_array is None:
        return None
    
    info = {
        'shape': img_array.shape,
        'dtype': str(img_array.dtype),
        'min': float(img_array.min()),
        'max': float(img_array.max()),
        'mean': float(img_array.mean()),
        'std': float(img_array.std())
    }
    
    if len(img_array.shape) == 3:
        info['channels'] = img_array.shape[2]
        info['height'] = img_array.shape[0]
        info['width'] = img_array.shape[1]
    else:
        info['channels'] = 1
        info['height'] = img_array.shape[0]
        info['width'] = img_array.shape[1]
    
    return info
