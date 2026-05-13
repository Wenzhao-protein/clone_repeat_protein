"""
交互式凝胶电泳图像分析演示脚本
用鼠标框选矩形区域，分析沿着长边的信号强度分布
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import minimum_filter
from matplotlib.widgets import RectangleSelector
from matplotlib.patches import Rectangle

# 导入scn_reader
sys.path.append('/home/wenzhao/github_repo/clone_repeat_protein/agarose_gel_analysis/input')
from scn_reader import read_scn_file_enhanced

class InteractiveGelAnalyzer:
    """交互式凝胶电泳图像分析器"""
    
    def __init__(self, image, filename):
        self.image = image
        self.filename = filename
        self.fig = None
        self.ax = None
        self.selector = None
        
    def analyze_panel(self, x1, y1, x2, y2):
        """分析选中的矩形区域"""
        # 确保坐标顺序正确
        x1, x2 = int(min(x1, x2)), int(max(x1, x2))
        y1, y2 = int(min(y1, y2)), int(max(y1, y2))
        
        # 提取panel区域
        panel = self.image[y1:y2, x1:x2]
        
        # 判断长边方向
        height, width = panel.shape
        is_vertical = height > width
        
        # 沿着长边计算平均信号强度
        if is_vertical:
            # 垂直方向（从上到下）
            profile = np.mean(panel, axis=1)  # 每行的平均值
            position = np.arange(len(profile))
            direction = "垂直方向 (上→下)"
        else:
            # 水平方向（从左到右）
            profile = np.mean(panel, axis=0)  # 每列的平均值
            position = np.arange(len(profile))
            direction = "水平方向 (左→右)"
        
        # 估算背景：使用滚动最小值
        window_size = max(10, len(profile) // 20)
        background = minimum_filter(profile, size=window_size, mode='nearest')
        
        # 去除背景
        corrected = profile - background
        corrected = np.maximum(corrected, 0)  # 确保非负
        
        return {
            'panel': panel,
            'position': position,
            'profile': profile,
            'background': background,
            'corrected': corrected,
            'direction': direction,
            'is_vertical': is_vertical,
            'coords': (x1, y1, x2, y2)
        }
    
    def plot_results(self, analysis):
        """绘制分析结果（4个子图）"""
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # 1. 原图 + 矩形框
        ax1 = axes[0, 0]
        ax1.imshow(self.image, cmap='gray')
        x1, y1, x2, y2 = analysis['coords']
        rect = Rectangle((x1, y1), x2-x1, y2-y1, 
                        fill=False, edgecolor='red', linewidth=2)
        ax1.add_patch(rect)
        ax1.set_title(f'原始图像\n{self.filename}', fontsize=10)
        ax1.axis('off')
        
        # 2. 提取的panel
        ax2 = axes[0, 1]
        ax2.imshow(analysis['panel'], cmap='gray')
        ax2.set_title(f'选中区域\n{analysis["direction"]}', fontsize=10)
        ax2.axis('off')
        
        # 3. 原始信号 + 背景
        ax3 = axes[1, 0]
        ax3.plot(analysis['position'], analysis['profile'], 
                label='原始信号', color='blue', linewidth=1.5)
        ax3.plot(analysis['position'], analysis['background'], 
                label='估算背景', color='red', linewidth=1.5, linestyle='--')
        ax3.set_xlabel('位置 (像素)')
        ax3.set_ylabel('信号强度')
        ax3.set_title('原始信号与背景', fontsize=10)
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. 去背景后的信号
        ax4 = axes[1, 1]
        ax4.plot(analysis['position'], analysis['corrected'], 
                color='green', linewidth=2)
        ax4.fill_between(analysis['position'], analysis['corrected'], 
                        alpha=0.3, color='green')
        ax4.set_xlabel('位置 (像素)')
        ax4.set_ylabel('校正后信号强度')
        ax4.set_title(f'去除背景后的信号\n峰值: {analysis["corrected"].max():.0f}, 平均: {analysis["corrected"].mean():.0f}', 
                     fontsize=10)
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def select_rectangle(self):
        """启动交互式矩形选择"""
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.ax.imshow(self.image, cmap='gray')
        self.ax.set_title(f'用鼠标框选矩形区域进行分析\n{self.filename}\n'
                         '提示: 拖动鼠标选择区域，释放后自动分析', 
                         fontsize=12)
        self.ax.axis('off')
        
        def onselect(eclick, erelease):
            """矩形选择回调函数"""
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
            
            print(f"\n选中区域: ({x1:.0f}, {y1:.0f}) → ({x2:.0f}, {y2:.0f})")
            print(f"区域大小: {abs(x2-x1):.0f} × {abs(y2-y1):.0f} 像素")
            
            # 分析选中区域
            analysis = self.analyze_panel(x1, y1, x2, y2)
            
            # 绘制结果
            result_fig = self.plot_results(analysis)
            plt.show()
        
        # 创建矩形选择器
        self.selector = RectangleSelector(
            self.ax, onselect,
            useblit=True,
            button=[1],  # 左键
            minspanx=5, minspany=5,
            spancoords='pixels',
            interactive=True,
            props=dict(facecolor='red', alpha=0.3, edgecolor='red', linewidth=2)
        )
        
        plt.show()
        return self.selector


if __name__ == "__main__":
    # 查找SCN文件
    import glob
    import os
    
    scn_files = sorted(glob.glob("*.scn"))
    if not scn_files:
        print("未找到SCN文件！")
        sys.exit(1)
    
    # 使用第一个文件作为演示
    test_file = scn_files[0]
    print(f"加载文件: {test_file}")
    
    # 读取图像数据
    img_data, img_info = read_scn_file_enhanced(test_file, prefer_mode='crop')
    print(f"图像尺寸: {img_data.shape}")
    print(f"图像信息: {img_info}")
    
    # 创建分析器并启动交互式选择
    print("\n启动交互式分析...")
    print("请在图像上拖动鼠标框选一个矩形区域")
    
    analyzer = InteractiveGelAnalyzer(img_data, os.path.basename(test_file))
    analyzer.select_rectangle()
