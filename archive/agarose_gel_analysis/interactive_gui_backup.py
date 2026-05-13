#!/usr/bin/env python3
"""
Interactive Gel Electrophoresis Image Analysis - Tkinter GUI Version

Features:
- Mouse drag to select rectangular region
- Real-time coordinate display
- Automatic analysis with 4-subplot results
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.ndimage import minimum_filter
import sys
import glob
import re

# Import scn_reader
from scn_reader import read_scn_file_enhanced, get_image_info


# ============================================================
# SANITIZATION UTILITY - CRITICAL FOR CLEAN GUI DISPLAY
# ============================================================

def sanitize_for_gui(text):
    """
    Sanitize text for GUI display - ONLY ASCII output.
    
    Removes:
    - Unicode escape sequences like \\uXXXX or \\UXXXXXXXX
    - Replacement characters (\\ufffd or �)
    - All non-ASCII characters
    - Control characters except newline
    
    Args:
        text: Input string (may contain non-ASCII or escape sequences)
        
    Returns:
        Clean ASCII-only string safe for GUI display
    """
    if text is None:
        return "N/A"
    
    if not isinstance(text, str):
        text = str(text)
    
    # Step 1: Remove literal unicode escape sequences like \\u51c8\\u5907
    # Pattern matches \\u followed by 4 hex digits or \\U followed by 8 hex digits
    text = re.sub(r'\\[uU][0-9a-fA-F]{4,8}', '', text)
    
    # Step 2: Remove replacement character (both literal and escaped forms)
    text = text.replace('\ufffd', '')
    text = text.replace('�', '')
    text = text.replace('\\ufffd', '')
    
    # Step 3: Convert any remaining non-ASCII to '?'
    # Keep only printable ASCII (0x20-0x7E) plus newline (0x0A)
    cleaned = []
    for char in text:
        code = ord(char)
        if 0x20 <= code <= 0x7E:  # Printable ASCII
            cleaned.append(char)
        elif code == 0x0A:  # Newline
            cleaned.append(char)
        elif code > 0x7E:  # Non-ASCII (e.g., Chinese)
            cleaned.append('?')
        # Ignore other control characters
    
    result = ''.join(cleaned).strip()
    
    # Step 4: If result is empty or only contains '?', return placeholder
    if not result or all(c == '?' for c in result):
        return "UNKNOWN"
    
    # Step 5: Replace common patterns with readable English
    result = result.replace('???', '')  # Remove consecutive question marks
    
    # Collapse multiple spaces
    result = re.sub(r'\s+', ' ', result)
    
    return result.strip() or "N/A"


def safe_decode_bytes(data):
    """
    Safely decode bytes to string, falling back through multiple encodings.
    
    Args:
        data: bytes object
        
    Returns:
        Decoded string (may still contain non-ASCII, use sanitize_for_gui after)
    """
    if isinstance(data, str):
        return data
    
    # Try UTF-8 first (strict mode)
    try:
        return data.decode('utf-8', errors='strict')
    except (UnicodeDecodeError, AttributeError):
        pass
    
    # Try UTF-8 with BOM
    try:
        return data.decode('utf-8-sig', errors='strict')
    except (UnicodeDecodeError, AttributeError):
        pass
    
    # Fall back to latin-1 (always succeeds but may produce garbage)
    try:
        return data.decode('latin-1', errors='strict')
    except Exception:
        pass
    
    # Last resort: return placeholder
    return "UNREADABLE_DATA"


class InteractiveGelAnalyzer:
    """Interactive Gel Electrophoresis Image Analyzer"""
    
    def __init__(self, image, filename):
        self.image = image
        self.filename = sanitize_for_gui(filename)
        self.rotation_angle = 0  # Track rotation angle
        
    def rotate_image(self, angle):
        """Rotate image by angle degrees"""
        from scipy.ndimage import rotate
        self.rotation_angle = (self.rotation_angle + angle) % 360
        self.image = rotate(self.image, angle, reshape=False, order=1)
        return self.image
        
    def analyze_panel(self, x1, y1, x2, y2, plasmid_coords=None, fragment_coords=None):
        """Analyze selected rectangular region with plasmid and fragment bands"""
        # Ensure correct coordinate order
        x1, x2 = int(min(x1, x2)), int(max(x1, x2))
        y1, y2 = int(min(y1, y2)), int(max(y1, y2))
        
        # Extract panel region
        panel = self.image[y1:y2, x1:x2]
        
        # Determine long edge direction
        height, width = panel.shape
        is_vertical = height > width
        
        # Calculate average signal intensity along long edge
        if is_vertical:
            profile = np.mean(panel, axis=1)
            position = np.arange(len(profile))
            direction = "Vertical (Top to Bottom)"
        else:
            profile = np.mean(panel, axis=0)
            position = np.arange(len(profile))
            direction = "Horizontal (Left to Right)"
        
        # Use LINEAR background estimation instead of minimum_filter
        from numpy.polynomial import Polynomial
        # Fit a linear polynomial (degree 1) to the profile
        p = Polynomial.fit(position, profile, deg=1)
        background = p(position)
        
        # Remove background
        corrected = profile - background
        corrected = np.maximum(corrected, 0)
        
        # Calculate plasmid and fragment intensities if provided
        plasmid_intensity = 0
        fragment_intensity = 0
        fragment_ratio = 0
        
        if plasmid_coords:
            px1, py1, px2, py2 = plasmid_coords
            px1, px2 = int(min(px1, px2)), int(max(px1, px2))
            py1, py2 = int(min(py1, py2)), int(max(py1, py2))
            plasmid_band = self.image[py1:py2, px1:px2]
            plasmid_intensity = np.sum(plasmid_band)
        
        if fragment_coords:
            fx1, fy1, fx2, fy2 = fragment_coords
            fx1, fx2 = int(min(fx1, fx2)), int(max(fx1, fx2))
            fy1, fy2 = int(min(fy1, fy2)), int(max(fy1, fy2))
            fragment_band = self.image[fy1:fy2, fx1:fx2]
            fragment_intensity = np.sum(fragment_band)
        
        # Calculate fragment ratio
        total_intensity = plasmid_intensity + fragment_intensity
        if total_intensity > 0:
            fragment_ratio = fragment_intensity / total_intensity
        
        return {
            'panel': panel,
            'position': position,
            'profile': profile,
            'background': background,
            'corrected': corrected,
            'direction': direction,
            'is_vertical': is_vertical,
            'coords': (x1, y1, x2, y2),
            'plasmid_intensity': plasmid_intensity,
            'fragment_intensity': fragment_intensity,
            'fragment_ratio': fragment_ratio,
            'plasmid_coords': plasmid_coords,
            'fragment_coords': fragment_coords
        }
    
    def plot_results(self, analysis):
        """Plot analysis results (4 subplots with A,B,C,D labels)"""
        fig, axes = plt.subplots(2, 2, figsize=(18, 14))
        
        # Larger font sizes
        title_fontsize = 16
        label_fontsize = 14
        legend_fontsize = 12
        subplot_label_fontsize = 24
        
        # 1. Original image + rectangle boxes (LABEL: A)
        ax1 = axes[0, 0]
        ax1.imshow(self.image, cmap='gray')
        x1, y1, x2, y2 = analysis['coords']
        
        # Draw panel rectangle (RED)
        from matplotlib.patches import Rectangle
        rect = Rectangle((x1, y1), x2-x1, y2-y1, 
                        fill=False, edgecolor='red', linewidth=3, label='Panel')
        ax1.add_patch(rect)
        
        # Draw plasmid band rectangle (GREEN) if available
        if analysis.get('plasmid_coords'):
            px1, py1, px2, py2 = analysis['plasmid_coords']
            prect = Rectangle((px1, py1), px2-px1, py2-py1,
                            fill=False, edgecolor='green', linewidth=3, label='Plasmid')
            ax1.add_patch(prect)
        
        # Draw fragment band rectangle (BLUE) if available
        if analysis.get('fragment_coords'):
            fx1, fy1, fx2, fy2 = analysis['fragment_coords']
            frect = Rectangle((fx1, fy1), fx2-fx1, fy2-fy1,
                            fill=False, edgecolor='blue', linewidth=3, label='Fragment')
            ax1.add_patch(frect)
        
        ax1.set_title(f'Original Image\n{self.filename}', fontsize=title_fontsize, pad=15)
        ax1.axis('off')
        ax1.legend(fontsize=legend_fontsize, loc='upper right')
        # Add subplot label A
        ax1.text(0.02, 0.98, 'A', transform=ax1.transAxes, fontsize=subplot_label_fontsize,
                fontweight='bold', va='top', ha='left', color='black',
                bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8))
        
        # 2. Extracted panel (LABEL: B)
        ax2 = axes[0, 1]
        ax2.imshow(analysis['panel'], cmap='gray')
        
        # Add intensity information as text
        info_text = f"{analysis['direction']}\n"
        if analysis['plasmid_intensity'] > 0 or analysis['fragment_intensity'] > 0:
            info_text += f"Plasmid Intensity: {analysis['plasmid_intensity']:.0f}\n"
            info_text += f"Fragment Intensity: {analysis['fragment_intensity']:.0f}\n"
            info_text += f"Fragment Ratio: {analysis['fragment_ratio']:.3f}"
        
        ax2.set_title(f'Selected Region\n{info_text}', fontsize=title_fontsize, pad=15)
        ax2.axis('off')
        # Add subplot label B
        ax2.text(0.02, 0.98, 'B', transform=ax2.transAxes, fontsize=subplot_label_fontsize,
                fontweight='bold', va='top', ha='left', color='black',
                bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8))
        
        # 3. Original signal + LINEAR background (LABEL: C)
        ax3 = axes[1, 0]
        ax3.plot(analysis['position'], analysis['profile'], 
                label='Raw Signal', color='blue', linewidth=2.5)
        ax3.plot(analysis['position'], analysis['background'], 
                label='Linear Background', color='red', linewidth=2.5, linestyle='--')
        ax3.set_xlabel('Position (pixels)', fontsize=label_fontsize)
        ax3.set_ylabel('Signal Intensity', fontsize=label_fontsize)
        ax3.set_title('Raw Signal vs Linear Background', fontsize=title_fontsize, pad=15)
        ax3.legend(fontsize=legend_fontsize)
        ax3.grid(True, alpha=0.3)
        ax3.tick_params(labelsize=12)
        # Add subplot label C
        ax3.text(0.02, 0.98, 'C', transform=ax3.transAxes, fontsize=subplot_label_fontsize,
                fontweight='bold', va='top', ha='left', color='black',
                bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8))
        
        # 4. Background-corrected signal (LABEL: D)
        ax4 = axes[1, 1]
        ax4.plot(analysis['position'], analysis['corrected'], 
                color='green', linewidth=3)
        ax4.fill_between(analysis['position'], analysis['corrected'], 
                        alpha=0.3, color='green')
        ax4.set_xlabel('Position (pixels)', fontsize=label_fontsize)
        ax4.set_ylabel('Corrected Signal Intensity', fontsize=label_fontsize)
        max_val = analysis['corrected'].max()
        mean_val = analysis['corrected'].mean()
        ax4.set_title(f'Background-Corrected Signal\nPeak: {max_val:.0f}, Mean: {mean_val:.0f}', 
                     fontsize=title_fontsize, pad=15)
        ax4.grid(True, alpha=0.3)
        ax4.tick_params(labelsize=12)
        # Add subplot label D
        ax4.text(0.02, 0.98, 'D', transform=ax4.transAxes, fontsize=subplot_label_fontsize,
                fontweight='bold', va='top', ha='left', color='black',
                bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8))
        
        plt.tight_layout(pad=3.0)
        return fig


class GelAnalysisGUI:
    """Gel Electrophoresis Image Analysis GUI Main Window"""
    
    def __init__(self, root, image_data, filename):
        self.root = root
        
        # Sanitize filename for window title (CRITICAL)
        safe_filename = sanitize_for_gui(filename)
        self.root.title(f"Interactive Gel Image Analysis - {safe_filename}")
        
        self.image_data = image_data
        self.filename = safe_filename
        self.analyzer = InteractiveGelAnalyzer(image_data, filename)
        
        # 矩形选择相关
        self.rect_start = None
        self.rect_end = None
        self.rect_id = None
        self.is_drawing = False
        
        # 创建UI
        self.create_widgets()
        
    def create_widgets(self):
        """Create GUI components"""
        # Top control panel
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Instruction labels (ALL ASCII)
        ttk.Label(control_frame, text="INSTRUCTIONS:", font=('Arial', 12, 'bold')).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
        ttk.Label(control_frame, text="1. Click and drag on image below to draw a rectangle").grid(row=1, column=0, columnspan=2, sticky=tk.W)
        ttk.Label(control_frame, text="2. Release mouse, then click 'Analyze Selection' button").grid(row=2, column=0, columnspan=2, sticky=tk.W)
        ttk.Label(control_frame, text="3. Or manually enter coordinates and click 'Analyze'").grid(row=3, column=0, columnspan=2, sticky=tk.W)
        
        # Coordinate input boxes
        ttk.Label(control_frame, text="Coordinates:", font=('Arial', 11, 'bold')).grid(row=4, column=0, columnspan=2, sticky=tk.W, pady=(10, 5))
        
        coord_frame = ttk.Frame(control_frame)
        coord_frame.grid(row=5, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(coord_frame, text="X1:").grid(row=0, column=0, padx=5)
        self.x1_var = tk.IntVar(value=100)
        self.x1_entry = ttk.Entry(coord_frame, textvariable=self.x1_var, width=8)
        self.x1_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(coord_frame, text="Y1:").grid(row=0, column=2, padx=5)
        self.y1_var = tk.IntVar(value=50)
        self.y1_entry = ttk.Entry(coord_frame, textvariable=self.y1_var, width=8)
        self.y1_entry.grid(row=0, column=3, padx=5)
        
        ttk.Label(coord_frame, text="X2:").grid(row=1, column=0, padx=5)
        self.x2_var = tk.IntVar(value=140)
        self.x2_entry = ttk.Entry(coord_frame, textvariable=self.x2_var, width=8)
        self.x2_entry.grid(row=1, column=1, padx=5)
        
        ttk.Label(coord_frame, text="Y2:").grid(row=1, column=2, padx=5)
        self.y2_var = tk.IntVar(value=150)
        self.y2_entry = ttk.Entry(coord_frame, textvariable=self.y2_var, width=8)
        self.y2_entry.grid(row=1, column=3, padx=5)
        
        # Buttons
        button_frame = ttk.Frame(control_frame)
        button_frame.grid(row=6, column=0, columnspan=2, pady=10)
        
        self.analyze_btn = ttk.Button(button_frame, text="Analyze Selection", command=self.analyze_selection)
        self.analyze_btn.grid(row=0, column=0, padx=5)
        
        self.clear_btn = ttk.Button(button_frame, text="Clear Selection", command=self.clear_selection)
        self.clear_btn.grid(row=0, column=1, padx=5)
        
        # Status label (ALWAYS SANITIZED)
        self.status_label = ttk.Label(control_frame, text="Ready", foreground="green")
        self.status_label.grid(row=7, column=0, columnspan=2, pady=5)
        
        # Image display area
        image_frame = ttk.Frame(self.root)
        image_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create image canvas
        self.fig = Figure(figsize=(10, 7))
        self.ax = self.fig.add_subplot(111)
        
        # Display image
        self.ax.imshow(self.image_data, cmap='gray')
        self.ax.set_title(f'{self.filename}\nImage size: {self.image_data.shape[1]} x {self.image_data.shape[0]} pixels', 
                         fontsize=12)
        self.ax.axis('off')
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=image_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # 绑定鼠标事件
        self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        
        # 配置网格权重
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
    def on_mouse_press(self, event):
        """Mouse press event"""
        if event.inaxes != self.ax:
            return
        
        self.rect_start = (event.xdata, event.ydata)
        self.is_drawing = True
        self.status_label.config(text="Drawing rectangle...", foreground="blue")
        
    def on_mouse_move(self, event):
        """Mouse move event"""
        if not self.is_drawing or event.inaxes != self.ax:
            return
        
        # 删除旧矩形
        if self.rect_id:
            self.rect_id.remove()
        
        # 绘制新矩形
        x1, y1 = self.rect_start
        x2, y2 = event.xdata, event.ydata
        
        from matplotlib.patches import Rectangle
        width = x2 - x1
        height = y2 - y1
        self.rect_id = Rectangle((x1, y1), width, height,
                                 fill=False, edgecolor='red', linewidth=2)
        self.ax.add_patch(self.rect_id)
        self.canvas.draw()
        
    def on_mouse_release(self, event):
        """Mouse release event"""
        if not self.is_drawing:
            return
        
        self.is_drawing = False
        
        if event.inaxes != self.ax or not self.rect_start:
            return
        
        self.rect_end = (event.xdata, event.ydata)
        
        # Update coordinate input boxes
        x1, y1 = self.rect_start
        x2, y2 = self.rect_end
        
        self.x1_var.set(int(min(x1, x2)))
        self.y1_var.set(int(min(y1, y2)))
        self.x2_var.set(int(max(x1, x2)))
        self.y2_var.set(int(max(y1, y2)))
        
        self.status_label.config(
            text=f"Selected region: ({int(min(x1, x2))}, {int(min(y1, y2))}) to ({int(max(x1, x2))}, {int(max(y1, y2))})",
            foreground="green"
        )
        
    def clear_selection(self):
        """Clear selection"""
        if self.rect_id:
            self.rect_id.remove()
            self.rect_id = None
            self.canvas.draw()
        
        self.rect_start = None
        self.rect_end = None
        self.status_label.config(text="Selection cleared", foreground="orange")
        
    def analyze_selection(self):
        """Analyze selected region"""
        try:
            x1 = self.x1_var.get()
            y1 = self.y1_var.get()
            x2 = self.x2_var.get()
            y2 = self.y2_var.get()
            
            # Validate coordinates
            if x2 <= x1 or y2 <= y1:
                messagebox.showerror("Error", "Bottom-right coordinate must be greater than top-left!")
                return
            
            if x2 - x1 < 5 or y2 - y1 < 5:
                messagebox.showerror("Error", "Rectangle too small, please select at least 5x5 pixels!")
                return
            
            self.status_label.config(text="Analyzing...", foreground="blue")
            self.root.update()
            
            # Analyze
            analysis = self.analyzer.analyze_panel(x1, y1, x2, y2)
            
            # Display results - create new window with results
            self.show_analysis_results(analysis)
            
            # Display statistics in status label (SANITIZED)
            peak = analysis['corrected'].max()
            mean = analysis['corrected'].mean()
            direction = sanitize_for_gui(analysis['direction'])
            self.status_label.config(
                text=f"Done! Direction: {direction}, Peak: {peak:.0f}, Mean: {mean:.0f}",
                foreground="green"
            )
            
        except Exception as e:
            error_msg = sanitize_for_gui(str(e))
            messagebox.showerror("Analysis Error", f"Error during analysis:\n{error_msg}")
            self.status_label.config(text="Analysis failed", foreground="red")
    
    def show_analysis_results(self, analysis):
        """Show analysis results in a new window"""
        # Create new top-level window
        result_window = tk.Toplevel(self.root)
        result_window.title("Analysis Results")
        result_window.geometry("1400x1000")
        
        # Create matplotlib figure
        fig = self.analyzer.plot_results(analysis)
        
        # Embed figure in Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=result_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add close button
        close_btn = ttk.Button(result_window, text="Close", 
                              command=result_window.destroy)
        close_btn.pack(pady=10)
        
        # Bring window to front
        result_window.lift()
        result_window.focus_force()


def main():
    """Main function"""
    print("=" * 80)
    print("Interactive Gel Electrophoresis Image Analysis - Tkinter GUI")
    print("=" * 80)
    
    # Find SCN files
    scn_files = sorted(glob.glob("*.scn"))
    
    if not scn_files:
        print("\nERROR: No SCN files found!")
        print("Please ensure you are running this program in a directory containing .scn files")
        sys.exit(1)
    
    print(f"\nFound {len(scn_files)} SCN files:")
    for i, f in enumerate(scn_files, 1):
        # Sanitize filename for console output
        safe_name = sanitize_for_gui(f)
        print(f"  {i}. {safe_name}")
    
    # Select file
    if len(scn_files) == 1:
        selected_file = scn_files[0]
        safe_name = sanitize_for_gui(selected_file)
        print(f"\nAuto-selected: {safe_name}")
    else:
        while True:
            try:
                choice = input(f"\nSelect file number (1-{len(scn_files)}): ")
                idx = int(choice) - 1
                if 0 <= idx < len(scn_files):
                    selected_file = scn_files[idx]
                    break
                else:
                    print("Invalid number, please try again")
            except ValueError:
                print("Please enter a number")
    
    # Read image
    safe_filename = sanitize_for_gui(selected_file)
    print(f"\nReading: {safe_filename}")
    img_data, dimensions, error = read_scn_file_enhanced(selected_file, prefer_mode='crop')
    
    if error:
        safe_error = sanitize_for_gui(error)
        print(f"ERROR: Read failed: {safe_error}")
        sys.exit(1)
    
    print(f"Success! Image size: {img_data.shape}")
    print("\nLaunching GUI window...")
    
    # Create GUI
    root = tk.Tk()
    app = GelAnalysisGUI(root, img_data, selected_file)
    
    print("GUI window opened")
    print("\nTips:")
    print("  - Drag mouse on image to select rectangular region")
    print("  - Or manually enter coordinates")
    print("  - Click 'Analyze Selection' button to view results")
    print("=" * 80)
    
    root.mainloop()


if __name__ == "__main__":
    main()
