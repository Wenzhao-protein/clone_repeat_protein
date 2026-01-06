#!/usr/bin/env python3
"""
Enhanced Interactive Gel Electrophoresis Image Analysis - Tkinter GUI

NEW FEATURES:
- Image rotation (±90° and ±1° buttons with cumulative tracking)
- Multiple selection modes: Panel (red), Plasmid Band (green), Fragment Band (blue)
- Linear background subtraction with percentile-based robust fitting
- Panel-relative band intensity calculations
- Fragment ratio computation
- 5 subplots with A-E labels (including original image without rectangles)
- Larger fonts (14-24pt) and window sizes (1800x1400, 2000x1600)
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use('TkAgg', force=True)  # Force TkAgg backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from scipy.ndimage import rotate as scipy_rotate
from numpy.polynomial import Polynomial
import sys
import os
import glob
import re
from pathlib import Path

# Import scn_reader
from scn_reader import read_scn_file, get_image_info


# ============================================================
# SANITIZATION UTILITY
# ============================================================

def sanitize_for_gui(text):
    """Sanitize text for GUI display - ONLY ASCII output."""
    if text is None:
        return "N/A"
    
    if not isinstance(text, str):
        text = str(text)
    
    # Remove unicode escape sequences
    text = re.sub(r'\\[uU][0-9a-fA-F]{4,8}', '', text)
    
    # Remove replacement character
    text = text.replace('\ufffd', '')
    text = text.replace('�', '')
    text = text.replace('\\ufffd', '')
    
    # Convert non-ASCII to '?'
    cleaned = []
    for char in text:
        code = ord(char)
        if 0x20 <= code <= 0x7E:
            cleaned.append(char)
        elif code == 0x0A:
            cleaned.append(char)
        elif code > 0x7E:
            cleaned.append('?')
    
    result = ''.join(cleaned).strip()
    
    if not result or all(c == '?' for c in result):
        return "UNKNOWN"
    
    result = result.replace('???', '')
    result = re.sub(r'\s+', ' ', result)
    
    return result.strip() or "N/A"


# ============================================================
# INTERACTIVE GEL ANALYZER CLASS
# ============================================================

class InteractiveGelAnalyzer:
    """Enhanced Interactive Gel Electrophoresis Image Analyzer"""
    
    def __init__(self, image, filename, filepath=None):
        self.original_image = image.copy()
        self.image = image
        self.filename = sanitize_for_gui(filename)
        self.filepath = filepath  # Store full file path for saving
        self.rotation_angle = 0
        
    def rotate_image(self, angle):
        """Rotate image by angle degrees (cumulative)"""
        self.rotation_angle = (self.rotation_angle + angle) % 360
        self.image = scipy_rotate(self.image, angle, reshape=False, order=1)
        return self.image
    
    def reset_rotation(self):
        """Reset image to original (no rotation)"""
        self.image = self.original_image.copy()
        self.rotation_angle = 0
        return self.image
        
    def analyze_panel(self, panel_coords, plasmid_coords=None, fragment_coords=None):
        """
        Analyze selected panel with optional plasmid and fragment bands.
        
        Returns dictionary with analysis results including:
        - Panel image, signal profile, background, corrected signal
        - Plasmid and fragment intensities (if coords provided)
        - Fragment ratio
        """
        x1, y1, x2, y2 = panel_coords
        x1, x2 = int(min(x1, x2)), int(max(x1, x2))
        y1, y2 = int(min(y1, y2)), int(max(y1, y2))
        
        # Extract panel region
        panel = self.image[y1:y2, x1:x2]
        height, width = panel.shape
        is_vertical = height > width
        
        # Calculate average signal along long edge
        if is_vertical:
            profile = np.mean(panel, axis=1)
            position = np.arange(len(profile))
            direction = "Vertical (Top to Bottom)"
        else:
            profile = np.mean(panel, axis=0)
            position = np.arange(len(profile))
            direction = "Horizontal (Left to Right)"
        
        # LINEAR background estimation with percentile-based robust fitting
        # Use lower percentile values to exclude signal peaks
        percentile_threshold = 80
        low_values_mask = profile < np.percentile(profile, percentile_threshold)
        
        if np.sum(low_values_mask) > len(profile) // 10:  # At least 10% points
            # Fit linear polynomial to low-intensity points
            p = Polynomial.fit(position[low_values_mask], profile[low_values_mask], deg=1)
        else:
            # Fallback: fit to all points
            p = Polynomial.fit(position, profile, deg=1)
        
        background = p(position)
        
        # Remove background
        corrected = profile - background
        corrected = np.maximum(corrected, 0)
        
        # Calculate band intensities if coords provided
        plasmid_intensity = 0
        fragment_intensity = 0
        plasmid_window = None
        fragment_window = None
        
        if plasmid_coords:
            # Map plasmid coordinates relative to panel
            px1, py1, px2, py2 = plasmid_coords
            px1, px2 = int(min(px1, px2)), int(max(px1, px2))
            py1, py2 = int(min(py1, py2)), int(max(py1, py2))
            
            # Calculate relative position within panel
            if is_vertical:
                rel_start = max(0, py1 - y1)
                rel_end = min(len(corrected), py2 - y1)
                plasmid_window = (rel_start, rel_end)
                if rel_start < rel_end:
                    plasmid_intensity = np.sum(corrected[rel_start:rel_end])
            else:
                rel_start = max(0, px1 - x1)
                rel_end = min(len(corrected), px2 - x1)
                plasmid_window = (rel_start, rel_end)
                if rel_start < rel_end:
                    plasmid_intensity = np.sum(corrected[rel_start:rel_end])
        
        if fragment_coords:
            # Map fragment coordinates relative to panel
            fx1, fy1, fx2, fy2 = fragment_coords
            fx1, fx2 = int(min(fx1, fx2)), int(max(fx1, fx2))
            fy1, fy2 = int(min(fy1, fy2)), int(max(fy1, fy2))
            
            # Calculate relative position within panel
            if is_vertical:
                rel_start = max(0, fy1 - y1)
                rel_end = min(len(corrected), fy2 - y1)
                fragment_window = (rel_start, rel_end)
                if rel_start < rel_end:
                    fragment_intensity = np.sum(corrected[rel_start:rel_end])
            else:
                rel_start = max(0, fx1 - x1)
                rel_end = min(len(corrected), fx2 - x1)
                fragment_window = (rel_start, rel_end)
                if rel_start < rel_end:
                    fragment_intensity = np.sum(corrected[rel_start:rel_end])
        
        # Calculate fragment ratio
        total_intensity = plasmid_intensity + fragment_intensity
        fragment_ratio = fragment_intensity / total_intensity if total_intensity > 0 else 0
        
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
            'fragment_coords': fragment_coords,
            'plasmid_window': plasmid_window,
            'fragment_window': fragment_window
        }
    
    def save_figure(self, fig, analysis):
        """Save figure to folder based on SCN filename and selection coordinates"""
        try:
            # Get SCN file path and create output folder
            scn_path = Path(self.filepath)
            scn_name = scn_path.stem  # filename without extension
            output_folder = scn_path.parent / scn_name
            
            # Create folder if not exists (不覆盖已有文件夹)
            output_folder.mkdir(exist_ok=True)
            
            # Get coordinates for filename
            x1, y1, x2, y2 = analysis['coords']
            filename_base = f"{int(x1)}_{int(x2)}"
            
            # Save as PNG (900 dpi)
            png_path = output_folder / f"{filename_base}.png"
            fig.savefig(png_path, dpi=900, bbox_inches='tight')
            print(f"Saved PNG: {png_path}")
            
            # Save as PDF
            pdf_path = output_folder / f"{filename_base}.pdf"
            fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
            print(f"Saved PDF: {pdf_path}")
            
        except Exception as e:
            print(f"Error saving figure: {e}")
    
    def plot_results(self, analysis):
        """Plot analysis results in 5 subplots with A-E labels
        Layout: A (row 0), B (row 1), C/D/E (row 2)
        使用add_axes明确定义每个subplot的位置
        """
        fig = plt.figure(figsize=(20, 16))
        
        # Font sizes
        title_fs = 16
        label_fs = 14
        legend_fs = 12
        sublabel_fs = 24
        
        # 明确定义每个subplot的位置 [left, bottom, width, height]
        # 所有值在0-1范围内，相对于figure
        
        # === SUBPLOT A: Original image (顶部，占整行) ===
        ax_a = fig.add_axes([0.08, 0.70, 0.84, 0.25])  # 更高
        ax_a.imshow(self.image, cmap='gray')
        ax_a.set_title(f'Original Image\n{self.filename}', fontsize=title_fs, pad=15)
        ax_a.axis('off')
        # Add label A
        ax_a.text(-0.03, 0.98, 'A', transform=ax_a.transAxes, fontsize=sublabel_fs,
                 fontweight='bold', va='top', ha='right', color='black')
        
        # === SUBPLOT B: Image with selections (中部，占整行) ===
        ax_b = fig.add_axes([0.08, 0.40, 0.84, 0.25])  # 更高
        ax_b.imshow(self.image, cmap='gray')
        x1, y1, x2, y2 = analysis['coords']
        
        from matplotlib.patches import Rectangle
        # Panel (RED)
        rect_panel = Rectangle((x1, y1), x2-x1, y2-y1,
                               fill=False, edgecolor='red', linewidth=3, label='Selected Panel')
        ax_b.add_patch(rect_panel)
        
        # Debug: print all relevant info
        print(f"\nDEBUG plot_results:")
        print(f"  Panel coords: ({x1}, {y1}) to ({x2}, {y2})")
        print(f"  is_vertical: {analysis['is_vertical']}")
        print(f"  plasmid_window: {analysis.get('plasmid_window')}")
        print(f"  fragment_window: {analysis.get('fragment_window')}")
        
        # Plasmid (GREEN) - 只框选实际参与计算的部分，用虚线
        if analysis.get('plasmid_window'):
            rel_start, rel_end = analysis['plasmid_window']
            print(f"  Plasmid: rel_start={rel_start}, rel_end={rel_end}")
            if rel_end > rel_start:
                if analysis['is_vertical']:
                    rect_plasmid = Rectangle((x1, y1 + rel_start), x2-x1, rel_end-rel_start,
                                             fill=False, edgecolor='green', linewidth=3, 
                                             linestyle='--', label='Plasmid')
                    print(f"    Drawing plasmid (vertical): ({x1}, {y1 + rel_start}), width={x2-x1}, height={rel_end-rel_start}")
                else:
                    rect_plasmid = Rectangle((x1 + rel_start, y1), rel_end-rel_start, y2-y1,
                                             fill=False, edgecolor='green', linewidth=3, 
                                             linestyle='--', label='Plasmid')
                    print(f"    Drawing plasmid (horizontal): ({x1 + rel_start}, {y1}), width={rel_end-rel_start}, height={y2-y1}")
                ax_b.add_patch(rect_plasmid)
            else:
                print(f"    Skipping plasmid (invalid window)")
        
        # Fragment (BLUE) - 只框选实际参与计算的部分，用虚线
        if analysis.get('fragment_window'):
            rel_start, rel_end = analysis['fragment_window']
            print(f"  Fragment: rel_start={rel_start}, rel_end={rel_end}")
            if rel_end > rel_start:
                if analysis['is_vertical']:
                    rect_fragment = Rectangle((x1, y1 + rel_start), x2-x1, rel_end-rel_start,
                                              fill=False, edgecolor='blue', linewidth=4, 
                                              linestyle='--', label='Fragment')
                    print(f"    Drawing fragment (vertical): ({x1}, {y1 + rel_start}), width={x2-x1}, height={rel_end-rel_start}")
                else:
                    rect_fragment = Rectangle((x1 + rel_start, y1), rel_end-rel_start, y2-y1,
                                              fill=False, edgecolor='blue', linewidth=4, 
                                              linestyle='--', label='Fragment')
                    print(f"    Drawing fragment (horizontal): ({x1 + rel_start}, {y1}), width={rel_end-rel_start}, height={y2-y1}")
                ax_b.add_patch(rect_fragment)
            else:
                print(f"    Skipping fragment (invalid window)")
        
        ax_b.set_title('Image with Selections', fontsize=title_fs, pad=15)
        ax_b.axis('off')
        ax_b.legend(fontsize=legend_fs, loc='lower right', framealpha=0.9)
        # Add label B
        ax_b.text(-0.03, 0.98, 'B', transform=ax_b.transAxes, fontsize=sublabel_fs,
                 fontweight='bold', va='top', ha='right', color='black')
        
        # === SUBPLOT C: Selected panel (底部左，调整位置缩小CD间距) ===
        ax_c = fig.add_axes([0.13, 0.05, 0.36, 0.28])
        ax_c.imshow(analysis['panel'], cmap='gray')
        ax_c.set_title(f'Selected Panel Region\n{analysis["direction"]}', fontsize=title_fs, pad=15)
        ax_c.axis('off')
        # Add label C
        ax_c.text(-0.08, 0.98, 'C', transform=ax_c.transAxes, fontsize=sublabel_fs,
                 fontweight='bold', va='top', ha='right', color='black')
        
        # === SUBPLOT D: Raw signal (底部中，位置不变) ===
        ax_d = fig.add_axes([0.50, 0.05, 0.18, 0.28])
        ax_d.plot(analysis['profile'], analysis['position'],
                 label='Raw Signal', color='blue', linewidth=2.5)
        ax_d.plot(analysis['background'], analysis['position'],
                 label='Linear Background', color='red', linewidth=2.5, linestyle='--')
        ax_d.set_ylabel('Position (pixels)', fontsize=label_fs)
        ax_d.set_xlabel('Signal Intensity', fontsize=label_fs)
        ax_d.set_title('Raw Signal vs Linear Background', fontsize=title_fs, pad=15)
        ax_d.legend(fontsize=legend_fs, loc='lower right')
        ax_d.tick_params(labelsize=12)
        ax_d.invert_yaxis()
        # Center axis labels
        ax_d.xaxis.set_label_coords(0.5, -0.08)
        ax_d.yaxis.set_label_coords(-0.12, 0.5)
        # Add label D
        ax_d.text(-0.15, 0.98, 'D', transform=ax_d.transAxes, fontsize=sublabel_fs,
                 fontweight='bold', va='top', ha='right', color='black')
        
        # === SUBPLOT E: Corrected signal (底部右，往左移，增加与D的间距) ===
        ax_e = fig.add_axes([0.74, 0.05, 0.18, 0.28])
        ax_e.plot(analysis['corrected'], analysis['position'],
                 color='green', linewidth=3, label='Corrected Signal')
        
        # Highlight plasmid and fragment windows
        if analysis.get('plasmid_window'):
            p_start, p_end = analysis['plasmid_window']
            ax_e.axhspan(p_start, p_end, color='green', alpha=0.2, label='Plasmid Region')
        
        if analysis.get('fragment_window'):
            f_start, f_end = analysis['fragment_window']
            ax_e.axhspan(f_start, f_end, color='blue', alpha=0.2, label='Fragment Region')
        
        ax_e.set_ylabel('Position (pixels)', fontsize=label_fs)
        ax_e.set_xlabel('Corrected Signal Intensity', fontsize=label_fs)
        
        # Build title with ratio as third line
        title_text = 'Background-Corrected Signal'
        
        if analysis['plasmid_intensity'] > 0 or analysis['fragment_intensity'] > 0:
            title_text += f'\nPlasmid: {analysis["plasmid_intensity"]:.0f}, Fragment: {analysis["fragment_intensity"]:.0f}'
        
        if analysis['fragment_ratio'] > 0:
            title_text += f'\nRatio: {analysis["fragment_ratio"]:.2f}'
        
        ax_e.set_title(title_text, fontsize=title_fs, pad=15)
        
        ax_e.legend(fontsize=legend_fs, loc='lower right')
        ax_e.tick_params(labelsize=12)
        ax_e.invert_yaxis()
        # Center axis labels
        ax_e.xaxis.set_label_coords(0.5, -0.08)
        ax_e.yaxis.set_label_coords(-0.12, 0.5)
        # Add label E
        ax_e.text(-0.15, 0.98, 'E', transform=ax_e.transAxes, fontsize=sublabel_fs,
                 fontweight='bold', va='top', ha='right', color='black')
        
        # Save figure if filename is available
        if hasattr(self, 'filepath') and self.filepath:
            self.save_figure(fig, analysis)
        
        return fig


# ============================================================
# GUI CLASS
# ============================================================

class GelAnalysisGUI:
    """Enhanced Gel Image Analysis GUI with multiple selection modes"""
    
    def __init__(self, root, image_data, filename, filepath=None):
        self.root = root
        
        safe_filename = sanitize_for_gui(filename)
        self.root.title(f"Enhanced Gel Analysis - {safe_filename}")
        
        # Adjusted window size to fit most screens
        self.root.geometry("1400x1000")
        
        self.image_data = image_data
        self.filename = safe_filename
        self.filepath = filepath  # Store full file path
        self.analyzer = InteractiveGelAnalyzer(image_data, filename, filepath)
        
        # Selection state
        self.rect_start = None
        self.rect_end = None
        self.is_drawing = False
        
        # Multiple rectangles for different band types
        self.selection_mode = tk.StringVar(value="panel")
        self.rectangles = {
            'panel': {'rect_id': None, 'coords': None, 'color': 'red'},
            'plasmid': {'rect_id': None, 'coords': None, 'color': 'green'},
            'fragment': {'rect_id': None, 'coords': None, 'color': 'blue'}
        }
        
        # Create UI
        self.create_widgets()
    
    def create_widgets(self):
        """Create enhanced GUI components"""
        # Main container
        main_container = ttk.Frame(self.root)
        main_container.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_container.columnconfigure(0, weight=1)
        main_container.rowconfigure(1, weight=1)
        
        # === TOP CONTROL PANEL ===
        control_frame = ttk.Frame(main_container, padding="15")
        control_frame.grid(row=0, column=0, sticky=(tk.W, tk.E))
        
        # Style configuration for larger fonts
        style = ttk.Style()
        style.configure('Large.TLabel', font=('Arial', 14))
        style.configure('Title.TLabel', font=('Arial', 16, 'bold'))
        style.configure('Large.TButton', font=('Arial', 13))
        style.configure('Large.TRadiobutton', font=('Arial', 13))
        
        # Title
        ttk.Label(control_frame, text="ENHANCED GEL IMAGE ANALYSIS",
                 style='Title.TLabel').grid(row=0, column=0, columnspan=6, sticky=tk.W, pady=(0,10))
        
        # Instructions
        instructions = [
            "1. Use rotation buttons to adjust image orientation (90° or 1° increments)",
            "2. Select mode: Panel (analysis region), Plasmid (reference), or Fragment (target)",
            "3. Click and drag on image to draw rectangles for each mode",
            "4. Click 'Analyze Panel' to process and view 5-subplot results"
        ]
        for i, text in enumerate(instructions, 1):
            ttk.Label(control_frame, text=text,
                     style='Large.TLabel').grid(row=i, column=0, columnspan=6, sticky=tk.W)
        
        # === ROTATION CONTROLS ===
        rotation_frame = ttk.LabelFrame(control_frame, text="Image Rotation", padding="10")
        rotation_frame.grid(row=5, column=0, columnspan=6, sticky=(tk.W, tk.E), pady=(15,10))
        
        ttk.Button(rotation_frame, text="◄ 90°",
                  command=lambda: self.rotate_image(-90),
                  style='Large.TButton').grid(row=0, column=0, padx=5)
        ttk.Button(rotation_frame, text="◄ 1°",
                  command=lambda: self.rotate_image(-1),
                  style='Large.TButton').grid(row=0, column=1, padx=5)
        ttk.Button(rotation_frame, text="Reset",
                  command=self.reset_rotation,
                  style='Large.TButton').grid(row=0, column=2, padx=5)
        ttk.Button(rotation_frame, text="1° ►",
                  command=lambda: self.rotate_image(1),
                  style='Large.TButton').grid(row=0, column=3, padx=5)
        ttk.Button(rotation_frame, text="90° ►",
                  command=lambda: self.rotate_image(90),
                  style='Large.TButton').grid(row=0, column=4, padx=5)
        
        # === SELECTION MODE ===
        mode_frame = ttk.LabelFrame(control_frame, text="Selection Mode", padding="10")
        mode_frame.grid(row=6, column=0, columnspan=6, sticky=(tk.W, tk.E), pady=10)
        
        ttk.Radiobutton(mode_frame, text="Panel (Red - Analysis Region)",
                       variable=self.selection_mode, value="panel",
                       style='Large.TRadiobutton').grid(row=0, column=0, padx=10, sticky=tk.W)
        ttk.Radiobutton(mode_frame, text="Plasmid Band (Green - Reference)",
                       variable=self.selection_mode, value="plasmid",
                       style='Large.TRadiobutton').grid(row=0, column=1, padx=10, sticky=tk.W)
        ttk.Radiobutton(mode_frame, text="Fragment Band (Blue - Target)",
                       variable=self.selection_mode, value="fragment",
                       style='Large.TRadiobutton').grid(row=0, column=2, padx=10, sticky=tk.W)
        
        # === ACTION BUTTONS ===
        button_frame = ttk.Frame(control_frame, padding="10")
        button_frame.grid(row=7, column=0, columnspan=6, pady=10)
        
        ttk.Button(button_frame, text="Analyze Panel",
                  command=self.analyze_selection,
                  style='Large.TButton').grid(row=0, column=0, padx=10)
        ttk.Button(button_frame, text="Clear All Selections",
                  command=self.clear_all_selections,
                  style='Large.TButton').grid(row=0, column=1, padx=10)
        ttk.Button(button_frame, text="Clear Current Mode",
                  command=self.clear_current_selection,
                  style='Large.TButton').grid(row=0, column=2, padx=10)
        
        # === STATUS LABEL ===
        self.status_label = ttk.Label(control_frame, text="Ready - Select Panel region to begin",
                                     style='Large.TLabel', foreground="green")
        self.status_label.grid(row=8, column=0, columnspan=6, pady=10)
        
        # === IMAGE DISPLAY AREA ===
        image_frame = ttk.Frame(main_container)
        image_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create matplotlib figure
        self.fig = Figure(figsize=(14, 10))
        self.ax = self.fig.add_subplot(111)
        
        # Display image
        self.update_image_display()
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=image_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Bind mouse events
        self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
    
    def update_image_display(self):
        """Update the displayed image"""
        self.ax.clear()
        self.ax.imshow(self.analyzer.image, cmap='gray')
        self.ax.set_title(
            f'{self.filename}\n' +
            f'Image: {self.analyzer.image.shape[1]} x {self.analyzer.image.shape[0]} px  |  ' +
            f'Rotation: {self.analyzer.rotation_angle}°',
            fontsize=16, pad=15
        )
        self.ax.axis('off')
        
        # Redraw all existing rectangles
        self.redraw_rectangles()
    
    def redraw_rectangles(self):
        """Redraw all existing rectangles on the image"""
        from matplotlib.patches import Rectangle
        for mode, data in self.rectangles.items():
            if data['coords']:
                x1, y1, x2, y2 = data['coords']
                rect = Rectangle((x1, y1), x2-x1, y2-y1,
                               fill=False, edgecolor=data['color'],
                               linewidth=4, label=mode.capitalize())
                self.ax.add_patch(rect)
        
        if any(data['coords'] for data in self.rectangles.values()):
            self.ax.legend(fontsize=14, loc='upper right', framealpha=0.9)
    
    def rotate_image(self, angle):
        """Rotate the image"""
        self.analyzer.rotate_image(angle)
        self.update_image_display()
        self.canvas.draw()
        self.status_label.config(
            text=f"Image rotated by {angle}° | Total rotation: {self.analyzer.rotation_angle}°",
            foreground="blue"
        )
    
    def reset_rotation(self):
        """Reset image rotation"""
        self.analyzer.reset_rotation()
        # Clear all rectangles since they're now invalid
        for mode in self.rectangles:
            self.rectangles[mode]['rect_id'] = None
            self.rectangles[mode]['coords'] = None
        self.update_image_display()
        self.canvas.draw()
        self.status_label.config(text="Rotation reset - Please reselect regions", foreground="orange")
    
    def on_mouse_press(self, event):
        """Handle mouse press event"""
        if event.inaxes != self.ax:
            return
        
        self.rect_start = (event.xdata, event.ydata)
        self.is_drawing = True
        mode = self.selection_mode.get()
        self.status_label.config(text=f"Drawing {mode} rectangle...", foreground="blue")
    
    def on_mouse_move(self, event):
        """Handle mouse move event"""
        if not self.is_drawing or event.inaxes != self.ax:
            return
        
        # Remove temporary rectangle if exists
        mode = self.selection_mode.get()
        if self.rectangles[mode]['rect_id']:
            self.rectangles[mode]['rect_id'].remove()
        
        # Draw new temporary rectangle
        x1, y1 = self.rect_start
        x2, y2 = event.xdata, event.ydata
        
        from matplotlib.patches import Rectangle
        color = self.rectangles[mode]['color']
        rect = Rectangle((x1, y1), x2-x1, y2-y1,
                        fill=False, edgecolor=color, linewidth=4)
        self.rectangles[mode]['rect_id'] = self.ax.add_patch(rect)
        self.canvas.draw()
    
    def on_mouse_release(self, event):
        """Handle mouse release event"""
        if not self.is_drawing:
            return
        
        self.is_drawing = False
        
        if event.inaxes != self.ax or not self.rect_start:
            return
        
        # Save coordinates
        mode = self.selection_mode.get()
        x1, y1 = self.rect_start
        x2, y2 = event.xdata, event.ydata
        
        coords = (int(min(x1, x2)), int(min(y1, y2)),
                 int(max(x1, x2)), int(max(y1, y2)))
        self.rectangles[mode]['coords'] = coords
        
        self.status_label.config(
            text=f"{mode.capitalize()} selected: ({coords[0]}, {coords[1]}) to ({coords[2]}, {coords[3]})",
            foreground="green"
        )
    
    def clear_current_selection(self):
        """Clear the current selection mode"""
        mode = self.selection_mode.get()
        if self.rectangles[mode]['rect_id']:
            self.rectangles[mode]['rect_id'].remove()
        self.rectangles[mode]['rect_id'] = None
        self.rectangles[mode]['coords'] = None
        self.update_image_display()
        self.canvas.draw()
        self.status_label.config(text=f"{mode.capitalize()} selection cleared", foreground="orange")
    
    def clear_all_selections(self):
        """Clear all selections"""
        for mode in self.rectangles:
            self.rectangles[mode]['rect_id'] = None
            self.rectangles[mode]['coords'] = None
        self.update_image_display()
        self.canvas.draw()
        self.status_label.config(text="All selections cleared", foreground="orange")
    
    def analyze_selection(self):
        """Analyze the selected panel"""
        try:
            # Check if panel is selected
            if not self.rectangles['panel']['coords']:
                messagebox.showerror("Error", "Please select a Panel region first!")
                return
            
            panel_coords = self.rectangles['panel']['coords']
            plasmid_coords = self.rectangles['plasmid']['coords']
            fragment_coords = self.rectangles['fragment']['coords']
            
            # Validate panel size
            x1, y1, x2, y2 = panel_coords
            if x2 - x1 < 5 or y2 - y1 < 5:
                messagebox.showerror("Error", "Panel region too small (minimum 5x5 pixels)!")
                return
            
            self.status_label.config(text="Analyzing...", foreground="blue")
            self.root.update()
            
            # Perform analysis
            analysis = self.analyzer.analyze_panel(panel_coords, plasmid_coords, fragment_coords)
            
            # Show results
            self.show_analysis_results(analysis)
            
            # Update status
            peak = analysis['corrected'].max()
            mean = analysis['corrected'].mean()
            frag_ratio = analysis['fragment_ratio']
            
            status_text = f"Analysis complete! Peak: {peak:.0f}, Mean: {mean:.0f}"
            if frag_ratio > 0:
                status_text += f", Fragment Ratio: {frag_ratio:.4f}"
            
            self.status_label.config(text=status_text, foreground="green")
        
        except Exception as e:
            error_msg = sanitize_for_gui(str(e))
            messagebox.showerror("Analysis Error", f"Error during analysis:\n{error_msg}")
            self.status_label.config(text="Analysis failed", foreground="red")
    
    def show_analysis_results(self, analysis):
        """Display analysis results in new window with 5 subplots"""
        result_window = tk.Toplevel(self.root)
        result_window.title("Analysis Results - 5 Subplots")
        result_window.geometry("2000x1600")
        
        # Create figure with 5 subplots
        fig = self.analyzer.plot_results(analysis)
        
        # Embed in window
        canvas = FigureCanvasTkAgg(fig, master=result_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Close button
        btn_frame = ttk.Frame(result_window, padding="10")
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X)
        
        ttk.Button(btn_frame, text="Close", command=result_window.destroy,
                  style='Large.TButton').pack(side=tk.RIGHT, padx=10)
        
        result_window.lift()
        result_window.focus_force()


# ============================================================
# MAIN FUNCTION
# ============================================================

def main():
    """Main function"""
    print("=" * 80)
    print("Enhanced Gel Image Analysis - Tkinter GUI v2.0")
    print("=" * 80)
    
    # Find SCN files
    scn_files = sorted(glob.glob("*.scn"))
    
    if not scn_files:
        print("\nERROR: No SCN files found!")
        print("Please run in directory containing .scn files")
        sys.exit(1)
    
    print(f"\nFound {len(scn_files)} SCN files:")
    for i, f in enumerate(scn_files, 1):
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
                choice = input(f"\nSelect file (1-{len(scn_files)}): ")
                idx = int(choice) - 1
                if 0 <= idx < len(scn_files):
                    selected_file = scn_files[idx]
                    break
                else:
                    print("Invalid number")
            except ValueError:
                print("Please enter a number")
    
    # Read image
    safe_filename = sanitize_for_gui(selected_file)
    print(f"\nReading: {safe_filename}")
    print(f"File path: {selected_file}")
    print(f"File exists: {os.path.exists(selected_file)}")
    
    try:
        img_data, dimensions, error = read_scn_file(selected_file, priority='crop')
        
        if error:
            print(f"ERROR: {error}")  # 不要sanitize错误信息
            print(f"Dimensions: {dimensions}")
            sys.exit(1)
        
        if img_data is None:
            print("ERROR: Image data is None but no error message")
            sys.exit(1)
        
        print(f"Success! Image shape: {img_data.shape}")
        print(f"Image dtype: {img_data.dtype}")
        print(f"Image range: [{img_data.min()}, {img_data.max()}]")
        
    except Exception as e:
        print(f"EXCEPTION during file reading: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    print("\nLaunching Enhanced GUI...")
    
    # Create GUI with full file path
    root = tk.Tk()
    app = GelAnalysisGUI(root, img_data, selected_file, filepath=os.path.abspath(selected_file))
    
    print("\n" + "=" * 80)
    print("GUI Features:")
    print("  - Image rotation: ±90° and ±1° buttons (cumulative)")
    print("  - Three selection modes: Panel (red), Plasmid (green), Fragment (blue)")
    print("  - Linear background subtraction with robust fitting")
    print("  - Panel-relative intensity calculations")
    print("  - 5-subplot analysis results (A-E labels)")
    print("  - Fragment ratio computation")
    print("  - Large fonts (14-24pt) and windows (1800x1400, 2000x1600)")
    print("=" * 80)
    
    root.mainloop()


if __name__ == "__main__":
    main()
