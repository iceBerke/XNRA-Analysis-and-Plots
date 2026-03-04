### XNRA Spectrum Viewer GUI - Enhanced Version ###
# Tkinter-based graphical interface for XNRA spectrum visualization
# VERSION 5: Improved Uncertainty Calculation tab with editable systematic σ and depth profiles

# Compatible with main v6 and xnra_uncertainty_v5

# Script developed by Berke Santos & Giuseppe Legrottaglie
# Last updated: 2025

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import json
from pathlib import Path

# Import plotting function
try:
    from xnra_main_v6 import plot_xnra_spectrum
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("Warning: xnra_main_v6 not found. Plotting disabled.")

# Import uncertainty calculation (v5 with full documentation)
try:
    from xnra_uncertainty_v5 import (
        parse_xnra_full,
        calculate_chi_squared_for_roi,
        calculate_uncertainty_for_element,
        create_depth_profile,
        DEFAULT_SYSTEMATIC_UNCERTAINTIES
    )
    UNCERTAINTY_AVAILABLE = True
except ImportError:
    # Fallback defaults if module not available
    DEFAULT_SYSTEMATIC_UNCERTAINTIES = {'RBS': 0.03, 'ERDA': 0.05, 'NRA': 0.05}
    UNCERTAINTY_AVAILABLE = False
    print("Warning: xnra_uncertainty_v5 not found. Uncertainty calculations disabled.")


class XNRAViewerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("XNRA Spectrum Viewer - v5")
        self.root.geometry("1200x750")
        
        # Store current filepath
        self.current_file = tk.StringVar()
        self.last_fig = None
        
        # Store ROIs for uncertainty calculation
        self.roi_list = []
        self.uncertainty_results = None
        self.uncertainty_results_df = None  # For depth profile export
        
        # Create main container
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
        # Create widgets
        self.create_file_section(main_frame)
        self.create_tabs(main_frame)
        self.create_action_buttons(main_frame)
        
        # Load default values
        self.load_defaults()
    
    def validate_number_input(self, var, is_int=False, allow_empty=False):
        """Validate number input when Enter is pressed"""
        value = str(var.get()).strip()
        
        if allow_empty and value == "":
            return True
        
        value = value.replace(',', '.')
        
        try:
            if is_int:
                num_value = int(float(value))
                var.set(num_value)
            else:
                num_value = float(value)
                var.set(num_value)
            return True
        except ValueError:
            messagebox.showerror(
                "Invalid Input",
                f"Invalid number format: '{value}'\n\n"
                "Please enter a valid number.\n"
                "Note: Use '.' as decimal separator (not ',')"
            )
            return False
    
    def validate_number_input_from_entry(self, entry, var, is_int=False, allow_empty=False):
        """Validate number input directly from Entry widget text"""
        value = entry.get().strip()
        
        if allow_empty and value == "":
            return True
        
        value = value.replace(',', '.')
        
        try:
            if is_int:
                num_value = int(float(value))
                var.set(num_value)
            else:
                num_value = float(value)
                var.set(num_value)
            return True
        except ValueError:
            messagebox.showerror(
                "Invalid Input",
                f"Invalid number format: '{value}'\n\n"
                "Please enter a valid number.\n"
                "Note: Use '.' as decimal separator (not ',')"
            )
            entry.delete(0, tk.END)
            entry.insert(0, str(var.get()))
            return False
    
    def bind_mousewheel(self, canvas):
        """Bind mouse wheel scrolling to canvas"""
        def on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        def on_enter(event):
            canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        def on_leave(event):
            canvas.unbind_all("<MouseWheel>")
        
        canvas.bind("<Enter>", on_enter)
        canvas.bind("<Leave>", on_leave)
    
    def create_file_section(self, parent):
        """File selection section"""
        file_frame = ttk.LabelFrame(parent, text="File Selection", padding="10")
        file_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        file_frame.columnconfigure(1, weight=1)
        
        ttk.Label(file_frame, text="XNRA File:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        
        file_entry = ttk.Entry(file_frame, textvariable=self.current_file, state='readonly')
        file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        
        ttk.Button(file_frame, text="Browse...", command=self.browse_file).grid(row=0, column=2, padx=(5, 0))
    
    def create_tabs(self, parent):
        """Create tabbed interface for settings"""
        self.notebook = ttk.Notebook(parent)
        self.notebook.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))
        
        # Create tabs
        self.create_general_tab()
        self.create_raw_data_tab()
        self.create_smoothed_data_tab()
        self.create_simulated_data_tab()
        self.create_appearance_tab()
        self.create_axes_tab()
        self.create_uncertainty_tab()  # Improved version!
    
    def create_general_tab(self):
        """General settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="General")
        
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.bind_mousewheel(canvas)
        
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=2)
        
        # STACK container
        stack = ttk.Frame(scrollable_frame)
        stack.grid(row=0, column=0, sticky="nsew", padx=(8, 0), pady=(0, 10))
        stack.columnconfigure(0, weight=1)

        # Title section
        title_frame = ttk.LabelFrame(stack, text="Plot Title", padding="10")
        title_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        
        ttk.Label(title_frame, text="Custom Title:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        self.custom_title = tk.StringVar(value="")
        ttk.Entry(title_frame, textvariable=self.custom_title, width=50).grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Label(title_frame, text="(empty = filename without extension)", font=('TkDefaultFont', 8, 'italic')).grid(row=1, column=0, sticky=tk.W, padx=(2, 0))

        # Figure settings
        fig_frame = ttk.LabelFrame(stack, text="Figure Settings", padding="10")
        fig_frame.grid(row=1, column=0, sticky="ew")
        
        self.fig_width = self.create_number_field(fig_frame, "Width (inches):", 0, 12.0)
        self.fig_height = self.create_number_field(fig_frame, "Height (inches):", 1, 6.0)
        self.fig_dpi = self.create_number_field(fig_frame, "DPI:", 2, 100, is_int=True)
        
        # Display options
        display_frame = ttk.LabelFrame(scrollable_frame, text="Display Options", padding="10")
        display_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
        # Legend section
        self.show_legend = tk.BooleanVar(value=True)
        ttk.Checkbutton(display_frame, text="Show Legend", variable=self.show_legend).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(display_frame, text="Legend Position:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.legend_position = tk.StringVar(value='right')
        ttk.Radiobutton(display_frame, text="Left", variable=self.legend_position, value='left',
                        command=self.update_legend_defaults).grid(row=2, column=0, sticky=tk.W, padx=(20, 0))
        ttk.Radiobutton(display_frame, text="Right", variable=self.legend_position, value='right',
                        command=self.update_legend_defaults).grid(row=3, column=0, sticky=tk.W, padx=(20, 0))
        
        ttk.Label(display_frame, text="Fine-tune Position:", font=('TkDefaultFont', 8, 'italic')).grid(row=4, column=0, sticky=tk.W, pady=(10, 2))
        self.legend_x_offset = self.create_number_field(display_frame, "X-offset:", 5, 0.98)
        self.legend_y_offset = self.create_number_field(display_frame, "Y-offset:", 6, 0.88)
        ttk.Label(display_frame, text="(0.0 = left/bottom, 1.0 = right/top)", font=('TkDefaultFont', 7, 'italic')).grid(row=7, column=0, sticky=tk.W)
        
        ttk.Separator(display_frame, orient='horizontal').grid(row=8, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        
        # Info box section
        self.show_info_box = tk.BooleanVar(value=True)
        ttk.Checkbutton(display_frame, text="Show Info Box", variable=self.show_info_box).grid(row=9, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(display_frame, text="Info Box Position:").grid(row=10, column=0, sticky=tk.W, pady=(5, 0))
        self.info_box_position = tk.StringVar(value='left')
        ttk.Radiobutton(display_frame, text="Left", variable=self.info_box_position, value='left',
                        command=self.update_info_box_defaults).grid(row=11, column=0, sticky=tk.W, padx=(20, 0))
        ttk.Radiobutton(display_frame, text="Right", variable=self.info_box_position, value='right',
                        command=self.update_info_box_defaults).grid(row=12, column=0, sticky=tk.W, padx=(20, 0))
        
        ttk.Label(display_frame, text="Fine-tune Position:", font=('TkDefaultFont', 8, 'italic')).grid(row=13, column=0, sticky=tk.W, pady=(10, 2))
        self.info_box_x_offset = self.create_number_field(display_frame, "X-offset:", 14, 0.05)
        self.info_box_y_offset = self.create_number_field(display_frame, "Y-offset:", 15, 0.85)
        ttk.Label(display_frame, text="(0.0 = left/bottom, 1.0 = right/top)", font=('TkDefaultFont', 7, 'italic')).grid(row=16, column=0, sticky=tk.W)
            
    def update_info_box_defaults(self):
        if self.info_box_position.get() == 'right':
            self.info_box_x_offset.set(0.95)
            self.info_box_y_offset.set(0.65)
        else:
            self.info_box_x_offset.set(0.05)
            self.info_box_y_offset.set(0.85)

    def update_legend_defaults(self):
        if self.legend_position.get() == 'right':
            self.legend_x_offset.set(0.98)
            self.legend_y_offset.set(0.88)
        else:
            self.legend_x_offset.set(0.02)
            self.legend_y_offset.set(0.88)

    def create_raw_data_tab(self):
        """Raw data settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Raw Data")
        
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.bind_mousewheel(canvas)
        
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        scrollable_frame.columnconfigure(2, weight=1)
        
        self.show_raw_data = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="Show Raw Data", variable=self.show_raw_data).grid(row=0, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))
        
        # Markers
        marker_frame = ttk.LabelFrame(scrollable_frame, text="Markers", padding="10")
        marker_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.show_raw_markers = tk.BooleanVar(value=True)
        ttk.Checkbutton(marker_frame, text="Show Markers", variable=self.show_raw_markers).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(marker_frame, text="Shape:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.raw_marker_shape = tk.StringVar(value='square')
        ttk.Combobox(marker_frame, textvariable=self.raw_marker_shape, values=['circle', 'square', 'triangle'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        ttk.Label(marker_frame, text="Fill:").grid(row=2, column=0, sticky=tk.W)
        self.raw_marker_fill = tk.StringVar(value='hollow')
        ttk.Combobox(marker_frame, textvariable=self.raw_marker_fill, values=['filled', 'hollow'], state='readonly', width=12).grid(row=2, column=1, sticky=tk.W, padx=(5, 0))
        
        self.raw_marker_color = self.create_color_field(marker_frame, "Color:", 3, 'black')
        self.raw_marker_size = self.create_number_field(marker_frame, "Size (points):", 4, 2.0)
        self.raw_marker_edge_width = self.create_number_field(marker_frame, "Edge Width:", 5, 1.0)
        
        # Marker Range
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.raw_marker_range_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.raw_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.raw_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.raw_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.raw_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 0.0)
        self.raw_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 200.0)
        
        # Line
        line_frame = ttk.LabelFrame(scrollable_frame, text="Line", padding="10")
        line_frame.grid(row=1, column=2, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.show_raw_line = tk.BooleanVar(value=False)
        ttk.Checkbutton(line_frame, text="Show Line", variable=self.show_raw_line).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        self.raw_line_color = self.create_color_field(line_frame, "Color:", 1, 'blue')
        self.raw_line_width = self.create_number_field(line_frame, "Width:", 2, 1.0)
        self.raw_line_alpha = self.create_number_field(line_frame, "Transparency:", 3, 1.0)
    
    def create_smoothed_data_tab(self):
        """Smoothed data settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Smoothed Data")
        
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.bind_mousewheel(canvas)
        
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        scrollable_frame.columnconfigure(2, weight=1)
        
        self.show_smoothed_data = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="Show Smoothed Data", variable=self.show_smoothed_data).grid(row=0, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))
        
        # Markers
        marker_frame = ttk.LabelFrame(scrollable_frame, text="Markers", padding="10")
        marker_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.show_smoothed_markers = tk.BooleanVar(value=True)
        ttk.Checkbutton(marker_frame, text="Show Markers", variable=self.show_smoothed_markers).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(marker_frame, text="Shape:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.smoothed_marker_shape = tk.StringVar(value='triangle')
        ttk.Combobox(marker_frame, textvariable=self.smoothed_marker_shape, values=['circle', 'square', 'triangle'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        ttk.Label(marker_frame, text="Fill:").grid(row=2, column=0, sticky=tk.W)
        self.smoothed_marker_fill = tk.StringVar(value='hollow')
        ttk.Combobox(marker_frame, textvariable=self.smoothed_marker_fill, values=['filled', 'hollow'], state='readonly', width=12).grid(row=2, column=1, sticky=tk.W, padx=(5, 0))
        
        self.smoothed_marker_color = self.create_color_field(marker_frame, "Color:", 3, 'orange')
        self.smoothed_marker_size = self.create_number_field(marker_frame, "Size (points):", 4, 3.0)
        self.smoothed_marker_edge_width = self.create_number_field(marker_frame, "Edge Width:", 5, 1.0)
        
        # Range
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.smoothed_marker_range_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.smoothed_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.smoothed_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.smoothed_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.smoothed_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 200.0)
        self.smoothed_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 400.0)
        
        # Line
        line_frame = ttk.LabelFrame(scrollable_frame, text="Line", padding="10")
        line_frame.grid(row=1, column=2, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.show_smoothed_line = tk.BooleanVar(value=True)
        ttk.Checkbutton(line_frame, text="Show Line", variable=self.show_smoothed_line).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        self.smoothed_line_color = self.create_color_field(line_frame, "Color:", 1, 'red')
        self.smoothed_line_width = self.create_number_field(line_frame, "Width:", 2, 1.0)
        self.smoothed_line_alpha = self.create_number_field(line_frame, "Transparency:", 3, 1.0)
    
    def create_simulated_data_tab(self):
        """Simulated data settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Simulated Data")
        
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.bind_mousewheel(canvas)
        
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        scrollable_frame.columnconfigure(2, weight=1)
        
        self.show_simulated_data = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="Show Simulated Data", variable=self.show_simulated_data).grid(row=0, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))
        
        # Markers
        marker_frame = ttk.LabelFrame(scrollable_frame, text="Markers", padding="10")
        marker_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.show_simulated_markers = tk.BooleanVar(value=False)
        ttk.Checkbutton(marker_frame, text="Show Markers", variable=self.show_simulated_markers).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(marker_frame, text="Shape:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.simulated_marker_shape = tk.StringVar(value='circle')
        ttk.Combobox(marker_frame, textvariable=self.simulated_marker_shape, values=['circle', 'square', 'triangle'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        ttk.Label(marker_frame, text="Fill:").grid(row=2, column=0, sticky=tk.W)
        self.simulated_marker_fill = tk.StringVar(value='filled')
        ttk.Combobox(marker_frame, textvariable=self.simulated_marker_fill, values=['filled', 'hollow'], state='readonly', width=12).grid(row=2, column=1, sticky=tk.W, padx=(5, 0))
        
        self.simulated_marker_color = self.create_color_field(marker_frame, "Color:", 3, 'green')
        self.simulated_marker_size = self.create_number_field(marker_frame, "Size (points):", 4, 2.0)
        self.simulated_marker_edge_width = self.create_number_field(marker_frame, "Edge Width:", 5, 1.0)
        
        # Range
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.simulated_marker_range_enabled = tk.BooleanVar(value=False)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.simulated_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.simulated_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.simulated_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.simulated_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 0.0)
        self.simulated_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 1000.0)
        
        # Line
        line_frame = ttk.LabelFrame(scrollable_frame, text="Line", padding="10")
        line_frame.grid(row=1, column=2, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.show_simulated_line = tk.BooleanVar(value=True)
        ttk.Checkbutton(line_frame, text="Show Line", variable=self.show_simulated_line).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        self.simulated_line_color = self.create_color_field(line_frame, "Color:", 1, 'green')
        self.simulated_line_width = self.create_number_field(line_frame, "Width:", 2, 2.0)
        self.simulated_line_alpha = self.create_number_field(line_frame, "Transparency:", 3, 0.8)
    
    def create_appearance_tab(self):
        """Plot appearance settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Appearance")
        
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.bind_mousewheel(canvas)
        
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        
        # Colors
        color_frame = ttk.LabelFrame(scrollable_frame, text="Colors", padding="10")
        color_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5), pady=(0, 10))
        
        self.background_color = self.create_color_field(color_frame, "Background:", 0, 'white')
        self.grid_color = self.create_color_field(color_frame, "Grid:", 1, 'white')
        
        self.show_grid = tk.BooleanVar(value=True)
        ttk.Checkbutton(color_frame, text="Show Grid", variable=self.show_grid).grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(5, 0))
        
        # Text sizes
        size_frame = ttk.LabelFrame(scrollable_frame, text="Text Sizes", padding="10")
        size_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
        self.title_size = self.create_number_field(size_frame, "Title:", 0, 16.0)
        self.axis_label_size = self.create_number_field(size_frame, "Axis Labels:", 1, 12.0)
        self.tick_label_size = self.create_number_field(size_frame, "Tick Labels:", 2, 10.0)
        self.legend_size = self.create_number_field(size_frame, "Legend:", 3, 10.0)
        self.info_box_size = self.create_number_field(size_frame, "Info Box:", 4, 10.0)
        
        # Text colors
        text_color_frame = ttk.LabelFrame(scrollable_frame, text="Text Colors", padding="10")
        text_color_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.title_color = self.create_color_field(text_color_frame, "Title:", 0, 'black')
        self.axis_label_color = self.create_color_field(text_color_frame, "Axis Labels:", 1, 'black')
        self.tick_label_color = self.create_color_field(text_color_frame, "Tick Labels:", 2, 'black')
        self.legend_text_color = self.create_color_field(text_color_frame, "Legend:", 3, 'black')
        self.info_box_text_color = self.create_color_field(text_color_frame, "Info Box:", 4, 'black')
        
        # Tick marks
        tick_frame = ttk.LabelFrame(scrollable_frame, text="Tick Marks", padding="10")
        tick_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.tick_length = self.create_number_field(tick_frame, "Length (points):", 0, 6.0)
        self.tick_width = self.create_number_field(tick_frame, "Width (points):", 1, 1.0)
        self.tick_color = self.create_color_field(tick_frame, "Color:", 2, 'black')
        
        ttk.Label(tick_frame, text="Direction:").grid(row=3, column=0, sticky=tk.W, pady=(5, 0))
        self.tick_direction = tk.StringVar(value='in')
        ttk.Combobox(tick_frame, textvariable=self.tick_direction, values=['in', 'out', 'inout'], state='readonly', width=12).grid(row=3, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
    
    def create_axes_tab(self):
        """Axes and limits settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Axes & Limits")
        
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.bind_mousewheel(canvas)
        
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        
        # X-axis mode
        mode_frame = ttk.LabelFrame(scrollable_frame, text="X-Axis Mode", padding="10")
        mode_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5), pady=(0, 10))
        
        ttk.Label(mode_frame, text="Display Mode:").grid(row=0, column=0, sticky=tk.W)
        self.x_axis_mode = tk.StringVar(value='energy')
        ttk.Combobox(mode_frame, textvariable=self.x_axis_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=0, column=1, sticky=tk.W, padx=(5, 0))
        
        # Channel limits
        channel_frame = ttk.LabelFrame(scrollable_frame, text="Channel Limits", padding="10")
        channel_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
        self.channel_min = self.create_optional_number_field(channel_frame, "Min:", 0)
        self.channel_max = self.create_optional_number_field(channel_frame, "Max:", 1)
        
        # Energy limits
        energy_frame = ttk.LabelFrame(scrollable_frame, text="Energy Limits (keV)", padding="10")
        energy_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.energy_min = self.create_optional_number_field(energy_frame, "Min:", 0)
        self.energy_max = self.create_optional_number_field(energy_frame, "Max:", 1)
        
        # Counts limits
        counts_frame = ttk.LabelFrame(scrollable_frame, text="Counts Limits (Y-axis)", padding="10")
        counts_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.counts_min = self.create_optional_number_field(counts_frame, "Min:", 0)
        self.counts_max = self.create_optional_number_field(counts_frame, "Max:", 1)
    
    # =========================================================================
    # IMPROVED UNCERTAINTY CALCULATION TAB
    # =========================================================================
    
    def create_uncertainty_tab(self):
        """
        Uncertainty calculation tab with:
        - Editable systematic uncertainty
        - ROI definition
        - Results display
        - Depth profile export
        """
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Uncertainty")
        
        # Main layout
        tab.columnconfigure(0, weight=1)
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(0, weight=1)
        
        # ===== LEFT PANEL: Settings =====
        left_panel = ttk.Frame(tab)
        left_panel.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        left_panel.rowconfigure(2, weight=1)
        left_panel.columnconfigure(0, weight=1)
        
        # Technique and Systematic Uncertainty frame
        tech_frame = ttk.LabelFrame(left_panel, text="Technique & Systematic Uncertainty", padding="10")
        tech_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        tech_frame.columnconfigure(1, weight=1)
        
        # Technique selection
        ttk.Label(tech_frame, text="Data Type:").grid(row=0, column=0, sticky=tk.W, padx=(0, 10))
        self.uncertainty_technique = tk.StringVar(value='RBS')
        tech_combo = ttk.Combobox(tech_frame, textvariable=self.uncertainty_technique, 
                                   values=['RBS', 'ERDA', 'NRA'], state='readonly', width=10)
        tech_combo.grid(row=0, column=1, sticky=tk.W)
        
        # Editable systematic uncertainty
        ttk.Label(tech_frame, text="Systematic σ (%):").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.systematic_uncertainty = tk.DoubleVar(value=3.0)
        sys_entry = ttk.Entry(tech_frame, textvariable=self.systematic_uncertainty, width=10)
        sys_entry.grid(row=1, column=1, sticky=tk.W, pady=(5, 0))
        
        # Info label explaining the systematic uncertainty
        sys_info = ttk.Label(tech_frame, 
                             text="Typical values: RBS ≈ 3% (stopping powers),\n"
                                  "ERDA ≈ 5%, NRA ≈ 5-10% (cross-sections)",
                             font=('TkDefaultFont', 8, 'italic'), foreground='gray')
        sys_info.grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(5, 0))
        
        # Update systematic uncertainty when technique changes
        def update_sys_unc_default(*args):
            tech = self.uncertainty_technique.get()
            default_sys = {'RBS': 3.0, 'ERDA': 5.0, 'NRA': 5.0}.get(tech, 5.0)
            self.systematic_uncertainty.set(default_sys)
        self.uncertainty_technique.trace('w', update_sys_unc_default)
        
        # Formula info frame
        formula_frame = ttk.LabelFrame(left_panel, text="Formulas Used", padding="10")
        formula_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        
        formulas_text = (
            "Statistical:  σ_stat = C / √(C × N_roi)\n"
            "    where C = concentration, N_roi = counts in ROI\n\n"
            "Systematic:  σ_sys = C × (σ_sys%/100)\n\n"
            "Combined:    σ_total = √(σ_stat² + σ_sys²)\n\n"
            "Chi-squared: χ²/ν = Σ[(raw-sim)²/raw] / (n-1)\n"
            "    Good fit: χ²/ν ≈ 1"
        )
        ttk.Label(formula_frame, text=formulas_text, font=('Consolas', 8), justify=tk.LEFT).pack(anchor=tk.W)
        
        # ROI definition section
        roi_frame = ttk.LabelFrame(left_panel, text="Regions of Interest (ROIs)", padding="10")
        roi_frame.grid(row=2, column=0, sticky="nsew", pady=(0, 10))
        roi_frame.columnconfigure(0, weight=1)
        roi_frame.rowconfigure(1, weight=1)
        
        # ROI list header
        header_frame = ttk.Frame(roi_frame)
        header_frame.grid(row=0, column=0, sticky="ew", pady=(0, 5))
        
        ttk.Label(header_frame, text="Name", width=12).pack(side=tk.LEFT, padx=2)
        ttk.Label(header_frame, text="Start Ch", width=8).pack(side=tk.LEFT, padx=2)
        ttk.Label(header_frame, text="End Ch", width=8).pack(side=tk.LEFT, padx=2)
        ttk.Label(header_frame, text="Elements (comma-sep)", width=20).pack(side=tk.LEFT, padx=2)
        
        # Scrollable ROI list
        roi_canvas = tk.Canvas(roi_frame, height=150)
        roi_scrollbar = ttk.Scrollbar(roi_frame, orient="vertical", command=roi_canvas.yview)
        self.roi_container = ttk.Frame(roi_canvas)
        
        self.roi_container.bind("<Configure>", lambda e: roi_canvas.configure(scrollregion=roi_canvas.bbox("all")))
        roi_canvas.create_window((0, 0), window=self.roi_container, anchor="nw")
        roi_canvas.configure(yscrollcommand=roi_scrollbar.set)
        
        roi_canvas.grid(row=1, column=0, sticky="nsew")
        roi_scrollbar.grid(row=1, column=1, sticky="ns")
        
        # ROI buttons
        roi_btn_frame = ttk.Frame(roi_frame)
        roi_btn_frame.grid(row=2, column=0, sticky="ew", pady=(10, 0))
        
        ttk.Button(roi_btn_frame, text="+ Add ROI", command=self.add_roi_row).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(roi_btn_frame, text="Load Defaults", command=self.load_default_rois).pack(side=tk.LEFT, padx=5)
        ttk.Button(roi_btn_frame, text="Clear All", command=self.clear_all_rois).pack(side=tk.LEFT, padx=5)
        
        # Action buttons
        calc_frame = ttk.Frame(left_panel)
        calc_frame.grid(row=3, column=0, sticky="ew")
        
        ttk.Button(calc_frame, text="Calculate Uncertainties", 
                   command=self.calculate_uncertainties).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(calc_frame, text="Export Raw Results", 
                   command=self.export_uncertainty_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(calc_frame, text="Export Depth Profile", 
                   command=self.export_depth_profile).pack(side=tk.LEFT, padx=5)
        
        # ===== RIGHT PANEL: Results =====
        right_panel = ttk.LabelFrame(tab, text="Results", padding="10")
        right_panel.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
        right_panel.rowconfigure(0, weight=1)
        right_panel.columnconfigure(0, weight=1)
        
        # Results treeview
        columns = ('Layer', 'Element', 'Conc (%)', '± Stat', '± Sys', '± Total', 'ROI', 'χ²/ν')
        self.results_tree = ttk.Treeview(right_panel, columns=columns, show='headings', height=15)
        
        for col in columns:
            self.results_tree.heading(col, text=col)
            width = 75 if col not in ('Element', 'ROI') else 55
            self.results_tree.column(col, width=width, anchor='center')
        
        # Scrollbars
        results_yscroll = ttk.Scrollbar(right_panel, orient="vertical", command=self.results_tree.yview)
        results_xscroll = ttk.Scrollbar(right_panel, orient="horizontal", command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=results_yscroll.set, xscrollcommand=results_xscroll.set)
        
        self.results_tree.grid(row=0, column=0, sticky="nsew")
        results_yscroll.grid(row=0, column=1, sticky="ns")
        results_xscroll.grid(row=1, column=0, sticky="ew")
        
        # Status label
        self.uncertainty_status = ttk.Label(right_panel, 
                                             text="Load a file and define ROIs to calculate uncertainties", 
                                             font=('TkDefaultFont', 9, 'italic'))
        self.uncertainty_status.grid(row=2, column=0, sticky="w", pady=(10, 0))
        
        # Load default ROIs
        self.load_default_rois()
    
    def add_roi_row(self, name="", start="", end="", elements=""):
        """Add a new ROI entry row"""
        row_frame = ttk.Frame(self.roi_container)
        row_frame.pack(fill=tk.X, pady=2)
        
        name_var = tk.StringVar(value=name)
        start_var = tk.StringVar(value=start)
        end_var = tk.StringVar(value=end)
        elements_var = tk.StringVar(value=elements)
        
        ttk.Entry(row_frame, textvariable=name_var, width=12).pack(side=tk.LEFT, padx=2)
        ttk.Entry(row_frame, textvariable=start_var, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Entry(row_frame, textvariable=end_var, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Entry(row_frame, textvariable=elements_var, width=20).pack(side=tk.LEFT, padx=2)
        
        def remove_row():
            row_frame.destroy()
            self.roi_list = [r for r in self.roi_list if r['frame'] != row_frame]
        
        ttk.Button(row_frame, text="×", width=3, command=remove_row).pack(side=tk.LEFT, padx=2)
        
        self.roi_list.append({
            'frame': row_frame,
            'name': name_var,
            'start': start_var,
            'end': end_var,
            'elements': elements_var
        })
    
    def load_default_rois(self):
        """Load default ROI definitions based on technique"""
        self.clear_all_rois()
        
        tech = self.uncertainty_technique.get()
        
        if tech == 'RBS':
            self.add_roi_row("Heavy", "700", "900", "Ba,Ca,Fe")
            self.add_roi_row("Medium", "380", "550", "Si,Al,Na,Mg")
            self.add_roi_row("Light", "100", "250", "O,C")
        elif tech == 'ERDA':
            self.add_roi_row("Hydrogen", "200", "400", "H")
        elif tech == 'NRA':
            self.add_roi_row("Carbon", "100", "250", "C")
            self.add_roi_row("Oxygen", "400", "600", "O")
            self.add_roi_row("Nitrogen", "250", "400", "N")
    
    def clear_all_rois(self):
        """Remove all ROI entries"""
        for roi in self.roi_list:
            roi['frame'].destroy()
        self.roi_list = []
    
    def get_rois_from_gui(self):
        """Extract ROI definitions from GUI"""
        rois = []
        for roi in self.roi_list:
            try:
                name = roi['name'].get().strip() or f"ROI_{len(rois)+1}"
                start = int(roi['start'].get())
                end = int(roi['end'].get())
                elements = [e.strip() for e in roi['elements'].get().split(',') if e.strip()]
                
                if elements and start < end:
                    rois.append({
                        'name': name,
                        'channels': (start, end),
                        'elements': elements
                    })
            except ValueError:
                continue
        return rois
    
    def calculate_uncertainties(self):
        """Calculate uncertainties for current file"""
        if not UNCERTAINTY_AVAILABLE:
            messagebox.showerror("Error", 
                "Uncertainty calculation module not available.\n\n"
                "Make sure xnra_uncertainty_v5.py is in the same directory.")
            return
        
        if not self.current_file.get():
            messagebox.showerror("Error", "Please select an XNRA file first!")
            return
        
        if not Path(self.current_file.get()).exists():
            messagebox.showerror("Error", "Selected file does not exist!")
            return
        
        rois = self.get_rois_from_gui()
        if not rois:
            messagebox.showerror("Error", "Please define at least one valid ROI!")
            return
        
        try:
            data = parse_xnra_full(self.current_file.get())
            technique = self.uncertainty_technique.get()
            
            # Get user-specified systematic uncertainty (as fraction)
            sys_unc = self.systematic_uncertainty.get() / 100.0
            
            # Clear previous results
            for item in self.results_tree.get_children():
                self.results_tree.delete(item)
            
            results = []
            
            for roi in rois:
                roi_name = roi['name']
                roi_start, roi_end = roi['channels']
                roi_elements = roi['elements']
                
                # Calculate chi-squared for this ROI
                if data['raw_counts'] is not None and data['simulated_counts'] is not None:
                    chi2_r, dof, roi_counts = calculate_chi_squared_for_roi(
                        data['raw_counts'],
                        data['simulated_counts'],
                        roi_start,
                        roi_end
                    )
                else:
                    chi2_r, dof, roi_counts = float('nan'), 0, 0
                
                for layer in data['layers']:
                    for element in roi_elements:
                        concentration = layer['concentrations'].get(element, 0)
                        
                        if concentration <= 0:
                            continue
                        
                        # Use user-specified systematic uncertainty
                        unc = calculate_uncertainty_for_element(
                            concentration, roi_counts, technique, sys_unc
                        )
                        
                        # Format chi² status
                        if chi2_r < 2:
                            chi2_status = f"{chi2_r:.2f} ✓"
                        elif chi2_r < 10:
                            chi2_status = f"{chi2_r:.2f} ⚠"
                        else:
                            chi2_status = f"{chi2_r:.1f} ✗"
                        
                        values = (
                            layer['layer_num'],
                            element,
                            f"{concentration*100:.2f}",
                            f"{unc['sigma_stat']*100:.3f}" if unc['sigma_stat'] == unc['sigma_stat'] else "N/A",
                            f"{unc['sigma_sys']*100:.3f}",
                            f"{unc['sigma_total']*100:.3f}",
                            roi_name,
                            chi2_status
                        )
                        self.results_tree.insert('', 'end', values=values)
                        
                        results.append({
                            'filename': data['filename'],
                            'layer': layer['layer_num'],
                            'thickness': layer['thickness'],
                            'element': element,
                            'concentration_percent': concentration * 100,
                            'sigma_stat_percent': unc['sigma_stat'] * 100 if unc['sigma_stat'] == unc['sigma_stat'] else None,
                            'sigma_sys_percent': unc['sigma_sys'] * 100,
                            'sigma_total_percent': unc['sigma_total'] * 100,
                            'roi_name': roi_name,
                            'roi_start': roi_start,
                            'roi_end': roi_end,
                            'roi_counts': roi_counts,
                            'chi2_reduced': chi2_r,
                            'technique': technique,
                            'systematic_used_percent': sys_unc * 100
                        })
            
            self.uncertainty_results = results
            
            # Create DataFrame for depth profile export
            import pandas as pd
            self.uncertainty_results_df = pd.DataFrame(results)
            
            self.uncertainty_status.config(
                text=f"✓ Calculated {len(results)} results from {len(rois)} ROI(s) | "
                     f"Systematic σ: {sys_unc*100:.1f}%"
            )
            
        except Exception as e:
            import traceback
            messagebox.showerror("Error", f"Calculation failed:\n{str(e)}\n\n{traceback.format_exc()}")
            self.uncertainty_status.config(text=f"✗ Error: {str(e)}")
    
    def export_uncertainty_results(self):
        """Export raw uncertainty results to Excel or CSV"""
        if not self.uncertainty_results:
            messagebox.showerror("Error", "No results to export. Calculate uncertainties first!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export Raw Uncertainty Results",
            defaultextension=".xlsx",
            filetypes=[
                ("Excel files", "*.xlsx"),
                ("CSV files", "*.csv"),
                ("All files", "*.*")
            ]
        )
        
        if not filename:
            return
        
        try:
            import pandas as pd
            df = pd.DataFrame(self.uncertainty_results)
            
            if filename.endswith('.xlsx'):
                df.to_excel(filename, index=False, sheet_name='Uncertainties')
            else:
                df.to_csv(filename, index=False)
            
            messagebox.showinfo("Success", f"Results exported to:\n{filename}")
        except ImportError:
            messagebox.showerror("Error", "pandas is required for export.\nInstall with: pip install pandas openpyxl")
        except Exception as e:
            messagebox.showerror("Error", f"Export failed:\n{str(e)}")
    
    def export_depth_profile(self):
        """Export depth profile table (concentrations with uncertainties)"""
        if self.uncertainty_results_df is None or len(self.uncertainty_results_df) == 0:
            messagebox.showerror("Error", "No results to export. Calculate uncertainties first!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export Depth Profile",
            defaultextension=".xlsx",
            filetypes=[
                ("Excel files", "*.xlsx"),
                ("CSV files", "*.csv"),
                ("All files", "*.*")
            ]
        )
        
        if not filename:
            return
        
        try:
            import pandas as pd
            
            # Create depth profile using the imported function
            df = self.uncertainty_results_df
            profile = create_depth_profile(df)
            
            if filename.endswith('.xlsx'):
                with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                    profile.to_excel(writer, sheet_name='Depth Profile', index=False)
                    df.to_excel(writer, sheet_name='Raw Data', index=False)
            else:
                profile.to_csv(filename, index=False)
            
            messagebox.showinfo("Success", f"Depth profile exported to:\n{filename}")
        except ImportError:
            messagebox.showerror("Error", "pandas is required for export.\nInstall with: pip install pandas openpyxl")
        except Exception as e:
            import traceback
            messagebox.showerror("Error", f"Export failed:\n{str(e)}\n\n{traceback.format_exc()}")
    
    # =========================================================================
    # ACTION BUTTONS AND HELPERS
    # =========================================================================
    
    def create_action_buttons(self, parent):
        """Action buttons section"""
        button_frame = ttk.Frame(parent)
        button_frame.grid(row=2, column=0, sticky=(tk.W, tk.E))
        
        ttk.Button(button_frame, text="Load Preset", command=self.load_preset).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(button_frame, text="Save Preset", command=self.save_preset).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Export Plot", command=self.export_plot).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Generate Plot", command=self.generate_plot, style='Accent.TButton').pack(side=tk.RIGHT, padx=(5, 0))
    
    def create_number_field(self, parent, label, row, default, is_int=False):
        """Create a labeled number entry field with validation"""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky=tk.W)
        
        var = tk.DoubleVar(value=default) if not is_int else tk.IntVar(value=int(default))
        entry = ttk.Entry(parent, textvariable=var, width=12)
        entry.grid(row=row, column=1, sticky=tk.W, padx=(5, 0))
        
        def on_enter(event):
            self.validate_number_input_from_entry(entry, var, is_int=is_int, allow_empty=False)
        entry.bind('<Return>', on_enter)
        
        return var
    
    def create_optional_number_field(self, parent, label, row):
        """Create a labeled optional number entry field"""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky=tk.W)
        
        var = tk.StringVar(value="")
        entry = ttk.Entry(parent, textvariable=var, width=12)
        entry.grid(row=row, column=1, sticky=tk.W, padx=(5, 0))
        
        ttk.Label(parent, text="(empty = auto)", font=('TkDefaultFont', 8, 'italic')).grid(row=row, column=2, sticky=tk.W, padx=(5, 0))
        
        def on_enter(event):
            self.validate_number_input_from_entry(entry, var, is_int=False, allow_empty=True)
        entry.bind('<Return>', on_enter)
        
        return var
    
    def create_color_field(self, parent, label, row, default):
        """Create a labeled color picker field"""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky=tk.W)
        
        var = tk.StringVar(value=default)
        frame = ttk.Frame(parent)
        frame.grid(row=row, column=1, sticky=tk.W, padx=(5, 0))
        
        entry = ttk.Entry(frame, textvariable=var, width=8)
        entry.pack(side=tk.LEFT, padx=(0, 3))
        
        color_preview = tk.Canvas(frame, width=18, height=18, bg=default, relief=tk.SUNKEN, borderwidth=1)
        color_preview.pack(side=tk.LEFT, padx=(0, 3))
        
        def pick_color():
            color = colorchooser.askcolor(initialcolor=var.get(), title=f"Choose {label}")
            if color[1]:
                var.set(color[1])
                color_preview.configure(bg=color[1])
        
        ttk.Button(frame, text="Pick", command=pick_color, width=5).pack(side=tk.LEFT)
        
        def update_preview(*args):
            try:
                color_preview.configure(bg=var.get())
            except:
                pass
        var.trace('w', update_preview)
        
        return var
    
    def browse_file(self):
        """Open file browser dialog"""
        filename = filedialog.askopenfilename(
            title="Select XNRA File",
            filetypes=[("XNRA files", "*.xnra"), ("All files", "*.*")]
        )
        if filename:
            self.current_file.set(filename)
    
    def get_optional_float(self, var):
        """Convert StringVar to float or None"""
        val = var.get().strip()
        if val == "":
            return None
        val = val.replace(',', '.')
        try:
            return float(val)
        except ValueError:
            return None
    
    def generate_plot(self):
        """Generate plot with current settings"""
        if not PLOTTING_AVAILABLE:
            messagebox.showerror("Error", "Plotting module not available.\n\nMake sure xnra_main_v6.py is in the same directory.")
            return
        
        if not self.current_file.get():
            messagebox.showerror("Error", "Please select an XNRA file first!")
            return
        
        if not Path(self.current_file.get()).exists():
            messagebox.showerror("Error", "Selected file does not exist!")
            return
        
        try:
            if self.custom_title.get().strip():
                plot_title = self.custom_title.get().strip()
            else:
                plot_title = Path(self.current_file.get()).stem
            
            params = {
                'filepath': self.current_file.get(),
                'title': plot_title,
                'figure_size': (self.fig_width.get(), self.fig_height.get()),
                'dpi': self.fig_dpi.get(),
                'show_legend': self.show_legend.get(),
                'show_info_box': self.show_info_box.get(),
                'legend_position': self.legend_position.get(),
                'legend_x_offset': self.legend_x_offset.get(),
                'legend_y_offset': self.legend_y_offset.get(),
                'info_box_position': self.info_box_position.get(),
                'info_box_x_offset': self.info_box_x_offset.get(),
                'info_box_y_offset': self.info_box_y_offset.get(),
                'show_raw_data': self.show_raw_data.get(),
                'show_raw_markers': self.show_raw_markers.get(),
                'raw_marker_shape': self.raw_marker_shape.get(),
                'raw_marker_fill': self.raw_marker_fill.get(),
                'raw_marker_color': self.raw_marker_color.get(),
                'raw_marker_size': self.raw_marker_size.get(),
                'raw_marker_edge_width': self.raw_marker_edge_width.get(),
                'raw_marker_range_enabled': self.raw_marker_range_enabled.get(),
                'raw_marker_range_mode': self.raw_marker_range_mode.get(),
                'raw_marker_range_min': self.raw_marker_range_min.get(),
                'raw_marker_range_max': self.raw_marker_range_max.get(),
                'show_raw_line': self.show_raw_line.get(),
                'raw_line_color': self.raw_line_color.get(),
                'raw_line_width': self.raw_line_width.get(),
                'raw_line_alpha': self.raw_line_alpha.get(),
                'show_smoothed_data': self.show_smoothed_data.get(),
                'show_smoothed_markers': self.show_smoothed_markers.get(),
                'smoothed_marker_shape': self.smoothed_marker_shape.get(),
                'smoothed_marker_fill': self.smoothed_marker_fill.get(),
                'smoothed_marker_color': self.smoothed_marker_color.get(),
                'smoothed_marker_size': self.smoothed_marker_size.get(),
                'smoothed_marker_edge_width': self.smoothed_marker_edge_width.get(),
                'smoothed_marker_range_enabled': self.smoothed_marker_range_enabled.get(),
                'smoothed_marker_range_mode': self.smoothed_marker_range_mode.get(),
                'smoothed_marker_range_min': self.smoothed_marker_range_min.get(),
                'smoothed_marker_range_max': self.smoothed_marker_range_max.get(),
                'show_smoothed_line': self.show_smoothed_line.get(),
                'smoothed_line_color': self.smoothed_line_color.get(),
                'smoothed_line_width': self.smoothed_line_width.get(),
                'smoothed_line_alpha': self.smoothed_line_alpha.get(),
                'show_simulated_data': self.show_simulated_data.get(),
                'show_simulated_markers': self.show_simulated_markers.get(),
                'simulated_marker_shape': self.simulated_marker_shape.get(),
                'simulated_marker_fill': self.simulated_marker_fill.get(),
                'simulated_marker_color': self.simulated_marker_color.get(),
                'simulated_marker_size': self.simulated_marker_size.get(),
                'simulated_marker_edge_width': self.simulated_marker_edge_width.get(),
                'simulated_marker_range_enabled': self.simulated_marker_range_enabled.get(),
                'simulated_marker_range_mode': self.simulated_marker_range_mode.get(),
                'simulated_marker_range_min': self.simulated_marker_range_min.get(),
                'simulated_marker_range_max': self.simulated_marker_range_max.get(),
                'show_simulated_line': self.show_simulated_line.get(),
                'simulated_line_color': self.simulated_line_color.get(),
                'simulated_line_width': self.simulated_line_width.get(),
                'simulated_line_alpha': self.simulated_line_alpha.get(),
                'background_color': self.background_color.get(),
                'show_grid': self.show_grid.get(),
                'grid_color': self.grid_color.get(),
                'title_size': self.title_size.get(),
                'axis_label_size': self.axis_label_size.get(),
                'tick_label_size': self.tick_label_size.get(),
                'legend_size': self.legend_size.get(),
                'info_box_size': self.info_box_size.get(),
                'title_color': self.title_color.get(),
                'axis_label_color': self.axis_label_color.get(),
                'tick_label_color': self.tick_label_color.get(),
                'legend_text_color': self.legend_text_color.get(),
                'info_box_text_color': self.info_box_text_color.get(),
                'tick_length': self.tick_length.get(),
                'tick_width': self.tick_width.get(),
                'tick_color': self.tick_color.get(),
                'tick_direction': self.tick_direction.get(),
                'x_axis_mode': self.x_axis_mode.get(),
                'channel_min': self.get_optional_float(self.channel_min),
                'channel_max': self.get_optional_float(self.channel_max),
                'energy_min': self.get_optional_float(self.energy_min),
                'energy_max': self.get_optional_float(self.energy_max),
                'counts_min': self.get_optional_float(self.counts_min),
                'counts_max': self.get_optional_float(self.counts_max),
                'show_plot': True
            }
            
            self.last_fig = plot_xnra_spectrum(**params)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate plot:\n{str(e)}")
    
    def save_preset(self):
        filename = filedialog.asksaveasfilename(
            title="Save Preset", defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if filename:
            try:
                preset = {'uncertainty_technique': self.uncertainty_technique.get(),
                          'systematic_uncertainty': self.systematic_uncertainty.get()}
                with open(filename, 'w') as f:
                    json.dump(preset, f, indent=4)
                messagebox.showinfo("Success", f"Preset saved!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save: {e}")
    
    def load_preset(self):
        filename = filedialog.askopenfilename(
            title="Load Preset", filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if filename:
            try:
                with open(filename, 'r') as f:
                    preset = json.load(f)
                if 'uncertainty_technique' in preset:
                    self.uncertainty_technique.set(preset['uncertainty_technique'])
                if 'systematic_uncertainty' in preset:
                    self.systematic_uncertainty.set(preset['systematic_uncertainty'])
                messagebox.showinfo("Success", "Preset loaded!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load: {e}")
    
    def export_plot(self):
        if self.last_fig is None:
            messagebox.showerror("Error", "Generate a plot first!")
            return
        filename = filedialog.asksaveasfilename(
            title="Export Plot", defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("All", "*.*")])
        if filename:
            try:
                self.last_fig.savefig(filename, dpi=self.fig_dpi.get(), bbox_inches='tight')
                messagebox.showinfo("Success", f"Plot exported!")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {e}")
    
    def load_defaults(self):
        pass


def main():
    root = tk.Tk()
    app = XNRAViewerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
