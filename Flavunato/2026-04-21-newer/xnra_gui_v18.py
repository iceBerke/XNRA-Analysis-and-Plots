### XNRA Spectrum Viewer GUI - Version 15 ###
# Tkinter-based graphical interface for XNRA spectrum visualization
# VERSION 15: 
#   - Fixed Q·Ω and σ display: now shows CONSTANT relative values (not scaled by concentration)
#   - Added isotope-specific uncertainty calculations
#   - Per-isotope cross-section selection
#   - NEW: Element/Iso Data tab for viewing/plotting element & isotope spectra
#   - Calibration extracted from .xnra file (no manual entry needed)
#   - Single file upload location (Element/Iso Data tab serves both visualization and uncertainty)
#   - Compatible with xnra_main_v6 and xnra_uncertainty_v9

# Script developed by Berke Santos & Giuseppe Legrottaglie
# Last updated: 2025

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import json
from pathlib import Path
import re

# Import plotting function
try:
    from xnra_main_v6 import plot_xnra_spectrum
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("Warning: xnra_main_v6 not found. Plotting disabled.")

# Import element spectrum functions (dedicated module)
try:
    from xnra_element_spectra_v1 import (
        parse_element_spectrum,
        get_element_spectrum_summary,
        get_available_elements,
        get_available_isotopes,
        calculate_element_counts_from_spectrum,
        plot_element_spectra,
        ELEMENT_COLORS,
    )
    ELEMENT_SPECTRUM_AVAILABLE = True
except ImportError:
    ELEMENT_SPECTRUM_AVAILABLE = False
    ELEMENT_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    print("Warning: xnra_element_spectra not found. Element spectrum features disabled.")

# Import uncertainty calculation (v9 with isotope support)
try:
    from xnra_uncertainty_v12 import (
        parse_xnra_full,
        calculate_chi_squared_for_roi,
        calculate_uncertainty_v9,
        create_depth_profile,
        interpret_chi2,
        CROSS_SECTION_UNCERTAINTIES,
        ISOTOPE_CROSS_SECTION_DEFAULTS,
        get_cross_section_uncertainty,
        get_isotope_cross_section_uncertainty,
    )
    UNCERTAINTY_AVAILABLE = True
except ImportError:
    # Fallback defaults if module not available
    CROSS_SECTION_UNCERTAINTIES = {
        'RBS: Mx(α,α)Mx': (0.00, 'Rutherford backscattering'),
        'NRA: ¹⁶O(d,α₀)¹⁴N': (0.02, 'Deuteron NRA for oxygen'),
        'NRA: ¹²C(d,p₀)¹³C': (0.05, 'Deuteron NRA for carbon'),
        'ERDA: ¹H(⁴He,¹H)⁴He': (0.02, 'Forward recoil for hydrogen'),
    }
    ISOTOPE_CROSS_SECTION_DEFAULTS = {}
    UNCERTAINTY_AVAILABLE = False
    ELEMENT_SPECTRUM_AVAILABLE = False
    print("Warning: xnra_uncertainty_v9 not found. Uncertainty calculations disabled.")


class XNRAViewerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("XNRA Spectrum Viewer - v15")
        self.root.geometry("1400x900")
        
        # Store current filepath
        self.current_file = tk.StringVar()
        self.last_fig = None
        
        # Store ROIs for uncertainty calculation
        self.roi_list = []
        self.uncertainty_results = None
        self.uncertainty_results_df = None
        
        # Store element spectrum data for simulation-guided A_A calculation
        self.element_spectrum_data = None
        self.element_spectrum_file = tk.StringVar()
        
        # Isotope-specific cross-section overrides
        self.isotope_xsec_overrides = {}
        self.isotope_override_widgets = []
        
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
        self.create_element_spectrum_tab()  # NEW: Element/Iso Data tab
        self.create_appearance_tab()
        self.create_axes_tab()
        self.create_uncertainty_tab_v9()
    
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
        
        stack = ttk.Frame(scrollable_frame)
        stack.grid(row=0, column=0, sticky="nsew", padx=(8, 0), pady=(0, 10))
        stack.columnconfigure(0, weight=1)

        title_frame = ttk.LabelFrame(stack, text="Plot Title", padding="10")
        title_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        
        ttk.Label(title_frame, text="Custom Title:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        self.custom_title = tk.StringVar(value="")
        ttk.Entry(title_frame, textvariable=self.custom_title, width=50).grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Label(title_frame, text="(empty = filename without extension)", font=('TkDefaultFont', 8, 'italic')).grid(row=1, column=0, sticky=tk.W, padx=(2, 0))

        fig_frame = ttk.LabelFrame(stack, text="Figure Settings", padding="10")
        fig_frame.grid(row=1, column=0, sticky="ew")
        
        self.fig_width = self.create_number_field(fig_frame, "Width (inches):", 0, 12.0)
        self.fig_height = self.create_number_field(fig_frame, "Height (inches):", 1, 6.0)
        self.fig_dpi = self.create_number_field(fig_frame, "DPI:", 2, 100, is_int=True)
        
        display_frame = ttk.LabelFrame(scrollable_frame, text="Display Options", padding="10")
        display_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
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
        
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.raw_marker_range_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.raw_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.raw_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.raw_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.raw_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 0.0)
        self.raw_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 200.0)
        
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
        
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.smoothed_marker_range_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.smoothed_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.smoothed_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.smoothed_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.smoothed_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 200.0)
        self.smoothed_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 400.0)
        
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
        
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.simulated_marker_range_enabled = tk.BooleanVar(value=False)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.simulated_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.simulated_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.simulated_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.simulated_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 0.0)
        self.simulated_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 1000.0)
        
        line_frame = ttk.LabelFrame(scrollable_frame, text="Line", padding="10")
        line_frame.grid(row=1, column=2, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.show_simulated_line = tk.BooleanVar(value=True)
        ttk.Checkbutton(line_frame, text="Show Line", variable=self.show_simulated_line).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        self.simulated_line_color = self.create_color_field(line_frame, "Color:", 1, 'green')
        self.simulated_line_width = self.create_number_field(line_frame, "Width:", 2, 2.0)
        self.simulated_line_alpha = self.create_number_field(line_frame, "Transparency:", 3, 0.8)
    
    def create_element_spectrum_tab(self):
        """Element/Iso Data tab - load and visualize element/isotope spectra"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Element/Iso Data")
        
        # Configure grid
        tab.columnconfigure(0, weight=1)
        tab.columnconfigure(1, weight=2)
        tab.rowconfigure(0, weight=1)
        
        left_panel = ttk.Frame(tab)
        left_panel.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        left_panel.columnconfigure(0, weight=1)
        
        # File upload section
        file_frame = ttk.LabelFrame(left_panel, text="Element/Iso Data File", padding="10")
        file_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        file_frame.columnconfigure(0, weight=1)
        
        info_text = ttk.Label(file_frame, 
            text="Load SIMNRA element spectrum file:\n"
                 "Plot → Show Element Spectra → File → Write Spectrum Data",
            font=('TkDefaultFont', 9), foreground='gray', justify=tk.LEFT)
        info_text.grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=(0, 10))
        
        btn_frame = ttk.Frame(file_frame)
        btn_frame.grid(row=1, column=0, sticky="ew")
        
        ttk.Button(btn_frame, text="Load File...", 
                   command=self.load_element_spectrum_for_viewer).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(btn_frame, text="Clear", 
                   command=self.clear_element_spectrum_viewer).pack(side=tk.LEFT)
        
        self.elem_viewer_file = tk.StringVar(value="No file loaded")
        ttk.Label(file_frame, textvariable=self.elem_viewer_file, 
                  font=('TkDefaultFont', 8, 'italic'), foreground='#666').grid(
            row=2, column=0, sticky=tk.W, pady=(5, 0))
        
        # Channel range info
        self.elem_channel_info = tk.StringVar(value="")
        ttk.Label(file_frame, textvariable=self.elem_channel_info,
                  font=('TkDefaultFont', 8)).grid(row=3, column=0, sticky=tk.W)
        
        elem_list_frame = ttk.LabelFrame(left_panel, text="Elements Detected", padding="10")
        elem_list_frame.grid(row=1, column=0, sticky="nsew", pady=(0, 10))
        elem_list_frame.columnconfigure(0, weight=1)
        elem_list_frame.rowconfigure(0, weight=1)
        
        elem_columns = ('Element', 'Total Counts')
        self.elem_tree = ttk.Treeview(elem_list_frame, columns=elem_columns, 
                                       show='headings', height=6)
        self.elem_tree.heading('Element', text='Element')
        self.elem_tree.heading('Total Counts', text='Total Counts')
        self.elem_tree.column('Element', width=80, anchor='center')
        self.elem_tree.column('Total Counts', width=100, anchor='center')
        
        elem_scroll = ttk.Scrollbar(elem_list_frame, orient="vertical", 
                                     command=self.elem_tree.yview)
        self.elem_tree.configure(yscrollcommand=elem_scroll.set)
        
        self.elem_tree.grid(row=0, column=0, sticky="nsew")
        elem_scroll.grid(row=0, column=1, sticky="ns")
        
        iso_list_frame = ttk.LabelFrame(left_panel, text="Isotopes Detected", padding="10")
        iso_list_frame.grid(row=2, column=0, sticky="nsew", pady=(0, 10))
        iso_list_frame.columnconfigure(0, weight=1)
        iso_list_frame.rowconfigure(0, weight=1)
        
        iso_columns = ('Isotope', 'Parent', 'Total Counts')
        self.iso_tree = ttk.Treeview(iso_list_frame, columns=iso_columns, 
                                      show='headings', height=6)
        self.iso_tree.heading('Isotope', text='Isotope')
        self.iso_tree.heading('Parent', text='Parent')
        self.iso_tree.heading('Total Counts', text='Total Counts')
        self.iso_tree.column('Isotope', width=70, anchor='center')
        self.iso_tree.column('Parent', width=60, anchor='center')
        self.iso_tree.column('Total Counts', width=90, anchor='center')
        
        iso_scroll = ttk.Scrollbar(iso_list_frame, orient="vertical", 
                                    command=self.iso_tree.yview)
        self.iso_tree.configure(yscrollcommand=iso_scroll.set)
        
        self.iso_tree.grid(row=0, column=0, sticky="nsew")
        iso_scroll.grid(row=0, column=1, sticky="ns")
        
        right_panel = ttk.LabelFrame(tab, text="Spectrum Visualization", padding="10")
        right_panel.grid(row=0, column=1, sticky="nsew")
        right_panel.columnconfigure(0, weight=1)
        right_panel.columnconfigure(1, weight=1)
        
        # Selection frame
        selection_frame = ttk.Frame(right_panel)
        selection_frame.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        # Elements selection
        elem_sel_frame = ttk.LabelFrame(selection_frame, text="Select Elements", padding="5")
        elem_sel_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))
        
        elem_btn_frame = ttk.Frame(elem_sel_frame)
        elem_btn_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Button(elem_btn_frame, text="All", width=5,
                   command=lambda: self.select_all_plot_items('elements')).pack(side=tk.LEFT, padx=1)
        ttk.Button(elem_btn_frame, text="None", width=5,
                   command=lambda: self.clear_all_plot_items('elements')).pack(side=tk.LEFT, padx=1)
        
        self.element_checkboxes_frame = ttk.Frame(elem_sel_frame)
        self.element_checkboxes_frame.pack(fill=tk.BOTH, expand=True)
        self.element_plot_vars = {}  # Will hold {element_name: BooleanVar}
        
        # Isotopes selection
        iso_sel_frame = ttk.LabelFrame(selection_frame, text="Select Isotopes", padding="5")
        iso_sel_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(5, 0))
        
        iso_btn_frame = ttk.Frame(iso_sel_frame)
        iso_btn_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Button(iso_btn_frame, text="All", width=5,
                   command=lambda: self.select_all_plot_items('isotopes')).pack(side=tk.LEFT, padx=1)
        ttk.Button(iso_btn_frame, text="None", width=5,
                   command=lambda: self.clear_all_plot_items('isotopes')).pack(side=tk.LEFT, padx=1)
        
        self.isotope_checkboxes_frame = ttk.Frame(iso_sel_frame)
        self.isotope_checkboxes_frame.pack(fill=tk.BOTH, expand=True)
        self.isotope_plot_vars = {}  # Will hold {isotope_name: BooleanVar}
        
        # Plot options
        options_frame = ttk.LabelFrame(right_panel, text="Plot Options", padding="10")
        options_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(0, 10))

        # Reference-spectra overlays (from the SIMNRA file's reserved columns)
        self.show_sim_overlay = tk.BooleanVar(value=False)
        self.show_raw_overlay = tk.BooleanVar(value=False)
        self.show_smoothed_overlay = tk.BooleanVar(value=False)
        ttk.Checkbutton(options_frame, text="Show Simulated",
                        variable=self.show_sim_overlay).grid(row=0, column=0, sticky=tk.W)
        ttk.Checkbutton(options_frame, text="Show Raw (Experimental)",
                        variable=self.show_raw_overlay).grid(row=0, column=1, sticky=tk.W, padx=(15, 0))
        ttk.Checkbutton(options_frame, text="Show Smoothed",
                        variable=self.show_smoothed_overlay).grid(row=0, column=2, sticky=tk.W, padx=(15, 0))

        # X-axis mode
        self.elem_x_axis_mode = tk.StringVar(value='channel')
        xaxis_frame = ttk.Frame(options_frame)
        xaxis_frame.grid(row=1, column=0, columnspan=4, sticky=tk.W, pady=(8, 0))
        ttk.Label(xaxis_frame, text="X-axis:").pack(side=tk.LEFT)
        ttk.Radiobutton(xaxis_frame, text="Channel", variable=self.elem_x_axis_mode,
                        value='channel').pack(side=tk.LEFT, padx=(5, 10))
        ttk.Radiobutton(xaxis_frame, text="Energy", variable=self.elem_x_axis_mode,
                        value='energy').pack(side=tk.LEFT)

        # Axis limits (empty = auto)
        limits_frame = ttk.LabelFrame(options_frame, text="Axis Limits (blank = auto)", padding="5")
        limits_frame.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(8, 0))

        self.elem_channel_min = tk.StringVar(value="")
        self.elem_channel_max = tk.StringVar(value="")
        self.elem_energy_min = tk.StringVar(value="")
        self.elem_energy_max = tk.StringVar(value="")
        self.elem_counts_min = tk.StringVar(value="")
        self.elem_counts_max = tk.StringVar(value="")

        def _lim_row(r, label, vmin, vmax, unit):
            ttk.Label(limits_frame, text=label).grid(row=r, column=0, sticky=tk.W)
            ttk.Label(limits_frame, text="min:").grid(row=r, column=1, sticky=tk.E, padx=(10, 2))
            ttk.Entry(limits_frame, textvariable=vmin, width=10).grid(row=r, column=2, sticky=tk.W)
            ttk.Label(limits_frame, text="max:").grid(row=r, column=3, sticky=tk.E, padx=(10, 2))
            ttk.Entry(limits_frame, textvariable=vmax, width=10).grid(row=r, column=4, sticky=tk.W)
            ttk.Label(limits_frame, text=unit).grid(row=r, column=5, sticky=tk.W, padx=(2, 0))

        _lim_row(0, "Channel:", self.elem_channel_min, self.elem_channel_max, "ch")
        _lim_row(1, "Energy:",  self.elem_energy_min,  self.elem_energy_max,  "keV")
        _lim_row(2, "Counts:",  self.elem_counts_min,  self.elem_counts_max,  "cts")

        # Internal cal vars (populated by browse_file, not user-editable)
        self.elem_cal_offset = tk.DoubleVar(value=0.0)
        self.elem_cal_gain = tk.DoubleVar(value=1.0)

        # Plot title
        title_frame = ttk.Frame(options_frame)
        title_frame.grid(row=3, column=0, columnspan=4, sticky=tk.W, pady=(8, 0))
        ttk.Label(title_frame, text="Title:").pack(side=tk.LEFT)
        self.elem_plot_title = tk.StringVar(value="Element Spectra")
        ttk.Entry(title_frame, textvariable=self.elem_plot_title, width=40).pack(side=tk.LEFT, padx=(5, 0))
        
        # Plot button
        plot_btn_frame = ttk.Frame(right_panel)
        plot_btn_frame.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        ttk.Button(plot_btn_frame, text="📊 Plot Element Spectra", 
                   command=self.plot_element_spectra_gui).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(plot_btn_frame, text="📊 Plot Isotope Spectra", 
                   command=self.plot_isotope_spectra_gui).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(plot_btn_frame, text="📊 Plot All Selected", 
                   command=self.plot_all_selected_spectra).pack(side=tk.LEFT)
        
        # Note at bottom
        note_label = ttk.Label(right_panel, 
            text="💡 Hover over lines to identify elements (requires: pip install mplcursors). Also used for A_A calculations.",
            font=('TkDefaultFont', 8, 'italic'), foreground='#666')
        note_label.grid(row=3, column=0, columnspan=2, sticky=tk.W, pady=(10, 0))
    
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
        
        color_frame = ttk.LabelFrame(scrollable_frame, text="Colors", padding="10")
        color_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5), pady=(0, 10))
        
        self.background_color = self.create_color_field(color_frame, "Background:", 0, 'white')
        self.grid_color = self.create_color_field(color_frame, "Grid:", 1, 'white')
        
        self.show_grid = tk.BooleanVar(value=True)
        ttk.Checkbutton(color_frame, text="Show Grid", variable=self.show_grid).grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(5, 0))
        
        size_frame = ttk.LabelFrame(scrollable_frame, text="Text Sizes", padding="10")
        size_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
        self.title_size = self.create_number_field(size_frame, "Title:", 0, 16.0)
        self.axis_label_size = self.create_number_field(size_frame, "Axis Labels:", 1, 12.0)
        self.tick_label_size = self.create_number_field(size_frame, "Tick Labels:", 2, 10.0)
        self.legend_size = self.create_number_field(size_frame, "Legend:", 3, 10.0)
        self.info_box_size = self.create_number_field(size_frame, "Info Box:", 4, 10.0)
        
        text_color_frame = ttk.LabelFrame(scrollable_frame, text="Text Colors", padding="10")
        text_color_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.title_color = self.create_color_field(text_color_frame, "Title:", 0, 'black')
        self.axis_label_color = self.create_color_field(text_color_frame, "Axis Labels:", 1, 'black')
        self.tick_label_color = self.create_color_field(text_color_frame, "Tick Labels:", 2, 'black')
        self.legend_text_color = self.create_color_field(text_color_frame, "Legend:", 3, 'black')
        self.info_box_text_color = self.create_color_field(text_color_frame, "Info Box:", 4, 'black')
        
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
        
        mode_frame = ttk.LabelFrame(scrollable_frame, text="X-Axis Mode", padding="10")
        mode_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5), pady=(0, 10))
        
        ttk.Label(mode_frame, text="Display Mode:").grid(row=0, column=0, sticky=tk.W)
        self.x_axis_mode = tk.StringVar(value='energy')
        ttk.Combobox(mode_frame, textvariable=self.x_axis_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=0, column=1, sticky=tk.W, padx=(5, 0))
        
        channel_frame = ttk.LabelFrame(scrollable_frame, text="Channel Limits", padding="10")
        channel_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
        self.channel_min = self.create_optional_number_field(channel_frame, "Min:", 0)
        self.channel_max = self.create_optional_number_field(channel_frame, "Max:", 1)
        
        energy_frame = ttk.LabelFrame(scrollable_frame, text="Energy Limits (keV)", padding="10")
        energy_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.energy_min = self.create_optional_number_field(energy_frame, "Min:", 0)
        self.energy_max = self.create_optional_number_field(energy_frame, "Max:", 1)
        
        counts_frame = ttk.LabelFrame(scrollable_frame, text="Counts Limits (Y-axis)", padding="10")
        counts_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.counts_min = self.create_optional_number_field(counts_frame, "Min:", 0)
        self.counts_max = self.create_optional_number_field(counts_frame, "Max:", 1)
    
    # UNCERTAINTY TAB v9 - French Expert Formulas + Isotope Support
    
    def create_uncertainty_tab_v9(self):
        """
        Uncertainty calculation tab implementing French IBA expert formulas
        with isotope support and FIXED Q·Ω display (constant relative values).
        """
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Uncertainty (v9)")
        
        tab.columnconfigure(0, weight=1)
        tab.columnconfigure(1, weight=1)
        tab.rowconfigure(0, weight=1)
        
        left_container = ttk.Frame(tab)
        left_container.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        left_container.rowconfigure(0, weight=1)
        left_container.columnconfigure(0, weight=1)
        
        left_canvas = tk.Canvas(left_container, highlightthickness=0)
        left_scrollbar = ttk.Scrollbar(left_container, orient="vertical", command=left_canvas.yview)
        left_panel = ttk.Frame(left_canvas)
        
        left_panel.bind("<Configure>", lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))
        left_canvas.create_window((0, 0), window=left_panel, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scrollbar.set)
        
        left_canvas.grid(row=0, column=0, sticky="nsew")
        left_scrollbar.grid(row=0, column=1, sticky="ns")
        
        def _on_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        def _bind_mousewheel(event):
            left_canvas.bind_all("<MouseWheel>", _on_mousewheel)
        def _unbind_mousewheel(event):
            left_canvas.unbind_all("<MouseWheel>")
        left_canvas.bind("<Enter>", _bind_mousewheel)
        left_canvas.bind("<Leave>", _unbind_mousewheel)
        
        def _configure_canvas(event):
            left_canvas.itemconfig(left_canvas.find_all()[0], width=event.width)
        left_canvas.bind("<Configure>", _configure_canvas)
        
        left_panel.columnconfigure(0, weight=1)
        
        inputs_frame = ttk.LabelFrame(left_panel, text="Uncertainty Components (CONSTANT)", padding="10")
        inputs_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        inputs_frame.columnconfigure(1, weight=1)
        
        # Explanation banner
        const_info = ttk.Label(inputs_frame, 
            text="⚠ Q·Ω and σ are RELATIVE uncertainties — constant for all elements/layers",
            font=('TkDefaultFont', 8, 'bold'), foreground='#0369A1')
        const_info.grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=(0, 10))
        
        # 1. Q·Ω uncertainty
        ttk.Label(inputs_frame, text="Q·Ω Uncertainty (%):", 
                  font=('TkDefaultFont', 9, 'bold')).grid(row=1, column=0, sticky=tk.W, padx=(0, 10))
        self.qomega_uncertainty = tk.DoubleVar(value=3.0)
        qomega_entry = ttk.Entry(inputs_frame, textvariable=self.qomega_uncertainty, width=10)
        qomega_entry.grid(row=1, column=1, sticky=tk.W)
        
        qomega_info = ttk.Label(inputs_frame, 
                                text="Detector/setup property — same for entire spectrum",
                                font=('TkDefaultFont', 8, 'italic'), foreground='gray')
        qomega_info.grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(0, 10))
        
        # 2. Cross-section reaction dropdown
        ttk.Label(inputs_frame, text="Default Reaction:", 
                  font=('TkDefaultFont', 9, 'bold')).grid(row=3, column=0, sticky=tk.W, padx=(0, 10))
        
        self.reaction_options = list(CROSS_SECTION_UNCERTAINTIES.keys()) + ['Custom...']
        self.selected_reaction = tk.StringVar(value='RBS: Mx(α,α)Mx')
        
        reaction_combo = ttk.Combobox(inputs_frame, textvariable=self.selected_reaction,
                                       values=self.reaction_options, state='readonly', width=25)
        reaction_combo.grid(row=3, column=1, sticky=tk.W)
        reaction_combo.bind('<<ComboboxSelected>>', self.on_reaction_changed)
        
        # Cross-section display
        self.xsec_frame = ttk.Frame(inputs_frame)
        self.xsec_frame.grid(row=4, column=0, columnspan=2, sticky=tk.W, pady=(5, 0))
        
        ttk.Label(self.xsec_frame, text="Cross-section σ (%):").pack(side=tk.LEFT)
        self.xsec_uncertainty = tk.DoubleVar(value=0.0)
        self.xsec_entry = ttk.Entry(self.xsec_frame, textvariable=self.xsec_uncertainty, 
                                     width=8, state='readonly')
        self.xsec_entry.pack(side=tk.LEFT, padx=(5, 10))
        
        self.xsec_desc = tk.StringVar(value="(Rutherford backscattering)")
        ttk.Label(self.xsec_frame, textvariable=self.xsec_desc, 
                  font=('TkDefaultFont', 8, 'italic'), foreground='gray').pack(side=tk.LEFT)
        
        ttk.Separator(inputs_frame, orient='horizontal').grid(row=5, column=0, columnspan=2, 
                                                               sticky='ew', pady=10)
        
        # 3. Statistical uncertainty info
        stat_label = ttk.Label(inputs_frame, text="Statistical (Poisson):", 
                               font=('TkDefaultFont', 9, 'bold'))
        stat_label.grid(row=6, column=0, sticky=tk.W)
        ttk.Label(inputs_frame, text="VARIES — depends on counts per element",
                  font=('TkDefaultFont', 8, 'italic'), foreground='#B45309').grid(row=6, column=1, sticky=tk.W)
        
        ttk.Label(inputs_frame, text="Formula: ΔA_A = 1/√N  (N = element counts in ROI)",
                  font=('TkDefaultFont', 8)).grid(row=7, column=0, columnspan=2, sticky=tk.W, pady=(2, 0))
        
        # Isotope mode has been removed from the UI. Keep variables defined
        # so downstream code (override helpers, calc loop) still works; they
        # just always report "disabled / no overrides".
        self.isotope_mode = tk.BooleanVar(value=False)
        self.isotope_override_container = ttk.Frame(left_panel)  # orphan, not gridded

        formula_frame = ttk.LabelFrame(left_panel, text="Adjustment Method", padding="10")
        formula_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))
        
        self.adjustment_method = tk.StringVar(value='none')
        
        standard_frame = ttk.Frame(formula_frame)
        standard_frame.pack(fill=tk.X, pady=(0, 3))
        
        ttk.Radiobutton(standard_frame, text="Standard (Eq. 3)", 
                        variable=self.adjustment_method, value='none',
                        command=self.update_formula_display).pack(side=tk.LEFT)
        
        ttk.Label(standard_frame, text="ΔC = √[(ΔQ·Ω)² + (Δσ)² + (1/√N)²]",
                  font=('Consolas', 8), foreground='gray').pack(side=tk.LEFT, padx=(10, 0))
        
        chi2_frame_radio = ttk.Frame(formula_frame)
        chi2_frame_radio.pack(fill=tk.X, pady=(3, 3))
        
        ttk.Radiobutton(chi2_frame_radio, text="χ²-Adjusted (Eq. 6)", 
                        variable=self.adjustment_method, value='chi2',
                        command=self.update_formula_display).pack(side=tk.LEFT)
        
        ttk.Label(chi2_frame_radio, text="... × √χ²/ν",
                  font=('Consolas', 8), foreground='gray').pack(side=tk.LEFT, padx=(10, 0))
        
        auc_frame_radio = ttk.Frame(formula_frame)
        auc_frame_radio.pack(fill=tk.X, pady=(3, 0))
        
        ttk.Radiobutton(auc_frame_radio, text="AUC-Adjusted", 
                        variable=self.adjustment_method, value='auc',
                        command=self.update_formula_display).pack(side=tk.LEFT)
        
        ttk.Label(auc_frame_radio, text="... × max(AUC, 1/AUC)",
                  font=('Consolas', 8), foreground='gray').pack(side=tk.LEFT, padx=(10, 0))
        
        self.formula_explanation = tk.StringVar(value="Using standard formula (no adjustment)")
        ttk.Label(formula_frame, textvariable=self.formula_explanation, 
                  font=('TkDefaultFont', 8, 'italic'), foreground='#666666').pack(anchor=tk.W, pady=(10, 0))

        # Layer selection: populated dynamically when an .xnra file is loaded
        layers_frame = ttk.LabelFrame(left_panel, text="Layers to Calculate", padding="10")
        layers_frame.grid(row=25, column=0, sticky="ew", pady=(0, 10))

        ttk.Label(layers_frame,
                  text="Select which sample layers to include in the calculation.",
                  font=('TkDefaultFont', 8, 'italic'), foreground='gray').pack(anchor=tk.W, pady=(0, 5))

        layers_btn_frame = ttk.Frame(layers_frame)
        layers_btn_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Button(layers_btn_frame, text="All", width=5,
                   command=lambda: self._set_all_layers(True)).pack(side=tk.LEFT, padx=1)
        ttk.Button(layers_btn_frame, text="None", width=5,
                   command=lambda: self._set_all_layers(False)).pack(side=tk.LEFT, padx=1)

        self.layers_checkbox_frame = ttk.Frame(layers_frame)
        self.layers_checkbox_frame.pack(fill=tk.X)
        self.layer_vars = {}  # {layer_num: BooleanVar}

        self._layers_placeholder = ttk.Label(self.layers_checkbox_frame,
            text="(load an .xnra file to populate)",
            font=('TkDefaultFont', 8, 'italic'), foreground='#999')
        self._layers_placeholder.pack(anchor=tk.W)

        elem_spec_frame = ttk.LabelFrame(left_panel, text="Element/Isotope Spectrum", padding="10")
        elem_spec_frame.grid(row=3, column=0, sticky="ew", pady=(0, 10))
        
        elem_spec_info = ttk.Label(elem_spec_frame, 
            text="💡 Load element spectrum in the 'Element/Iso Data' tab.\n"
                 "Required for accurate A_A (element counts) calculation.",
            font=('TkDefaultFont', 8), foreground='gray', justify=tk.LEFT)
        elem_spec_info.pack(anchor=tk.W, pady=(0, 5))
        
        self.elem_spec_status = tk.StringVar(value="⚠ Not loaded — go to Element/Iso Data tab")
        self.elem_spec_status_label = ttk.Label(elem_spec_frame, textvariable=self.elem_spec_status,
                  font=('TkDefaultFont', 8, 'italic'), foreground='#B45309')
        self.elem_spec_status_label.pack(anchor=tk.W)
        
        roi_frame = ttk.LabelFrame(left_panel, text="Regions of Interest (ROIs)", padding="10")
        roi_frame.grid(row=4, column=0, sticky="nsew", pady=(0, 10))
        roi_frame.columnconfigure(0, weight=1)
        roi_frame.rowconfigure(1, weight=1)
        
        header_frame = ttk.Frame(roi_frame)
        header_frame.grid(row=0, column=0, sticky="ew", pady=(0, 5))
        
        ttk.Label(header_frame, text="Name", width=12).pack(side=tk.LEFT, padx=2)
        ttk.Label(header_frame, text="Start Ch", width=8).pack(side=tk.LEFT, padx=2)
        ttk.Label(header_frame, text="End Ch", width=8).pack(side=tk.LEFT, padx=2)
        ttk.Label(header_frame, text="Elements/Isotopes", width=22).pack(side=tk.LEFT, padx=2)
        
        roi_canvas = tk.Canvas(roi_frame, height=100)
        roi_scrollbar = ttk.Scrollbar(roi_frame, orient="vertical", command=roi_canvas.yview)
        self.roi_container = ttk.Frame(roi_canvas)
        
        self.roi_container.bind("<Configure>", lambda e: roi_canvas.configure(scrollregion=roi_canvas.bbox("all")))
        roi_canvas.create_window((0, 0), window=self.roi_container, anchor="nw")
        roi_canvas.configure(yscrollcommand=roi_scrollbar.set)
        
        roi_canvas.grid(row=1, column=0, sticky="nsew")
        roi_scrollbar.grid(row=1, column=1, sticky="ns")
        
        roi_btn_frame = ttk.Frame(roi_frame)
        roi_btn_frame.grid(row=2, column=0, sticky="ew", pady=(10, 0))
        
        ttk.Button(roi_btn_frame, text="+ Add ROI", command=self.add_roi_row).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(roi_btn_frame, text="Load Defaults", command=self.load_default_rois).pack(side=tk.LEFT, padx=5)
        ttk.Button(roi_btn_frame, text="Clear All", command=self.clear_all_rois).pack(side=tk.LEFT, padx=5)
        
        calc_frame = ttk.Frame(left_panel)
        calc_frame.grid(row=5, column=0, sticky="ew")
        
        ttk.Button(calc_frame, text="Calculate Uncertainties", 
                   command=self.calculate_uncertainties_v9).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(calc_frame, text="Export Results", 
                   command=self.export_uncertainty_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(calc_frame, text="Export Depth Profile", 
                   command=self.export_depth_profile).pack(side=tk.LEFT, padx=5)
        
        right_panel = ttk.LabelFrame(tab, text="Results", padding="10")
        right_panel.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
        right_panel.rowconfigure(0, weight=1)
        right_panel.columnconfigure(0, weight=1)
        
        # Results treeview with CLARIFIED column headers
        # Q·Ω(rel) and σ(rel) indicate these are CONSTANT relative values
        columns = ('Layer', 'Element', 'Conc(%)', 'Stat(rel)', 'Q·Ω(rel)', 'σ(rel)', '±Total', 'ROI', 'A_A', 'f_A', 'χ²/ν', 'Status', 'AUC')
        self.results_tree = ttk.Treeview(right_panel, columns=columns, show='headings', height=20)
        
        col_widths = {
            'Layer': 42, 'Element': 60, 'Conc(%)': 55, 
            'Stat(rel)': 75, 'Q·Ω(rel)': 75, 'σ(rel)': 70, '±Total': 75,
            'ROI': 55, 'A_A': 55, 'f_A': 45, 'χ²/ν': 45, 'Status': 55, 'AUC': 45
        }
        
        for col in columns:
            self.results_tree.heading(col, text=col)
            self.results_tree.column(col, width=col_widths.get(col, 55), anchor='center')
        
        results_yscroll = ttk.Scrollbar(right_panel, orient="vertical", command=self.results_tree.yview)
        results_xscroll = ttk.Scrollbar(right_panel, orient="horizontal", command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=results_yscroll.set, xscrollcommand=results_xscroll.set)
        
        self.results_tree.grid(row=0, column=0, sticky="nsew")
        results_yscroll.grid(row=0, column=1, sticky="ns")
        results_xscroll.grid(row=1, column=0, sticky="ew")
        
        # Legend explaining columns
        legend_text = "Stat(rel), Q·Ω(rel) and σ(rel) are RELATIVE uncertainties (%). ±Total is ABSOLUTE (% of sample)."
        legend_label = ttk.Label(right_panel, text=legend_text, font=('TkDefaultFont', 8, 'italic'), foreground='#0369A1')
        legend_label.grid(row=2, column=0, sticky="w", pady=(5, 0))
        
        self.uncertainty_status = ttk.Label(right_panel, 
                                             text="Load a file and define ROIs to calculate uncertainties", 
                                             font=('TkDefaultFont', 9, 'italic'))
        self.uncertainty_status.grid(row=3, column=0, sticky="w", pady=(5, 0))
        
        self.load_default_rois()
    
    def add_isotope_override_row(self, isotope="", uncertainty=""):
        """Add a row for isotope-specific cross-section override"""
        row_frame = ttk.Frame(self.isotope_override_container)
        row_frame.pack(fill=tk.X, pady=2)
        
        isotope_var = tk.StringVar(value=isotope)
        unc_var = tk.StringVar(value=uncertainty)
        
        ttk.Label(row_frame, text="Isotope:").pack(side=tk.LEFT, padx=(0, 2))
        ttk.Entry(row_frame, textvariable=isotope_var, width=8).pack(side=tk.LEFT, padx=2)
        
        ttk.Label(row_frame, text="σ (%):").pack(side=tk.LEFT, padx=(10, 2))
        ttk.Entry(row_frame, textvariable=unc_var, width=6).pack(side=tk.LEFT, padx=2)
        
        def remove_row():
            row_frame.destroy()
            self.isotope_override_widgets = [w for w in self.isotope_override_widgets if w['frame'] != row_frame]
        
        ttk.Button(row_frame, text="×", width=3, command=remove_row).pack(side=tk.LEFT, padx=5)
        
        self.isotope_override_widgets.append({
            'frame': row_frame,
            'isotope': isotope_var,
            'uncertainty': unc_var
        })
    
    def clear_isotope_overrides(self):
        """Clear all isotope override rows"""
        for widget in self.isotope_override_widgets:
            widget['frame'].destroy()
        self.isotope_override_widgets = []
    
    def get_isotope_overrides(self):
        """Get isotope cross-section overrides from GUI"""
        overrides = {}
        for widget in self.isotope_override_widgets:
            try:
                isotope = widget['isotope'].get().strip()
                unc_str = widget['uncertainty'].get().strip()
                if isotope and unc_str:
                    unc = float(unc_str) / 100.0  # Convert from % to fraction
                    overrides[isotope] = unc
            except ValueError:
                continue
        return overrides
    
    def update_formula_display(self):
        """Update the formula explanation text when selection changes"""
        method = self.adjustment_method.get()
        
        explanations = {
            'none': "Using standard formula (no adjustment)",
            'chi2': (
                "Using χ²-adjusted formula: statistical uncertainty scaled by √χ²/ν\n"
                "(Recommended when χ²/ν > 1)"
            ),
            'auc': (
                "Using AUC-adjusted formula: statistical uncertainty scaled by max(AUC, 1/AUC)\n"
                "(Alternative method based on area ratio)"
            )
        }
        self.formula_explanation.set(explanations.get(method, explanations['none']))
    
    def on_reaction_changed(self, event=None):
        """Update cross-section uncertainty when reaction selection changes"""
        reaction = self.selected_reaction.get()
        
        if reaction == 'Custom...':
            self.xsec_entry.config(state='normal')
            self.xsec_desc.set("(enter custom value)")
        elif reaction in CROSS_SECTION_UNCERTAINTIES:
            unc, desc = CROSS_SECTION_UNCERTAINTIES[reaction]
            self.xsec_uncertainty.set(unc * 100)
            self.xsec_entry.config(state='readonly')
            self.xsec_desc.set(f"({desc})")
        else:
            self.xsec_uncertainty.set(5.0)
            self.xsec_entry.config(state='readonly')
            self.xsec_desc.set("(unknown reaction)")
    
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
        ttk.Entry(row_frame, textvariable=elements_var, width=22).pack(side=tk.LEFT, padx=2)
        
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
        """Load default ROI definitions based on selected reaction"""
        self.clear_all_rois()
        
        reaction = self.selected_reaction.get()
        
        if 'RBS' in reaction:
            self.add_roi_row("Heavy", "700", "900", "Ba,Ca,Fe")
            self.add_roi_row("Medium", "380", "550", "Si,Al,Na,Mg")
            self.add_roi_row("Light", "100", "250", "O,C")
        elif 'ERDA' in reaction or '¹H' in reaction:
            self.add_roi_row("Hydrogen", "200", "400", "H")
        elif 'NRA' in reaction:
            if '¹⁶O' in reaction:
                self.add_roi_row("Oxygen", "400", "600", "O,16O")
            elif '¹²C' in reaction:
                self.add_roi_row("Carbon", "100", "250", "C,12C,13C")
            elif '¹⁵N' in reaction:
                self.add_roi_row("Nitrogen", "250", "400", "N,14N,15N")
            else:
                self.add_roi_row("ROI1", "100", "300", "C,O,N")
        else:
            self.add_roi_row("Heavy", "700", "900", "Ba,Ca,Fe")
            self.add_roi_row("Light", "100", "300", "O,C")
    
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
    
    
    def load_element_spectrum_for_viewer(self):
        """Load element spectrum for the Element/Iso Data viewer tab"""
        filename = filedialog.askopenfilename(
            title="Load Element/Iso Data",
            filetypes=[
                ("Text files", "*.txt"),
                ("Data files", "*.dat"),
                ("All files", "*.*")
            ]
        )
        
        if not filename:
            return
        
        try:
            self.element_spectrum_data = parse_element_spectrum(filename)
            self.element_spectrum_file.set(filename)
            
            # Update viewer tab display
            self.elem_viewer_file.set(Path(filename).name)
            
            # Get summary
            summary = get_element_spectrum_summary(self.element_spectrum_data)
            
            # Update channel info
            self.elem_channel_info.set(
                f"Channels: {summary['channel_min']} - {summary['channel_max']} "
                f"({summary['n_channels']} points)"
            )
            
            # Populate elements table
            self.elem_tree.delete(*self.elem_tree.get_children())
            for elem, counts in summary['element_totals'].items():
                self.elem_tree.insert('', 'end', values=(elem, f"{counts:,}"))
            
            # Populate isotopes table
            self.iso_tree.delete(*self.iso_tree.get_children())
            for iso, counts in summary['isotope_totals'].items():
                parent = summary['isotope_to_element'].get(iso, '—')
                self.iso_tree.insert('', 'end', values=(iso, parent, f"{counts:,}"))
            
            # Create checkboxes for element selection
            self.update_element_checkboxes(summary['elements'])
            self.update_isotope_checkboxes(summary['isotopes'])
            
            # Also update the Uncertainty tab status
            status_text = f"✓ Loaded: {Path(filename).name}\n"
            status_text += f"   {len(summary['elements'])} elements, {len(summary['isotopes'])} isotopes\n"
            status_text += f"   {summary['n_channels']} channels"
            self.elem_spec_status.set(status_text)
            self.elem_spec_status_label.config(foreground='#059669')
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load element spectrum:\n{str(e)}")
            self.element_spectrum_data = None
            self.elem_viewer_file.set("Error loading file")
            self.elem_channel_info.set("")
    
    def clear_element_spectrum_viewer(self):
        """Clear element spectrum from viewer tab"""
        self.element_spectrum_data = None
        self.element_spectrum_file.set("")
        self.elem_viewer_file.set("No file loaded")
        self.elem_channel_info.set("")
        
        # Clear tables
        self.elem_tree.delete(*self.elem_tree.get_children())
        self.iso_tree.delete(*self.iso_tree.get_children())
        
        # Clear checkboxes
        for widget in self.element_checkboxes_frame.winfo_children():
            widget.destroy()
        for widget in self.isotope_checkboxes_frame.winfo_children():
            widget.destroy()
        self.element_plot_vars = {}
        self.isotope_plot_vars = {}
        
        # Update uncertainty tab
        self.elem_spec_status.set("⚠ Not loaded (using approximate method)")
        self.elem_spec_status_label.config(foreground='#B45309')
    
    def update_element_checkboxes(self, elements):
        """Create checkboxes for element selection"""
        # Clear existing
        for widget in self.element_checkboxes_frame.winfo_children():
            widget.destroy()
        self.element_plot_vars = {}
        
        # Create new checkboxes in a grid
        for i, elem in enumerate(elements):
            var = tk.BooleanVar(value=True)
            self.element_plot_vars[elem] = var
            cb = ttk.Checkbutton(self.element_checkboxes_frame, text=elem, variable=var)
            cb.grid(row=i // 4, column=i % 4, sticky=tk.W, padx=2, pady=1)
    
    def update_isotope_checkboxes(self, isotopes):
        """Create checkboxes for isotope selection"""
        # Clear existing
        for widget in self.isotope_checkboxes_frame.winfo_children():
            widget.destroy()
        self.isotope_plot_vars = {}
        
        # Create new checkboxes in a grid
        for i, iso in enumerate(isotopes):
            var = tk.BooleanVar(value=False)  # Isotopes unchecked by default
            self.isotope_plot_vars[iso] = var
            cb = ttk.Checkbutton(self.isotope_checkboxes_frame, text=iso, variable=var)
            cb.grid(row=i // 4, column=i % 4, sticky=tk.W, padx=2, pady=1)
    
    def select_all_plot_items(self, item_type):
        """Select all elements or isotopes for plotting"""
        if item_type == 'elements':
            for var in self.element_plot_vars.values():
                var.set(True)
        elif item_type == 'isotopes':
            for var in self.isotope_plot_vars.values():
                var.set(True)
    
    def clear_all_plot_items(self, item_type):
        """Clear all element or isotope selections"""
        if item_type == 'elements':
            for var in self.element_plot_vars.values():
                var.set(False)
        elif item_type == 'isotopes':
            for var in self.isotope_plot_vars.values():
                var.set(False)
    
    def get_selected_elements(self):
        """Get list of selected elements for plotting"""
        return [elem for elem, var in self.element_plot_vars.items() if var.get()]
    
    def get_selected_isotopes(self):
        """Get list of selected isotopes for plotting"""
        return [iso for iso, var in self.isotope_plot_vars.items() if var.get()]
    
    def _get_elem_limits(self):
        """Return (ch_min, ch_max, en_min, en_max, cts_min, cts_max) with empty→None."""
        def _f(var):
            s = var.get().strip().replace(',', '.')
            if not s:
                return None
            try:
                return float(s)
            except ValueError:
                return None
        return (_f(self.elem_channel_min), _f(self.elem_channel_max),
                _f(self.elem_energy_min), _f(self.elem_energy_max),
                _f(self.elem_counts_min), _f(self.elem_counts_max))

    def plot_element_spectra_gui(self):
        """Plot selected element spectra"""
        if self.element_spectrum_data is None:
            messagebox.showwarning("Warning", "Please load an element spectrum file first!")
            return

        if self.elem_x_axis_mode.get() == 'energy' and not self.current_file.get():
            messagebox.showwarning("Warning",
                "Please load an .xnra file (top of window) to get energy calibration,\n"
                "or switch to Channel mode.")
            return

        elements = self.get_selected_elements()
        if not elements and not (self.show_sim_overlay.get() or self.show_raw_overlay.get()
                                  or self.show_smoothed_overlay.get()):
            messagebox.showwarning("Warning", "Please select at least one element or enable an overlay!")
            return

        ch_min, ch_max, en_min, en_max, cts_min, cts_max = self._get_elem_limits()
        try:
            fig = plot_element_spectra(
                self.element_spectrum_data,
                elements_to_plot=elements,
                isotopes_to_plot=None,
                show_simulated=self.show_sim_overlay.get(),
                show_raw=self.show_raw_overlay.get(),
                show_smoothed=self.show_smoothed_overlay.get(),
                x_axis_mode=self.elem_x_axis_mode.get(),
                cal_offset=self.elem_cal_offset.get(),
                cal_gain=self.elem_cal_gain.get(),
                title=self.elem_plot_title.get(),
                channel_min=ch_min, channel_max=ch_max,
                energy_min=en_min, energy_max=en_max,
                counts_min=cts_min, counts_max=cts_max,
                show_plot=True
            )
        except Exception as e:
            messagebox.showerror("Error", f"Failed to plot element spectra:\n{str(e)}")

    def plot_isotope_spectra_gui(self):
        """Plot selected isotope spectra"""
        if self.element_spectrum_data is None:
            messagebox.showwarning("Warning", "Please load an element spectrum file first!")
            return

        if self.elem_x_axis_mode.get() == 'energy' and not self.current_file.get():
            messagebox.showwarning("Warning",
                "Please load an .xnra file (top of window) to get energy calibration,\n"
                "or switch to Channel mode.")
            return

        isotopes = self.get_selected_isotopes()
        if not isotopes and not (self.show_sim_overlay.get() or self.show_raw_overlay.get()
                                  or self.show_smoothed_overlay.get()):
            messagebox.showwarning("Warning", "Please select at least one isotope or enable an overlay!")
            return

        ch_min, ch_max, en_min, en_max, cts_min, cts_max = self._get_elem_limits()
        try:
            fig = plot_element_spectra(
                self.element_spectrum_data,
                elements_to_plot=None,
                isotopes_to_plot=isotopes,
                show_simulated=self.show_sim_overlay.get(),
                show_raw=self.show_raw_overlay.get(),
                show_smoothed=self.show_smoothed_overlay.get(),
                x_axis_mode=self.elem_x_axis_mode.get(),
                cal_offset=self.elem_cal_offset.get(),
                cal_gain=self.elem_cal_gain.get(),
                title=self.elem_plot_title.get() or "Isotope Spectra",
                channel_min=ch_min, channel_max=ch_max,
                energy_min=en_min, energy_max=en_max,
                counts_min=cts_min, counts_max=cts_max,
                show_plot=True
            )
        except Exception as e:
            messagebox.showerror("Error", f"Failed to plot isotope spectra:\n{str(e)}")

    def plot_all_selected_spectra(self):
        """Plot all selected elements and isotopes together"""
        if self.element_spectrum_data is None:
            messagebox.showwarning("Warning", "Please load an element spectrum file first!")
            return

        if self.elem_x_axis_mode.get() == 'energy' and not self.current_file.get():
            messagebox.showwarning("Warning",
                "Please load an .xnra file (top of window) to get energy calibration,\n"
                "or switch to Channel mode.")
            return

        elements = self.get_selected_elements()
        isotopes = self.get_selected_isotopes()

        if not elements and not isotopes and not (
                self.show_sim_overlay.get() or self.show_raw_overlay.get()
                or self.show_smoothed_overlay.get()):
            messagebox.showwarning("Warning", "Please select at least one element/isotope or enable an overlay!")
            return

        ch_min, ch_max, en_min, en_max, cts_min, cts_max = self._get_elem_limits()
        try:
            fig = plot_element_spectra(
                self.element_spectrum_data,
                elements_to_plot=elements if elements else None,
                isotopes_to_plot=isotopes if isotopes else None,
                show_simulated=self.show_sim_overlay.get(),
                show_raw=self.show_raw_overlay.get(),
                show_smoothed=self.show_smoothed_overlay.get(),
                x_axis_mode=self.elem_x_axis_mode.get(),
                cal_offset=self.elem_cal_offset.get(),
                cal_gain=self.elem_cal_gain.get(),
                title=self.elem_plot_title.get(),
                channel_min=ch_min, channel_max=ch_max,
                energy_min=en_min, energy_max=en_max,
                counts_min=cts_min, counts_max=cts_max,
                show_plot=True
            )
        except Exception as e:
            messagebox.showerror("Error", f"Failed to plot spectra:\n{str(e)}")
    
    def calculate_uncertainties_v9(self):
        """Calculate uncertainties using v9 formulas with isotope support"""
        if not UNCERTAINTY_AVAILABLE:
            messagebox.showerror("Error", 
                "Uncertainty calculation module not available.\n\n"
                "Make sure xnra_uncertainty_v9.py is in the same directory.")
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
            import numpy as np
            
            data = parse_xnra_full(self.current_file.get())
            
            qomega_unc = self.qomega_uncertainty.get() / 100.0
            xsec_unc = self.xsec_uncertainty.get() / 100.0
            reaction = self.selected_reaction.get()
            adj_method = self.adjustment_method.get()
            isotope_mode = self.isotope_mode.get()
            isotope_overrides = self.get_isotope_overrides()
            
            for item in self.results_tree.get_children():
                self.results_tree.delete(item)
            
            results = []
            
            for roi in rois:
                roi_name = roi['name']
                roi_start, roi_end = roi['channels']
                roi_elements = roi['elements']
                
                if data['raw_counts'] is not None and data['simulated_counts'] is not None:
                    chi2_r, dof, roi_counts, auc_ratio = calculate_chi_squared_for_roi(
                        data['raw_counts'],
                        data['simulated_counts'],
                        roi_start,
                        roi_end
                    )
                    chi2_symbol, chi2_status_text = interpret_chi2(chi2_r)
                else:
                    chi2_r, dof, roi_counts, auc_ratio = float('nan'), 0, 0, float('nan')
                    chi2_symbol, chi2_status_text = "?", "N/A"
                
                selected_layers = self._get_selected_layer_nums()
                for layer in data['layers']:
                    if selected_layers and layer.get('layer_num') not in selected_layers:
                        continue
                    for element_or_isotope in roi_elements:
                        # Determine if this is an isotope
                        is_isotope = element_or_isotope[0].isdigit() if element_or_isotope else False
                        
                        if isotope_mode and is_isotope:
                            concentration = layer.get('isotope_concentrations', {}).get(element_or_isotope, 0)
                            # Get isotope-specific cross-section
                            iso_xsec_unc, iso_reaction = get_isotope_cross_section_uncertainty(
                                element_or_isotope, isotope_overrides
                            )
                        else:
                            base_element = re.sub(r'^\d+', '', element_or_isotope)
                            concentration = layer['concentrations'].get(base_element, 0)
                            iso_xsec_unc = xsec_unc
                            iso_reaction = reaction
                        
                        if concentration <= 0:
                            continue
                        
                        # Calculate element counts
                        element_fraction = None
                        if self.element_spectrum_data is not None and data['raw_counts'] is not None:
                            element_counts, element_fraction, total_exp = calculate_element_counts_from_spectrum(
                                self.element_spectrum_data,
                                data['raw_counts'],
                                element_or_isotope,
                                roi_start,
                                roi_end,
                                use_isotope=(isotope_mode and is_isotope)
                            )
                            aa_method = 'spectrum'
                        else:
                            element_counts = int(roi_counts * concentration)
                            element_fraction = concentration
                            aa_method = 'approximate'
                        
                        unc = calculate_uncertainty_v9(
                            concentration, 
                            element_counts, 
                            qomega_unc, 
                            iso_xsec_unc,
                            chi2_reduced=chi2_r,
                            auc_ratio=auc_ratio,
                            adjustment_method=adj_method
                        )
                        
                        # Format displays
                        if np.isnan(chi2_r):
                            chi2_display = "N/A"
                        else:
                            chi2_display = f"{chi2_r:.2f}"
                        
                        status_display = f"{chi2_symbol} {chi2_status_text}"
                        
                        if np.isnan(auc_ratio):
                            auc_display = "N/A"
                        else:
                            auc_display = f"{auc_ratio:.3f}"
                        
                        def fmt_unc_abs(val):
                            """Format ABSOLUTE uncertainty (concentration-scaled), scientific 2dp"""
                            if val is None or (isinstance(val, float) and np.isnan(val)):
                                return "N/A"
                            return f"{val*100:.2e}"

                        def fmt_unc_rel(val):
                            """Format RELATIVE uncertainty (constant), scientific 2dp"""
                            if val is None or (isinstance(val, float) and np.isnan(val)):
                                return "N/A"
                            return f"{val*100:.2e}"
                        
                        aa_display = f"{element_counts:,}" if element_counts > 0 else "N/A"
                        
                        if element_fraction is not None:
                            fa_display = f"{element_fraction:.3f}"
                        else:
                            fa_display = "N/A"
                        
                        # KEY FIX: Display RELATIVE (constant) values for Q·Ω and σ
                        values = (
                            layer['layer_num'],
                            element_or_isotope,
                            f"{concentration*100:.2f}",
                            fmt_unc_rel(unc['sigma_stat_rel']),   # Stat: RELATIVE
                            fmt_unc_rel(unc['sigma_qomega_rel']), # Q·Ω: RELATIVE (CONSTANT)
                            fmt_unc_rel(unc['sigma_xsec_rel']),   # σ: RELATIVE (CONSTANT)
                            fmt_unc_abs(unc['sigma_total']),      # Total: ABSOLUTE (varies)
                            roi_name,
                            aa_display,
                            fa_display,
                            chi2_display,
                            status_display,
                            auc_display
                        )
                        self.results_tree.insert('', 'end', values=values)
                        
                        results.append({
                            'filename': data['filename'],
                            'reaction': iso_reaction,
                            'is_isotope': is_isotope,
                            'layer': layer['layer_num'],
                            'thickness': layer['thickness'],
                            'element': element_or_isotope,
                            'concentration_percent': concentration * 100,
                            # Relative uncertainties (CONSTANT for Q·Ω and σ)
                            'sigma_stat_rel_percent': unc['sigma_stat_rel'] * 100 if not np.isnan(unc['sigma_stat_rel']) else None,
                            'sigma_qomega_rel_percent': unc['sigma_qomega_rel'] * 100,  # CONSTANT
                            'sigma_xsec_rel_percent': unc['sigma_xsec_rel'] * 100,      # CONSTANT per reaction
                            'sigma_total_rel_percent': unc['sigma_total_rel'] * 100,
                            # Absolute uncertainties (concentration-scaled)
                            'sigma_stat_percent': unc['sigma_stat'] * 100 if not np.isnan(unc['sigma_stat']) else None,
                            'sigma_qomega_percent': unc['sigma_qomega'] * 100,
                            'sigma_xsec_percent': unc['sigma_xsec'] * 100,
                            'sigma_total_percent': unc['sigma_total'] * 100,
                            # All versions
                            'sigma_stat_standard_percent': unc['sigma_stat_standard'] * 100 if not np.isnan(unc['sigma_stat_standard']) else None,
                            'sigma_stat_chi2_adjusted_percent': unc['sigma_stat_chi2_adjusted'] * 100 if not np.isnan(unc['sigma_stat_chi2_adjusted']) else None,
                            'sigma_stat_auc_adjusted_percent': unc['sigma_stat_auc_adjusted'] * 100 if not np.isnan(unc['sigma_stat_auc_adjusted']) else None,
                            'sigma_total_standard_percent': unc['sigma_total_standard'] * 100,
                            'sigma_total_chi2_adjusted_percent': unc['sigma_total_chi2_adjusted'] * 100,
                            'sigma_total_auc_adjusted_percent': unc['sigma_total_auc_adjusted'] * 100,
                            'roi_name': roi_name,
                            'roi_start': roi_start,
                            'roi_end': roi_end,
                            'roi_counts': roi_counts,
                            'element_counts': element_counts,
                            'element_fraction': element_fraction,
                            'aa_method': aa_method,
                            'chi2_reduced': chi2_r,
                            'chi2_status': chi2_status_text,
                            'chi2_factor': unc['chi2_factor'],
                            'auc_ratio': auc_ratio,
                            'auc_factor': unc['auc_factor'],
                            'adjustment_method': adj_method,
                            'adjustment_factor': unc['adjustment_factor'],
                            'qomega_input_percent': qomega_unc * 100,
                            'xsec_input_percent': iso_xsec_unc * 100,
                        })
            
            self.uncertainty_results = results
            
            import pandas as pd
            self.uncertainty_results_df = pd.DataFrame(results)
            
            aa_methods_used = set(r['aa_method'] for r in results) if results else set()
            if 'spectrum' in aa_methods_used:
                aa_method_str = "A_A: spectrum ✓"
            else:
                aa_method_str = "A_A: approximate ⚠"
            
            iso_str = " | Isotope mode ✓" if isotope_mode else ""
            
            method_names = {'none': 'Standard', 'chi2': 'χ²-adjusted', 'auc': 'AUC-adjusted'}
            self.uncertainty_status.config(
                text=f"✓ {len(results)} results | {method_names.get(adj_method, adj_method)} | "
                     f"Q·Ω: {qomega_unc*100:.1f}% (const) | {aa_method_str}{iso_str}"
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
            title="Export Uncertainty Results",
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
                df.to_excel(filename, index=False, sheet_name='Uncertainties_v9')
            else:
                df.to_csv(filename, index=False)
            
            messagebox.showinfo("Success", f"Results exported to:\n{filename}")
        except ImportError:
            messagebox.showerror("Error", "pandas is required for export.\nInstall with: pip install pandas openpyxl")
        except Exception as e:
            messagebox.showerror("Error", f"Export failed:\n{str(e)}")
    
    def export_depth_profile(self):
        """Export depth profile table"""
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
            
            df = self.uncertainty_results_df
            adj_method = self.adjustment_method.get()
            profile = create_depth_profile(df, adjustment_method=adj_method)
            
            method_names = {'none': 'Standard (none)', 'chi2': 'χ²-adjusted', 'auc': 'AUC-adjusted'}
            method_formulas = {
                'none': 'Eq. 3 (Standard)',
                'chi2': 'Eq. 6 (χ²-adjusted)',
                'auc': 'AUC-adjusted (×max(AUC,1/AUC))'
            }
            
            if filename.endswith('.xlsx'):
                with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                    profile.to_excel(writer, sheet_name='Depth Profile', index=False)
                    df.to_excel(writer, sheet_name='Raw Data', index=False)
                    
                    summary = pd.DataFrame([{
                        'Parameter': 'Q·Ω Uncertainty (%)',
                        'Value': self.qomega_uncertainty.get(),
                        'Note': 'CONSTANT for all elements'
                    }, {
                        'Parameter': 'Cross-section Uncertainty (%)',
                        'Value': self.xsec_uncertainty.get(),
                        'Note': 'Per-reaction (may vary by isotope)'
                    }, {
                        'Parameter': 'Reaction',
                        'Value': self.selected_reaction.get()
                    }, {
                        'Parameter': 'Adjustment Method',
                        'Value': method_names.get(adj_method, adj_method)
                    }, {
                        'Parameter': 'Formula',
                        'Value': method_formulas.get(adj_method, 'Standard')
                    }, {
                        'Parameter': 'Isotope Mode',
                        'Value': 'Enabled' if self.isotope_mode.get() else 'Disabled'
                    }])
                    summary.to_excel(writer, sheet_name='Parameters', index=False)
            else:
                profile.to_csv(filename, index=False)
            
            messagebox.showinfo("Success", f"Depth profile exported to:\n{filename}")
        except ImportError:
            messagebox.showerror("Error", "pandas is required for export.\nInstall with: pip install pandas openpyxl")
        except Exception as e:
            import traceback
            messagebox.showerror("Error", f"Export failed:\n{str(e)}\n\n{traceback.format_exc()}")
    
    
    def create_action_buttons(self, parent):
        """Action buttons section"""
        button_frame = ttk.Frame(parent)
        button_frame.grid(row=2, column=0, sticky=(tk.W, tk.E))

        ttk.Button(button_frame, text="Export Plot", command=self.export_plot).pack(side=tk.LEFT, padx=(0, 5))
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
        var.trace_add('write', update_preview)
        
        return var
    
    def browse_file(self):
        """Open file browser dialog"""
        filename = filedialog.askopenfilename(
            title="Select XNRA File",
            filetypes=[("XNRA files", "*.xnra"), ("All files", "*.*")]
        )
        if filename:
            self.current_file.set(filename)

            try:
                data = parse_xnra_full(filename)
                self.elem_cal_offset.set(data.get('cal_offset', 0.0))
                self.elem_cal_gain.set(data.get('cal_gain', 1.0))
                self._populate_layer_checkboxes(data.get('layers', []))
            except Exception:
                self.elem_cal_offset.set(0.0)
                self.elem_cal_gain.set(1.0)
                self._populate_layer_checkboxes([])

    def _populate_layer_checkboxes(self, layers):
        """Rebuild the layer-selection checkbox list. All layers checked by default."""
        for widget in self.layers_checkbox_frame.winfo_children():
            widget.destroy()
        self.layer_vars = {}

        if not layers:
            ttk.Label(self.layers_checkbox_frame,
                      text="(no layers found in file)",
                      font=('TkDefaultFont', 8, 'italic'), foreground='#999'
                      ).pack(anchor=tk.W)
            return

        for layer in layers:
            n = layer.get('layer_num')
            if n is None:
                continue
            var = tk.BooleanVar(value=True)
            self.layer_vars[n] = var
            thick = layer.get('thickness')
            label = f"Layer {n}"
            if thick is not None:
                label += f"  ({thick:g} ×10¹⁵ at/cm²)"
            ttk.Checkbutton(self.layers_checkbox_frame, text=label, variable=var
                            ).pack(anchor=tk.W)

    def _set_all_layers(self, value):
        for var in self.layer_vars.values():
            var.set(bool(value))

    def _get_selected_layer_nums(self):
        """Return set of currently-selected layer numbers; empty means 'no filter'."""
        return {n for n, var in self.layer_vars.items() if var.get()}
    
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
        """Initialize default values"""
        self.on_reaction_changed()
        self.update_formula_display()


def main():
    root = tk.Tk()
    app = XNRAViewerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
