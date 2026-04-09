### XNRA Spectrum Viewer GUI - Enhanced Version ###
# Tkinter-based graphical interface for XNRA spectrum visualization
# VERSION 2: Added input validation and Enter key functionality

# Script developed by Berke Santos & Giuseppe Legrottaglie
# Last updated: 16/12/2025

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser
import json
from pathlib import Path
from xnra_main_v3 import plot_xnra_spectrum

class XNRAViewerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("XNRA Spectrum Viewer")
        self.root.geometry("1200x700")  # Wider for horizontal layout
        
        # Store current filepath
        self.current_file = tk.StringVar()
        self.last_fig = None  # Store last generated figure for export
        
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
        """
        Validate number input when Enter is pressed
        Returns True if valid, False otherwise
        """
        # Get value and convert to string (handles DoubleVar, IntVar, StringVar)
        value = str(var.get()).strip()
        
        # Handle empty input for optional fields
        if allow_empty and value == "":
            return True
        
        # Replace comma with period (common European decimal separator)
        value = value.replace(',', '.')
        
        try:
            if is_int:
                # Convert to int
                num_value = int(float(value))  # float first to handle "5.0" -> 5
                var.set(num_value)
            else:
                # Convert to float
                num_value = float(value)
                var.set(num_value)
            
            return True
            
        except ValueError:
            # Invalid input - show error
            messagebox.showerror(
                "Invalid Input",
                f"Invalid number format: '{value}'\n\n"
                "Please enter a valid number.\n"
                "Note: Use '.' as decimal separator (not ',')"
            )
            return False
    
    def validate_number_input_from_entry(self, entry, var, is_int=False, allow_empty=False):
        """
        Validate number input directly from Entry widget text
        Returns True if valid, False otherwise
        """
        # Get value directly from Entry widget (not from variable)
        value = entry.get().strip()
        
        # Handle empty input for optional fields
        if allow_empty and value == "":
            return True
        
        # Replace comma with period (common European decimal separator)
        value = value.replace(',', '.')
        
        try:
            if is_int:
                # Convert to int
                num_value = int(float(value))  # float first to handle "5.0" -> 5
                var.set(num_value)
            else:
                # Convert to float
                num_value = float(value)
                var.set(num_value)
            
            return True
            
        except ValueError:
            # Invalid input - show error and restore previous valid value
            messagebox.showerror(
                "Invalid Input",
                f"Invalid number format: '{value}'\n\n"
                "Please enter a valid number.\n"
                "Note: Use '.' as decimal separator (not ',')"
            )
            # Force Entry to show the current valid value from variable
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
        self.create_appearance_tab()
        self.create_axes_tab()
    
    def create_general_tab(self):
        """General settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="General")
        
        # Create scrollable frame
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Enable mouse wheel scrolling
        self.bind_mousewheel(canvas)
        
        # Configure columns for horizontal layout
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        
        # Figure settings (LEFT)
        fig_frame = ttk.LabelFrame(scrollable_frame, text="Figure Settings", padding="10")
        fig_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5), pady=(0, 10))
        
        self.fig_width = self.create_number_field(fig_frame, "Width (inches):", 0, 12.0)
        self.fig_height = self.create_number_field(fig_frame, "Height (inches):", 1, 6.0)
        self.fig_dpi = self.create_number_field(fig_frame, "DPI:", 2, 100, is_int=True)
        
        # Display options (RIGHT)
        display_frame = ttk.LabelFrame(scrollable_frame, text="Display Options", padding="10")
        display_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
        self.show_legend = tk.BooleanVar(value=True)
        ttk.Checkbutton(display_frame, text="Show Legend", variable=self.show_legend).grid(row=0, column=0, sticky=tk.W)
        
        self.show_info_box = tk.BooleanVar(value=True)
        ttk.Checkbutton(display_frame, text="Show Info Box", variable=self.show_info_box).grid(row=1, column=0, sticky=tk.W)
        
        ttk.Label(display_frame, text="Info Box Position:").grid(row=2, column=0, sticky=tk.W, pady=(10, 0))
        self.info_box_position = tk.StringVar(value='left')
        ttk.Radiobutton(display_frame, text="Left", variable=self.info_box_position, value='left').grid(row=3, column=0, sticky=tk.W, padx=(20, 0))
        ttk.Radiobutton(display_frame, text="Right", variable=self.info_box_position, value='right').grid(row=4, column=0, sticky=tk.W, padx=(20, 0))
        
        # Title Settings (full width, below Figure and Display)
        title_frame = ttk.LabelFrame(scrollable_frame, text="Plot Title", padding="10")
        title_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        title_frame.columnconfigure(1, weight=1)
        
        ttk.Label(title_frame, text="Custom Title:").grid(row=0, column=0, sticky=tk.W, padx=(0, 5))
        self.custom_title = tk.StringVar(value="")
        title_entry = ttk.Entry(title_frame, textvariable=self.custom_title)
        title_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5)
        ttk.Label(title_frame, text="(empty = filename without extension)", font=('TkDefaultFont', 8, 'italic')).grid(row=0, column=2, sticky=tk.W, padx=(5, 0))
    
    def create_raw_data_tab(self):
        """Raw data settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Raw Data")
        
        # Create scrollable frame
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Enable mouse wheel scrolling
        self.bind_mousewheel(canvas)
        
        # Configure columns for horizontal layout
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        scrollable_frame.columnconfigure(2, weight=1)
        
        # Master toggle (full width)
        self.show_raw_data = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="Show Raw Data", variable=self.show_raw_data).grid(row=0, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))
        
        # Markers (LEFT)
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
        
        # Marker Range (MIDDLE)
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.raw_marker_range_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.raw_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.raw_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.raw_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.raw_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 0.0)
        self.raw_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 200.0)
        
        # Line (RIGHT)
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
        
        # Create scrollable frame
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Enable mouse wheel scrolling
        self.bind_mousewheel(canvas)
        
        # Configure columns for horizontal layout
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        scrollable_frame.columnconfigure(2, weight=1)
        
        # Master toggle (full width)
        self.show_smoothed_data = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="Show Smoothed Data", variable=self.show_smoothed_data).grid(row=0, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))
        
        # Markers (LEFT)
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
        
        # Marker Range (MIDDLE)
        range_frame = ttk.LabelFrame(scrollable_frame, text="Marker Range", padding="10")
        range_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=5)
        
        self.smoothed_marker_range_enabled = tk.BooleanVar(value=True)
        ttk.Checkbutton(range_frame, text="Enable Range Limit", variable=self.smoothed_marker_range_enabled).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        ttk.Label(range_frame, text="Mode:").grid(row=1, column=0, sticky=tk.W, pady=(5, 0))
        self.smoothed_marker_range_mode = tk.StringVar(value='channel')
        ttk.Combobox(range_frame, textvariable=self.smoothed_marker_range_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=1, column=1, sticky=tk.W, padx=(5, 0), pady=(5, 0))
        
        self.smoothed_marker_range_min = self.create_number_field(range_frame, "Min:", 2, 200.0)
        self.smoothed_marker_range_max = self.create_number_field(range_frame, "Max:", 3, 400.0)
        
        # Line (RIGHT)
        line_frame = ttk.LabelFrame(scrollable_frame, text="Line", padding="10")
        line_frame.grid(row=1, column=2, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.show_smoothed_line = tk.BooleanVar(value=True)
        ttk.Checkbutton(line_frame, text="Show Line", variable=self.show_smoothed_line).grid(row=0, column=0, columnspan=2, sticky=tk.W)
        
        self.smoothed_line_color = self.create_color_field(line_frame, "Color:", 1, 'red')
        self.smoothed_line_width = self.create_number_field(line_frame, "Width:", 2, 1.0)
        self.smoothed_line_alpha = self.create_number_field(line_frame, "Transparency:", 3, 1.0)
    
    def create_appearance_tab(self):
        """Plot appearance settings tab"""
        tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(tab, text="Appearance")
        
        # Create scrollable frame
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Enable mouse wheel scrolling
        self.bind_mousewheel(canvas)
        
        # Configure columns for horizontal layout
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        
        # ROW 1: Colors (LEFT) + Text Sizes (RIGHT)
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
        
        # ROW 2: Text Colors (LEFT) + Tick Marks (RIGHT)
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
        
        # Create scrollable frame
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Enable mouse wheel scrolling
        self.bind_mousewheel(canvas)
        
        # Configure columns for horizontal layout
        scrollable_frame.columnconfigure(0, weight=1)
        scrollable_frame.columnconfigure(1, weight=1)
        
        # ROW 1: X-axis mode (LEFT) + Channel Limits (RIGHT)
        mode_frame = ttk.LabelFrame(scrollable_frame, text="X-Axis Mode", padding="10")
        mode_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5), pady=(0, 10))
        
        ttk.Label(mode_frame, text="Display Mode:").grid(row=0, column=0, sticky=tk.W)
        self.x_axis_mode = tk.StringVar(value='energy')
        ttk.Combobox(mode_frame, textvariable=self.x_axis_mode, values=['channel', 'energy'], state='readonly', width=12).grid(row=0, column=1, sticky=tk.W, padx=(5, 0))
        
        channel_frame = ttk.LabelFrame(scrollable_frame, text="Channel Limits", padding="10")
        channel_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0), pady=(0, 10))
        
        self.channel_min = self.create_optional_number_field(channel_frame, "Min:", 0)
        self.channel_max = self.create_optional_number_field(channel_frame, "Max:", 1)
        
        # ROW 2: Energy Limits (LEFT) + Counts Limits (RIGHT)
        energy_frame = ttk.LabelFrame(scrollable_frame, text="Energy Limits (keV)", padding="10")
        energy_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), padx=(0, 5))
        
        self.energy_min = self.create_optional_number_field(energy_frame, "Min:", 0)
        self.energy_max = self.create_optional_number_field(energy_frame, "Max:", 1)
        
        counts_frame = ttk.LabelFrame(scrollable_frame, text="Counts Limits (Y-axis)", padding="10")
        counts_frame.grid(row=1, column=1, sticky=(tk.W, tk.E, tk.N), padx=(5, 0))
        
        self.counts_min = self.create_optional_number_field(counts_frame, "Min:", 0)
        self.counts_max = self.create_optional_number_field(counts_frame, "Max:", 1)
    
    def create_action_buttons(self, parent):
        """Action buttons section"""
        button_frame = ttk.Frame(parent)
        button_frame.grid(row=2, column=0, sticky=(tk.W, tk.E))
        
        ttk.Button(button_frame, text="Load Preset", command=self.load_preset).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(button_frame, text="Save Preset", command=self.save_preset).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Export Plot", command=self.export_plot).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Generate Plot", command=self.generate_plot, style='Accent.TButton').pack(side=tk.RIGHT, padx=(5, 0))
    
    # Helper methods for creating fields
    def create_number_field(self, parent, label, row, default, is_int=False):
        """Create a labeled number entry field with validation"""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky=tk.W)
        
        var = tk.DoubleVar(value=default) if not is_int else tk.IntVar(value=int(default))
        entry = ttk.Entry(parent, textvariable=var, width=12)
        entry.grid(row=row, column=1, sticky=tk.W, padx=(5, 0))
        
        # Bind Enter key to validation
        def on_enter(event):
            self.validate_number_input_from_entry(entry, var, is_int=is_int, allow_empty=False)
        
        entry.bind('<Return>', on_enter)
        
        return var
    
    def create_optional_number_field(self, parent, label, row):
        """Create a labeled optional number entry field (can be empty for None)"""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky=tk.W)
        
        var = tk.StringVar(value="")
        entry = ttk.Entry(parent, textvariable=var, width=12)
        entry.grid(row=row, column=1, sticky=tk.W, padx=(5, 0))
        
        ttk.Label(parent, text="(empty = auto)", font=('TkDefaultFont', 8, 'italic')).grid(row=row, column=2, sticky=tk.W, padx=(5, 0))
        
        # Bind Enter key to validation
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
        
        # Update preview when entry changes
        def update_preview(*args):
            try:
                color_preview.configure(bg=var.get())
            except:
                pass
        var.trace('w', update_preview)
        
        return var
    
    # Action methods
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
        
        # Replace comma with period
        val = val.replace(',', '.')
        
        try:
            return float(val)
        except ValueError:
            return None
    
    def generate_plot(self):
        """Generate plot with current settings"""
        if not self.current_file.get():
            messagebox.showerror("Error", "Please select an XNRA file first!")
            return
        
        if not Path(self.current_file.get()).exists():
            messagebox.showerror("Error", "Selected file does not exist!")
            return
        
        try:
            # Determine plot title
            if self.custom_title.get().strip():
                plot_title = self.custom_title.get().strip()
            else:
                # Use filename without extension
                plot_title = Path(self.current_file.get()).stem
            
            # Gather all parameters
            params = {
                'filepath': self.current_file.get(),
                'title': plot_title,  # Add custom title
                # Figure settings
                'figure_size': (self.fig_width.get(), self.fig_height.get()),
                'dpi': self.fig_dpi.get(),
                'show_legend': self.show_legend.get(),
                'show_info_box': self.show_info_box.get(),
                'info_box_position': self.info_box_position.get(),
                'info_box_x_offset': self.info_box_x_offset.get(),  
                'info_box_y_offset': self.info_box_y_offset.get(), 
                # Raw data
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
                # Smoothed data
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
                # Appearance
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
                # Axes
                'x_axis_mode': self.x_axis_mode.get(),
                'channel_min': self.get_optional_float(self.channel_min),
                'channel_max': self.get_optional_float(self.channel_max),
                'energy_min': self.get_optional_float(self.energy_min),
                'energy_max': self.get_optional_float(self.energy_max),
                'counts_min': self.get_optional_float(self.counts_min),
                'counts_max': self.get_optional_float(self.counts_max),
                # Display 
                'show_plot': True
            }
            
            # Generate plot
            self.last_fig = plot_xnra_spectrum(**params)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate plot:\n{str(e)}")
    
    def save_preset(self):
        """Save current settings to JSON file"""
        filename = filedialog.asksaveasfilename(
            title="Save Preset",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if not filename:
            return
        
        try:
            preset = self.get_current_settings()
            with open(filename, 'w') as f:
                json.dump(preset, f, indent=4)
            messagebox.showinfo("Success", f"Preset saved to:\n{filename}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save preset:\n{str(e)}")
    
    def load_preset(self):
        """Load settings from JSON file"""
        filename = filedialog.askopenfilename(
            title="Load Preset",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if not filename:
            return
        
        try:
            with open(filename, 'r') as f:
                preset = json.load(f)
            self.apply_settings(preset)
            messagebox.showinfo("Success", f"Preset loaded from:\n{filename}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load preset:\n{str(e)}")
    
    def export_plot(self):
        """Export last generated plot to file"""
        if self.last_fig is None:
            messagebox.showerror("Error", "Please generate a plot first!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export Plot",
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("JPEG files", "*.jpg"),
                ("TIFF files", "*.tiff"),
                ("PDF files", "*.pdf"),
                ("All files", "*.*")
            ]
        )
        
        if not filename:
            return
        
        try:
            self.last_fig.savefig(filename, dpi=self.fig_dpi.get(), bbox_inches='tight')
            messagebox.showinfo("Success", f"Plot exported to:\n{filename}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to export plot:\n{str(e)}")
    
    def get_current_settings(self):
        """Get all current settings as dictionary"""
        return {
            # Figure
            'figure_width': self.fig_width.get(),
            'figure_height': self.fig_height.get(),
            'dpi': self.fig_dpi.get(),
            'custom_title': self.custom_title.get(),  # Add custom title
            'show_legend': self.show_legend.get(),
            'show_info_box': self.show_info_box.get(),
            'info_box_position': self.info_box_position.get(),
            # Raw data
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
            # Smoothed data
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
            # Appearance
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
            # Axes
            'x_axis_mode': self.x_axis_mode.get(),
            'channel_min': self.channel_min.get(),
            'channel_max': self.channel_max.get(),
            'energy_min': self.energy_min.get(),
            'energy_max': self.energy_max.get(),
            'counts_min': self.counts_min.get(),
            'counts_max': self.counts_max.get()
        }
    
    def apply_settings(self, settings):
        """Apply settings from dictionary"""
        # Figure
        if 'figure_width' in settings: self.fig_width.set(settings['figure_width'])
        if 'figure_height' in settings: self.fig_height.set(settings['figure_height'])
        if 'dpi' in settings: self.fig_dpi.set(settings['dpi'])
        if 'custom_title' in settings: self.custom_title.set(settings['custom_title'])  # Add custom title
        if 'show_legend' in settings: self.show_legend.set(settings['show_legend'])
        if 'show_info_box' in settings: self.show_info_box.set(settings['show_info_box'])
        if 'info_box_position' in settings: self.info_box_position.set(settings['info_box_position'])
        # Raw data
        if 'show_raw_data' in settings: self.show_raw_data.set(settings['show_raw_data'])
        if 'show_raw_markers' in settings: self.show_raw_markers.set(settings['show_raw_markers'])
        if 'raw_marker_shape' in settings: self.raw_marker_shape.set(settings['raw_marker_shape'])
        if 'raw_marker_fill' in settings: self.raw_marker_fill.set(settings['raw_marker_fill'])
        if 'raw_marker_color' in settings: self.raw_marker_color.set(settings['raw_marker_color'])
        if 'raw_marker_size' in settings: self.raw_marker_size.set(settings['raw_marker_size'])
        if 'raw_marker_edge_width' in settings: self.raw_marker_edge_width.set(settings['raw_marker_edge_width'])
        if 'raw_marker_range_enabled' in settings: self.raw_marker_range_enabled.set(settings['raw_marker_range_enabled'])
        if 'raw_marker_range_mode' in settings: self.raw_marker_range_mode.set(settings['raw_marker_range_mode'])
        if 'raw_marker_range_min' in settings: self.raw_marker_range_min.set(settings['raw_marker_range_min'])
        if 'raw_marker_range_max' in settings: self.raw_marker_range_max.set(settings['raw_marker_range_max'])
        if 'show_raw_line' in settings: self.show_raw_line.set(settings['show_raw_line'])
        if 'raw_line_color' in settings: self.raw_line_color.set(settings['raw_line_color'])
        if 'raw_line_width' in settings: self.raw_line_width.set(settings['raw_line_width'])
        if 'raw_line_alpha' in settings: self.raw_line_alpha.set(settings['raw_line_alpha'])
        # Smoothed data
        if 'show_smoothed_data' in settings: self.show_smoothed_data.set(settings['show_smoothed_data'])
        if 'show_smoothed_markers' in settings: self.show_smoothed_markers.set(settings['show_smoothed_markers'])
        if 'smoothed_marker_shape' in settings: self.smoothed_marker_shape.set(settings['smoothed_marker_shape'])
        if 'smoothed_marker_fill' in settings: self.smoothed_marker_fill.set(settings['smoothed_marker_fill'])
        if 'smoothed_marker_color' in settings: self.smoothed_marker_color.set(settings['smoothed_marker_color'])
        if 'smoothed_marker_size' in settings: self.smoothed_marker_size.set(settings['smoothed_marker_size'])
        if 'smoothed_marker_edge_width' in settings: self.smoothed_marker_edge_width.set(settings['smoothed_marker_edge_width'])
        if 'smoothed_marker_range_enabled' in settings: self.smoothed_marker_range_enabled.set(settings['smoothed_marker_range_enabled'])
        if 'smoothed_marker_range_mode' in settings: self.smoothed_marker_range_mode.set(settings['smoothed_marker_range_mode'])
        if 'smoothed_marker_range_min' in settings: self.smoothed_marker_range_min.set(settings['smoothed_marker_range_min'])
        if 'smoothed_marker_range_max' in settings: self.smoothed_marker_range_max.set(settings['smoothed_marker_range_max'])
        if 'show_smoothed_line' in settings: self.show_smoothed_line.set(settings['show_smoothed_line'])
        if 'smoothed_line_color' in settings: self.smoothed_line_color.set(settings['smoothed_line_color'])
        if 'smoothed_line_width' in settings: self.smoothed_line_width.set(settings['smoothed_line_width'])
        if 'smoothed_line_alpha' in settings: self.smoothed_line_alpha.set(settings['smoothed_line_alpha'])
        # Appearance
        if 'background_color' in settings: self.background_color.set(settings['background_color'])
        if 'show_grid' in settings: self.show_grid.set(settings['show_grid'])
        if 'grid_color' in settings: self.grid_color.set(settings['grid_color'])
        if 'title_size' in settings: self.title_size.set(settings['title_size'])
        if 'axis_label_size' in settings: self.axis_label_size.set(settings['axis_label_size'])
        if 'tick_label_size' in settings: self.tick_label_size.set(settings['tick_label_size'])
        if 'legend_size' in settings: self.legend_size.set(settings['legend_size'])
        if 'info_box_size' in settings: self.info_box_size.set(settings['info_box_size'])
        if 'title_color' in settings: self.title_color.set(settings['title_color'])
        if 'axis_label_color' in settings: self.axis_label_color.set(settings['axis_label_color'])
        if 'tick_label_color' in settings: self.tick_label_color.set(settings['tick_label_color'])
        if 'legend_text_color' in settings: self.legend_text_color.set(settings['legend_text_color'])
        if 'info_box_text_color' in settings: self.info_box_text_color.set(settings['info_box_text_color'])
        if 'tick_length' in settings: self.tick_length.set(settings['tick_length'])
        if 'tick_width' in settings: self.tick_width.set(settings['tick_width'])
        if 'tick_color' in settings: self.tick_color.set(settings['tick_color'])
        if 'tick_direction' in settings: self.tick_direction.set(settings['tick_direction'])
        # Axes
        if 'x_axis_mode' in settings: self.x_axis_mode.set(settings['x_axis_mode'])
        if 'channel_min' in settings: self.channel_min.set(settings['channel_min'])
        if 'channel_max' in settings: self.channel_max.set(settings['channel_max'])
        if 'energy_min' in settings: self.energy_min.set(settings['energy_min'])
        if 'energy_max' in settings: self.energy_max.set(settings['energy_max'])
        if 'counts_min' in settings: self.counts_min.set(settings['counts_min'])
        if 'counts_max' in settings: self.counts_max.set(settings['counts_max'])
    
    def load_defaults(self):
        """Load default values (already set in creation)"""
        pass


def main():
    root = tk.Tk()
    app = XNRAViewerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
