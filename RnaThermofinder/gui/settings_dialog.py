import tkinter as tk
from tkinter import ttk, messagebox


class SettingsDialog:
    """Dialog for configuring analysis parameters"""

    def __init__(self, parent, current_settings=None):
        self.result = None
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Analysis Settings")
        self.dialog.geometry("600x700")  # ✨ Made taller
        self.dialog.resizable(False, False)

        # Make dialog modal
        self.dialog.transient(parent)
        self.dialog.grab_set()

        # Default settings
        self.settings = current_settings or {
            # Hairpin settings
            'au_min': 50, 'au_max': 60,
            'gc_min': 0, 'gc_max': 30,
            'gu_min': 15, 'gu_max': 25,
            'mfe_25_min': -17, 'mfe_25_max': -10,
            'mfe_37_min': -13, 'mfe_37_max': -6,
            'mfe_42_min': -7, 'mfe_42_max': -2,

            # ✨ NEW: Original sequence settings
            'orig_mfe_25_min': -30, 'orig_mfe_25_max': -10,
            'orig_mfe_37_min': -25, 'orig_mfe_37_max': -5,
            'orig_mfe_42_min': -20, 'orig_mfe_42_max': -2,
            'orig_au_min': 0, 'orig_au_max': 100,
            'orig_gc_min': 0, 'orig_gc_max': 100,
            'orig_gu_min': 0, 'orig_gu_max': 100,
        }

        self._create_widgets()

        # Center dialog on parent
        self.dialog.update_idletasks()
        x = parent.winfo_x() + (parent.winfo_width() // 2) - (self.dialog.winfo_width() // 2)
        y = parent.winfo_y() + (parent.winfo_height() // 2) - (self.dialog.winfo_height() // 2)
        self.dialog.geometry(f"+{x}+{y}")

    def _create_widgets(self):
        """Create all dialog widgets"""
        main_frame = ttk.Frame(self.dialog, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Title
        title_label = ttk.Label(
            main_frame,
            text="Configure Analysis Ranges",
            font=("Arial", 12, "bold")
        )
        title_label.pack(pady=(0, 20))

        # Create notebook for organized tabs
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.BOTH, expand=True, pady=(0, 20))

        # Base Pair Composition Tab
        bp_frame = ttk.Frame(notebook, padding="15")
        notebook.add(bp_frame, text="Base Pair Composition")
        self._create_bp_settings(bp_frame)

        # MFE Tab
        mfe_frame = ttk.Frame(notebook, padding="15")
        notebook.add(mfe_frame, text="MFE at Temperatures")
        self._create_mfe_settings(mfe_frame)

        # Button frame
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))

        ttk.Button(
            button_frame,
            text="Reset to Defaults",
            command=self._reset_defaults
        ).pack(side=tk.LEFT, padx=5)

        ttk.Button(
            button_frame,
            text="Cancel",
            command=self.dialog.destroy
        ).pack(side=tk.RIGHT, padx=5)

        ttk.Button(
            button_frame,
            text="Save",
            command=self._save_settings
        ).pack(side=tk.RIGHT, padx=5)

    def _create_range_input(self, parent, label_text, min_key, max_key, row):
        """Create a min-max range input"""
        ttk.Label(parent, text=label_text, font=("Arial", 10)).grid(
            row=row, column=0, sticky=tk.W, pady=8
        )

        # Min input
        ttk.Label(parent, text="Min:").grid(row=row, column=1, sticky=tk.E, padx=(10, 5))
        min_var = tk.DoubleVar(value=self.settings.get(min_key, 0))
        min_entry = ttk.Entry(parent, textvariable=min_var, width=10)
        min_entry.grid(row=row, column=2, sticky=tk.W)

        # Max input
        ttk.Label(parent, text="Max:").grid(row=row, column=3, sticky=tk.E, padx=(15, 5))
        max_var = tk.DoubleVar(value=self.settings.get(max_key, 100))
        max_entry = ttk.Entry(parent, textvariable=max_var, width=10)
        max_entry.grid(row=row, column=4, sticky=tk.W)

        return min_var, max_var

    def _create_bp_settings(self, parent):
        """Create base pair composition settings"""
        # ✨ Hairpin Section
        ttk.Label(
            parent,
            text="Hairpin Base Pair Percentages:",
            font=("Arial", 10, "bold")
        ).grid(row=0, column=0, columnspan=5, sticky=tk.W, pady=(0, 10))

        # Hairpin AU pairs
        self.au_min_var, self.au_max_var = self._create_range_input(
            parent, "Hairpin AU Pairs (%)", 'au_min', 'au_max', 1
        )

        # Hairpin GC pairs
        self.gc_min_var, self.gc_max_var = self._create_range_input(
            parent, "Hairpin GC Pairs (%)", 'gc_min', 'gc_max', 2
        )

        # Hairpin GU pairs
        self.gu_min_var, self.gu_max_var = self._create_range_input(
            parent, "Hairpin GU Pairs (%)", 'gu_min', 'gu_max', 3
        )

        # Info label
        info_label = ttk.Label(
            parent,
            text="💡 Typical RNA thermometer hairpin ranges:\nAU: 50-60%, GC: 0-30%, GU: 15-25%",
            font=("Arial", 8),
            foreground="gray"
        )
        info_label.grid(row=4, column=0, columnspan=5, sticky=tk.W, pady=(10, 0))

        # ✨ Separator
        ttk.Separator(parent, orient=tk.HORIZONTAL).grid(
            row=5, column=0, columnspan=5, sticky=(tk.W, tk.E), pady=20
        )

        # ✨ NEW: Original Sequence Section
        ttk.Label(
            parent,
            text="Original Sequence Base Pair Percentages:",
            font=("Arial", 10, "bold")
        ).grid(row=6, column=0, columnspan=5, sticky=tk.W, pady=(0, 10))

        # Original AU pairs
        self.orig_au_min_var, self.orig_au_max_var = self._create_range_input(
            parent, "Original AU Pairs (%)", 'orig_au_min', 'orig_au_max', 7
        )

        # Original GC pairs
        self.orig_gc_min_var, self.orig_gc_max_var = self._create_range_input(
            parent, "Original GC Pairs (%)", 'orig_gc_min', 'orig_gc_max', 8
        )

        # Original GU pairs
        self.orig_gu_min_var, self.orig_gu_max_var = self._create_range_input(
            parent, "Original GU Pairs (%)", 'orig_gu_min', 'orig_gu_max', 9
        )

        # Info label
        info_label2 = ttk.Label(
            parent,
            text="💡 Wide ranges (0-100) effectively disable filtering",
            font=("Arial", 8),
            foreground="gray"
        )
        info_label2.grid(row=10, column=0, columnspan=5, sticky=tk.W, pady=(10, 0))

    def _create_mfe_settings(self, parent):
        """Create MFE temperature settings"""
        # ✨ Hairpin Section
        ttk.Label(
            parent,
            text="Hairpin MFE Ranges (kcal/mol):",
            font=("Arial", 10, "bold")
        ).grid(row=0, column=0, columnspan=5, sticky=tk.W, pady=(0, 10))

        # Hairpin MFE at 25°C
        self.mfe_25_min_var, self.mfe_25_max_var = self._create_range_input(
            parent, "Hairpin MFE at 25°C", 'mfe_25_min', 'mfe_25_max', 1
        )

        # Hairpin MFE at 37°C
        self.mfe_37_min_var, self.mfe_37_max_var = self._create_range_input(
            parent, "Hairpin MFE at 37°C", 'mfe_37_min', 'mfe_37_max', 2
        )

        # Hairpin MFE at 42°C
        self.mfe_42_min_var, self.mfe_42_max_var = self._create_range_input(
            parent, "Hairpin MFE at 42°C", 'mfe_42_min', 'mfe_42_max', 3
        )

        # Info label
        info_label = ttk.Label(
            parent,
            text="💡 Typical RNA thermometer hairpin MFE ranges:\n"
                 "25°C: -17 to -10, 37°C: -13 to -6, 42°C: -7 to -2",
            font=("Arial", 8),
            foreground="gray"
        )
        info_label.grid(row=4, column=0, columnspan=5, sticky=tk.W, pady=(10, 0))

        # ✨ Separator
        ttk.Separator(parent, orient=tk.HORIZONTAL).grid(
            row=5, column=0, columnspan=5, sticky=(tk.W, tk.E), pady=20
        )

        # ✨ NEW: Original Sequence Section
        ttk.Label(
            parent,
            text="Original Sequence MFE Ranges (kcal/mol):",
            font=("Arial", 10, "bold")
        ).grid(row=6, column=0, columnspan=5, sticky=tk.W, pady=(0, 10))

        # Original MFE at 25°C
        self.orig_mfe_25_min_var, self.orig_mfe_25_max_var = self._create_range_input(
            parent, "Original MFE at 25°C", 'orig_mfe_25_min', 'orig_mfe_25_max', 7
        )

        # Original MFE at 37°C
        self.orig_mfe_37_min_var, self.orig_mfe_37_max_var = self._create_range_input(
            parent, "Original MFE at 37°C", 'orig_mfe_37_min', 'orig_mfe_37_max', 8
        )

        # Original MFE at 42°C
        self.orig_mfe_42_min_var, self.orig_mfe_42_max_var = self._create_range_input(
            parent, "Original MFE at 42°C", 'orig_mfe_42_min', 'orig_mfe_42_max', 9
        )

        # Info label
        info_label2 = ttk.Label(
            parent,
            text="💡 Original sequence MFE is typically more negative than hairpin alone",
            font=("Arial", 8),
            foreground="gray"
        )
        info_label2.grid(row=10, column=0, columnspan=5, sticky=tk.W, pady=(10, 0))

    def _validate_settings(self):
        """Validate that min < max for all ranges"""
        checks = [
            # Hairpin settings
            (self.au_min_var.get(), self.au_max_var.get(), "Hairpin AU%"),
            (self.gc_min_var.get(), self.gc_max_var.get(), "Hairpin GC%"),
            (self.gu_min_var.get(), self.gu_max_var.get(), "Hairpin GU%"),
            (self.mfe_25_min_var.get(), self.mfe_25_max_var.get(), "Hairpin MFE at 25°C"),
            (self.mfe_37_min_var.get(), self.mfe_37_max_var.get(), "Hairpin MFE at 37°C"),
            (self.mfe_42_min_var.get(), self.mfe_42_max_var.get(), "Hairpin MFE at 42°C"),

            # ✨ NEW: Original sequence settings
            (self.orig_mfe_25_min_var.get(), self.orig_mfe_25_max_var.get(), "Original MFE at 25°C"),
            (self.orig_mfe_37_min_var.get(), self.orig_mfe_37_max_var.get(), "Original MFE at 37°C"),
            (self.orig_mfe_42_min_var.get(), self.orig_mfe_42_max_var.get(), "Original MFE at 42°C"),
            (self.orig_au_min_var.get(), self.orig_au_max_var.get(), "Original AU%"),
            (self.orig_gc_min_var.get(), self.orig_gc_max_var.get(), "Original GC%"),
            (self.orig_gu_min_var.get(), self.orig_gu_max_var.get(), "Original GU%"),
        ]

        for min_val, max_val, name in checks:
            if min_val >= max_val:
                messagebox.showerror(
                    "Invalid Range",
                    f"{name}: Minimum must be less than maximum"
                )
                return False

        return True

    def _save_settings(self):
        """Save settings and close dialog"""
        if not self._validate_settings():
            return

        self.result = {
            # Hairpin settings
            'au_min': self.au_min_var.get(),
            'au_max': self.au_max_var.get(),
            'gc_min': self.gc_min_var.get(),
            'gc_max': self.gc_max_var.get(),
            'gu_min': self.gu_min_var.get(),
            'gu_max': self.gu_max_var.get(),
            'mfe_25_min': self.mfe_25_min_var.get(),
            'mfe_25_max': self.mfe_25_max_var.get(),
            'mfe_37_min': self.mfe_37_min_var.get(),
            'mfe_37_max': self.mfe_37_max_var.get(),
            'mfe_42_min': self.mfe_42_min_var.get(),
            'mfe_42_max': self.mfe_42_max_var.get(),

            # ✨ NEW: Original sequence settings
            'orig_mfe_25_min': self.orig_mfe_25_min_var.get(),
            'orig_mfe_25_max': self.orig_mfe_25_max_var.get(),
            'orig_mfe_37_min': self.orig_mfe_37_min_var.get(),
            'orig_mfe_37_max': self.orig_mfe_37_max_var.get(),
            'orig_mfe_42_min': self.orig_mfe_42_min_var.get(),
            'orig_mfe_42_max': self.orig_mfe_42_max_var.get(),
            'orig_au_min': self.orig_au_min_var.get(),
            'orig_au_max': self.orig_au_max_var.get(),
            'orig_gc_min': self.orig_gc_min_var.get(),
            'orig_gc_max': self.orig_gc_max_var.get(),
            'orig_gu_min': self.orig_gu_min_var.get(),
            'orig_gu_max': self.orig_gu_max_var.get(),
        }
        self.dialog.destroy()

    def _reset_defaults(self):
        """Reset all settings to default values"""
        defaults = {
            # Hairpin defaults
            'au_min': 50, 'au_max': 60,
            'gc_min': 0, 'gc_max': 30,
            'gu_min': 15, 'gu_max': 25,
            'mfe_25_min': -17, 'mfe_25_max': -10,
            'mfe_37_min': -13, 'mfe_37_max': -6,
            'mfe_42_min': -7, 'mfe_42_max': -2,

            # ✨ NEW: Original sequence defaults
            'orig_mfe_25_min': -30, 'orig_mfe_25_max': -10,
            'orig_mfe_37_min': -25, 'orig_mfe_37_max': -5,
            'orig_mfe_42_min': -20, 'orig_mfe_42_max': -2,
            'orig_au_min': 0, 'orig_au_max': 100,
            'orig_gc_min': 0, 'orig_gc_max': 100,
            'orig_gu_min': 0, 'orig_gu_max': 100,
        }

        # Update hairpin variables
        self.au_min_var.set(defaults['au_min'])
        self.au_max_var.set(defaults['au_max'])
        self.gc_min_var.set(defaults['gc_min'])
        self.gc_max_var.set(defaults['gc_max'])
        self.gu_min_var.set(defaults['gu_min'])
        self.gu_max_var.set(defaults['gu_max'])
        self.mfe_25_min_var.set(defaults['mfe_25_min'])
        self.mfe_25_max_var.set(defaults['mfe_25_max'])
        self.mfe_37_min_var.set(defaults['mfe_37_min'])
        self.mfe_37_max_var.set(defaults['mfe_37_max'])
        self.mfe_42_min_var.set(defaults['mfe_42_min'])
        self.mfe_42_max_var.set(defaults['mfe_42_max'])

        # ✨ NEW: Update original sequence variables
        self.orig_mfe_25_min_var.set(defaults['orig_mfe_25_min'])
        self.orig_mfe_25_max_var.set(defaults['orig_mfe_25_max'])
        self.orig_mfe_37_min_var.set(defaults['orig_mfe_37_min'])
        self.orig_mfe_37_max_var.set(defaults['orig_mfe_37_max'])
        self.orig_mfe_42_min_var.set(defaults['orig_mfe_42_min'])
        self.orig_mfe_42_max_var.set(defaults['orig_mfe_42_max'])
        self.orig_au_min_var.set(defaults['orig_au_min'])
        self.orig_au_max_var.set(defaults['orig_au_max'])
        self.orig_gc_min_var.set(defaults['orig_gc_min'])
        self.orig_gc_max_var.set(defaults['orig_gc_max'])
        self.orig_gu_min_var.set(defaults['orig_gu_min'])
        self.orig_gu_max_var.set(defaults['orig_gu_max'])

    def show(self):
        """Show dialog and wait for result"""
        self.dialog.wait_window()
        return self.result