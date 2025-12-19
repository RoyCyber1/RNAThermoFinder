import tkinter as tk
from tkinter import ttk, messagebox
from typing import Optional
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from settings_manager import SettingsManager


class SequenceSettingsDialog:
    """Dialog for configuring sequence preprocessing options"""

    def __init__(self, parent, settings_manager: SettingsManager):
        self.parent = parent
        self.settings_manager = settings_manager
        self.dialog = None

        # Variables for settings
        self.append_enabled_var = None
        self.append_sequence_var = None
        self.append_position_var = None
        self.preview_var = None

    def show(self):
        """Display the sequence settings dialog"""
        self.dialog = tk.Toplevel(self.parent)
        self.dialog.title("Sequence Processing Settings")
        self.dialog.geometry("700x500")
        self.dialog.resizable(False, False)

        # Make dialog modal
        self.dialog.transient(self.parent)
        self.dialog.grab_set()

        self._create_widgets()
        self._load_current_settings()
        self._update_preview()

        # Center the dialog
        self.dialog.update_idletasks()
        x = (self.dialog.winfo_screenwidth() // 2) - (self.dialog.winfo_width() // 2)
        y = (self.dialog.winfo_screenheight() // 2) - (self.dialog.winfo_height() // 2)
        self.dialog.geometry(f"+{x}+{y}")

    def _create_widgets(self):
        """Create dialog widgets"""
        main_frame = ttk.Frame(self.dialog, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Title
        title_label = ttk.Label(
            main_frame,
            text="🧬 Sequence Preprocessing Options",
            font=("Arial", 14, "bold")
        )
        title_label.pack(pady=(0, 20))

        # Append Sequence Section
        append_frame = ttk.LabelFrame(main_frame, text="Append Sequence", padding="15")
        append_frame.pack(fill=tk.X, pady=(0, 15))

        # Enable/Disable checkbox
        self.append_enabled_var = tk.BooleanVar()
        enable_check = ttk.Checkbutton(
            append_frame,
            text="Enable sequence appending",
            variable=self.append_enabled_var,
            command=self._on_enable_changed
        )
        enable_check.pack(anchor=tk.W, pady=(0, 10))

        # Sequence input
        seq_input_frame = ttk.Frame(append_frame)
        seq_input_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(seq_input_frame, text="Sequence to append:").pack(side=tk.LEFT, padx=(0, 10))

        self.append_sequence_var = tk.StringVar()
        seq_entry = ttk.Entry(
            seq_input_frame,
            textvariable=self.append_sequence_var,
            width=20
        )
        seq_entry.pack(side=tk.LEFT, padx=(0, 10))

        # Validate button
        validate_btn = ttk.Button(
            seq_input_frame,
            text="Validate",
            command=self._validate_sequence,
            width=10
        )
        validate_btn.pack(side=tk.LEFT)

        # Trace changes to update preview
        self.append_sequence_var.trace_add("write", lambda *args: self._update_preview())

        # Position selection
        position_frame = ttk.Frame(append_frame)
        position_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Label(position_frame, text="Append position:").pack(side=tk.LEFT, padx=(0, 10))

        self.append_position_var = tk.StringVar()

        start_radio = ttk.Radiobutton(
            position_frame,
            text="Start (5' end)",
            variable=self.append_position_var,
            value="start",
            command=self._update_preview
        )
        start_radio.pack(side=tk.LEFT, padx=(0, 15))

        end_radio = ttk.Radiobutton(
            position_frame,
            text="End (3' end)",
            variable=self.append_position_var,
            value="end",
            command=self._update_preview
        )
        end_radio.pack(side=tk.LEFT)

        # Info text
        info_label = ttk.Label(
            append_frame,
            text="💡 Tip: AUG is commonly appended to study RNA thermometer hairpin melting effects",
            font=("Arial", 9),
            foreground="gray",
            wraplength=550
        )
        info_label.pack(anchor=tk.W, pady=(5, 0))

        # Preview Section
        preview_frame = ttk.LabelFrame(main_frame, text="Preview", padding="15")
        preview_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        preview_label = ttk.Label(
            preview_frame,
            text="Example transformation:",
            font=("Arial", 10, "bold")
        )
        preview_label.pack(anchor=tk.W, pady=(0, 10))

        self.preview_var = tk.StringVar()
        preview_text = ttk.Label(
            preview_frame,
            textvariable=self.preview_var,
            font=("Courier", 10),
            foreground="blue",
            wraplength=550,
            justify=tk.LEFT
        )
        preview_text.pack(anchor=tk.W)

        # Buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X)

        save_btn = ttk.Button(
            button_frame,
            text="💾 Save Settings",
            command=self._save_settings,
            width=20
        )
        save_btn.pack(side=tk.LEFT, padx=(0, 10))

        cancel_btn = ttk.Button(
            button_frame,
            text="Cancel",
            command=self.dialog.destroy,
            width=15
        )
        cancel_btn.pack(side=tk.LEFT)

    def _load_current_settings(self):
        """Load current settings from manager"""
        seq_settings = self.settings_manager.settings.get("sequence_processing", {})

        self.append_enabled_var.set(seq_settings.get("append_sequence_enabled", False))
        self.append_sequence_var.set(seq_settings.get("append_sequence", "AUG"))
        self.append_position_var.set(seq_settings.get("append_position", "end"))

    def _on_enable_changed(self):
        """Handle enable/disable checkbox change"""
        self._update_preview()

    def _validate_sequence(self):
        """Validate that sequence contains only valid RNA nucleotides"""
        sequence = self.append_sequence_var.get().upper()

        if not sequence:
            messagebox.showwarning(
                "Empty Sequence",
                "Please enter a sequence to append."
            )
            return False

        valid_chars = set("ACGU")
        invalid_chars = set(sequence) - valid_chars

        if invalid_chars:
            messagebox.showerror(
                "Invalid Sequence",
                f"Sequence contains invalid characters: {', '.join(sorted(invalid_chars))}\n\n"
                "Only A, C, G, U are allowed for RNA sequences."
            )
            return False
        else:
            messagebox.showinfo(
                "Valid Sequence",
                f"✓ Sequence '{sequence}' is valid!"
            )
            # Auto-uppercase
            self.append_sequence_var.set(sequence)
            return True

    def _update_preview(self):
        """Update the preview text"""
        if not self.append_enabled_var.get():
            self.preview_var.set(
                "Sequence appending is disabled.\n\nOriginal: AUGCGAUUCGAGCUAG\nResult:   AUGCGAUUCGAGCUAG")
            return

        append_seq = self.append_sequence_var.get().upper()
        position = self.append_position_var.get()

        example_seq = "AUGCGAUUCGAGCUAG"

        if position == "start":
            result_seq = append_seq + example_seq
        else:
            result_seq = example_seq + append_seq

        preview_text = f"Original: {example_seq}\nResult:   {result_seq}"

        if append_seq:
            # Highlight the appended part
            if position == "start":
                preview_text += f"\n\nAppended '{append_seq}' at 5' end (start)"
            else:
                preview_text += f"\n\nAppended '{append_seq}' at 3' end (end)"

        self.preview_var.set(preview_text)

    def _save_settings(self):
        """Save settings and close dialog"""
        # Validate sequence if enabled
        if self.append_enabled_var.get():
            if not self._validate_sequence():
                return

        # Update settings
        self.settings_manager.settings["sequence_processing"] = {
            "append_sequence_enabled": self.append_enabled_var.get(),
            "append_sequence": self.append_sequence_var.get().upper(),
            "append_position": self.append_position_var.get()
        }

        # Save to file
        self.settings_manager.save_settings()

        messagebox.showinfo(
            "Settings Saved",
            "Sequence processing settings have been saved successfully!"
        )

        self.dialog.destroy()