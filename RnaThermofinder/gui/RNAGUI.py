import os
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, ttk, messagebox, scrolledtext
from .settings_dialog import SettingsDialog

# Import from core - use relative imports
#from ..core import FastaParse
#from ..core import HairpinAnalysis

from RnaThermofinder.core import FastaParse
from RnaThermofinder.core import HairpinAnalysis



class RNAThermoFinderGUI:
    """Main GUI application for RNA Thermometer Finder"""

    def __init__(self, root):
        self.root = root
        self.root.title("RNA Thermometer Finder")
        self.root.geometry("900x500")

        # State variables
        self.sequences = []
        self.results = []
        self.analysis_settings = {
            'au_min': 50, 'au_max': 60,
            'gc_min': 0, 'gc_max': 30,
            'gu_min': 15, 'gu_max': 25,
            'mfe_25_min': -17, 'mfe_25_max': -10,
            'mfe_37_min': -13, 'mfe_37_max': -6,
            'mfe_42_min': -7, 'mfe_42_max': -2,
        }
        self.status_var= tk.StringVar(value="Ready")

        # Set output directory (use absolute path)
        project_root = Path(__file__).parent.parent.parent
        self.output_dir = project_root / "Data" / "Outputs"
        self.output_dir.mkdir(parents=True, exist_ok=True)


        # Initialize UI
        self._create_widgets()
        self._create_menu()

    def open_settings(self):
        """Open settings dialog"""
        dialog = SettingsDialog(self.root, self.analysis_settings)
        result = dialog.show()
        if result:
            self.analysis_settings = result
            self.log("✓ Settings updated")

    def _create_widgets(self):
        """Create and layout all widgets"""
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights for resizing
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)

        # File selection section
        ttk.Label(main_frame, text="Input File:", font=("Arial", 10, "bold")).grid(
            row=0, column=0, sticky=tk.W, pady=(0, 10)
        )
        self.file_path_var = tk.StringVar()
        ttk.Entry(main_frame, textvariable=self.file_path_var, width=60).grid(
            row=0, column=1, sticky=(tk.W, tk.E), padx=5
        )
        ttk.Button(main_frame, text="Browse", command=self.browse_file).grid(
            row=0, column=2, padx=(5, 0)
        )

        # Results display with scrollbar
        ttk.Label(main_frame, text="Analysis Output:", font=("Arial", 10, "bold")).grid(
            row=1, column=0, sticky=tk.NW, pady=(10, 5)
        )

        # Use ScrolledText for automatic scrollbar
        self.results_text = scrolledtext.ScrolledText(
            main_frame,
            height=30,
            width=90,
            wrap=tk.WORD,
            font=("Courier", 9)
        )
        self.results_text.grid(
            row=2, column=0, columnspan=3,
            sticky=(tk.W, tk.E, tk.N, tk.S),
            pady=(0, 10)
        )

        # Button frame
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=3, column=0, columnspan=3, pady=10)

        # Buttons
        self.analyze_btn = ttk.Button(
            button_frame,
            text="🧬 Analyze",
            command=self.run_analysis
        )
        self.analyze_btn.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            button_frame,
            text="🗑️ Clear",
            command=self.clear_output
        ).pack(side=tk.LEFT, padx=5)

        self.export_btn = ttk.Button(
            button_frame,
            text="💾 Export",
            command=self.export_results,
            state=tk.DISABLED
        )
        self.export_btn.pack(side=tk.LEFT, padx=5)

        ttk.Button(
            button_frame,
            text="⚙️ Settings",
            command=self.open_settings
        ).pack(side=tk.LEFT, padx=5)

        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode='indeterminate')
        self.progress.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(10, 0))

        # Status bar
        status_bar = ttk.Label(
            self.root,
            textvariable=self.status_var,
            relief=tk.SUNKEN,
            anchor=tk.W
        )
        status_bar.grid(row=1, column=0, sticky=(tk.W, tk.E))
    def _create_menu(self):
        """Create menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Open", command=self.browse_file)
        file_menu.add_command(label="Export Results", command=self.export_results)
        file_menu.add_command(label="Settings", command=self.open_settings)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)


    def browse_file(self):
        """Open file dialog for sequence selection"""
        filename = filedialog.askopenfilename(
            title="Select sequence file",
            filetypes=[
                ("FASTA files", "*.fasta *.fa"),
                ("CSV files", "*.csv"),
                ("TSV files", "*.tsv"),
                ("All files", "*.*")
            ]
        )
        if filename:
            self.file_path_var.set(filename)

    def log(self, message):
        """Add message to output text widget"""
        self.results_text.insert(tk.END, message + "\n")
        self.results_text.see(tk.END)
        self.root.update_idletasks()

    def browse_output(self):
        """Open directory dialog for output selection"""
        directory = filedialog.askdirectory(title="Select Output Directory")
        if directory:
            self.output_dir = Path(directory)
            self.output_var.set(str(directory))
            self.log(f"Output directory: {directory}")

    def _update_log(self, message):
        """Internal method to update log (runs in main thread)"""
        self.results_text.insert(tk.END, message + "\n")
        self.results_text.see(tk.END)
        self.root.update_idletasks()

    def clear_output(self):
        """Clear the output text"""
        self.results_text.delete(1.0, tk.END)
        self.results = []
        self.sequences = []
        self.status_var.set("Ready")
        self.export_btn.config(state=tk.DISABLED)

    def run_analysis(self):
        """Execute RNA thermometer analysis in a separate thread"""
        file_path = self.file_path_var.get()

        if not file_path:
            messagebox.showwarning("No File", "Please select a FASTA file first")
            return

        if not Path(file_path).exists():
            messagebox.showerror("Error", "Selected file does not exist")
            return

        # Disable button and start progress
        self.analyze_btn.config(state=tk.DISABLED)
        self.progress.start(10)
        self.clear_output()
        self.status_var.set("Loading sequences...")

        # Run in separate thread to keep GUI responsive
        thread = threading.Thread(target=self._perform_analysis, args=(file_path,))
        thread.daemon = True
        thread.start()

    def _perform_analysis(self, file_path):
        """Perform the actual RNA analysis (runs in separate thread)"""
        try:
            self.status_var.set("Parsing FASTA file...")

            # # Parse FASTA file (convert to RNA)
            # self.sequences = FastaParse.read_fasta(file_path, convert_to_rna=True)
            # self.log(f"Loaded {len(self.sequences)} sequences\n")

            # Detect file type by extension
            file_path_lower = file_path.lower()
            if file_path_lower.endswith((".fa", ".fasta")):
                # Parse FASTA
                self.sequences = FastaParse.read_fasta(file_path, convert_to_rna=True, validate=True)
            elif file_path_lower.endswith(".csv")or file_path_lower.endswith(".tsv"):
                # Parse CSV or TSV
                self.sequences = FastaParse.read_csv_tsv_sequences(file_path, skip_rows=33, seq_col=10, convert_to_rna=True)
            else:
                raise ValueError(f"Unsupported file type: {file_path}")


            self.status_var.set(f"Analyzing {len(self.sequences)} sequences...")

            # Run analysis with log callback (NO PARENTHESES!)

            # Control how many sequences processed
            #max_sequences = 200  # process only 100 sequences
            #self.sequences = self.sequences[:max_sequences]

            self.results = HairpinAnalysis.calculate_results_final(
                self.sequences,
                self.output_dir,
                self.analysis_settings,
                self.log  # ← Pass function reference, not self.log()
            )

            self.status_var.set(f"✅ Analysis complete! Processed {len(self.sequences)} sequences")
            self.export_btn.config(state=tk.NORMAL)

        except Exception as e:
            self.status_var.set("❌ Error occurred")
            self.log(f"\n❌ ERROR: {str(e)}")

            # Show detailed error
            import traceback
            error_details = traceback.format_exc()
            self.log(error_details)
            messagebox.showerror("Analysis Error", f"An error occurred:\n{str(e)}")

        finally:
            # Re-enable button and stop progress
            self.root.after(0, lambda: self.analyze_btn.config(state=tk.NORMAL))
            self.root.after(0, self.progress.stop)




    def _display_results(self):
        """Display results in text widget"""
        self.results_text.delete(1.0, tk.END)
        if self.results:
            self.log(f"Displaying {len(self.results)} results")
            for result in self.results:
                self.log(str(result))
        else:
            self.log("No results to display")

    def _open_folder(self, folder_path: Path):
        """Open folder in file explorer (cross-platform)"""
        import sys
        import subprocess

        try:
            if sys.platform == 'win32':  # Windows
                os.startfile(folder_path)
            elif sys.platform == 'darwin':  # macOS
                subprocess.call(['open', folder_path])
            else:  # Linux
                subprocess.call(['xdg-open', folder_path])
        except Exception as e:
            self.log(f"Could not open folder: {str(e)}")


    def export_results(self):
        """Export results to user-selected location"""
        if not self.results:
            messagebox.showwarning("No Results", "Run analysis first")
            return

        try:
            # Ask user where to save the file
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            default_filename = f"rna_results_{timestamp}.csv"

            output_file = filedialog.asksaveasfilename(
                title="Save Results As",
                defaultextension=".csv",
                filetypes=[
                    ("CSV files", "*.csv"),
                    ("All files", "*.*")
                ],
                initialfile=default_filename,
                initialdir=str(Path.home() / "Downloads")  # Start in Downloads
            )

            # User cancelled
            if not output_file:
                return

            # Copy the CSV from output directory to selected location
            source_csv = self.output_dir / "rna_results.csv"
            if source_csv.exists():
                import shutil
                shutil.copy2(source_csv, output_file)

                self.log(f"\n✅ Results exported to: {output_file}")
                messagebox.showinfo(
                    "Export Success",
                    f"Results exported successfully!\n\n{Path(output_file).name}"
                )

                # Ask if user wants to open the folder
                if messagebox.askyesno("Open Folder?", "Would you like to open the folder containing the file?"):
                    self._open_folder(Path(output_file).parent)
            else:
                messagebox.showerror("Error", "Results file not found. Please run analysis first.")

        except Exception as e:
            self.log(f"❌ Export failed: {str(e)}")
            messagebox.showerror("Export Error", f"Failed to export:\n{str(e)}")


def main():
    """Entry point for GUI application"""
    root = tk.Tk()
    app = RNAThermoFinderGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()

