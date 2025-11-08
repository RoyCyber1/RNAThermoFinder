# User Guide

## GUI Mode

### Starting the Application
```bash
python main.py
```

### Analyzing Sequences

1. Click **Browse** to select a FASTA file
2. Click **Analyze** to start analysis
3. Monitor progress in the output window
4. Results are saved automatically to `Data/Outputs/`

### Exporting Results

- Click **Export ▼** for export options
- Choose **Export to Downloads** for quick save
- Choose **Export to...** to select custom location

## Command-Line Mode
```python
from pathlib import Path
from RnaThermofinder.core import FastaParse, analysis

# Load sequences
sequences = FastaParse.read_fasta("input.fasta")

# Run analysis
output_dir = Path("results")
results = analysis.calculate_results_final(
    sequences,
    output_dir,
    progress_callback=print
)

print(f"Analyzed {len(results)} sequences")
```

## Understanding Results

[Explanation of CSV columns, quality scores, etc.]