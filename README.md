# RNA Thermometer Finder v2.0

A Python tool for identifying and analyzing RNA thermometers in bacterial sequences with advanced filtering and quality scoring.

## 🆕 What's New in v2.0

### Major Features
- **Original Sequence Quality Score**: New 0-6 metric for complete sequence evaluation
- **Enhanced Filtering**: MFE and composition range checks for original sequences
- **Smart Performance**: Conditional calculations based on selected CSV columns
- **Dual Quality Metrics**: Separate scores for hairpin and original sequence

### New Analysis Settings
- Original sequence MFE ranges (25°C, 37°C, 42°C)
- Original sequence composition ranges (AU%, GC%, GU%)

### New CSV Columns (7 total)
- `Original_MFE_25C_InRange` / `37C` / `42C`
- `Original_AU%_InRange` / `GC%` / `GU%`
- `Quality_Score_Original` (0-6)

---

## Features

- 🧬 Detect terminal hairpin structures
- 🌡️ Analyze temperature-dependent MFE (25°C, 37°C, 42°C)
- 📊 Calculate base pair composition (AU%, GC%, GU%)
- 🔍 Identify ribosome binding sites
- 📈 Dual quality scoring (hairpin + original sequence, 0-6 each)
- ⚙️ Customizable filtering ranges
- 📋 Flexible CSV output configuration
- 🖥️ User-friendly GUI

---

## Quick Start

### For Users (No Python Required)

**Download the application:**
- [📥 v2.0.0 - macOS Download](../../releases/tag/v2.0.0)

**Run:**
1. Unzip the downloaded file
2. Double-click `RNAThermoFinder.app`
3. If macOS blocks it: Right-click → Open → Open

### For Developers

**Prerequisites:**
- Python 3.8+
- ViennaRNA package

**Installation:**
```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/RNAThermoFinder.git
cd RNAThermoFinder

# Install ViennaRNA
brew install viennarna  # macOS
# OR
sudo apt-get install viennarna  # Linux

# Install Python dependencies
pip install -r requirements.txt

# Run
python main.py
```

---

## Usage

### Basic Workflow

1. **Open File**: Click "Browse" to select a FASTA file
2. **Configure** (Optional):
   - Click "⚙️ Analysis Settings" to set filter ranges
   - Click "📊 CSV Output" to choose columns
3. **Analyze**: Click "🧬 Analyze" to run analysis
4. **Export**: Results auto-save to `Data/Outputs/rna_results.csv`

### Analysis Settings

**Hairpin Filters** (Tab 1):
- AU%, GC%, GU% ranges
- MFE ranges at 25°C, 37°C, 42°C

**Original Sequence Filters** (Tab 1 - Bottom):
- Original AU%, GC%, GU% ranges
- Original MFE ranges at all temperatures

### CSV Output Configuration

**Presets:**
- **Hairpin Preset**: Focus on terminal hairpin analysis
- **Preset 2**: Includes original sequence data + range checks

**Custom**: Select individual columns across 11 categories

---

## Output

### CSV Columns

**Basic Info:**
- Name, Complete Sequence, Complete Structure

**Original Sequence:**
- MFE at 25/37/42°C
- AU%, GC%, GU% composition
- Range check results (In Range / Not in Range)
- Quality Score (0-6)

**Terminal Hairpin:**
- Hairpin sequence & structure
- MFE at 25/37/42°C
- AU%, GC%, GU% composition
- Range check results
- Quality Score (0-6)

**RBS Analysis:**
- RBS sequence, structure, paired %

### Quality Scores

**Hairpin Quality (0-6):** Based on 6 hairpin criteria
**Original Quality (0-6):** Based on 6 original sequence criteria

*Higher scores = more criteria met*

---

## Performance Optimization

v2.0 includes smart calculation logic:
- Original sequence MFE only calculated if needed
- Original composition only calculated if needed
- Range checks only run when columns enabled

**Result:** Faster analysis when using hairpin-only mode

---

## Building from Source
```bash
# Install build tools
pip install pyinstaller

# Build application
python build_app.py

# Output:
# macOS: dist/RNAThermoFinder.app
# Windows: dist/RNAThermoFinder.exe
# Linux: dist/RNAThermoFinder
```

---

## Upgrading from v1.0

**Breaking Changes:** None - fully backward compatible

**New Features:**
- Old CSV files will work
- Old settings will work
- New columns disabled by default

**To use new features:**
1. Open "Analysis Settings"
2. Set original sequence ranges
3. Open "CSV Output Settings"
4. Enable desired original sequence columns

---

## Example Use Cases

### Use Case 1: Quick Hairpin Screening
```
Settings: Hairpin Preset (default)
Output: Hairpin-focused columns only
Speed: ⚡ Fast (minimal calculations)
```

### Use Case 2: Comprehensive Analysis
```
Settings: Preset 2
Output: All original + hairpin data
Speed: 🐌 Slower (full calculations)
Use: When you need complete sequence context
```

### Use Case 3: Custom Research
```
Settings: Custom column selection
Output: Only columns you need
Speed: ⚡ Optimized (calculates what you select)
```

---

## Citation

If you use this tool in your research, please cite:
```
Vaknin, R. (2024). RNA Thermometer Finder v2.0: 
Advanced filtering and quality scoring for RNA thermometer discovery.
GitHub: https://github.com/RoyCyber1/RNAThermoFinder
```

---

## Changelog

### v2.0.0 (December 2024)
**Added:**
- Original sequence quality scoring (0-6)
- Original sequence MFE range filters
- Original sequence composition range filters
- 7 new CSV output columns
- Conditional calculation system
- Auto-enable logic for settings

**Improved:**
- Performance when using hairpin-only mode
- Settings dialog organization
- CSV output flexibility

**Fixed:**
- Quality score display name consistency
- Fallback CSV header naming

### v1.0.0 (Initial Release)
- Basic hairpin detection
- Temperature-dependent MFE analysis
- RBS identification
- Quality scoring

---

## License

MIT License - See LICENSE file

## Author

Roy Vaknin  
Email: roycyber13@gmail.com  
GitHub:https://github.com/RoyCyber1/RNAThermoFinder

---

## Acknowledgments

- ViennaRNA package for RNA folding algorithms
- Dr. Abdelsayed for experimental validation data
- SCREAM team for testing and feedback