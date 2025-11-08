cat > README.md << 'EOF'
# RNA Thermometer Finder

A Python tool for identifying and analyzing RNA thermometers in bacterial sequences.

## Features

- 🧬 Detect terminal hairpin structures
- 🌡️ Analyze temperature-dependent MFE (25°C, 37°C, 42°C)
- 📊 Calculate base pair composition (AU%, GC%, GU%)
- 🔍 Identify ribosome binding sites
- 📈 Quality scoring and ranking (0-6 criteria)
- 🖥️ User-friendly GUI

## Quick Start

### For Users (No Python Required)

**Download the application:**
- [📥 macOS Download](../../releases/latest)

**Run:**
1. Unzip the downloaded file
2. Double-click `RNAThermoFinder.app`
3. If macOS blocks it: Right-click → Open

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

# Install Python dependencies
pip install -r requirements.txt

# Run
python main.py
```

## Usage

1. Click "Browse" to select a FASTA file
2. Click "Analyze" to run analysis
3. Results are saved to `Data/Outputs/rna_results.csv`
4. Best candidates (highest quality scores) appear at the top

## Output

Results include:
- Sequence information
- Hairpin structure details
- MFE at multiple temperatures
- Base pair composition
- RBS accessibility
- Quality score (0-6)

## Building from Source
```bash
# Install build tools
pip install pyinstaller

# Build application
python build_app.py

# Package for distribution
cd dist
zip -r ../RNAThermoFinder-v1.0-macOS.zip RNAThermoFinder.app
```

## License

MIT License - See LICENSE file

## Author

Roy Vaknin
EOF