# Installation Guide

## System Requirements

- **Operating System**: macOS, Linux, or Windows
- **Python**: 3.8 or higher
- **RAM**: 4GB minimum, 8GB recommended
- **CPU**: 2+ cores recommended

## Step-by-Step Installation

### 1. Install Python

Download from [python.org](https://www.python.org/downloads/)

### 2. Install ViennaRNA

#### macOS
```bash
brew install viennarna
```

#### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install viennarna
```

#### Windows
Download and install from [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)

### 3. Install RNA Thermometer Finder
```bash
git clone https://github.com/yourusername/RNAThermoFinder.git
cd RNAThermoFinder
pip install -r requirements.txt
pip install -e .
```

### 4. Verify Installation
```bash
python main.py
```

## Troubleshooting

[Common issues and solutions]