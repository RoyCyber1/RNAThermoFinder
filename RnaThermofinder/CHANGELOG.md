# Changelog

All notable changes to RNA Thermometer Finder will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-12-18

### Added
- Original sequence quality scoring system (0-6 metric)
- Original sequence MFE range filters (25°C, 37°C, 42°C)
- Original sequence composition range filters (AU%, GC%, GU%)
- 7 new CSV output columns:
  - `Original_MFE_25C_InRange`
  - `Original_MFE_37C_InRange`
  - `Original_MFE_42C_InRange`
  - `Original_AU%_InRange`
  - `Original_GC%_InRange`
  - `Original_GU%_InRange`
  - `Quality_Score_Original`
- Conditional calculation system (only calculates what's needed)
- Auto-enable logic in CSV settings dialog
- "Complete Sequence Range Checks" group in CSV Output Settings
- Original sequence filter inputs in Analysis Settings

### Changed
- Analysis Settings dialog now includes original sequence ranges (bottom of tabs)
- CSV Output Settings expanded to 800px height
- Quality score naming: `Quality_Score` → `Quality_Score_Hairpin`
- Settings manager includes `calculate_original_range_checks` setting
- Performance optimization: calculations based on enabled columns

### Fixed
- Fallback CSV headers now include all original sequence columns
- Column mapping consistency in analysis_helpers.py
- Quality score display names in settings_manager.py

## [1.0.0] - 2024-11-10

### Added
- Initial release
- Terminal hairpin detection
- Temperature-dependent MFE analysis (25°C, 37°C, 42°C)
- Base pair composition analysis (AU%, GC%, GU%)
- RBS identification
- Quality scoring (0-6 for hairpin)
- GUI interface
- CSV export
- Analysis settings dialog
- macOS .app bundle distribution