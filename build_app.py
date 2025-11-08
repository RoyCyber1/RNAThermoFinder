"""
Build script for creating standalone executables
"""

import PyInstaller.__main__
import sys
import platform
import os


def build_app():
    """Build standalone application"""

    system = platform.system()

    # Find main.py location
    if os.path.exists('main.py'):
        main_script = 'main.py'
    elif os.path.exists('RnaThermofinder/main.py'):
        main_script = 'RnaThermofinder/main.py'
    else:
        print("❌ Error: Cannot find main.py")
        sys.exit(1)

    print(f"✓ Found entry point: {main_script}")

    # Base options for all platforms
    base_options = [
        main_script,
        '--name=RNAThermoFinder',
        '--clean',
        '--noconfirm',

        # Hidden imports
        '--hidden-import=RNA',
        '--hidden-import=tkinter',
        '--hidden-import=pandas',
        '--hidden-import=numpy',
        '--hidden-import=Bio',
        '--hidden-import=PIL',
        '--collect-all=RnaThermofinder',

        # Exclude unnecessary
        '--exclude-module=matplotlib',
        '--exclude-module=IPython',
        '--exclude-module=pytest',
        '--exclude-module=jupyter',
    ]

    # Platform-specific options
    if system == "Darwin":  # macOS
        print("Building for macOS...")
        options = base_options + [
            '--windowed',  # Creates .app bundle
            '--onedir',  # Recommended for macOS (not onefile)
            '--osx-bundle-identifier=com.yourname.rnathermofinder',
        ]
        output_msg = "open dist/RNAThermoFinder.app"

    elif system == "Windows":
        print("Building for Windows...")
        options = base_options + [
            '--windowed',  # No console
            '--onefile',  # Single .exe (works fine on Windows)
            # '--icon=icon.ico',  # Add if you have an icon
        ]
        output_msg = "dist\\RNAThermoFinder.exe"

    else:  # Linux
        print("Building for Linux...")
        options = base_options + [
            '--onefile',  # Single executable
        ]
        output_msg = "./dist/RNAThermoFinder"

    print(f"This may take several minutes...")

    try:
        PyInstaller.__main__.run(options)
        print("\n" + "=" * 60)
        print("✅ Build complete!")
        print("=" * 60)
        print(f"\nTo run: {output_msg}")

        if system == "Darwin":
            print("\nYour .app is in: dist/RNAThermoFinder.app")
            print("To distribute: Compress it to a .zip file")
        elif system == "Windows":
            print("\nYour .exe is in: dist/RNAThermoFinder.exe")
        else:
            print("\nYour executable is in: dist/RNAThermoFinder")

    except Exception as e:
        print(f"\n❌ Build failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    build_app()