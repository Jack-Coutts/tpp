#!/usr/bin/env python3
"""
Windows build script for TPP Plotter
Run this on Windows to create TPP Plotter.exe
"""

import subprocess
import sys


def build_windows():
    """Build Windows executable without console"""
    print("Building TPP Plotter for Windows...")

    args = [
        "uv",
        "run",
        "pyinstaller",
        "--noconsole",  # No console window
        "--onefile",  # Single executable file
        "--name=TPP Plotter",
        "--add-data=src/styles.qss;.",  # Include QSS file (Windows uses ; as separator)
        "src/main.py",
    ]

    try:
        subprocess.run(args, check=True)
        print("✅ Windows build successful!")
        print("Output: dist/TPP Plotter.exe")
    except subprocess.CalledProcessError as e:
        print(f"❌ Build failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    build_windows()
