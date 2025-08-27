#!/usr/bin/env python3
"""
macOS build script for TPP Plotter
Run this on macOS to create TPP Plotter.app
"""

import subprocess
import sys


def build_mac():
    """Build macOS app bundle without console"""
    print("Building TPP Plotter for macOS...")

    args = [
        "uv",
        "run",
        "pyinstaller",
        "--noconsole",  # No console window, creates .app bundle
        "--onedir",  # Directory bundle (recommended for macOS)
        "--name=TPP Plotter",
        "--add-data=src/styles.qss:.",  # Include QSS file (macOS/Linux uses : as separator)
        "src/main.py",
    ]

    try:
        subprocess.run(args, check=True)
        print("✅ macOS build successful!")
        print("Output: dist/TPP Plotter.app")
        print("Users can drag this to Applications folder")
    except subprocess.CalledProcessError as e:
        print(f"❌ Build failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    build_mac()
