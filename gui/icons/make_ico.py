#!/usr/bin/env python3
"""Generate gui/icons/icon.ico from icon.png (tauri-build requires an .ico on
Windows; the repo otherwise only carries the macOS .icns).

The Windows SDK rc.exe rejects PNG-compressed ICO entries ("old DIB",
RC2176), so every entry is written as a classic 32-bpp BMP DIB. Sizes follow
the Tauri recommendation; icon.png must be at least 256x256.

Usage:  python make_ico.py
"""
from pathlib import Path

from PIL import Image

here = Path(__file__).parent
SIZES = (16, 24, 32, 48, 64, 128, 256)

src = Image.open(here / "icon.png").convert("RGBA")
if min(src.size) < max(SIZES):
    raise SystemExit("icon.png must be at least %dx%d" % (max(SIZES), max(SIZES)))
src.save(here / "icon.ico", format="ICO", sizes=[(s, s) for s in SIZES])
print("wrote", here / "icon.ico", (here / "icon.ico").stat().st_size, "bytes")

