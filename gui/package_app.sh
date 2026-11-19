#!/bin/sh
# Package ChuanFlow GUI as a macOS .app bundle.
#
# Tauri 2 embeds the frontendDist assets into the binary at compile time, so
# the bundle needs only the executable + Info.plist + icon. The 1D engine
# executables are bundled into Contents/Resources/bin, where the GUI finds
# them through the `default_engines` command (resolved relative to the
# running binary, so Finder-launched CWD=/ does not matter).
#
# Usage:  gui/package_app.sh            (expects ../build to hold the engines)
# Output: gui/dist/ChuanFlow GUI.app    (ad-hoc signed, x86_64)
# Bundle: steady + transient + coupled-1d3d engines in Resources/bin, the
#         coupling template case in Resources/cases (default_engines finds
#         all four relative to the running binary).
set -eu

cd "$(dirname "$0")"
ROOT=$(pwd)/..

echo "== cargo build --release (chflow-gui) =="
cargo build --release

APP="dist/ChuanFlow GUI.app"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources/bin" "$APP/Contents/Resources/cases"

cp target/release/chflow-gui "$APP/Contents/MacOS/chflow-gui"
cp "$ROOT"/build/src/1d/apps/chflow-1d-steady \
   "$ROOT"/build/src/1d/apps/chflow-1d-transient \
   "$ROOT"/build/src/coupling/apps/chflow-coupled-1d3d \
   "$APP/Contents/Resources/bin/"
cp "$ROOT"/tests/cases/coupled1d3d/case.dict \
   "$APP/Contents/Resources/cases/coupled1d3d.dict"

cat > "$APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
	<key>CFBundleDevelopmentRegion</key>
	<string>zh_CN</string>
	<key>CFBundleExecutable</key>
	<string>chflow-gui</string>
	<key>CFBundleIconFile</key>
	<string>chflow</string>
	<key>CFBundleIdentifier</key>
	<string>com.chuanflow.gui</string>
	<key>CFBundleInfoDictionaryVersion</key>
	<string>6.0</string>
	<key>CFBundleName</key>
	<string>ChuanFlow GUI</string>
	<key>CFBundleDisplayName</key>
	<string>ChuanFlow GUI</string>
	<key>CFBundlePackageType</key>
	<string>APPL</string>
	<key>CFBundleShortVersionString</key>
	<string>0.1.0</string>
	<key>CFBundleVersion</key>
	<string>0.1.0</string>
	<key>LSMinimumSystemVersion</key>
	<string>12.0</string>
	<key>NSHighResolutionCapable</key>
	<true/>
	<key>NSPrincipalClass</key>
	<string>NSApplication</string>
</dict>
</plist>
PLIST

# icon: png -> icns (best effort; the app runs fine without it)
if command -v sips >/dev/null 2>&1; then
    sips -s format icns icons/icon.png --out "$APP/Contents/Resources/chflow.icns" >/dev/null 2>&1 || true
fi

echo "== ad-hoc codesign =="
codesign -s - --force --deep "$APP"

echo "== packaged: $APP =="
"$APP/Contents/Resources/bin/chflow-1d-steady" \
    "$ROOT"/tests/cases/steadyAllBelowThreshold/case.dict >/dev/null && \
    echo "== bundled solver smoke run: ok =="
"$APP/Contents/Resources/bin/chflow-coupled-1d3d" \
    "$APP/Contents/Resources/cases/coupled1d3d.dict" >/dev/null && \
    echo "== bundled coupled solver smoke run: ok =="
