#!/bin/sh
# Package ChuanFlow GUI as a macOS .app bundle (and optionally a DMG).
#
# Tauri 2 embeds the frontendDist assets into the binary at compile time, so
# the bundle needs only the executable + Info.plist + icon. The 1D engine
# executables are bundled into Contents/Resources/bin, where the GUI finds
# them through the `default_engines` command (resolved relative to the
# running binary, so Finder-launched CWD=/ does not matter).
#
# Usage:  gui/package_app.sh [--universal] [--dmg]
#         (expects ../build to hold the engines)
# Output: gui/dist/ChuanFlow GUI.app    (ad-hoc signed, x86_64)
#         gui/dist/ChuanFlow-GUI.dmg    (only with --dmg)
#
# --universal: when arm64 slices of the GUI and the engines exist
#   (target/aarch64-apple-darwin/release/chflow-gui and build-arm64/ trees),
#   lipo-merge them into fat binaries. NOTE: this Intel host can cross-compile
#   an arm64 slice but cannot run or verify it (no Rosetta on Intel), so the
#   default build stays x86_64-only and the flag is opt-in. An Apple Silicon
#   machine produces and verifies both slices natively.
#
# Bundle: steady + transient + coupled-1d3d + two-fluid engines in Resources/bin,
#         the coupled template case and the editor template dictionaries in
#         Resources/cases (case_templates / default_engines find all of them
#         relative to the running binary).
set -eu

UNIVERSAL=0
MAKE_DMG=0
for arg in "$@"; do
    case "$arg" in
        --universal) UNIVERSAL=1 ;;
        --dmg) MAKE_DMG=1 ;;
        *) echo "unknown option: $arg (want --universal / --dmg)" >&2; exit 2 ;;
    esac
done

cd "$(dirname "$0")"
ROOT=$(pwd)/..

echo "== cargo build --release (chflow-gui) =="
cargo build --release

APP="dist/ChuanFlow GUI.app"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources/bin" "$APP/Contents/Resources/cases"

# Merge arm64+x86_64 slices when both exist (--universal); otherwise ship the
# native single-architecture build untouched. lipo refuses mixed-arch sets,
# which is exactly the guard we want.
merge_universal() {
    dst="$1"; arm="$2"; x86="$3"
    if [ "$UNIVERSAL" = 1 ] && [ -f "$arm" ] && [ -f "$x86" ]; then
        lipo -create -output "$dst" "$arm" "$x86"
        echo "   universal: $(lipo -archs "$dst")"
    else
        cp "$x86" "$dst"
    fi
}

merge_universal "$APP/Contents/MacOS/chflow-gui" \
    target/aarch64-apple-darwin/release/chflow-gui \
    target/release/chflow-gui
for exe in chflow-1d-steady chflow-1d-transient chflow-1d-twofluid; do
    merge_universal "$APP/Contents/Resources/bin/$exe" \
        "$ROOT/build-arm64/src/1d/apps/$exe" \
        "$ROOT/build/src/1d/apps/$exe"
done
merge_universal "$APP/Contents/Resources/bin/chflow-coupled-1d3d" \
    "$ROOT/build-arm64/src/coupling/apps/chflow-coupled-1d3d" \
    "$ROOT/build/src/coupling/apps/chflow-coupled-1d3d"

cp "$ROOT"/tests/cases/coupled1d3d/case.dict \
   "$APP/Contents/Resources/cases/coupled1d3d.dict"
for tpl in steadyNetwork steadyEquipmentValve steadyEquipmentPump \
           steadyNetworkMultiphase transientBlowdown twoFluidTerrainSlug; do
    cp "$ROOT/tests/cases/$tpl/case.dict" "$APP/Contents/Resources/cases/$tpl.dict"
done

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
"$APP/Contents/Resources/bin/chflow-1d-twofluid" \
    "$APP/Contents/Resources/cases/twoFluidTerrainSlug.dict" >/dev/null && \
    echo "== bundled two-fluid slug smoke run: ok =="
# GUI binary itself must at least report its architecture
lipo -archs "$APP/Contents/MacOS/chflow-gui"

if [ "$MAKE_DMG" = 1 ]; then
    echo "== hdiutil: ChuanFlow-GUI.dmg =="
    cat > "dist/安装说明-拖入 Applications.txt" <<'NOTE'
ChuanFlow — 自研多相流仿真平台（1D 稳态管网 / 瞬态多相 / 3D FVM / 1D↔3D 耦合）。

安装：把这个窗口里的「ChuanFlow GUI.app」拖进「应用程序」文件夹。
求解器随包（Contents/Resources/bin），无需另行安装。

签名说明：本包为 ad-hoc 签名（无 Apple Developer ID）。首次打开若被
Gatekeeper 拦截，请右键 app → 打开，或在系统设置 → 隐私与安全性中放行。
NOTE
    rm -f dist/ChuanFlow-GUI.dmg /tmp/cf_dmg_rw.dmg
    # This host's hdiutil refuses -srcfolder ("功能未实现"), so build the
    # image the long way: blank HFS+ image -> mount -> copy -> convert UDZO.
    hdiutil create -size 512m -fs 'Journaled HFS+' -volname "ChuanFlow GUI" \
        /tmp/cf_dmg_rw.dmg >/dev/null
    hdiutil attach /tmp/cf_dmg_rw.dmg -nobrowse -mountpoint /tmp/cf_dmg_mount >/dev/null
    cp -R "$APP" /tmp/cf_dmg_mount/
    cp "dist/安装说明-拖入 Applications.txt" /tmp/cf_dmg_mount/
    hdiutil detach /tmp/cf_dmg_mount >/dev/null
    hdiutil convert /tmp/cf_dmg_rw.dmg -format UDZO -o dist/ChuanFlow-GUI.dmg >/dev/null
    rm -f /tmp/cf_dmg_rw.dmg
    echo "== dmg written: dist/ChuanFlow-GUI.dmg =="
    # verify the image mounts and contains the app
    hdiutil attach dist/ChuanFlow-GUI.dmg -readonly -nobrowse -mountpoint /tmp/cf_dmg_verify >/dev/null
    test -d "/tmp/cf_dmg_verify/ChuanFlow GUI.app" && echo "== dmg mount verify: app present =="
    hdiutil detach /tmp/cf_dmg_verify >/dev/null
fi
