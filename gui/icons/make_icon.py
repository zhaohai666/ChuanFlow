#!/usr/bin/env python3
"""
Generate the ChuanFlow application icon: a 1024x1024 PNG plus a macOS
multi-resolution .iconset / .icns.

Design: a dark industrial squircle carrying a glowing 1D pipeline (green ->
cyan stream with discrete fluid parcels) through an orange bidirectional
coupling hub into an isometric 3D CFD block. Mirrors the product identity:
1D network + flow assurance + 3D CFD, coupled.

Every size is rendered independently at 4x supersampling. Small sizes use a
zoom factor (graphics occupy more of the tile) and thicker strokes so the
pipe, hub and cube stay legible down to 16pt. No fonts or external assets -
fully reproducible from stdlib + Pillow.

Usage:  python3 make_icon.py            # writes icon.png, chflow.iconset, chflow.icns
"""

from PIL import Image, ImageDraw, ImageFilter
import math
import os
import shutil
import subprocess

HERE = os.path.dirname(__file__)
S = 4  # supersampling factor for every tile


def render(size, zoom=1.0, boost=1.0, brighten=0):
    """Render one icon tile at `size` pts (4x supersampled internally)."""
    N = size * S

    def px(v):
        # v expressed relative to a 1024 canvas
        return max(1, int(round(v * S * size / 1024.0 * boost)))

    def Z(p):
        # zoom all *graphics* about the centre; background squircle stays put
        return (0.5 + (p[0] - 0.5) * zoom, 0.5 + (p[1] - 0.5) * zoom)

    def tv(p):
        p = Z(p)
        return (p[0] * N, p[1] * N)

    # ---- background squircle: deep blue gradient --------------------------
    grad = Image.new("RGBA", (1, N))
    for y in range(N):
        t = y / max(1, N - 1)
        top, bottom = (17 + brighten, 38 + brighten, 68 + brighten), (8 + brighten, 18 + brighten, 34 + brighten)
        grad.putpixel((0, y), tuple(int(top[i] + (bottom[i] - top[i]) * t) for i in range(3)) + (255,))
    grad = grad.resize((N, N), Image.Resampling.BILINEAR)
    mask = Image.new("L", (N, N), 0)
    ImageDraw.Draw(mask).rounded_rectangle(
        [px(8), px(8), N - px(8), N - px(8)], radius=px(224), fill=255
    )
    bg = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    bg.paste(grad, (0, 0), mask)

    # soft radial glow behind the graphics
    glow = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    gd = ImageDraw.Draw(glow)
    a, b = Z((0.12, 0.18)), Z((0.90, 0.88))
    gd.ellipse([a[0] * N, a[1] * N, b[0] * N, b[1] * N], fill=(38, 116, 190, 90))
    glow = glow.filter(ImageFilter.GaussianBlur(px(90)))
    bg.alpha_composite(glow)

    # ---- isometric 3D CFD block -------------------------------------------
    cx, cy = Z((0.700, 0.545))
    a, b, h = 0.205 * zoom, 0.112 * zoom, 0.168 * zoom
    top = [(cx, cy - b), (cx + a, cy), (cx, cy + b), (cx - a, cy)]
    bot = [(x, y + h) for x, y in top]
    vt = lambda p: (p[0] * N, p[1] * N)

    side = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    sd = ImageDraw.Draw(side)
    sd.polygon([vt(top[1]), vt(top[2]), vt(bot[2]), vt(bot[1])], fill=(14, 62, 86, 200))
    sd.polygon([vt(top[2]), vt(top[3]), vt(bot[3]), vt(bot[2])], fill=(10, 46, 70, 220))
    bg.alpha_composite(side)

    def lerp(p, q, f):
        return (p[0] + (q[0] - p[0]) * f, p[1] + (q[1] - p[1]) * f)

    grid = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    gridd = ImageDraw.Draw(grid)
    grid_col = (72 + brighten, 164 + brighten, 220 + brighten, 135)
    for k in range(1, 5):
        f = k / 5
        gridd.line([vt(lerp(top[3], top[0], f)), vt(lerp(top[2], top[1], f))], fill=grid_col, width=px(1.8))
        gridd.line([vt(lerp(top[0], top[1], f)), vt(lerp(top[3], top[2], f))], fill=grid_col, width=px(1.8))
        gridd.line([vt(lerp(top[1], top[2], f)), vt(lerp(bot[1], bot[2], f))], fill=grid_col, width=px(1.8))
        gridd.line([vt(lerp(top[2], top[3], f)), vt(lerp(bot[2], bot[3], f))], fill=grid_col, width=px(1.8))
    bg.alpha_composite(grid)

    ed = ImageDraw.Draw(bg)
    edge = (255, 159, 67, 255)
    for i in range(4):
        ed.line([vt(top[i]), vt(top[(i + 1) % 4])], fill=edge, width=px(5))
    for i in (1, 2, 3):
        ed.line([vt(top[i]), vt(bot[i])], fill=edge, width=px(4))
    for i in range(4):
        ed.line([vt(bot[i]), vt(bot[(i + 1) % 4])], fill=(214, 120, 45, 220), width=px(3))

    # ---- 1D pipeline -------------------------------------------------------
    ctrl = [Z((0.075, 0.79)), Z((0.30, 0.80)), Z((0.31, 0.40)), Z((0.605, 0.44))]
    pipe = []
    for i in range(201):
        t = i / 200
        u = 1 - t
        pipe.append((
            N * (u**3 * ctrl[0][0] + 3*u*u*t * ctrl[1][0] + 3*u*t*t * ctrl[2][0] + t**3 * ctrl[3][0]),
            N * (u**3 * ctrl[0][1] + 3*u*u*t * ctrl[1][1] + 3*u*t*t * ctrl[2][1] + t**3 * ctrl[3][1]),
        ))

    def ribbon(pts, half_width):
        left, right = [], []
        for i, (x, y) in enumerate(pts):
            x0, y0 = pts[max(0, i - 1)]
            x1, y1 = pts[min(len(pts) - 1, i + 1)]
            dx, dy = x1 - x0, y1 - y0
            d = math.hypot(dx, dy) or 1
            nx, ny = -dy / d * half_width, dx / d * half_width
            left.append((x + nx, y + ny))
            right.append((x - nx, y - ny))
        return left + right[::-1]

    halo = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    ImageDraw.Draw(halo).polygon(ribbon(pipe, px(45)), fill=(72, 205, 190, 120))
    halo = halo.filter(ImageFilter.GaussianBlur(px(22)))
    bg.alpha_composite(halo)

    pd = ImageDraw.Draw(bg)
    pd.polygon(ribbon(pipe, px(38)), fill=(15, 54, 72, 255))
    pd.polygon(ribbon(pipe, px(31)), fill=(34, 112, 138, 255))
    pd.polygon(ribbon(pipe, px(24)), fill=(9, 38, 58, 255))
    for i in range(0, len(pipe) - 2, 2):
        f = i / (len(pipe) - 1)
        col = (int(62 + 10 * f) + brighten, int(191 + 14 * f) + brighten, int(111 + 121 * f) + brighten, 255)
        pd.line([pipe[i], pipe[i + 2]], fill=col, width=px(22))

    def glow_dot(cx, cy, r, color, alpha=255, blur=0):
        layer = Image.new("RGBA", (N, N), (0, 0, 0, 0))
        ImageDraw.Draw(layer).ellipse([cx-r, cy-r, cx+r, cy+r], fill=color + (alpha,))
        if blur:
            layer = layer.filter(ImageFilter.GaussianBlur(blur))
        bg.alpha_composite(layer)

    for idx, f in enumerate((0.16, 0.38, 0.60, 0.82)):
        x, y = pipe[int(f * (len(pipe) - 1))]
        r = px(15 if idx % 2 else 12)
        glow_dot(x, y, r, (230, 255, 245), 255, px(8))
        glow_dot(x, y, int(r * 0.55), (255, 255, 255), 255)

    # ---- coupling hub ------------------------------------------------------
    hx, hy = pipe[-1]
    glow_dot(hx, hy, px(78), (255, 159, 67), 150, px(26))
    hd = ImageDraw.Draw(bg)
    rr = px(52)
    hexpts = [
        (hx + rr * math.cos(math.radians(60*k - 30)),
         hy + rr * math.sin(math.radians(60*k - 30)))
        for k in range(6)
    ]
    hd.polygon(hexpts, fill=(18, 40, 62, 255), outline=(255, 183, 94, 255), width=px(5))
    hd.line([(hx - px(22), hy), (hx + px(22), hy)], fill=(255, 214, 140, 255), width=px(5))
    hd.polygon([(hx + px(28), hy), (hx + px(14), hy - px(10)), (hx + px(14), hy + px(10))],
               fill=(255, 214, 140, 255))
    hd.polygon([(hx - px(28), hy), (hx - px(14), hy - px(10)), (hx - px(14), hy + px(10))],
               fill=(255, 214, 140, 255))

    # ---- flow-assurance waveform ------------------------------------------
    wd = ImageDraw.Draw(bg)
    wave = [
        (px(x), px(868) + int(math.sin(x / 105) * px(9)) * (1 if 190 < x < 830 else 0))
        for x in range(180, 845, 6)
    ]
    wd.line(wave, fill=(78, 161, 255, 130), width=px(4), joint="curve")

    return bg.resize((size, size), Image.Resampling.LANCZOS)


def main():
    # master PNG
    icon = render(1024)
    icon_path = os.path.join(HERE, "icon.png")
    icon.save(icon_path, "PNG")
    print(f"wrote {icon_path} ({os.path.getsize(icon_path)} bytes)")

    # macOS iconset: small tiles get zoomed + brightened for legibility
    setdir = os.path.join(HERE, "chflow.iconset")
    if os.path.isdir(setdir):
        shutil.rmtree(setdir)
    os.makedirs(setdir)
    tiles = [
        (16, 1.40, 1.12, 18), (32, 1.28, 1.10, 12), (64, 1.15, 1.06, 6),
        (128, 1.06, 1.03, 2), (256, 1.0, 1.0, 0), (512, 1.0, 1.0, 0),
    ]
    for size, zoom, boost, brighten in tiles:
        for scale in (1, 2):
            s = size * scale
            name = f"icon_{size}x{size}" + (f"@{scale}x" if scale == 2 else "") + ".png"
            render(s, zoom if scale == 1 else max(1.0, zoom - 0.06),
                   boost if scale == 1 else max(1.0, boost - 0.03),
                   brighten).save(os.path.join(setdir, name), "PNG")
    icns = os.path.join(HERE, "chflow.icns")
    if shutil.which("iconutil"):
        subprocess.run(["iconutil", "-c", "icns", setdir, "-o", icns], check=True)
        print(f"wrote {icns} ({os.path.getsize(icns)} bytes)")


if __name__ == "__main__":
    main()
