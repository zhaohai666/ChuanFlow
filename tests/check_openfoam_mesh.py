#!/usr/bin/env python3
"""OpenFOAM ASCII polyMesh 导入的独立验收脚本（纯标准库，不信任被测内核）。

用法：
    python3 check_openfoam_mesh.py <openfoam-mesh 可执行文件> <case.dict>

脚本自己做三件事：
  1. 用独立的 Python 解析器读取同一批 points/faces/owner/neighbour/boundary，
     重算拓扑、patch、面积矢量、体积（散度定理 Σ_f Sf·Cf/3）与包围盒；
  2. 运行被测 app，逐项核对它打印的数字；
  3. 跑 gated twin（错误 format / 缺失 mesh 必须失败，大小写容错必须仍成功）
     与三个拒收探针（坏 faces / owner 越界 / 缺 boundary 必须失败）。

任何一项不满足即以非零退出码结束，并打印失败原因。
"""

import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"

fails = []
notes = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


def close(a, b, rel=1e-9, abs_=1e-12):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


# ---------------------------------------------------------------------------
# 独立的 OpenFOAM ASCII 列表解析器
# ---------------------------------------------------------------------------

def strip_comments(text):
    out = []
    i = 0
    while i < len(text):
        if text[i:i + 2] == "//":
            j = text.find("\n", i)
            if j < 0:
                break
            i = j
        elif text[i:i + 2] == "/*":
            j = text.find("*/", i + 2)
            if j < 0:
                raise RuntimeError("未闭合的 /* 注释")
            i = j + 2
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def tokenize(text):
    """把文本切成单词/括号/花括号/分号，丢掉空白与注释。"""
    text = strip_comments(text)
    toks = []
    i = 0
    while i < len(text):
        c = text[i]
        if c.isspace():
            i += 1
        elif c in "(){};":
            toks.append(c)
            i += 1
        else:
            j = i
            while j < len(text) and not text[j].isspace() and text[j] not in "(){};":
                j += 1
            toks.append(text[i:j])
            i = j
    return toks


def matching(toks, open_idx, op, cl):
    depth = 0
    for j in range(open_idx, len(toks)):
        if toks[j] == op:
            depth += 1
        elif toks[j] == cl:
            depth -= 1
            if depth == 0:
                return j
    raise RuntimeError("括号不匹配")


def find_list(toks, label):
    """定位文件最外层的 `N ( ... )`，跳过 FoamFile 头里的花括号内容。"""
    depth = 0
    for i in range(len(toks) - 1):
        t = toks[i]
        if t == "{":
            depth += 1
        elif t == "}":
            depth -= 1
        elif depth == 0 and re.fullmatch(r"[+-]?\d+", t) and toks[i + 1] == "(":
            close = matching(toks, i + 1, "(", ")")
            return int(t), toks[i + 2:close]
    raise RuntimeError(label + " 中没有 N ( ... ) 列表")


def read_list(path, label):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return find_list(tokenize(fh.read()), label)


class Cursor:
    def __init__(self, toks, where):
        self.toks = toks
        self.i = 0
        self.where = where

    def take(self):
        if self.i >= len(self.toks):
            raise RuntimeError(self.where + " 列表被截断")
        t = self.toks[self.i]
        self.i += 1
        return t

    def expect(self, tok):
        got = self.take()
        if got != tok:
            raise RuntimeError(self.where + " 期望 " + tok + "，得到 " + got)

    def integer(self):
        t = self.take()
        if not re.fullmatch(r"[+-]?\d+", t):
            raise RuntimeError(self.where + " 期望整数，得到 " + t)
        return int(t)

    def number(self):
        t = self.take()
        if not re.fullmatch(NUM, t):
            raise RuntimeError(self.where + " 期望数字，得到 " + t)
        return float(t)

    def done(self):
        return self.i >= len(self.toks)


def unquote(s):
    if len(s) >= 2 and s[0] in "\"'" and s[-1] == s[0]:
        return s[1:-1]
    return s


def parse_points(path):
    n, toks = read_list(path, "points")
    c = Cursor(toks, "points")
    pts = []
    for _ in range(n):
        c.expect("(")
        p = (c.number(), c.number(), c.number())
        c.expect(")")
        pts.append(p)
    if not c.done():
        raise RuntimeError("points 有多余条目")
    return pts


def parse_faces(path):
    n, toks = read_list(path, "faces")
    c = Cursor(toks, "faces")
    faces = []
    for _ in range(n):
        count = c.integer()
        if count < 3:
            raise RuntimeError("面顶点数小于 3")
        c.expect("(")
        face = [c.integer() for _ in range(count)]
        c.expect(")")
        faces.append(face)
    if not c.done():
        raise RuntimeError("faces 有多余条目")
    return faces


def parse_labels(path, what):
    n, toks = read_list(path, what)
    c = Cursor(toks, what)
    vals = [c.integer() for _ in range(n)]
    if not c.done():
        raise RuntimeError(what + " 有多余条目")
    return vals


def parse_boundary(path):
    n, toks = read_list(path, "boundary")
    c = Cursor(toks, "boundary")
    patches = []
    for _ in range(n):
        name = unquote(c.take())
        c.expect("{")
        # 逐 token 宽松扫描，只认三个关键字，忽略 inGroups 等附加项
        block = []
        depth = 0
        k = c.i
        while k < len(toks):
            t = toks[k]
            if t == "}" and depth == 0:
                break
            if t == "{":
                depth += 1
            elif t == "}":
                depth -= 1
            block.append(t)
            k += 1
        info = {"type": None, "nfaces": None, "startface": None}
        for j, t in enumerate(block[:-1]):
            key = t.lower()
            nxt = block[j + 1]
            if key == "type":
                info[key] = unquote(nxt)
            elif key in ("nfaces", "startface") and re.fullmatch(r"[+-]?\d+", nxt):
                info[key] = int(nxt)
        if any(v is None for v in info.values()):
            raise RuntimeError("patch " + name + " 缺少 type/nFaces/startFace")
        patches.append({"name": name, "type": info["type"],
                        "nFaces": info["nfaces"], "startFace": info["startface"]})
        c.i = k + 1
    if not c.done():
        raise RuntimeError("boundary 有多余条目")
    return patches


def map_type(raw):
    key = unquote(raw).lower()
    if key in ("wall", "67", "patch"):
        return "wall"
    if key in ("empty", "symmetryplane", "symmetry"):
        return "symmetry"
    return unquote(raw)


# ---------------------------------------------------------------------------
# 独立几何重算：不调用被测程序的任何代码
# ---------------------------------------------------------------------------

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def mag(a):
    return math.sqrt(dot(a, a))


def recompute(meshdir, overrides=None):
    """返回独立重算得到的期望值字典。"""
    pts = parse_points(os.path.join(meshdir, "points"))
    faces = parse_faces(os.path.join(meshdir, "faces"))
    owner = parse_labels(os.path.join(meshdir, "owner"), "owner")
    neighbour = parse_labels(os.path.join(meshdir, "neighbour"), "neighbour")
    patches = parse_boundary(os.path.join(meshdir, "boundary"))

    n_points = len(pts)
    n_faces = len(faces)
    n_internal = len(neighbour)
    n_cells = max(max(owner, default=0), max(neighbour, default=0)) + 1

    # 面心与原始面积矢量（扇形三角剖分）
    cf = []
    sf = []
    for face in faces:
        if any(v < 0 or v >= n_points for v in face):
            raise RuntimeError("面引用了不存在的点")
        c = (0.0, 0.0, 0.0)
        for v in face:
            c = add(c, pts[v])
        c = scale(c, 1.0 / len(face))
        cf.append(c)

        v0 = pts[face[0]]
        n = (0.0, 0.0, 0.0)
        for k in range(1, len(face) - 1):
            n = add(n, cross(sub(pts[face[k]], v0), sub(pts[face[k + 1]], v0)))
        sf.append(scale(n, 0.5))

    # 单元中心的一阶估计（与内核相同思路，但在 Python 里独立实现）
    est = [(0.0, 0.0, 0.0) for _ in range(n_cells)]
    cnt = [0] * n_cells
    for f in range(n_faces):
        o = owner[f]
        est[o] = add(est[o], cf[f])
        cnt[o] += 1
        if f < n_internal:
            nb = neighbour[f]
            est[nb] = add(est[nb], cf[f])
            cnt[nb] += 1
    for c in range(n_cells):
        if cnt[c]:
            est[c] = scale(est[c], 1.0 / cnt[c])

    # 定向：内部面 owner->neighbour，边界面朝外
    for f in range(n_faces):
        o = owner[f]
        if f < n_internal:
            reference = sub(est[neighbour[f]], est[o])
        else:
            reference = sub(cf[f], est[o])
        if dot(sf[f], reference) < 0:
            sf[f] = scale(sf[f], -1.0)

    areas = [mag(s) for s in sf]

    external_volume = sum(dot(sf[f], cf[f]) for f in range(n_internal, n_faces)) / 3.0
    external_area = sum(areas[f] for f in range(n_internal, n_faces))

    cell_volumes = [0.0] * n_cells
    for f in range(n_faces):
        contrib = dot(sf[f], cf[f]) / 3.0
        cell_volumes[owner[f]] += contrib
        if f < n_internal:
            cell_volumes[neighbour[f]] -= contrib
    cell_volumes = [abs(v) for v in cell_volumes]

    lo = tuple(min(p[d] for p in pts) for d in range(3))
    hi = tuple(max(p[d] for p in pts) for d in range(3))

    out_patches = []
    for p in patches:
        t = p["type"]
        if overrides and p["name"] in overrides:
            t = overrides[p["name"]]
        out_patches.append({"name": p["name"], "type": map_type(t),
                            "start": p["startFace"], "size": p["nFaces"]})

    return {
        "nPoints": n_points, "nCells": n_cells,
        "nFaces": n_faces, "nInternalFaces": n_internal,
        "nBoundaryFaces": n_faces - n_internal,
        "patches": out_patches,
        "cellVolumes": cell_volumes,
        "cellVolumeSum": sum(cell_volumes),
        "externalVolume": external_volume,
        "externalArea": external_area,
        "areas": areas,
        "min": lo, "max": hi,
    }


# ---------------------------------------------------------------------------
# app 输出解析
# ---------------------------------------------------------------------------

PREFIX = r"(?:\[chflow\]\s*)?openfoam-mesh:\s*"


def one_int(text, key):
    m = re.search(PREFIX + re.escape(key) + r"\s+(-?\d+)\s*;", text)
    if not m:
        raise RuntimeError("app 输出缺少 " + key)
    return int(m.group(1))


def one_float(text, key):
    m = re.search(PREFIX + re.escape(key) + r"\s+(" + NUM + r")\s*;", text)
    if not m:
        raise RuntimeError("app 输出缺少 " + key)
    return float(m.group(1))


def one_vec(text, key):
    m = re.search(PREFIX + re.escape(key) + r"\s*\(?\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)?\s*;", text)
    if not m:
        raise RuntimeError("app 输出缺少 " + key)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def all_patches(text):
    pat = PREFIX + r"patch\s+(\S+)\s+type\s+(\S+)\s+start\s+(\d+)\s+size\s+(\d+)\s*;"
    return [(m.group(1), m.group(2), int(m.group(3)), int(m.group(4)))
            for m in re.finditer(pat, text)]


def indexed_values(text, key):
    pat = PREFIX + re.escape(key) + r"\s+(\d+)\s+(" + NUM + r")\s*;"
    return {int(m.group(1)): float(m.group(2)) for m in re.finditer(pat, text)}


def run(exe, case):
    proc = subprocess.run([exe, case], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=60)
    return proc.returncode, proc.stdout


def meshdir_from_case(case):
    """从 case.dict 里读出 directory（仅用于定位被测文件），不做几何。"""
    with open(case, "r", encoding="utf-8") as fh:
        text = strip_comments(fh.read())
    m = re.search(r"\bdirectory\s+([^\s;]+)\s*;", text)
    if not m:
        m = re.search(r"\bfile\s+([^\s;]+)\s*;", text)
    if not m:
        raise RuntimeError("case.dict 没有 directory/file")
    path = m.group(1).strip("\"'")
    if not os.path.isabs(path):
        path = os.path.join(os.path.dirname(os.path.abspath(case)), path)
    return os.path.abspath(path)


def compare_app(exe, case, expected, tag):
    rc, out = run(exe, case)
    if rc != 0:
        check(False, tag + ": app 返回 " + str(rc) + "\n" + out)
        return None
    for key in ("nPoints", "nCells", "nFaces", "nInternalFaces", "nBoundaryFaces",
                "cellVolumeSum"):
        got = one_int(out, key) if key != "cellVolumeSum" else one_float(out, key)
        want = expected[key]
        check((got == want) if isinstance(want, int) else close(got, want),
              tag + ": " + key + " 应为 " + str(want) + "，实为 " + str(got))

    got_patches = all_patches(out)
    want_patches = [(p["name"], p["type"], p["start"], p["size"])
                    for p in expected["patches"]]
    check(got_patches == want_patches,
          tag + ": patch 列表不一致\n  期望 " + str(want_patches) + "\n  实为 " + str(got_patches))

    got_lo, got_hi = one_vec(out, "minX"), one_vec(out, "maxX")
    check(all(close(a, b) for a, b in zip(got_lo, expected["min"])),
          tag + ": minX 不一致 " + str(got_lo) + " vs " + str(expected["min"]))
    check(all(close(a, b) for a, b in zip(got_hi, expected["max"])),
          tag + ": maxX 不一致 " + str(got_hi) + " vs " + str(expected["max"]))

    cell_vols = indexed_values(out, "cellVolume")
    # 注意 cellVolumeSum 也匹配 cellVolume 的前缀？正则要求整数下标，故不会误匹配。
    check(len(cell_vols) == expected["nCells"],
          tag + ": cellVolume 条数应为 " + str(expected["nCells"]))
    for c, want in enumerate(expected["cellVolumes"]):
        check(c in cell_vols and close(cell_vols[c], want),
              tag + ": cellVolume[" + str(c) + "] 应为 " + str(want)
              + "，实为 " + str(cell_vols.get(c)))

    face_areas = indexed_values(out, "faceArea")
    check(len(face_areas) == expected["nFaces"],
          tag + ": faceArea 条数应为 " + str(expected["nFaces"]))
    for f, want in enumerate(expected["areas"]):
        check(f in face_areas and close(face_areas[f], want),
              tag + ": faceArea[" + str(f) + "] 应为 " + str(want)
              + "，实为 " + str(face_areas.get(f)))
    return out


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print("usage: check_openfoam_mesh.py <openfoam-mesh exe> <case.dict>",
              file=sys.stderr)
        return 2
    exe, case = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    case_dir = os.path.dirname(case)
    meshdir = meshdir_from_case(case)

    # 1) 独立重算 + 与 app 逐项核对
    expected = recompute(meshdir)
    compare_app(exe, case, expected, "主用例")

    # 闭式几何锚：与解析实现和被测程序都无关的手算常数
    check(expected["nCells"] == 2, "闭式锚: nCells 应为 2")
    check(expected["nFaces"] == 11, "闭式锚: nFaces 应为 11")
    check(expected["nInternalFaces"] == 1, "闭式锚: nInternalFaces 应为 1")
    check(close(expected["cellVolumeSum"], 2.0), "闭式锚: 总体积应为 2")
    check(close(expected["externalVolume"], 2.0), "闭式锚: 散度定理体积应为 2")
    check(close(expected["externalArea"], 10.0), "闭式锚: 总外表面积应为 10")
    check(all(close(a, 1.0) for a in expected["areas"]), "闭式锚: 每个面面积应为 1")
    check(expected["min"] == (0.0, 0.0, 0.0), "闭式锚: min 应为 (0,0,0)")
    check(expected["max"] == (2.0, 1.0, 1.0), "闭式锚: max 应为 (2,1,1)")
    want_types = {"left": "wall", "right": "wall", "bottomY": "wall",
                  "topY": "wall", "frontZ": "symmetry", "backZ": "symmetry"}
    for p in expected["patches"]:
        check(p["type"] == want_types.get(p["name"]),
              "闭式锚: patch " + p["name"] + " 类型映射错误: " + p["type"])

    with tempfile.TemporaryDirectory(prefix="cf_openfoam_check_") as td:
        root = os.path.join(td, "case")
        shutil.copytree(case_dir, root)

        # 2) gated twin A：format 换成 gmsh 必须失败
        bad_format = os.path.join(td, "bad_format.dict")
        with open(bad_format, "w", encoding="utf-8") as fh:
            fh.write("mesh { format gmsh; directory \"" + meshdir + "\"; }\n")
        rc, out = run(exe, bad_format)
        check(rc != 0, "gated twin: format=gmsh 必须失败，实际 rc=" + str(rc))

        # 3) gated twin B：删掉 mesh 块必须失败
        no_mesh = os.path.join(td, "no_mesh.dict")
        with open(no_mesh, "w", encoding="utf-8") as fh:
            fh.write("// no mesh block at all\n")
        rc, out = run(exe, no_mesh)
        check(rc != 0, "gated twin: 没有 mesh 块必须失败，实际 rc=" + str(rc))

        # 4) 正向容错：format 值大写 + polyMesh 文件插入多余空白，结果须一致
        robust_dir = os.path.join(root, "polyMeshRobust")
        shutil.copytree(os.path.join(root, "polyMesh"), robust_dir)
        for name in ("points", "faces", "owner", "neighbour", "boundary"):
            fp = os.path.join(robust_dir, name)
            with open(fp, "r", encoding="utf-8") as fh:
                body = fh.read()
            with open(fp, "w", encoding="utf-8") as fh:
                fh.write(body.replace("\n", "\n\t\n"))
        robust_case = os.path.join(td, "robust.dict")
        with open(robust_case, "w", encoding="utf-8") as fh:
            fh.write("mesh { FORMAT OpenFOAM; directory \"" + robust_dir + "\"; }\n")
        rc, out = run(exe, robust_case)
        check(rc == 0, "正向容错: 大写 format + 多余空白应成功，实际 rc=" + str(rc) + "\n" + out)
        if rc == 0:
            check(one_int(out, "nFaces") == expected["nFaces"]
                  and close(one_float(out, "cellVolumeSum"), expected["cellVolumeSum"]),
                  "正向容错: 关键统计应与主用例一致")

        # 5) 拒收探针：坏 faces（顶点数 2）
        probe_dir = os.path.join(td, "bad_faces")
        shutil.copytree(root, probe_dir)
        fp = os.path.join(probe_dir, "polyMesh", "faces")
        with open(fp, "r", encoding="utf-8") as fh:
            body = fh.read()
        check("4(4 5 7 6)" in body, "探针准备: 找不到预期的 faces 行")
        with open(fp, "w", encoding="utf-8") as fh:
            fh.write(body.replace("4(4 5 7 6)", "2(4 5)", 1))
        probe_case = os.path.join(td, "bad_faces.dict")
        with open(probe_case, "w", encoding="utf-8") as fh:
            fh.write("mesh { format openfoam; directory \""
                     + os.path.join(probe_dir, "polyMesh") + "\"; }\n")
        rc, out = run(exe, probe_case)
        check(rc != 0, "拒收探针: 坏 faces（顶点数 2）必须失败，实际 rc=" + str(rc))

        # 6) 拒收探针：owner 越界（第一个 owner 改成 999）
        probe_dir = os.path.join(td, "bad_owner")
        shutil.copytree(root, probe_dir)
        fp = os.path.join(probe_dir, "polyMesh", "owner")
        with open(fp, "r", encoding="utf-8") as fh:
            body = fh.read()
        m = re.search(r"(\b11\s*\(\s*)", body)
        check(m is not None, "探针准备: 找不到 owner 列表开头")
        if m:
            pos = m.end()
            idx = body.find("0", pos)
            check(idx >= 0, "探针准备: 找不到第一个 owner 值")
            if idx >= 0:
                body = body[:idx] + "999" + body[idx + 1:]
                with open(fp, "w", encoding="utf-8") as fh:
                    fh.write(body)
        probe_case = os.path.join(td, "bad_owner.dict")
        with open(probe_case, "w", encoding="utf-8") as fh:
            fh.write("mesh { format openfoam; directory \""
                     + os.path.join(probe_dir, "polyMesh") + "\"; }\n")
        rc, out = run(exe, probe_case)
        check(rc != 0, "拒收探针: owner 越界必须失败，实际 rc=" + str(rc))

        # 7) 拒收探针：缺 boundary 文件
        probe_dir = os.path.join(td, "no_boundary")
        shutil.copytree(root, probe_dir)
        os.remove(os.path.join(probe_dir, "polyMesh", "boundary"))
        probe_case = os.path.join(td, "no_boundary.dict")
        with open(probe_case, "w", encoding="utf-8") as fh:
            fh.write("mesh { format openfoam; directory \""
                     + os.path.join(probe_dir, "polyMesh") + "\"; }\n")
        rc, out = run(exe, probe_case)
        check(rc != 0, "拒收探针: 缺 boundary 必须失败，实际 rc=" + str(rc))

    if fails:
        print("check_openfoam_mesh: FAIL (" + str(len(fails)) + " 项)")
        for f in fails:
            print("  - " + f)
        return 1
    print("check_openfoam_mesh: PASS (独立重算/闭式锚/gated twin/拒收探针全部通过)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
