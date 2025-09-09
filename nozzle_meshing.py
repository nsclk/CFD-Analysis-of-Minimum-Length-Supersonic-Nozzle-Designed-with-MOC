import gmsh
import numpy as np

TABLE = """\
Point_Number,X,Y
1,-3.0000000000,1.0000000000
2,-2.9387755102,0.9999773219
3,-2.8775510204,0.9998185859
4,-2.8163265306,0.9993878263
5,-2.7551020408,0.9985493778
6,-2.6938775510,0.9971682226
7,-2.6326530612,0.9951104740
8,-2.5714285714,0.9922439974
9,-2.5102040816,0.9884391656
10,-2.4489795918,0.9835697453
11,-2.3877551020,0.9775139105
12,-2.3265306122,0.9701553776
13,-2.2653061224,0.9613846550
14,-2.2040816327,0.9511003979
15,-2.1428571429,0.9392108564
16,-2.0816326531,0.9256354064
17,-2.0204081633,0.9103061449
18,-1.9591836735,0.8931695328
19,-1.8979591837,0.8741880649
20,-1.8367346939,0.8533419420
21,-1.7755102041,0.8306307189
22,-1.7142857143,0.8060748979
23,-1.6530612245,0.7797174328
24,-1.5918367347,0.7516251076
25,-1.5306122449,0.7218897461
26,-1.4693877551,0.6906292078
27,-1.4081632653,0.6579881195
28,-1.3469387755,0.6241382883
29,-1.2857142857,0.5892787353
30,-1.2244897959,0.5536352869
31,-1.1632653061,0.5174596529
32,-1.1020408163,0.4810279172
33,-1.0408163265,0.4446383603
34,-0.9795918367,0.4086085269
35,-0.9183673469,0.3732714474
36,-0.8571428571,0.3389709131
37,-0.7959183673,0.3060557028
38,-0.7346938776,0.2748726476
39,-0.6734693878,0.2457584161
40,-0.6122448980,0.2190298965
41,-0.5510204082,0.1949730420
42,-0.4897959184,0.1738300402
43,-0.4285714286,0.1557846609
44,-0.3673469388,0.1409456259
45,-0.3061224490,0.1293278394
46,-0.2448979592,0.1208313081
47,-0.1836734694,0.1152175709
48,-0.1224489796,0.1120834526
49,-0.0612244898,0.1108319434
50,0.0000000000,0.1106400000
51,0.0000000000,0.1106400000
52,0.1524513002,0.1828869208
53,0.1936315611,0.2019578210
54,0.2180445661,0.2128724224
55,0.2396190682,0.2222897536
56,0.2598187301,0.2308949236
57,0.2793203144,0.2389995349
58,0.2985010091,0.2467724227
59,0.3175982740,0.2543154098
60,0.3367755452,0.2616944569
61,0.3561534849,0.2689544719
62,0.3758264655,0.2761271086
63,0.3958719709,0.2832351987
64,0.4163562997,0.2902954228
65,0.4373382033,0.2973199969
66,0.4588713083,0.3043177768
67,0.4810057905,0.3112950030
68,0.5037895713,0.3182558143
69,0.5272691999,0.3252026082
70,0.5514905233,0.3321362963
71,0.5764992103,0.3390564852
72,0.6023411726,0.3459616044
73,0.6290629136,0.3528489948
74,0.6567118250,0.3597149668
75,0.6853364455,0.3665548356
76,0.7149866938,0.3733629372
77,0.7457140818,0.3801326299
78,0.7775719156,0.3868562826
79,0.8106154876,0.3935252518
80,0.8449022643,0.4001298491
81,0.8804920723,0.4066592993
82,0.9174472841,0.4131016898
83,0.9558330076,0.4194439124
84,0.9957172792,0.4256715964
85,1.0371712624,0.4317690338
86,1.0802694541,0.4377190966
87,1.1250898995,0.4435031450
88,1.1717144164,0.4491009275
89,1.2202288308,0.4544904714
90,1.2707232243,0.4596479636
91,1.3232921954,0.4645476214
92,1.3780351354,0.4691615517
93,1.4350565197,0.4734595996
94,1.4944662175,0.4774091831
95,1.5563798190,0.4809751150
96,1.6209189844,0.4841194107
97,1.6882118129,0.4868010795
98,1.7583932371,0.4889758993
99,1.8316054413,0.4905961733
100,1.9079983082,0.4916104669
101,1.9877298939,0.4919633228
"""

# ------------------------------
# Manual control (minimal patch)
# ------------------------------
MANUAL_SPLITS = False        # True => use manual x-splits / total block count
NBLK_TOTAL    = None         # e.g. 6  (used only if MANUAL_SPLITS=True and X_SPLITS is empty)
X_SPLITS      = []           # interior x-positions; e.g. [-1.5, 0.0, 0.54, 1.067, 1.602]

# Per-block point-count overrides (optional). Keys are block indices: 0..(num_blocks-1)
NX_OVERRIDE   = {0: 45, 1: 45, 2: 45, 3: 60, 4: 60, 5: 45, 6: 45, 7: 45}
NY_OVERRIDE   = {0: 45,  1: 45,  2: 45,  3: 45,  4: 45,  5: 45, 6: 45, 7: 45}

# ------------------------------
# Original globals (unchanged)
# ------------------------------
NY_GLOBAL  = 90     # across height (shared if not overridden)
NBLK_LEFT  = 4      # blocks before throat (auto mode)
NBLK_RIGHT = 4      # blocks after throat (auto mode)
SHOW_GUI   = True
OUTFILE    = "nozzle_xy_conformal.su2"

def parse_xy(table: str) -> np.ndarray:
    """
    Parse either whitespace-separated or CSV lines.
    Accepts rows like: 'idx, x, y' or 'x, y' or 'idx x y'.
    Skips headers automatically. Removes consecutive duplicates.
    """
    pts = []
    for line in table.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        # choose separator
        sep = ',' if ',' in line else None
        parts = [p.strip() for p in line.split(sep) if p.strip() != '']
        # need at least x,y
        if len(parts) < 2:
            continue
        # try last two fields as x,y
        try:
            x = float(parts[-2])
            y = float(parts[-1])
        except ValueError:
            # header or bad row
            continue
        pts.append((x, y))
    if not pts:
        raise ValueError("No valid (x,y) rows parsed from TABLE.")
    arr = np.array(pts, dtype=float)
    # drop consecutive duplicates
    out = [arr[0]]
    for q in arr[1:]:
        if not np.allclose(q, out[-1], atol=1e-14):
            out.append(q)
    return np.array(out)

def find_throat_idx(xy): 
    return int(np.argmin(xy[:,1]))

def indices_by_equal_x(x, i0, i1, nseg):
    xa, xb = x[i0], x[i1]
    cuts = np.linspace(xa, xb, nseg+1)
    idxs = [i0]
    for c in cuts[1:-1]:
        j_local = int(np.argmin(np.abs(x[i0:i1+1] - c)))
        j = i0 + j_local
        j = max(j, idxs[-1] + 1)
        j = min(j, i1 - 1)
        idxs.append(j)
    idxs.append(i1)
    return idxs

def blocks_by_custom_splits(x, xcuts):
    """xcuts: interior x positions (sorted). Returns [(i0,i1), ...]."""
    idxs = [0]
    for c in sorted(xcuts):
        j = int(np.argmin(np.abs(x - c)))
        j = max(j, idxs[-1] + 1)
        j = min(j, len(x) - 2)
        idxs.append(j)
    idxs.append(len(x) - 1)
    return [(idxs[k], idxs[k+1]) for k in range(len(idxs)-1)]

def block_width_and_yavg(xy, i0, i1):
    width = abs(xy[i1,0] - xy[i0,0])
    yavg  = float(np.mean(xy[i0:i1+1,1]))
    return width, max(yavg, 1e-9)

def main():
    xy = parse_xy(TABLE)
    x  = xy[:,0]; y = xy[:,1]
    i_th = find_throat_idx(xy)

    # -----------------------------
    # Block selection (manual/auto)
    # -----------------------------
    if MANUAL_SPLITS:
        if X_SPLITS:
            blocks = blocks_by_custom_splits(x, X_SPLITS)
        elif isinstance(NBLK_TOTAL, int) and NBLK_TOTAL >= 1:
            idxs = indices_by_equal_x(x, 0, len(xy)-1, NBLK_TOTAL)
            blocks = [(idxs[k], idxs[k+1]) for k in range(len(idxs)-1)]
        else:
            left_idxs  = indices_by_equal_x(x, 0,    i_th,      NBLK_LEFT)
            right_idxs = indices_by_equal_x(x, i_th, len(xy)-1, NBLK_RIGHT)
            blocks = [(left_idxs[k], left_idxs[k+1]) for k in range(len(left_idxs)-1)]
            blocks += [(right_idxs[k], right_idxs[k+1]) for k in range(len(right_idxs)-1)]
    else:
        left_idxs  = indices_by_equal_x(x, 0,    i_th,      NBLK_LEFT)
        right_idxs = indices_by_equal_x(x, i_th, len(xy)-1, NBLK_RIGHT)
        blocks = [(left_idxs[k], left_idxs[k+1]) for k in range(len(left_idxs)-1)]
        blocks += [(right_idxs[k], right_idxs[k+1]) for k in range(len(right_idxs)-1)]

    # per-block Nx (auto suggestion for near-square cells)
    NY = max(3, int(NY_GLOBAL))
    Nx_list = []
    for (i0,i1) in blocks:
        w, yavg = block_width_and_yavg(xy, i0, i1)
        nx = int(round(w * (NY-1)/yavg)) + 1
        Nx_list.append(max(4, nx))

    gmsh.initialize()
    gmsh.model.add("nozzle_xy_conformal")
    geo  = gmsh.model.geo
    mesh = gmsh.model.mesh

    # points
    h_len = 0.02 * (x.max() - x.min() if x.max()>x.min() else 1.0)
    wall_pts = []
    tag = 1000
    for xi, yi in xy:
        geo.addPoint(float(xi), float(yi), 0.0, h_len, tag)
        wall_pts.append(tag); tag += 1

    # cache axis bottom points and vertical interface lines so they are REUSED
    axis_pts = {}           # key: rounded x -> point tag
    vlines   = {}           # key: wall index i -> line tag (bottom->top)

    def get_axis_pt(xcoord):
        key = round(float(xcoord), 12)
        nonlocal tag
        if key not in axis_pts:
            geo.addPoint(float(xcoord), 0.0, 0.0, h_len, tag)
            axis_pts[key] = tag; tag += 1
        return axis_pts[key]

    def get_vertical_line(i):
        """Create ONCE and reuse: bottom->top line at wall index i."""
        nonlocal tag
        if i not in vlines:
            bot = get_axis_pt(xy[i,0])
            top = wall_pts[i]
            vlines[i] = geo.addLine(bot, top)  # bottom->top
        return vlines[i]

    # walls per block
    def add_wall_spline(i0,i1):
        return geo.addSpline(wall_pts[i0:i1+1])

    surfaces = []
    wall_curves = []
    axis_segments = []
    inlet_line = outlet_line = None

    for b,(i0,i1) in enumerate(blocks):
        vL = get_vertical_line(i0)          # bottom->top
        vR = get_vertical_line(i1)          # bottom->top
        wC = add_wall_spline(i0,i1)
        wall_curves.append(wC)

        # axis segment (right->left for CCW)
        aR = get_axis_pt(xy[i1,0])
        aL = get_axis_pt(xy[i0,0])
        aX = geo.addLine(aR, aL)
        axis_segments.append(aX)

        if b == 0: inlet_line  = vL
        if b == len(blocks)-1: outlet_line = vR

        loop = geo.addCurveLoop([ +vL, +wC, -vR, +aX ])
        srf  = geo.addPlaneSurface([loop])
        surfaces.append((srf, vL, wC, vR, aX))

    geo.synchronize()

    # transfinite & recombine
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Smoothing", 10)
    gmsh.option.setNumber("Mesh.Algorithm", 5)

    for (idx,(srf, vL, wC, vR, aX)) in enumerate(surfaces):
        # apply overrides if provided, otherwise use auto Nx and global NY
        nx = int(NX_OVERRIDE.get(idx, Nx_list[idx]))
        ny = int(NY_OVERRIDE.get(idx, NY))

        mesh.setTransfiniteCurve(wC, nx, meshType="Progression", coef=1.0)
        mesh.setTransfiniteCurve(aX, nx, meshType="Progression", coef=1.0)
        mesh.setTransfiniteCurve(vL, ny, meshType="Progression", coef=1.0)
        mesh.setTransfiniteCurve(vR, ny, meshType="Progression", coef=1.0)
        mesh.setTransfiniteSurface(srf)

    # physical groups: ONLY outer boundaries
    dim = 1
    gmsh.model.addPhysicalGroup(dim, wall_curves, 1001); gmsh.model.setPhysicalName(dim, 1001, "WALL")
    gmsh.model.addPhysicalGroup(dim, [inlet_line], 1002); gmsh.model.setPhysicalName(dim, 1002, "INLET")
    gmsh.model.addPhysicalGroup(dim, [outlet_line],1003); gmsh.model.setPhysicalName(dim, 1003, "OUTLET")
    gmsh.model.addPhysicalGroup(dim, axis_segments,1004); gmsh.model.setPhysicalName(dim, 1004, "SYMMETRY")

    gmsh.model.addPhysicalGroup(2, [s for (s,_,_,_,_) in surfaces], 2001)
    gmsh.model.setPhysicalName(2, 2001, "DOMAIN")

    mesh.generate(2)
    gmsh.write(OUTFILE)

    # sanity prints
    xl = xy[blocks[0][0],0]; xr = xy[blocks[-1][1],0]
    print(f"[INFO] INLET at x≈{xl:.6f}, OUTLET at x≈{xr:.6f}")
    print(f"[INFO] blocks={len(blocks)}, Nx/NY per block:")
    for i,(i0,i1) in enumerate(blocks):
        print(f"  - block {i}: x∈[{x[i0]:.6f},{x[i1]:.6f}]  Nx={int(NX_OVERRIDE.get(i, Nx_list[i]))}  Ny={int(NY_OVERRIDE.get(i, NY))}")
    print(f"[INFO] Wrote: {OUTFILE}")

    if SHOW_GUI: gmsh.fltk.run()
    gmsh.finalize()

if __name__ == "__main__":
    main()
