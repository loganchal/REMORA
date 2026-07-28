"""Compare the REMORA V5 port harness output against the atlas Python
reference (stretching_vstretching5 / compute_zlevels, ported from zlevs_flex.m
and QC'd against the Moana hindcast grid)."""
import sys
import numpy as np

# Load only the two reference functions (module import has heavy side-effects)
import ast
src_path = "/Users/loganchalmers/Documents/MATLAB/atlas/data_tools/run_shdr_pilot.py"
tree = ast.parse(open(src_path).read())
wanted = [n for n in tree.body if isinstance(n, ast.FunctionDef)
          and n.name in ("stretching_vstretching5", "compute_zlevels")]
ns = {"np": np}
exec(compile(ast.Module(body=wanted, type_ignores=[]), src_path, "exec"), ns)
stretching_vstretching5 = ns["stretching_vstretching5"]
compute_zlevels = ns["compute_zlevels"]

N, theta_s, theta_b, hc = 50, 6.0, 2.0, 250.0

# parse harness output
sw = np.zeros(N + 1); cw = np.zeros(N + 1)
sr = np.zeros(N); cr = np.zeros(N)
zr = {}; zw = {}
for line in open(sys.argv[1]):
    p = line.split()
    if p[0] == "W":
        sw[int(p[1])], cw[int(p[1])] = float(p[2]), float(p[3])
    elif p[0] == "R":
        sr[int(p[1])], cr[int(p[1])] = float(p[2]), float(p[3])
    elif p[0] in ("ZR", "ZW"):
        zr.setdefault((float(p[2]), float(p[3]), p[0]), {})[int(p[1])] = float(p[4])

s_r_ref, C_r_ref = stretching_vstretching5(theta_s, theta_b, hc, N, kgrid=0)
s_w_ref, C_w_ref = stretching_vstretching5(theta_s, theta_b, hc, N, kgrid=1)

def rel(a, b):
    return np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-30))

print(f"s_r  max rel diff: {rel(sr, s_r_ref):.3e}")
print(f"C_r  max rel diff: {rel(cr, C_r_ref):.3e}")
print(f"s_w  max rel diff: {rel(sw, s_w_ref):.3e}")
print(f"C_w  max rel diff: {rel(cw, C_w_ref):.3e}")

worst = 0.0
for (h, zeta, kind), zk in zr.items():
    grid_h = np.full((1, 1), h)
    grid_z = np.full((1, 1), zeta)
    kg = 0 if kind == "ZR" else 1
    ref = compute_zlevels(grid_h, grid_z, theta_s, theta_b, hc, N, kgrid=kg)[:, 0, 0]
    got = np.array([zk[k] for k in sorted(zk)])
    d = np.max(np.abs(got - ref) / np.maximum(np.abs(ref), 1e-12))
    worst = max(worst, d)
    print(f"{kind} h={h:6g} zeta={zeta:5g}: max rel diff {d:.3e}")
print(f"WORST depth rel diff: {worst:.3e}")
print("PASS" if worst < 1e-13 and rel(cw, C_w_ref) < 1e-13 else "FAIL")
