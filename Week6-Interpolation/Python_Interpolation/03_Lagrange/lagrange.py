"""
Nội suy Lagrange (Lagrange Interpolation)
Nhập: input.txt hoặc file đề thi Midterm (.txt, .csv, .xlsx)
Xuất: output.txt
"""
import os, sys, argparse
import numpy as np
import pandas as pd

try:
    from pps_io import parse_universal_input, select_window_for_method
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pps_io import parse_universal_input, select_window_for_method

def lagrange(points):
    x = [p[0] for p in points]; y = [p[1] for p in points]; n = len(points)
    L_polys, P_coeffs = [], np.zeros(n, dtype=float)
    step_data = []
    for i in range(n):
        roots = [x[j] for j in range(n) if j != i]
        denom = 1.0
        for r in roots: denom *= (x[i] - r)
        poly = np.array([1.0])
        for r in roots: poly = np.convolve(poly, [-r, 1.0])
        basis = poly / denom
        contrib = y[i] * basis
        P_coeffs += contrib
        L_polys.append(basis)
        step_data.append({"i": i, "x_i": x[i], "y_i": y[i], "Mẫu số": round(denom, 8),
                          "Đa thức cơ sở L_i(x)": [round(v, 8) for v in basis.tolist()]})
    return pd.DataFrame(step_data), pd.DataFrame({"Bậc (x^k)": np.arange(n), "Hệ số": P_coeffs}), P_coeffs

def poly_str(coeffs):
    terms = []
    for deg, c in enumerate(coeffs):
        if abs(c) < 1e-12: continue
        sign = "+" if c >= 0 and len(terms) > 0 else ("-" if c < 0 and len(terms) > 0 else "")
        val = abs(c) if len(terms) > 0 else c
        term = f"{val:.8g}" if deg == 0 else (f"{val:.8g}*x" if deg == 1 else f"{val:.8g}*x^{deg}")
        terms.append(f"{sign} {term}" if len(terms) > 0 else term)
    return " ".join(terms) if terms else "0"

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Lagrange")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức")
    p.add_argument("-n", "--n_points", type=int, help="Số mốc nội suy")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(1.0, 3.0), (2.0, 5.0), (3.0, 7.0), (4.0, 9.0)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, _, msg = select_window_for_method(points, eval_x, "lagrange", n_pts_cfg, degree_cfg)

    steps_df, coeff_df, coeffs = lagrange(sub_points)
    out = ["="*75, "KẾT QUẢ NỘI SUY LAGRANGE", "="*75]
    if msg: out.append(msg)
    out.extend([
        f"Số mốc tính toán: {len(sub_points)}", "
1. Các bước tính toán và đa thức cơ sở L_i(x):",
        steps_df.to_string(index=False), "
2. Bảng hệ số đa thức P(x):", coeff_df.to_string(index=False),
        f"
3. Đa thức nội suy P(x):
P(x) = {poly_str(coeffs)}"
    ])
    if eval_x is not None:
        val = sum(c * (eval_x ** deg) for deg, c in enumerate(coeffs))
        out.append(f"
4. Giá trị tại x = {eval_x}:
P({eval_x}) = {val:.10f}")
    content = "
".join(out)
    print("
" + content)
    with open(args.output_file, "w", encoding="utf-8") as f: f.write(content + "
")
    print(f"
-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
