"""
Nội suy Newton tỷ sai phân (mốc bất kỳ)
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

def divided_diff(points, condition=1):
    x = [p[0] for p in points]; y = [p[1] for p in points]
    m = len(points)
    table = np.full((m, m), np.nan, dtype=float)
    table[:, 0] = y.copy()
    for j in range(1, m):
        for i in range(m - j):
            table[i, j] = (table[i + 1, j - 1] - table[i, j - 1]) / (x[i + j] - x[i])
    data = {"x_i": x, "y_i": table[:, 0]}
    for j in range(1, m): data[f"Cấp {j}"] = [table[i, j] if i < m - j else np.nan for i in range(m)]
    res = [table[:m - j, j][0] if condition == 1 else table[:m - j, j][-1] for j in range(m)]
    return pd.DataFrame(data), res

def newton_div(points, condition=1):
    n = len(points) - 1
    d_table, d_coeffs = divided_diff(points, condition)
    x_nodes = [p[0] for p in points]
    nodes = list(reversed(x_nodes)) if condition == 0 else x_nodes
    sdata, N_total = [], np.zeros(n + 1, dtype=float)
    B_poly = np.array([1.0])
    for i in range(n + 1):
        Di = d_coeffs[i]
        if i > 0: B_poly = np.convolve(B_poly, [-nodes[i - 1], 1.0])
        Ni = Di * B_poly
        N_total[:len(Ni)] += Ni
        sdata.append({"i": i, "D_i": Di, "B_i(x) (0->cao)": [round(v, 8) for v in B_poly.tolist()],
                      "N_i(x)": [round(v, 8) for v in Ni.tolist()]})
    return d_table, pd.DataFrame(sdata), pd.DataFrame({"Bậc (x^k)": np.arange(n + 1), "Hệ số": N_total}), N_total

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
    p = argparse.ArgumentParser(description="Nội suy Newton tỷ sai phân")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("--type", "-t", choices=["forward", "backward", "tien", "lui"])
    p.add_argument("-x", "--eval_x", type=float)
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức")
    p.add_argument("-n", "--n_points", type=int, help="Số mốc nội suy")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(0.87, 2.4), (1.24, 2.2), (2.99, 2.0), (3.67, 1.8), (4.23, 1.6)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, _, msg = select_window_for_method(points, eval_x, "newton_divided", n_pts_cfg, degree_cfg)

    ctype = args.type if args.type else params.get("type", "forward")
    condition = 0 if ("back" in ctype or "lui" in ctype) else 1
    tname = "TIẾN (Forward)" if condition == 1 else "LÙI (Backward)"

    d_table, s_df, c_df, coeffs = newton_div(sub_points, condition)
    out = ["="*75, f"KẾT QUẢ NỘI SUY NEWTON TỶ SAI PHÂN - {tname}", "="*75]
    if msg: out.append(msg)
    out.extend([
        f"Số mốc tính toán: {len(sub_points)}", "
1. Bảng tỷ sai phân:",
        d_table.to_string(index=True), "
2. Các bước dựng đa thức:", s_df.to_string(index=False),
        "
3. Bảng hệ số P(x):", c_df.to_string(index=False), f"
4. Đa thức nội suy P(x):
P(x) = {poly_str(coeffs)}"
    ])
    if eval_x is not None:
        val = sum(c * (eval_x ** deg) for deg, c in enumerate(coeffs))
        out.append(f"
5. Giá trị tại x = {eval_x}:
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
