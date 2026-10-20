"""
Nội suy Newton sai phân hữu hạn (mốc cách đều)
Nhập: input.txt hoặc file đề thi Midterm (.txt, .csv, .xlsx)
Xuất: output.txt
"""
import os, sys, math, argparse
import numpy as np
import pandas as pd

try:
    from pps_io import parse_universal_input, select_window_for_method
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pps_io import parse_universal_input, select_window_for_method

def even_diff(points, condition=1):
    xv = [p[0] for p in points]; yv = [p[1] for p in points]; n = len(yv)
    tbl = np.full((n, n), np.nan); tbl[:, 0] = yv
    for j in range(1, n):
        for i in range(n - j): tbl[i, j] = tbl[i + 1, j - 1] - tbl[i, j - 1]
    data = {"x_i": xv, "y_i": yv}
    for j in range(1, n): data[f"Cấp {j}"] = tbl[:, j]
    coeffs = tbl[0, :].tolist() if condition == 1 else [tbl[n - 1 - j, j] for j in range(n)]
    return pd.DataFrame(data), coeffs

def newton_fg(points, condition=1):
    n = len(points) - 1
    d_tbl, d_coeffs = even_diff(points, condition)
    sdata, N_total = [], np.zeros(n + 1, dtype=float)
    B_prev = np.array([1.0])
    for i in range(n + 1):
        Di = d_coeffs[i] / math.factorial(i)
        if i == 0: Bc = np.array([1.0])
        else:
            k = i - 1
            Bc = (np.concatenate(([0.0], B_prev)) - np.concatenate((k * B_prev, [0.0]))) if condition == 1 \
                 else (np.concatenate(([0.0], B_prev)) + np.concatenate((k * B_prev, [0.0])))
        B_prev = Bc.copy()
        Ni = Di * Bc
        N_total[:len(Ni)] += Ni
        sdata.append({"i": i, "C_i": d_coeffs[i], "D_i": Di, "B_i(t) (0->cao)": [round(v, 8) for v in Bc.tolist()],
                      "N_i(t)": [round(v, 8) for v in Ni.tolist()]})
    return d_tbl, pd.DataFrame(sdata), pd.DataFrame({"Bậc (t^k)": np.arange(n + 1), "Hệ số": N_total}), N_total

def poly_str(coeffs):
    terms = []
    for deg, c in enumerate(coeffs):
        if abs(c) < 1e-12: continue
        sign = "+" if c >= 0 and len(terms) > 0 else ("-" if c < 0 and len(terms) > 0 else "")
        val = abs(c) if len(terms) > 0 else c
        term = f"{val:.8g}" if deg == 0 else (f"{val:.8g}*t" if deg == 1 else f"{val:.8g}*t^{deg}")
        terms.append(f"{sign} {term}" if len(terms) > 0 else term)
    return " ".join(terms) if terms else "0"

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Newton mốc cách đều")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("--type", "-t", choices=["forward", "backward", "tien", "lui"], help="Loại nội suy (tiến/lùi)")
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức (ví dụ 4)")
    p.add_argument("-n", "--n_points", type=int, help="Số mốc nội suy (ví dụ 5)")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(1.4, 0.9523), (1.5, 0.9661), (1.6, 0.9753), (1.7, 0.9838), (1.8, 0.9891), (1.9, 0.9928), (2.0, 0.9)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, _, msg = select_window_for_method(points, eval_x, "newton_finite", n_pts_cfg, degree_cfg)

    ctype = args.type if args.type else params.get("type", "forward")
    condition = 0 if ("back" in ctype or "lui" in ctype) else 1
    tname = "TIẾN (Forward)" if condition == 1 else "LÙI (Backward)"
    h = sub_points[1][0] - sub_points[0][0]
    x_base = sub_points[0][0] if condition == 1 else sub_points[-1][0]

    d_tbl, s_df, c_df, coeffs = newton_fg(sub_points, condition)
    out = ["="*75, f"KẾT QUẢ NỘI SUY NEWTON MỐC CÁCH ĐỀU - {tname}", "="*75]
    if msg: out.append(msg)
    out.extend([
        f"Số mốc tính toán: {len(sub_points)}, Bước h: {h:.6f}, Mốc gốc x0: {x_base}",
        f"Chuẩn hóa: t = (x - {x_base}) / {h:.6f}", "
1. Bảng sai phân hữu hạn:",
        d_tbl.to_string(index=True), "
2. Các bước dựng đa thức:", s_df.to_string(index=False),
        "
3. Bảng hệ số P(t):", c_df.to_string(index=False), f"
4. Đa thức nội suy P(t):
P(t) = {poly_str(coeffs)}"
    ])
    if eval_x is not None:
        t_val = (eval_x - x_base) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(coeffs))
        out.append(f"
5. Giá trị tại x = {eval_x}:
t = {t_val:.6f}
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
