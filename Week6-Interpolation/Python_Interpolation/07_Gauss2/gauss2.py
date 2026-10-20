"""
Nội suy Gauss II (Gauss nghịch)
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

def even_diff(points):
    xv = [p[0] for p in points]; yv = [p[1] for p in points]; n = len(yv)
    tbl = np.full((n, n), np.nan); tbl[:, 0] = yv
    for j in range(1, n):
        for i in range(n - j): tbl[i, j] = tbl[i + 1, j - 1] - tbl[i, j - 1]
    data = {"x_i": xv, "y_i": yv}
    for j in range(1, n):
        col = np.full(n, np.nan); col[:(n - j)] = tbl[:(n - j), j]; data[f"Cấp {j}"] = col
    return pd.DataFrame(data)

def gauss2(points, x0_idx):
    n = len(points) - 1
    if (n + 1) % 2 == 0: raise ValueError("Gauss II yêu cầu số mốc lẻ.")
    df_tbl = even_diff(points)
    C_c = [df_tbl[f"Cấp {i}" if i > 0 else "y_i"].iloc[x0_idx - (i + 1) // 2] for i in range(n + 1)]
    sdata, N_total = [], np.zeros(n + 1, dtype=float)
    Bc = np.array([1.0])
    for i in range(n + 1):
        Di = C_c[i] / math.factorial(i)
        if i == 1: Bc = np.array([0.0, 1.0])
        elif i > 1:
            k = (i - 1) // 2 if i % 2 == 1 else i // 2
            Bc = np.convolve(Bc, [k, 1.0]) if i % 2 == 0 else np.convolve(Bc, [-k, 1.0])
        Ni = Di * Bc
        N_total[:len(Ni)] += Ni
        sdata.append({"i": i, "C_i": C_c[i], "D_i": Di, "B_i(t) (0->cao)": [round(v, 8) for v in Bc.tolist()],
                      "N_i(t)": [round(v, 8) for v in Ni.tolist()]})
    return df_tbl, pd.DataFrame(sdata), pd.DataFrame({"Bậc (t^k)": np.arange(n + 1), "Hệ số": N_total}), N_total

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
    p = argparse.ArgumentParser(description="Nội suy Gauss II (Gauss nghịch)")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("--x0_index", type=int, help="Chỉ số mốc trung tâm x0")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức (ví dụ 4)")
    p.add_argument("-n", "--n_points", type=int, help="Số mốc nội suy (ví dụ 5)")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(9.2, 9.4341319), (9.3, 9.4307764), (9.4, 9.4261142), (9.5, 9.4211191),
                  (9.6, 9.4170553), (9.7, 9.4147476), (9.8, 9.4152900), (9.9, 9.4196762), (10.0, 9.4288617)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, x0_idx, msg = select_window_for_method(points, eval_x, "gauss2", n_pts_cfg, degree_cfg)
    if args.x0_index is not None:
        x0_idx = args.x0_index
    elif "x0_index" in params and len(points) == len(sub_points):
        x0_idx = params["x0_index"]

    if len(sub_points) % 2 == 0:
        sub_points = sub_points[:-1]
        x0_idx = min(x0_idx, len(sub_points) - 1)

    h = sub_points[1][0] - sub_points[0][0]
    x0_val = sub_points[x0_idx][0]

    df_tbl, s_df, c_df, coeffs = gauss2(sub_points, x0_idx)
    out = ["="*75, "KẾT QUẢ NỘI SUY GAUSS II (GAUSS NGHỊCH)", "="*75]
    if msg: out.append(msg)
    out.extend([
        f"Số mốc tính toán: {len(sub_points)}, Bước h: {h:.6f}, Mốc x0 (idx {x0_idx}): {x0_val}",
        f"Chuẩn hóa: t = (x - {x0_val}) / {h:.6f}", "\n1. Bảng sai phân hữu hạn:",
        df_tbl.to_string(index=True), "\n2. Các bước dựng đa thức:", s_df.to_string(index=False),
        "\n3. Bảng hệ số P(t):", c_df.to_string(index=False), f"\n4. Đa thức nội suy P(t):\nP(t) = {poly_str(coeffs)}"
    ])
    if eval_x is not None:
        t_val = (eval_x - x0_val) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(coeffs))
        out.append(f"\n5. Giá trị tại x = {eval_x}:\nt = {t_val:.6f}\nP({eval_x}) = {val:.10f}")
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
