"""
generate_clean_code.py - Ghi code trực tiếp cho các phương pháp, đảm bảo chuẩn 100% cú pháp
"""
import os, shutil

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def write_method_files(folder, filename_stem, code_content):
    f_main = os.path.join(BASE_DIR, folder, "main.py")
    f_method = os.path.join(BASE_DIR, folder, f"{filename_stem}.py")
    with open(f_main, "w", encoding="utf-8") as f:
        f.write(code_content)
    with open(f_method, "w", encoding="utf-8") as f:
        f.write(code_content)
    print(f"Đã ghi thành công: {folder}/main.py và {folder}/{filename_stem}.py")

# ==============================================================================
# 07_Gauss2
# ==============================================================================
gauss2_text = r'''"""
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
    content = "\n".join(out)
    print("\n" + content)
    with open(args.output_file, "w", encoding="utf-8") as f: f.write(content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
'''
write_method_files("07_Gauss2", "gauss2", gauss2_text)

# ==============================================================================
# 08_Stirling
# ==============================================================================
with open(os.path.join(BASE_DIR, "08_Stirling", "stirling.py"), "r", encoding="utf-8") as f:
    stirling_text = f.read()
write_method_files("08_Stirling", "stirling", stirling_text)

# ==============================================================================
# 09_Bessel
# ==============================================================================
with open(os.path.join(BASE_DIR, "09_Bessel", "bessel.py"), "r", encoding="utf-8") as f:
    bessel_text = f.read()
write_method_files("09_Bessel", "bessel", bessel_text)

# ==============================================================================
# 05_Newton_Finite
# ==============================================================================
newton_finite_text = r'''"""
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
        f"Chuẩn hóa: t = (x - {x_base}) / {h:.6f}", "\n1. Bảng sai phân hữu hạn:",
        d_tbl.to_string(index=True), "\n2. Các bước dựng đa thức:", s_df.to_string(index=False),
        "\n3. Bảng hệ số P(t):", c_df.to_string(index=False), f"\n4. Đa thức nội suy P(t):\nP(t) = {poly_str(coeffs)}"
    ])
    if eval_x is not None:
        t_val = (eval_x - x_base) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(coeffs))
        out.append(f"\n5. Giá trị tại x = {eval_x}:\nt = {t_val:.6f}\nP({eval_x}) = {val:.10f}")
    content = "\n".join(out)
    print("\n" + content)
    with open(args.output_file, "w", encoding="utf-8") as f: f.write(content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
'''
write_method_files("05_Newton_Finite", "newton_finite", newton_finite_text)

# ==============================================================================
# 04_Newton_Divided
# ==============================================================================
newton_divided_text = r'''"""
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
        f"Số mốc tính toán: {len(sub_points)}", "\n1. Bảng tỷ sai phân:",
        d_table.to_string(index=True), "\n2. Các bước dựng đa thức:", s_df.to_string(index=False),
        "\n3. Bảng hệ số P(x):", c_df.to_string(index=False), f"\n4. Đa thức nội suy P(x):\nP(x) = {poly_str(coeffs)}"
    ])
    if eval_x is not None:
        val = sum(c * (eval_x ** deg) for deg, c in enumerate(coeffs))
        out.append(f"\n5. Giá trị tại x = {eval_x}:\nP({eval_x}) = {val:.10f}")
    content = "\n".join(out)
    print("\n" + content)
    with open(args.output_file, "w", encoding="utf-8") as f: f.write(content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
'''
write_method_files("04_Newton_Divided", "newton_divided", newton_divided_text)

# ==============================================================================
# 03_Lagrange
# ==============================================================================
lagrange_text = r'''"""
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
        f"Số mốc tính toán: {len(sub_points)}", "\n1. Các bước tính toán và đa thức cơ sở L_i(x):",
        steps_df.to_string(index=False), "\n2. Bảng hệ số đa thức P(x):", coeff_df.to_string(index=False),
        f"\n3. Đa thức nội suy P(x):\nP(x) = {poly_str(coeffs)}"
    ])
    if eval_x is not None:
        val = sum(c * (eval_x ** deg) for deg, c in enumerate(coeffs))
        out.append(f"\n4. Giá trị tại x = {eval_x}:\nP({eval_x}) = {val:.10f}")
    content = "\n".join(out)
    print("\n" + content)
    with open(args.output_file, "w", encoding="utf-8") as f: f.write(content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
'''
write_method_files("03_Lagrange", "lagrange", lagrange_text)

print("Tất cả code phương pháp đã được cập nhật hoàn chỉnh!")
