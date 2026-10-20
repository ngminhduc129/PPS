"""
solve_midterm.py - All-in-One Midterm Solver for Numerical Analysis Coursework (Week 6)
Công cụ giải nhanh đề thi giữa kỳ: Tự động đọc file đề thi (.txt, .csv, .xlsx),
tự trích xuất mốc (x_i, y_i), điểm cần tính x* (hoặc x_bar), y* (y_bar),
tự cắt cửa sổ điểm tối ưu và tính toán kết quả theo phương pháp được chọn.

Ví dụ sử dụng:
    python solve_midterm.py --file ../../Midterm/20251GKtest1/20251GKtest1.csv --method stirling
    python solve_midterm.py --file ../../Midterm/20241GKtest1/20241GKtest1.txt --method gauss1 --degree 4
    python solve_midterm.py --file ../../Midterm/20251GKtest2/20251GKtest2.csv --method bessel
    python solve_midterm.py --file ../../Midterm/20251GKtest1/20251GKtest1.csv --compare
"""
import os, sys, argparse, math
import numpy as np
import pandas as pd

from pps_io import parse_universal_input, select_window_for_method

# --- CÁC HÀM THUẬT TOÁN NỘI SUY ---

def calc_even_diff(points):
    """Bảng sai phân hữu hạn mốc cách đều"""
    xv = [p[0] for p in points]; yv = [p[1] for p in points]; n = len(yv)
    tbl = np.full((n, n), np.nan); tbl[:, 0] = yv
    for j in range(1, n):
        for i in range(n - j): tbl[i, j] = tbl[i + 1, j - 1] - tbl[i, j - 1]
    data = {"x_i": xv, "y_i": yv}
    for j in range(1, n):
        col = np.full(n, np.nan); col[:(n - j)] = tbl[:(n - j), j]; data[f"Cấp {j}"] = col
    return pd.DataFrame(data)

def poly_to_str(coeffs, var="t"):
    terms = []
    for deg, c in enumerate(coeffs):
        if abs(c) < 1e-12: continue
        sign = "+" if c >= 0 and len(terms) > 0 else ("-" if c < 0 and len(terms) > 0 else "")
        val = abs(c) if len(terms) > 0 else c
        term = f"{val:.8g}" if deg == 0 else (f"{val:.8g}*{var}" if deg == 1 else f"{val:.8g}*{var}^{deg}")
        terms.append(f"{sign} {term}" if len(terms) > 0 else term)
    return " ".join(terms) if terms else "0"

def solve_gauss1(points, x0_idx, eval_x=None):
    n = len(points) - 1
    if (n + 1) % 2 == 0: n -= 1; points = points[:-1]
    h = points[1][0] - points[0][0]; x0_val = points[x0_idx][0]
    df_tbl = calc_even_diff(points)
    C_c = [df_tbl[f"Cấp {i}" if i > 0 else "y_i"].iloc[x0_idx - i // 2] for i in range(n + 1)]
    sdata, N_total = [], np.zeros(n + 1, dtype=float)
    Bc = np.array([1.0])
    for i in range(n + 1):
        Di = C_c[i] / math.factorial(i)
        if i == 1: Bc = np.array([0.0, 1.0])
        elif i > 1:
            k = i // 2 if i % 2 == 0 else (i - 1) // 2
            Bc = np.convolve(Bc, [-k, 1.0]) if i % 2 == 0 else np.convolve(Bc, [k, 1.0])
        Ni = Di * Bc
        N_total[:len(Ni)] += Ni
        sdata.append({"i": i, "C_i": C_c[i], "D_i": Di, "B_i(t)": [round(v, 6) for v in Bc.tolist()], "N_i(t)": [round(v, 6) for v in Ni.tolist()]})
    c_df = pd.DataFrame({"Bậc (t^k)": np.arange(n + 1), "Hệ số": N_total})
    val = None
    if eval_x is not None:
        t_val = (eval_x - x0_val) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(N_total))
    return {"name": "Gauss I", "table": df_tbl, "steps": pd.DataFrame(sdata), "coeffs_df": c_df, "coeffs": N_total,
            "x0": x0_val, "h": h, "poly_str": poly_to_str(N_total, "t"), "val": val, "var": "t"}

def solve_gauss2(points, x0_idx, eval_x=None):
    n = len(points) - 1
    if (n + 1) % 2 == 0: n -= 1; points = points[:-1]
    h = points[1][0] - points[0][0]; x0_val = points[x0_idx][0]
    df_tbl = calc_even_diff(points)
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
        sdata.append({"i": i, "C_i": C_c[i], "D_i": Di, "B_i(t)": [round(v, 6) for v in Bc.tolist()], "N_i(t)": [round(v, 6) for v in Ni.tolist()]})
    c_df = pd.DataFrame({"Bậc (t^k)": np.arange(n + 1), "Hệ số": N_total})
    val = None
    if eval_x is not None:
        t_val = (eval_x - x0_val) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(N_total))
    return {"name": "Gauss II", "table": df_tbl, "steps": pd.DataFrame(sdata), "coeffs_df": c_df, "coeffs": N_total,
            "x0": x0_val, "h": h, "poly_str": poly_to_str(N_total, "t"), "val": val, "var": "t"}

def solve_stirling(points, x0_idx, eval_x=None):
    n = len(points) - 1
    if (n + 1) % 2 == 0: n -= 1; points = points[:-1]
    h = points[1][0] - points[0][0]; x0_val = points[x0_idx][0]
    df_tbl = calc_even_diff(points)
    C_c = []
    for i in range(n + 1):
        col = f"Cấp {i}" if i > 0 else "y_i"
        if i == 0: C_c.append(df_tbl[col].iloc[x0_idx])
        elif i % 2 == 1:
            k = (i + 1) // 2
            C_c.append((df_tbl[col].iloc[x0_idx - k] + df_tbl[col].iloc[x0_idx - k + 1]) / 2.0)
        else:
            k = i // 2
            C_c.append(df_tbl[col].iloc[x0_idx - k])
    B_list = [np.array([1.0]), np.array([0.0, 1.0])]
    for i in range(2, n + 1):
        if i % 2 == 0: Bi = np.convolve(B_list[i - 1], [0.0, 1.0])
        else:
            k = (i - 1) // 2
            Bi = np.convolve(B_list[i - 2], [-float(k**2), 0.0, 1.0])
        B_list.append(Bi)
    sdata, N_total = [], np.zeros(n + 1, dtype=float)
    for i in range(n + 1):
        Di = C_c[i] / math.factorial(i)
        Ni = Di * B_list[i]
        N_total[:len(Ni)] += Ni
        sdata.append({"i": i, "C_i": C_c[i], "D_i": Di, "B_i(t)": [round(v, 6) for v in B_list[i].tolist()], "N_i(t)": [round(v, 6) for v in Ni.tolist()]})
    c_df = pd.DataFrame({"Bậc (t^k)": np.arange(n + 1), "Hệ số": N_total})
    val = None
    if eval_x is not None:
        t_val = (eval_x - x0_val) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(N_total))
    return {"name": "Stirling", "table": df_tbl, "steps": pd.DataFrame(sdata), "coeffs_df": c_df, "coeffs": N_total,
            "x0": x0_val, "h": h, "poly_str": poly_to_str(N_total, "t"), "val": val, "var": "t"}

def solve_bessel(points, x0_idx, eval_x=None):
    n = len(points) - 1
    if (n + 1) % 2 != 0: n -= 1; points = points[:-1]
    h = points[1][0] - points[0][0]; x0_val = points[x0_idx][0]
    df_tbl = calc_even_diff(points)
    C_c = []
    for i in range(n + 1):
        col = f"Cấp {i}" if i > 0 else "y_i"
        if i % 2 == 0:
            k = i // 2
            C_c.append((df_tbl[col].iloc[x0_idx - k] + df_tbl[col].iloc[x0_idx - k + 1]) / 2.0)
        else:
            k = (i - 1) // 2
            C_c.append(df_tbl[col].iloc[x0_idx - k])
    B_list = [np.array([1.0]), np.array([-0.5, 1.0])]
    for i in range(2, n + 1):
        if i % 2 == 0:
            k = i // 2
            nfact = np.convolve(np.array([-float(k), 1.0]), np.array([float(k - 1), 1.0]))
            Bi = np.convolve(B_list[i - 2], nfact)
        else:
            Bi = np.convolve(B_list[i - 1], [-0.5, 1.0])
        B_list.append(Bi)
    sdata, N_total = [], np.zeros(n + 1, dtype=float)
    for i in range(n + 1):
        Di = C_c[i] / math.factorial(i)
        Ni = Di * B_list[i]
        N_total[:len(Ni)] += Ni
        sdata.append({"i": i, "C_i": C_c[i], "D_i": Di, "B_i(t)": [round(v, 6) for v in B_list[i].tolist()], "N_i(t)": [round(v, 6) for v in Ni.tolist()]})
    c_df = pd.DataFrame({"Bậc (t^k)": np.arange(n + 1), "Hệ số": N_total})
    val = None
    if eval_x is not None:
        t_val = (eval_x - x0_val) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(N_total))
    return {"name": "Bessel", "table": df_tbl, "steps": pd.DataFrame(sdata), "coeffs_df": c_df, "coeffs": N_total,
            "x0": x0_val, "h": h, "poly_str": poly_to_str(N_total, "t"), "val": val, "var": "t"}

def solve_newton_finite(points, is_forward=True, eval_x=None):
    h = points[1][0] - points[0][0]; n = len(points) - 1
    df_tbl = calc_even_diff(points)
    x0_val = points[0][0] if is_forward else points[-1][0]
    coeffs_diff = [df_tbl[f"Cấp {i}" if i > 0 else "y_i"].iloc[0] for i in range(n + 1)] if is_forward else \
                  [df_tbl[f"Cấp {i}" if i > 0 else "y_i"].dropna().iloc[-1] for i in range(n + 1)]
    sdata, P_t = [], np.zeros(n + 1, dtype=float)
    B_t = np.array([1.0])
    for i in range(n + 1):
        Di = coeffs_diff[i] / math.factorial(i)
        if i == 1: B_t = np.array([0.0, 1.0])
        elif i > 1:
            term = np.array([-(i - 1), 1.0]) if is_forward else np.array([(i - 1), 1.0])
            B_t = np.convolve(B_t, term)
        Ni = Di * B_t
        P_t[:len(Ni)] += Ni
        sdata.append({"i": i, "Sai phân": coeffs_diff[i], "D_i": Di, "B_i(t)": [round(v, 6) for v in B_t.tolist()]})
    val = None
    if eval_x is not None:
        t_val = (eval_x - x0_val) / h
        val = sum(c * (t_val ** deg) for deg, c in enumerate(P_t))
    method_name = "Newton Tiến (Finite Forward)" if is_forward else "Newton Lùi (Finite Backward)"
    return {"name": method_name, "table": df_tbl, "steps": pd.DataFrame(sdata), "coeffs_df": pd.DataFrame({"Bậc (t^k)": np.arange(n+1), "Hệ số": P_t}),
            "coeffs": P_t, "x0": x0_val, "h": h, "poly_str": poly_to_str(P_t, "t"), "val": val, "var": "t"}

def solve_lagrange(points, eval_x=None):
    n = len(points)
    xs = [p[0] for p in points]; ys = [p[1] for p in points]
    val = None
    if eval_x is not None:
        total = 0.0
        for i in range(n):
            term = ys[i]
            for j in range(n):
                if i != j:
                    term *= (eval_x - xs[j]) / (xs[i] - xs[j])
            total += term
        val = total
    return {"name": "Lagrange", "table": pd.DataFrame({"x_i": xs, "y_i": ys}), "steps": pd.DataFrame(),
            "coeffs_df": pd.DataFrame(), "coeffs": None, "x0": xs[0], "h": None, "poly_str": "Đa thức Lagrange P_n(x)", "val": val, "var": "x"}

# --- MAIN CONTROLLER ---

def main():
    parser = argparse.ArgumentParser(description="All-in-One Midterm Numerical Analysis Solver")
    parser.add_argument("--file", "-f", default=None, help="Đường dẫn tới file đề thi (.txt, .csv, .xlsx)")
    parser.add_argument("--method", "-m", default="stirling",
                        choices=["gauss1", "gauss2", "stirling", "bessel", "newton_forward", "newton_backward", "newton_divided", "lagrange"],
                        help="Phương pháp nội suy cần giải")
    parser.add_argument("--degree", "-d", type=int, default=None, help="Bậc đa thức (ví dụ 4)")
    parser.add_argument("--n_points", "-n", type=int, default=None, help="Số mốc nội suy (ví dụ 5)")
    parser.add_argument("--eval_x", "-x", type=float, default=None, help="Điểm cần tính x* (tự động nhận từ file nếu có)")
    parser.add_argument("--all", action="store_true", help="Sử dụng toàn bộ mốc dữ liệu (không cắt cửa sổ)")
    parser.add_argument("--compare", action="store_true", help="Chạy và so sánh kết quả tất cả các phương pháp")
    parser.add_argument("--out", "-o", default="midterm_result.txt", help="File xuất kết quả")
    args = parser.parse_args()

    # Tìm file mặc định nếu không truyền
    target_file = args.file
    if not target_file:
        s_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(os.path.dirname(s_dir))
        candidates = [
            "input.txt",
            os.path.join(s_dir, "input.txt"),
            os.path.join(repo_root, "Midterm", "20251GKtest1", "20251GKtest1.csv"),
            os.path.join(repo_root, "Midterm", "20241GKtest1", "20241GKtest1.csv"),
            "../../Midterm/20251GKtest1/20251GKtest1.csv",
            "Midterm/20251GKtest1/20251GKtest1.csv"
        ]
        for candidate in candidates:
            if os.path.exists(candidate):
                target_file = candidate
                break

    if not target_file or not os.path.exists(target_file):
        print("[Lỗi] Không tìm thấy file dữ liệu đầu vào. Hãy dùng tham số --file <đường_dẫn>")
        return

    print("=" * 80)
    print(f" ĐỌC VÀ PHÂN TÍCH ĐỀ THI: {os.path.abspath(target_file)}")
    print("=" * 80)

    points, params = parse_universal_input(target_file)
    if not points:
        print("[Lỗi] Không tìm thấy mốc dữ liệu (x_i, y_i) trong file!")
        return

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    eval_y = params.get("eval_y", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    print(f"-> Tổng số mốc dữ liệu: {len(points)} điểm")
    print(f"-> Điểm cần tính x* (eval_x): {eval_x}")
    if eval_y is not None:
        print(f"-> Điểm cần tính ngược y* (eval_y): {eval_y}")

    # Chế độ COMPARE: So sánh tất cả các phương pháp
    methods_to_run = ["gauss1", "gauss2", "stirling", "bessel", "newton_forward", "lagrange"] if args.compare else [args.method]

    output_lines = [
        "=" * 80,
        f"KẾT QUẢ GIẢI ĐỀ THI GIỮA KỲ: {os.path.basename(target_file)}",
        "=" * 80,
        f"Tổng số mốc dữ liệu: {len(points)}",
        f"Điểm cần tính x*: {eval_x}",
    ]
    if eval_y is not None:
        output_lines.append(f"Điểm cần tính ngược y*: {eval_y}")

    for m in methods_to_run:
        map_method = "bessel" if m == "bessel" else ("gauss1" if m == "gauss1" else ("gauss2" if m == "gauss2" else ("stirling" if m == "stirling" else "newton_finite")))
        sub_pts, x0_idx, msg = select_window_for_method(points, eval_x, map_method, n_pts_cfg, degree_cfg)
        
        output_lines.append("\n" + "-" * 75)
        output_lines.append(f"PHƯƠNG PHÁP: {m.upper()}")
        output_lines.append(msg)
        output_lines.append("-" * 75)

        res = None
        if m == "gauss1":
            res = solve_gauss1(sub_pts, x0_idx, eval_x)
        elif m == "gauss2":
            res = solve_gauss2(sub_pts, x0_idx, eval_x)
        elif m == "stirling":
            res = solve_stirling(sub_pts, x0_idx, eval_x)
        elif m == "bessel":
            res = solve_bessel(sub_pts, x0_idx, eval_x)
        elif m == "newton_forward":
            res = solve_newton_finite(sub_pts, is_forward=True, eval_x=eval_x)
        elif m == "newton_backward":
            res = solve_newton_finite(sub_pts, is_forward=False, eval_x=eval_x)
        elif m == "lagrange":
            res = solve_lagrange(sub_pts, eval_x=eval_x)

        if res:
            output_lines.append(f"1. Bảng số liệu & Sai phân ({len(sub_pts)} mốc):")
            output_lines.append(res["table"].to_string(index=True))
            if not res["steps"].empty:
                output_lines.append(f"\n2. Bảng các bước tính {res['name']}:")
                output_lines.append(res["steps"].to_string(index=False))
            if not res["coeffs_df"].empty and "Hệ số" in res["coeffs_df"].columns:
                output_lines.append(f"\n3. Bảng hệ số đa thức P({res['var']}):")
                output_lines.append(res["coeffs_df"].to_string(index=False))
            output_lines.append(f"\n4. Dạng đa thức:\nP({res['var']}) = {res['poly_str']}")
            if res["val"] is not None:
                output_lines.append(f"\n5. KẾT QUẢ NỘI SUY TẠI x* = {eval_x}:")
                output_lines.append(f"   >>> P({eval_x}) = {res['val']:.10f} <<<")

    final_text = "\n".join(output_lines)
    print("\n" + final_text)

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(final_text + "\n")
    print(f"\n[Thành công] Đã lưu kết quả chi tiết ra: '{os.path.abspath(args.out)}'")

if __name__ == "__main__":
    main()
