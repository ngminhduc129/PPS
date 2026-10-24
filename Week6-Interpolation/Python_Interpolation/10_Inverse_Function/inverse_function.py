"""
Nội suy Hàm Ngược (Inverse Function Interpolation)
Đổi vai trò biến số x và y để xây dựng đa thức nội suy ngược x = P(y), từ đó tính nghiệm x* = P(y*).
1. BẢNG PHÂN TÍCH CÁC KHOẢNG ĐƠN ĐIỆU CỦA HÀM SỐ
2. BẢNG MỐC NỘI SUY NGƯỢC (HOÁN ĐỔI VAI TRÒ x VÀ y)
3. BẢNG TỶ SAI PHÂN THEO BIẾN y (DIVIDED DIFFERENCES)
4. ĐA THỨC NỘI SUY NGƯỢC x = P(y) (DẠNG NEWTON & CHÍNH TẮC)
5. KẾT QUẢ TÍNH NGHIỆM x* = P(y*) VÀ ĐẠO HÀM TẠI y*

Nhập: input.txt hoặc file đề thi Midterm (.txt, .csv, .xlsx)
Xuất: output.txt và in ra màn hình.
"""
import os, sys, math, argparse
import numpy as np
import pandas as pd

try:
    from pps_io import parse_universal_input, select_window_for_method
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pps_io import parse_universal_input, select_window_for_method

def find_monotonic_intervals(points):
    """
    Tìm tất cả các khoảng đơn điệu liên tục (tăng hoặc giảm ngặt) của bảng số liệu (x_i, y_i).
    """
    if len(points) < 2:
        return pd.DataFrame(columns=["STT", "Loại đơn điệu", "Điểm đầu (x, y)", "Điểm cuối (x, y)", "Miền giá trị y"])

    intervals = []
    start_pt = points[0]
    curr_dir = np.sign(points[1][1] - points[0][1])
    if curr_dir == 0:
        curr_dir = 1

    for i in range(1, len(points) - 1):
        diff = points[i + 1][1] - points[i][1]
        new_dir = np.sign(diff)
        if new_dir != 0 and new_dir != curr_dir:
            int_type = "Đồng biến (Tăng)" if curr_dir > 0 else "Nghịch biến (Giảm)"
            y_min = min(start_pt[1], points[i][1])
            y_max = max(start_pt[1], points[i][1])
            intervals.append({
                "STT": len(intervals) + 1,
                "Loại đơn điệu": int_type,
                "Điểm đầu (x, y)": f"({start_pt[0]:.4f}, {start_pt[1]:.6f})",
                "Điểm cuối (x, y)": f"({points[i][0]:.4f}, {points[i][1]:.6f})",
                "Miền giá trị y": f"[{y_min:.6f}, {y_max:.6f}]",
                "start_idx": points.index(start_pt),
                "end_idx": i
            })
            start_pt = points[i]
            curr_dir = new_dir

    int_type = "Đồng biến (Tăng)" if curr_dir > 0 else "Nghịch biến (Giảm)"
    y_min = min(start_pt[1], points[-1][1])
    y_max = max(start_pt[1], points[-1][1])
    intervals.append({
        "STT": len(intervals) + 1,
        "Loại đơn điệu": int_type,
        "Điểm đầu (x, y)": f"({start_pt[0]:.4f}, {start_pt[1]:.6f})",
        "Điểm cuối (x, y)": f"({points[-1][0]:.4f}, {points[-1][1]:.6f})",
        "Miền giá trị y": f"[{y_min:.6f}, {y_max:.6f}]",
        "start_idx": points.index(start_pt),
        "end_idx": len(points) - 1
    })
    return pd.DataFrame(intervals)

def divided_diff_inverse(y_arr, x_arr):
    """
    Tính bảng tỷ sai phân với biến độc lập là y và giá trị hàm là x.
    """
    n = len(y_arr)
    tbl = np.full((n, n), np.nan, dtype=float)
    tbl[:, 0] = x_arr

    for j in range(1, n):
        for i in range(n - j):
            denom = y_arr[i + j] - y_arr[i]
            if abs(denom) < 1e-14:
                tbl[i, j] = 0.0
            else:
                tbl[i, j] = (tbl[i + 1, j - 1] - tbl[i, j - 1]) / denom

    data = {"y_i (Biến)": y_arr, "x_i (Hàm)": x_arr}
    for j in range(1, n):
        col = np.full(n, np.nan)
        col[:(n - j)] = tbl[:(n - j), j]
        data[f"Tỷ SP cấp {j}"] = col
    df_tbl = pd.DataFrame(data)
    d_coeffs = [tbl[0, j] for j in range(n)]
    return df_tbl, d_coeffs

def build_inverse_polynomial(y_arr, d_coeffs):
    """
    Khai triển đa thức Newton thành dạng chính tắc P(y) = c0 + c1*y + c2*y^2 + ...
    """
    n = len(y_arr)
    p_total = np.zeros(n, dtype=float)
    basis = np.array([1.0])
    steps = []

    for i in range(n):
        di = d_coeffs[i]
        term = di * basis
        p_total[:len(term)] += term
        steps.append({
            "i": i,
            "Tỷ sai phân D_i": di,
            "Cơ sở omega_i(y) (0->cao)": [round(v, 8) for v in basis.tolist()],
            "Số hạng N_i(y)": [round(v, 8) for v in term.tolist()]
        })
        if i < n - 1:
            basis = np.convolve(basis, [-y_arr[i], 1.0])

    df_steps = pd.DataFrame(steps)
    df_coeffs = pd.DataFrame({"Bậc (y^k)": np.arange(n), "Hệ số": p_total})
    return df_steps, df_coeffs, p_total

def horner_eval_and_derivatives(coeffs, c):
    """
    Tính giá trị đa thức P(c) và đạo hàm các cấp bằng lược đồ Horner.
    coeffs: [c0, c1, ..., cn] theo bậc tăng dần
    """
    deg = len(coeffs) - 1
    # Chuyển về dạng bậc cao xuống thấp cho Horner
    a = list(reversed(coeffs))
    
    # Tính giá trị P(c)
    val = 0.0
    for coef in a:
        val = val * c + coef
        
    # Tính đạo hàm các cấp
    curr_poly = np.poly1d(a)
    derivs = []
    for k in range(1, min(deg + 1, 4)):
        curr_poly = np.polyder(curr_poly)
        derivs.append((k, curr_poly(c)))
        
    return val, derivs

def poly_str(coeffs, var="y"):
    terms = []
    for deg, c in enumerate(coeffs):
        if abs(c) < 1e-12:
            continue
        sign = "+" if c >= 0 and len(terms) > 0 else ("-" if c < 0 and len(terms) > 0 else "")
        val = abs(c) if len(terms) > 0 else c
        term = f"{val:.8g}" if deg == 0 else (f"{val:.8g}*{var}" if deg == 1 else f"{val:.8g}*{var}^{deg}")
        terms.append(f"{sign} {term}" if len(terms) > 0 else term)
    return " ".join(terms) if terms else "0"

def solve_inverse_function(points, eval_y=None, degree=None, n_points=None):
    """
    Thực hiện giải bài toán nội suy ngược bằng phương pháp hàm ngược.
    """
    total_pts = len(points)
    df_intervals = find_monotonic_intervals(points)
    
    # 1. Xác định khoảng đơn điệu phù hợp với eval_y
    selected_interval = None
    if eval_y is not None and not df_intervals.empty:
        for _, row in df_intervals.iterrows():
            s_idx = row["start_idx"]
            e_idx = row["end_idx"]
            y_min = min(points[s_idx][1], points[e_idx][1])
            y_max = max(points[s_idx][1], points[e_idx][1])
            if y_min <= eval_y <= y_max:
                selected_interval = (s_idx, e_idx, row["Loại đơn điệu"])
                break

    # Nếu không nằm trọn khoảng nào, lấy khoảng có mốc y gần eval_y nhất
    if selected_interval is None and eval_y is not None:
        diffs = [abs(p[1] - eval_y) for p in points]
        best_k = int(np.argmin(diffs))
        for _, row in df_intervals.iterrows():
            if row["start_idx"] <= best_k <= row["end_idx"]:
                selected_interval = (row["start_idx"], row["end_idx"], row["Loại đơn điệu"])
                break

    if selected_interval is None:
        sub_all = points
        msg_monotone = "Sử dụng toàn bộ mốc dữ liệu có sẵn."
    else:
        s_idx, e_idx, m_type = selected_interval
        sub_all = points[s_idx:e_idx + 1]
        msg_monotone = f"Đã chọn khoảng đơn điệu {m_type} từ mốc #{s_idx} ({points[s_idx][0]:.4f}) đến #{e_idx} ({points[e_idx][0]:.4f})."

    # 2. Cắt cửa sổ điểm lân cận quanh eval_y trong khoảng đơn điệu
    target_n = degree + 1 if degree is not None else (int(n_points) if n_points not in [None, "all"] else min(len(sub_all), 5))
    target_n = min(target_n, len(sub_all))

    if eval_y is not None and len(sub_all) > target_n:
        diffs_y = [abs(p[1] - eval_y) for p in sub_all]
        k_local = int(np.argmin(diffs_y))
        st = max(0, min(k_local - target_n // 2, len(sub_all) - target_n))
        chosen_points = sub_all[st:st + target_n]
    else:
        chosen_points = sub_all[:target_n]

    # 3. Hoán đổi vai trò biến số: x thành hàm, y thành biến
    y_raw = np.array([p[1] for p in chosen_points], dtype=float)
    x_raw = np.array([p[0] for p in chosen_points], dtype=float)

    # Đảm bảo y_arr tăng dần để lập bảng tỷ sai phân chuẩn xác
    sort_idx = np.argsort(y_raw)
    y_sorted = y_raw[sort_idx]
    x_sorted = x_raw[sort_idx]

    # 4. Xây dựng bảng tỷ sai phân và đa thức
    df_diff, d_coeffs = divided_diff_inverse(y_sorted, x_sorted)
    df_steps, df_coeffs, poly_c = build_inverse_polynomial(y_sorted, d_coeffs)
    poly_representation = poly_str(poly_c, var="y")

    # 5. Tính giá trị nghiệm tại eval_y
    val_x = None
    derivs = []
    if eval_y is not None:
        val_x, derivs = horner_eval_and_derivatives(poly_c, eval_y)

    return {
        "df_intervals": df_intervals,
        "msg_monotone": msg_monotone,
        "chosen_points": chosen_points,
        "df_diff": df_diff,
        "df_steps": df_steps,
        "df_coeffs": df_coeffs,
        "poly_coeffs": poly_c,
        "poly_str": poly_representation,
        "eval_y": eval_y,
        "val_x": val_x,
        "derivs": derivs
    }

def format_inverse_report(res):
    lines = [
        "=" * 80,
        "BÁO CÁO PHƯƠNG PHÁP NỘI SUY HÀM NGƯỢC (INVERSE FUNCTION INTERPOLATION)",
        "=" * 80,
        f"Nguyên lý: Hoán đổi vai trò biến số, xây dựng đa thức nội suy x = P(y).",
        f"Thông tin khoảng đơn điệu: {res['msg_monotone']}",
        "\n1. BẢNG CÁC KHOẢNG ĐƠN ĐIỆU CỦA HÀM SỐ TRÊN TẬP DỮ LIỆU:",
        res["df_intervals"].to_string(index=False) if not res["df_intervals"].empty else "Không đủ điểm để phân tích.",
        "\n2. CÁC MỐC NỘI SUY ĐƯỢC CHỌN VÀ HOÁN ĐỔI VAI TRÒ (y_i -> BIẾN, x_i -> HÀM):"
    ]
    df_pts = pd.DataFrame(res["chosen_points"], columns=["x_i (Gốc)", "y_i (Gốc)"])
    lines.append(df_pts.to_string(index=True))

    lines.append("\n3. BẢNG TỶ SAI PHÂN VỚI BIẾN y (DIVIDED DIFFERENCES TABLE):")
    lines.append(res["df_diff"].to_string(index=False))

    lines.append("\n4. ĐA THỨC NỘI SUY NGƯỢC x = P(y):")
    lines.append(res["df_coeffs"].to_string(index=False))
    lines.append(f"\n-> Dạng chính tắc đa thức nội suy ngược:\n   P(y) = {res['poly_str']}")

    lines.append("\n5. KẾT QUẢ TÍNH TOÁN NGHIỆM TẠI ĐIỂM y* (eval_y):")
    if res["eval_y"] is not None:
        lines.append(f"   * Giá trị cần tìm nghiệm: y* = {res['eval_y']}")
        lines.append(f"   * Nghiệm x* tìm được:     x* = P({res['eval_y']}) = {res['val_x']:.10f}")
        if res["derivs"]:
            lines.append("   * Đạo hàm của hàm ngược x = f^(-1)(y) tại y*:")
            for order, dval in res["derivs"]:
                lines.append(f"     + Đạo hàm cấp {order}: (x)^{{({order})}}({res['eval_y']}) = {dval:.10f}")
    else:
        lines.append("   [Lưu ý] Không có giá trị y* (eval_y) trong input. Chỉ xây dựng đa thức tổng quát.")

    lines.append("=" * 80)
    return "\n".join(lines)

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy hàm ngược - Hoán đổi vai trò biến số")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-y", "--eval_y", type=float, help="Điểm cần tính ngược y*")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức nội suy ngược")
    p.add_argument("-n", "--n_points", help="Số mốc nội suy (ví dụ 4, 5 hoặc 'all')")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        # Dữ liệu mặc định kiểm thử
        points = [(1.3, 2.3170624), (1.4, 2.4401799), (1.5, 2.5529043), (1.6, 2.6551084), (1.7, 2.7467702)]

    eval_y = args.eval_y if args.eval_y is not None else params.get("eval_y", None)
    degree_cfg = args.degree or params.get("degree", None)
    n_pts_cfg = args.n_points or params.get("n_points", None)

    res = solve_inverse_function(points, eval_y=eval_y, degree=degree_cfg, n_points=n_pts_cfg)
    final_output = format_inverse_report(res)

    print("\n" + final_output)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_output + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
