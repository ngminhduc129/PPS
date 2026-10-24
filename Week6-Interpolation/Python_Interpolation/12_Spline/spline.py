"""
Nội suy Ghép trơn (Spline Interpolation)
Bao gồm đầy đủ 3 loại ghép trơn kinh điển:
1. Linear Spline (Ghép trơn bậc 1)
2. Quadratic Spline (Ghép trơn bậc 2 với điều kiện biên đạo hàm z0)
3. Cubic Spline (Ghép trơn bậc 3 tự nhiên / Natural Cubic Spline giải qua hệ 3 đường chéo Thomas)

Xuất kết quả đúng chuẩn:
1. BẢNG MỐC DỮ LIỆU VÀ KHOẢNG BƯỚC h_i
2. BẢNG HỆ SỐ GÓC / MOMENT TẠI CÁC NÚT (z_i hoặc M_i)
3. BẢNG HỆ SỐ CỦA CÁC ĐOẠN SPLINE S_i(x)
4. DẠNG GIẢI TÍCH CÁC HÀM TỪNG ĐOẠN S_i(x) THEO BIẾN x
5. KẾT QUẢ NỘI SUY VÀ CÁC CẤP ĐẠO HÀM TẠI ĐIỂM x* (eval_x)

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

# ==============================================================================
# 1. LINEAR SPLINE (BẬC 1)
# ==============================================================================
def solve_linear_spline(points, eval_x=None):
    """
    Xây dựng Spline tuyến tính bậc 1: S_i(x) = a_1*x + a_0 trên [x_i, x_{i+1}]
    """
    n = len(points) - 1
    segments = []
    for i in range(n):
        x_i, y_i = points[i]
        x_ip1, y_ip1 = points[i + 1]
        h_i = x_ip1 - x_i
        if h_i <= 0:
            raise ValueError(f"Dữ liệu x không tăng nghiêm ngặt tại chỉ số {i}.")
        a_1 = (y_ip1 - y_i) / h_i
        a_0 = y_i - a_1 * x_i
        formula = f"{a_1:.8g}*x + {a_0:.8g}" if a_0 >= 0 else f"{a_1:.8g}*x - {abs(a_0):.8g}"
        segments.append({
            "Đoạn i": i,
            "Khoảng [x_i, x_{i+1}]": f"[{x_i:.4f}, {x_ip1:.4f}]",
            "x_i": x_i, "x_ip1": x_ip1, "h_i": h_i,
            "a_0": a_0, "a_1": a_1,
            "Công thức S_i(x)": formula
        })
    df_segments = pd.DataFrame(segments)

    eval_res = None
    if eval_x is not None:
        seg_match = None
        for seg in segments:
            if seg["x_i"] <= eval_x <= (seg["x_ip1"] + 1e-12):
                seg_match = seg
                break
        if seg_match is None:
            # Lấy đoạn gần nhất nếu ở biên
            seg_match = segments[0] if eval_x < segments[0]["x_i"] else segments[-1]
        y_val = seg_match["a_1"] * eval_x + seg_match["a_0"]
        eval_res = {
            "eval_x": eval_x,
            "segment_idx": seg_match["Đoạn i"],
            "interval": seg_match["Khoảng [x_i, x_{i+1}]"],
            "val": y_val,
            "deriv1": seg_match["a_1"],
            "deriv2": 0.0
        }

    return {"type": "linear", "points": points, "df_segments": df_segments, "eval_res": eval_res}

# ==============================================================================
# 2. QUADRATIC SPLINE (BẬC 2)
# ==============================================================================
def solve_quadratic_spline(points, z0_condition=0.0, eval_x=None):
    """
    Xây dựng Spline bậc 2: S_i(x) = a_2*x^2 + a_1*x + a_0 trên [x_i, x_{i+1}]
    với điều kiện biên hệ số góc z_0 = S'(x_0).
    """
    n = len(points) - 1
    x = np.array([p[0] for p in points], dtype=float)
    y = np.array([p[1] for p in points], dtype=float)
    h = np.diff(x)
    z = np.zeros(n + 1, dtype=float)

    z[0] = z0_condition
    z_steps = [{"i": 0, "x_i": x[0], "y_i": y[0], "h_i": np.nan, "z_i = S'(x_i)": z[0], "Ghi chú": "Điều kiện biên (z0)"}]

    for i in range(n):
        z[i + 1] = (2.0 * (y[i + 1] - y[i]) / h[i]) - z[i]
        z_steps.append({
            "i": i + 1,
            "x_i": x[i + 1],
            "y_i": y[i + 1],
            "h_i": h[i],
            "z_i = S'(x_i)": z[i + 1],
            "Ghi chú": f"2*({y[i+1]:.4f} - {y[i]:.4f})/{h[i]:.4f} - {z[i]:.4f}"
        })
    df_z = pd.DataFrame(z_steps)

    segments = []
    for i in range(n):
        a_2 = (z[i + 1] - z[i]) / (2.0 * h[i])
        a_1 = z[i] - 2.0 * a_2 * x[i]
        a_0 = y[i] - a_2 * (x[i]**2) - a_1 * x[i]
        terms = [f"{a_2:.8g}*x^2", f"{'+ ' if a_1>=0 else '- '}{abs(a_1):.8g}*x", f"{'+ ' if a_0>=0 else '- '}{abs(a_0):.8g}"]
        formula = " ".join(terms)
        segments.append({
            "Đoạn i": i,
            "Khoảng [x_i, x_{i+1}]": f"[{x[i]:.4f}, {x[i+1]:.4f}]",
            "x_i": x[i], "x_ip1": x[i + 1], "h_i": h[i],
            "a_0": a_0, "a_1": a_1, "a_2": a_2,
            "Công thức S_i(x)": formula
        })
    df_segments = pd.DataFrame(segments)

    eval_res = None
    if eval_x is not None:
        seg_match = None
        for seg in segments:
            if seg["x_i"] <= eval_x <= (seg["x_ip1"] + 1e-12):
                seg_match = seg
                break
        if seg_match is None:
            seg_match = segments[0] if eval_x < segments[0]["x_i"] else segments[-1]
        y_val = seg_match["a_2"] * (eval_x**2) + seg_match["a_1"] * eval_x + seg_match["a_0"]
        d1 = 2.0 * seg_match["a_2"] * eval_x + seg_match["a_1"]
        d2 = 2.0 * seg_match["a_2"]
        eval_res = {
            "eval_x": eval_x,
            "segment_idx": seg_match["Đoạn i"],
            "interval": seg_match["Khoảng [x_i, x_{i+1}]"],
            "val": y_val,
            "deriv1": d1,
            "deriv2": d2
        }

    return {"type": "quadratic", "points": points, "df_steps": df_z, "df_segments": df_segments, "eval_res": eval_res}

# ==============================================================================
# 3. CUBIC SPLINE (BẬC 3 TỰ NHIÊN - NATURAL CUBIC SPLINE)
# ==============================================================================
def solve_cubic_spline(points, eval_x=None):
    """
    Xây dựng Natural Cubic Spline (Ghép trơn bậc 3 tự nhiên M_0 = M_n = 0)
    S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
    Giải hệ moment M_i = S''(x_i) bằng thuật toán 3 đường chéo Thomas.
    """
    n = len(points) - 1
    x = np.array([p[0] for p in points], dtype=float)
    y = np.array([p[1] for p in points], dtype=float)
    h = np.diff(x)

    if n < 2:
        raise ValueError("Cần ít nhất 3 mốc để xây dựng Cubic Spline.")

    # Ma trận 3 đường chéo cho các nút trong i = 1, ..., n-1
    # h_{i-1} * M_{i-1} + 2*(h_{i-1} + h_i) * M_i + h_i * M_{i+1} = 6 * ( (y_{i+1} - y_i)/h_i - (y_i - y_{i-1})/h_{i-1} )
    num_internal = n - 1
    A = np.zeros((num_internal, num_internal), dtype=float)
    B = np.zeros(num_internal, dtype=float)

    for i in range(1, n):
        row = i - 1
        h_prev = h[i - 1]
        h_curr = h[i]
        delta_curr = (y[i + 1] - y[i]) / h_curr
        delta_prev = (y[i] - y[i - 1]) / h_prev

        if row > 0:
            A[row, row - 1] = h_prev
        A[row, row] = 2.0 * (h_prev + h_curr)
        if row < num_internal - 1:
            A[row, row + 1] = h_curr

        B[row] = 6.0 * (delta_curr - delta_prev)

    # Giải hệ Thomas
    M_internal = np.linalg.solve(A, B)
    M = np.zeros(n + 1, dtype=float)
    M[1:n] = M_internal
    M[0] = 0.0  # Natural spline
    M[n] = 0.0  # Natural spline

    df_moments = pd.DataFrame({
        "i": np.arange(n + 1),
        "x_i": x,
        "y_i": y,
        "Moment M_i = S''(x_i)": M,
        "Điều kiện": ["Biên tự nhiên (M0=0)"] + ["Nút trong"] * (n - 1) + ["Biên tự nhiên (Mn=0)"]
    })

    segments = []
    for i in range(n):
        x_i = x[i]
        x_ip1 = x[i + 1]
        h_i = h[i]

        # Hệ số theo dạng delta x = (x - x_i):
        # S_i(x) = a_i + b_i*(x - x_i) + c_i*(x - x_i)^2 + d_i*(x - x_i)^3
        a_i = y[i]
        c_i = M[i] / 2.0
        d_i = (M[i + 1] - M[i]) / (6.0 * h_i)
        b_i = (y[i + 1] - y[i]) / h_i - (h_i / 6.0) * (2.0 * M[i] + M[i + 1])

        # Khai triển về dạng chính tắc S_i(x) = A3*x^3 + A2*x^2 + A1*x + A0
        # (x - x_i) = x - x0
        p_basis = np.array([-x_i, 1.0])  # [-x_i, 1]
        p_cube = np.convolve(np.convolve(p_basis, p_basis), p_basis) * d_i
        p_sq = np.convolve(p_basis, p_basis) * c_i
        p_lin = p_basis * b_i

        poly_tot = np.zeros(4, dtype=float)
        poly_tot[0] += a_i
        poly_tot[:len(p_lin)] += p_lin
        poly_tot[:len(p_sq)] += p_sq
        poly_tot[:len(p_cube)] += p_cube

        A0, A1, A2, A3 = poly_tot[0], poly_tot[1], poly_tot[2], poly_tot[3]
        formula_std = f"{A3:.8g}*x^3 + {A2:.8g}*x^2 + {A1:.8g}*x + {A0:.8g}"
        formula_local = f"{a_i:.8g} + {b_i:.8g}*(x - {x_i:.4f}) + {c_i:.8g}*(x - {x_i:.4f})^2 + {d_i:.8g}*(x - {x_i:.4f})^3"

        segments.append({
            "Đoạn i": i,
            "Khoảng [x_i, x_{i+1}]": f"[{x_i:.4f}, {x_ip1:.4f}]",
            "x_i": x_i, "x_ip1": x_ip1, "h_i": h_i,
            "a_i": a_i, "b_i": b_i, "c_i": c_i, "d_i": d_i,
            "A3": A3, "A2": A2, "A1": A1, "A0": A0,
            "Dạng theo (x - x_i)": formula_local,
            "Dạng khai triển x": formula_std
        })

    df_segments = pd.DataFrame(segments)

    eval_res = None
    if eval_x is not None:
        seg_match = None
        for seg in segments:
            if seg["x_i"] <= eval_x <= (seg["x_ip1"] + 1e-12):
                seg_match = seg
                break
        if seg_match is None:
            seg_match = segments[0] if eval_x < segments[0]["x_i"] else segments[-1]
        dx = eval_x - seg_match["x_i"]
        y_val = seg_match["a_i"] + seg_match["b_i"] * dx + seg_match["c_i"] * (dx**2) + seg_match["d_i"] * (dx**3)
        d1 = seg_match["b_i"] + 2.0 * seg_match["c_i"] * dx + 3.0 * seg_match["d_i"] * (dx**2)
        d2 = 2.0 * seg_match["c_i"] + 6.0 * seg_match["d_i"] * dx
        eval_res = {
            "eval_x": eval_x,
            "segment_idx": seg_match["Đoạn i"],
            "interval": seg_match["Khoảng [x_i, x_{i+1}]"],
            "val": y_val,
            "deriv1": d1,
            "deriv2": d2
        }

    return {"type": "cubic", "points": points, "df_steps": df_moments, "df_segments": df_segments, "eval_res": eval_res}

# ==============================================================================
# BÁO CÁO FORMAT ĐẸP
# ==============================================================================
def format_spline_report(res):
    stype = res["type"].upper()
    lines = [
        "=" * 80,
        f"BÁO CÁO NỘI SUY GHÉP TRƠN {stype} SPLINE INTERPOLATION",
        "=" * 80,
        f"Tổng số mốc dữ liệu: {len(res['points'])} điểm ({len(res['df_segments'])} đoạn ghép trơn).",
        "\n1. BẢNG MỐC DỮ LIỆU ĐẦU VÀO:"
    ]
    df_pts = pd.DataFrame(res["points"], columns=["x_i", "y_i"])
    lines.append(df_pts.to_string(index=True))

    if "df_steps" in res:
        title_step = "BẢNG HỆ SỐ GÓC z_i TẠI CÁC NÚT (S'(x_i)):" if res["type"] == "quadratic" else "BẢNG MOMENT M_i = S''(x_i) TẠI CÁC NÚT (HỆ 3 ĐƯỜNG CHÉO):"
        lines.append(f"\n2. {title_step}")
        lines.append(res["df_steps"].to_string(index=False))

    lines.append("\n3. BẢNG TỔNG HỢP CÁC ĐOẠN GHÉP TRƠN S_i(x):")
    if res["type"] == "linear":
        cols = ["Đoạn i", "Khoảng [x_i, x_{i+1}]", "a_0", "a_1", "Công thức S_i(x)"]
    elif res["type"] == "quadratic":
        cols = ["Đoạn i", "Khoảng [x_i, x_{i+1}]", "a_0", "a_1", "a_2", "Công thức S_i(x)"]
    else:
        cols = ["Đoạn i", "Khoảng [x_i, x_{i+1}]", "a_i", "b_i", "c_i", "d_i", "Dạng khai triển x"]
    lines.append(res["df_segments"][cols].to_string(index=False))

    lines.append("\n4. CHI TIẾT CÔNG THỨC TỪNG ĐOẠN:")
    for _, r in res["df_segments"].iterrows():
        i = r["Đoạn i"]
        inv = r["Khoảng [x_i, x_{i+1}]"]
        if res["type"] == "cubic":
            lines.append(f"   * Đoạn {i} trên {inv}:")
            lines.append(f"     + Dạng cục bộ:  S_{i}(x) = {r['Dạng theo (x - x_i)']}")
            lines.append(f"     + Khai triển x: S_{i}(x) = {r['Dạng khai triển x']}")
        else:
            lines.append(f"   * Đoạn {i} trên {inv}: S_{i}(x) = {r['Công thức S_i(x)']}")

    lines.append("\n5. KẾT QUẢ NỘI SUY TẠI ĐIỂM x* (eval_x):")
    if res["eval_res"] is not None:
        er = res["eval_res"]
        lines.append(f"   * Điểm cần tính:      x* = {er['eval_x']}")
        lines.append(f"   * Thuộc đoạn nội suy: {er['interval']} (Đoạn #{er['segment_idx']})")
        lines.append(f"   >>> GIÁ TRỊ NỘI SUY:  S({er['eval_x']}) = {er['val']:.10f} <<<")
        lines.append(f"   * Đạo hàm cấp 1:      S'({er['eval_x']}) = {er['deriv1']:.10f}")
        lines.append(f"   * Đạo hàm cấp 2:      S''({er['eval_x']}) = {er['deriv2']:.10f}")
    else:
        lines.append("   [Lưu ý] Không có điểm x* (eval_x) trong input. Đã xây dựng toàn bộ biểu thức Spline.")

    lines.append("=" * 80)
    return "\n".join(lines)

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Ghép trơn Spline (Linear, Quadratic, Cubic)")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("--type", choices=["linear", "quadratic", "cubic", "all"], default="cubic", help="Loại spline")
    p.add_argument("--z0", type=float, default=0.0, help="Điều kiện biên z0 = S'(x0) cho Quadratic Spline")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(1.4, 2.4347052), (1.5, 2.5473627), (1.6, 2.6495552), (1.7, 2.7412610),
                  (1.8, 2.8225627), (1.9, 2.8936474), (2.0, 2.9548000), (2.1, 3.0064203), (2.2, 3.0489801)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", 1.75)
    types_to_run = ["linear", "quadratic", "cubic"] if args.type == "all" else [args.type]

    output_reports = []
    for st in types_to_run:
        if st == "linear":
            res = solve_linear_spline(points, eval_x=eval_x)
        elif st == "quadratic":
            res = solve_quadratic_spline(points, z0_condition=args.z0, eval_x=eval_x)
        else:
            res = solve_cubic_spline(points, eval_x=eval_x)
        output_reports.append(format_spline_report(res))

    full_output = "\n\n".join(output_reports)
    print("\n" + full_output)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(full_output + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
