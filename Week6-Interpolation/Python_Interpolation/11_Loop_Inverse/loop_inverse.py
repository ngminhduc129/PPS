"""
Phương pháp Lặp Nội suy Ngược (Inverse Interpolation by Iteration)
Dùng đa thức nội suy sai phân hữu hạn (Newton tiến / lùi) để đưa về dạng lặp t = phi(t), từ đó tìm nghiệm x* sao cho f(x*) = y*.
1. BẢNG MỐC DỮ LIỆU ĐƯỢC CHỌN VÀ BƯỚC LƯỚI h
2. BẢNG SAI PHÂN HỮU HẠN TIẾN/LÙI
3. CÔNG THỨC HÀM LẶP phi(t) VÀ ĐIỂM XẤP XỈ ĐẦU t^(0)
4. BẢNG CHI TIẾT CÁC BƯỚC LẶP (k, t_k, x_k, diff) VÀ ĐÁNH GIÁ HỘI TỤ
5. KẾT QUẢ NGHIỆM x* VÀ KIỂM TRA LẠI P(x*) ≈ y*

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

def build_finite_diff_table(points):
    """
    Xây dựng bảng sai phân hữu hạn tiến Δ cho các mốc cách đều.
    """
    xv = [p[0] for p in points]
    yv = [p[1] for p in points]
    n = len(yv)
    tbl = np.full((n, n), np.nan, dtype=float)
    tbl[:, 0] = yv

    for j in range(1, n):
        for i in range(n - j):
            tbl[i, j] = tbl[i + 1, j - 1] - tbl[i, j - 1]

    data = {"x_i": xv, "y_i": yv}
    for j in range(1, n):
        col = np.full(n, np.nan)
        col[:(n - j)] = tbl[:(n - j), j]
        data[f"Δ^{j} y"] = col
    df_tbl = pd.DataFrame(data)
    return df_tbl, tbl

def eval_newton_forward_t(t, y0, d_coeffs, h=1.0):
    """
    Tính giá trị đa thức Newton tiến tại tham số t.
    d_coeffs = [Δy0, Δ^2 y0, Δ^3 y0, ...]
    """
    val = y0
    prod = 1.0
    for i, dy in enumerate(d_coeffs, 1):
        prod *= (t - (i - 1))
        val += (dy / math.factorial(i)) * prod
    return val

def eval_newton_backward_t(t, yn, d_coeffs, h=1.0):
    """
    Tính giá trị đa thức Newton lùi tại tham số t (t = (x - x_n)/h <= 0).
    d_coeffs = [Δy_{n-1}, Δ^2 y_{n-2}, ...]
    """
    val = yn
    prod = 1.0
    for i, dy in enumerate(d_coeffs, 1):
        prod *= (t + (i - 1))
        val += (dy / math.factorial(i)) * prod
    return val

def solve_loop_inverse(points, eval_y, is_forward=True, tol=1e-8, max_iter=30):
    """
    Thực hiện phương pháp lặp nội suy ngược.
    is_forward: True nếu dùng Newton tiến (thích hợp khi y* gần đầu), False nếu dùng Newton lùi.
    """
    n_pts = len(points)
    if n_pts < 2:
        raise ValueError("Cần ít nhất 2 mốc dữ liệu để lặp nội suy ngược.")

    h = points[1][0] - points[0][0]
    df_diff, tbl = build_finite_diff_table(points)

    if is_forward:
        # Newton tiến: x = x_0 + t*h
        x0 = points[0][0]
        y0 = points[0][1]
        c1 = tbl[0, 1]  # Δy_0
        if abs(c1) < 1e-14:
            raise ValueError("Sai phân bậc 1 Δy_0 gần bằng 0, không thể chia để lập hàm lặp.")

        # Các hệ số sai phân bậc cao C_i = Δ^i y_0
        higher_diffs = [tbl[0, i] for i in range(2, n_pts)]

        # Điểm xuất phát t^(0)
        t0 = (eval_y - y0) / c1

        def phi(t_val):
            correction = 0.0
            prod = 1.0
            for i_idx, delta_i in enumerate(higher_diffs, 2):
                prod = 1.0
                for j in range(i_idx):
                    prod *= (t_val - j)
                fact_i = math.factorial(i_idx)
                correction += (delta_i / (fact_i * c1)) * prod
            return t0 - correction

        # Chuỗi giải thích công thức lặp
        phi_formula_parts = [f"t^(0) = (y* - y_0) / Δy_0 = ({eval_y} - {y0:.6f}) / {c1:.6f} = {t0:.8f}"]
        corr_strs = []
        for i_idx, delta_i in enumerate(higher_diffs, 2):
            fact_i = math.factorial(i_idx)
            coef_str = f"({delta_i:.6f} / ({fact_i} * {c1:.6f}))"
            poly_term = "*".join([f"(t - {j})" for j in range(i_idx)])
            corr_strs.append(f"{coef_str} * {poly_term}")
        phi_formula = f"phi(t) = {t0:.8f} - [ " + " + ".join(corr_strs) + " ]" if corr_strs else f"phi(t) = {t0:.8f}"

    else:
        # Newton lùi: x = x_n + t*h (với t <= 0)
        xn = points[-1][0]
        yn = points[-1][1]
        c1 = tbl[n_pts - 2, 1]  # Δy_{n-1}
        if abs(c1) < 1e-14:
            raise ValueError("Sai phân bậc 1 Δy_{n-1} gần bằng 0, không thể chia để lập hàm lặp.")

        higher_diffs = [tbl[n_pts - 1 - i, i] for i in range(2, n_pts)]
        t0 = (eval_y - yn) / c1

        def phi(t_val):
            correction = 0.0
            for i_idx, delta_i in enumerate(higher_diffs, 2):
                prod = 1.0
                for j in range(i_idx):
                    prod *= (t_val + j)
                fact_i = math.factorial(i_idx)
                correction += (delta_i / (fact_i * c1)) * prod
            return t0 - correction

        phi_formula_parts = [f"t^(0) = (y* - y_n) / Δy_{{n-1}} = ({eval_y} - {yn:.6f}) / {c1:.6f} = {t0:.8f}"]
        corr_strs = []
        for i_idx, delta_i in enumerate(higher_diffs, 2):
            fact_i = math.factorial(i_idx)
            coef_str = f"({delta_i:.6f} / ({fact_i} * {c1:.6f}))"
            poly_term = "*".join([f"(t + {j})" for j in range(i_idx)])
            corr_strs.append(f"{coef_str} * {poly_term}")
        phi_formula = f"phi(t) = {t0:.8f} - [ " + " + ".join(corr_strs) + " ]" if corr_strs else f"phi(t) = {t0:.8f}"

    # Thực hiện các bước lặp
    t_curr = t0
    history = []
    base_x = points[0][0] if is_forward else points[-1][0]
    x_curr = base_x + t_curr * h
    history.append({
        "Bước k": 0,
        "t^(k)": t_curr,
        "x^(k)": x_curr,
        "|t^(k) - t^(k-1)|": np.nan,
        "|x^(k) - x^(k-1)|": np.nan
    })

    converged = False
    for it in range(1, max_iter + 1):
        t_prev = t_curr
        t_curr = phi(t_prev)
        diff_t = abs(t_curr - t_prev)
        x_curr = base_x + t_curr * h
        diff_x = diff_t * h

        history.append({
            "Bước k": it,
            "t^(k)": t_curr,
            "x^(k)": x_curr,
            "|t^(k) - t^(k-1)|": diff_t,
            "|x^(k) - x^(k-1)|": diff_x
        })

        if diff_t < tol:
            converged = True
            break

    df_history = pd.DataFrame(history)
    x_final = base_x + t_curr * h

    # Tính kiểm tra lại giá trị P(x*)
    if is_forward:
        diff_list = [tbl[0, j] for j in range(1, n_pts)]
        check_y = eval_newton_forward_t(t_curr, points[0][1], diff_list, h)
    else:
        diff_list = [tbl[n_pts - 1 - j, j] for j in range(1, n_pts)]
        check_y = eval_newton_backward_t(t_curr, points[-1][1], diff_list, h)

    res_error = abs(check_y - eval_y)

    return {
        "points": points,
        "h": h,
        "eval_y": eval_y,
        "is_forward": is_forward,
        "df_diff": df_diff,
        "t0": t0,
        "phi_formula": phi_formula,
        "df_history": df_history,
        "converged": converged,
        "t_final": t_curr,
        "x_final": x_final,
        "check_y": check_y,
        "res_error": res_error,
        "iterations": len(history) - 1
    }

def format_loop_report(res):
    f_type = "Newton Tiến (Forward)" if res["is_forward"] else "Newton Lùi (Backward)"
    lines = [
        "=" * 80,
        "BÁO CÁO PHƯƠNG PHÁP LẶP NỘI SUY NGƯỢC (INVERSE INTERPOLATION LOOP)",
        "=" * 80,
        f"Phương pháp sử dụng: Lặp theo đa thức {f_type}",
        f"Bước lưới mốc:      h = {res['h']:.8g}",
        f"Giá trị mục tiêu:   y* = {res['eval_y']}",
        "\n1. BẢNG MỐC DỮ LIỆU ĐƯỢC CHỌN ĐỂ LẶP NỘI SUY:"
    ]
    df_pts = pd.DataFrame(res["points"], columns=["x_i", "y_i"])
    lines.append(df_pts.to_string(index=True))

    lines.append("\n2. BẢNG SAI PHÂN HỮU HẠN:")
    lines.append(res["df_diff"].to_string(index=True))

    lines.append("\n3. THIẾT LẬP CÔNG THỨC LẶP:")
    lines.append(f"   * Giá trị xấp xỉ ban đầu: t^(0) = {res['t0']:.10f}")
    lines.append(f"   * Công thức hàm lặp:      {res['phi_formula']}")

    lines.append("\n4. BẢNG QUÁ TRÌNH LẶP TÌM NGHIỆM:")
    lines.append(res["df_history"].to_string(index=False))

    status_str = f"HỘI TỤ sau {res['iterations']} bước lặp" if res["converged"] else f"CHƯA HỘI TỤ sau {res['iterations']} bước"
    lines.append(f"\n-> Tình trạng: {status_str}")

    lines.append("\n5. KẾT QUẢ NGHIỆM VÀ ĐÁNH GIÁ ĐỘ CHÍNH XÁC:")
    lines.append(f"   * Giá trị tham số hội tụ: t* = {res['t_final']:.10f}")
    lines.append(f"   * NGHIỆM x* TÌM ĐƯỢC:     x* = {res['x_final']:.10f}")
    lines.append(f"   * Thử lại: P(x*) = {res['check_y']:.10f}")
    lines.append(f"   * Sai số dư |P(x*) - y*| = {res['res_error']:.2e}")
    lines.append("=" * 80)
    return "\n".join(lines)

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Phương pháp Lặp nội suy ngược")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-y", "--eval_y", type=float, help="Điểm cần tính ngược y*")
    p.add_argument("-n", "--n_points", help="Số mốc lân cận dùng để lặp (mặc định 4)")
    p.add_argument("--type", choices=["forward", "backward", "auto"], default="auto", help="Hướng lặp sai phân")
    p.add_argument("--tol", type=float, default=1e-8, help="Sai số hội tụ (tolerance)")
    p.add_argument("--max_iter", type=int, default=30, help="Số lần lặp tối đa")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(2.3, -1.9774), (2.4, -1.7483), (2.5, -1.5366), (2.6, -1.3411), (2.7, -1.1607)]

    eval_y = args.eval_y if args.eval_y is not None else params.get("eval_y", -1.43)
    n_pts_cfg = args.n_points or params.get("n_points", 4)

    # Cắt cửa sổ lân cận quanh eval_y
    sub_points, _, msg = select_window_for_method(points, eval_x=None, method="loop_inverse", n_points=n_pts_cfg, eval_y=eval_y)

    # Tự động chọn forward hay backward nếu để auto
    if args.type == "auto":
        # Thử cả 2 hướng và ưu tiên kết quả hội tụ
        res_fwd = None
        res_bwd = None
        try:
            res_fwd = solve_loop_inverse(sub_points, eval_y, is_forward=True, tol=args.tol, max_iter=args.max_iter)
        except Exception:
            pass
        try:
            res_bwd = solve_loop_inverse(sub_points, eval_y, is_forward=False, tol=args.tol, max_iter=args.max_iter)
        except Exception:
            pass

        if res_fwd and res_fwd["converged"]:
            res = res_fwd
        elif res_bwd and res_bwd["converged"]:
            res = res_bwd
        elif res_fwd:
            res = res_fwd
        else:
            res = res_bwd
    else:
        is_fwd = (args.type == "forward")
        res = solve_loop_inverse(sub_points, eval_y, is_forward=is_fwd, tol=args.tol, max_iter=args.max_iter)
    final_output = format_loop_report(res)

    print("\n" + msg)
    print(final_output)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(msg + "\n" + final_output + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
