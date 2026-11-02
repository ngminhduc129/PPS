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
    from pps_io import format_report_output, parse_universal_input, select_window_for_method
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pps_io import format_report_output, parse_universal_input, select_window_for_method

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

def build_divided_difference_coefficients(points):
    """Return Newton divided-difference coefficients in the given node order."""
    x_values = [point[0] for point in points]
    differences = [point[1] for point in points]
    coefficients = [differences[0]]

    for order in range(1, len(points)):
        differences = [
            (differences[i + 1] - differences[i])
            / (x_values[i + order] - x_values[i])
            for i in range(len(differences) - 1)
        ]
        coefficients.append(differences[0])

    return coefficients

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
    coefficient_points = points if is_forward else list(reversed(points))
    divided_difference_coefficients = build_divided_difference_coefficients(
        coefficient_points
    )

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
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            t_curr = phi(t_prev)
        if not math.isfinite(t_curr):
            history.append({
                "Bước k": it,
                "t^(k)": t_curr,
                "x^(k)": np.nan,
                "|t^(k) - t^(k-1)|": np.nan,
                "|x^(k) - x^(k-1)|": np.nan
            })
            break
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
        "divided_difference_nodes": [point[0] for point in coefficient_points],
        "divided_difference_coefficients": divided_difference_coefficients,
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

def _bracketing_run(points, bracket_index):
    direction = np.sign(points[bracket_index + 1][1] - points[bracket_index][1])
    start = bracket_index
    end = bracket_index + 1

    while start > 0:
        previous_direction = np.sign(points[start][1] - points[start - 1][1])
        if previous_direction != direction:
            break
        start -= 1

    while end < len(points) - 1:
        next_direction = np.sign(points[end + 1][1] - points[end][1])
        if next_direction != direction:
            break
        end += 1

    return start, end

def _select_bracket_stencil(run_points, local_bracket_index, n_points):
    if n_points == "all":
        return run_points

    count = min(max(2, int(n_points)), len(run_points))
    start = max(0, min(
        local_bracket_index - (count - 2) // 2,
        len(run_points) - count
    ))
    return run_points[start:start + count]

def solve_all_loop_inverse(
    points, eval_y, n_points=4, direction="auto", tol=1e-8, max_iter=100
):
    """Solve each data-bracketed root independently on its monotone branch."""
    if len(points) < 2:
        raise ValueError("Cần ít nhất 2 mốc dữ liệu để lặp nội suy ngược.")

    x_steps = [points[i + 1][0] - points[i][0] for i in range(len(points) - 1)]
    h = x_steps[0]
    spacing_tol = 1e-10 * max(1.0, abs(h))
    if h == 0 or any(abs(step - h) > spacing_tol for step in x_steps[1:]):
        raise ValueError("Phương pháp lặp sai phân hữu hạn yêu cầu các mốc X cách đều.")

    y_scale = max(1.0, max(abs(point[1]) for point in points), abs(eval_y))
    y_tol = 1e-12 * y_scale
    solutions = []
    exact_indices = [
        i for i, (_, y_value) in enumerate(points)
        if abs(y_value - eval_y) <= y_tol
    ]

    for index in exact_indices:
        solutions.append({
            "kind": "exact",
            "bracket": (points[index], points[index]),
            "root": points[index][0],
            "result": None,
            "fallback": False
        })

    for bracket_index in range(len(points) - 1):
        y_left = points[bracket_index][1]
        y_right = points[bracket_index + 1][1]
        if not ((y_left - eval_y) * (y_right - eval_y) < 0):
            continue

        run_start, run_end = _bracketing_run(points, bracket_index)
        run_points = points[run_start:run_end + 1]
        stencil = _select_bracket_stencil(
            run_points, bracket_index - run_start, n_points
        )
        x_left, x_right = sorted((
            points[bracket_index][0],
            points[bracket_index + 1][0]
        ))

        direction_results = []
        directions = (
            (True, False) if direction == "auto"
            else (direction == "forward",)
        )
        for is_forward in directions:
            try:
                result = solve_loop_inverse(
                    stencil, eval_y, is_forward=is_forward,
                    tol=tol, max_iter=max_iter
                )
            except (ValueError, OverflowError, ZeroDivisionError) as error:
                direction_results.append({
                    "is_forward": is_forward,
                    "result": None,
                    "error": str(error),
                    "valid": False
                })
                continue

            root = float(result["x_final"])
            valid = (
                result["converged"]
                and math.isfinite(root)
                and x_left - spacing_tol <= root <= x_right + spacing_tol
                and math.isfinite(float(result["res_error"]))
                and result["res_error"] <= max(10 * tol, 1e-10) * y_scale
            )
            direction_results.append({
                "is_forward": is_forward,
                "result": result,
                "error": None,
                "valid": valid
            })

        solutions.append({
            "kind": "bracket",
            "bracket": (points[bracket_index], points[bracket_index + 1]),
            "run_indices": (run_start, run_end),
            "direction_results": direction_results
        })

    solutions.sort(key=lambda item: item["bracket"][0][0])
    return solutions

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

    coeffs = res["divided_difference_coefficients"]
    coeff_vector = ", ".join(f"{coefficient:.10g}" for coefficient in coeffs)
    nodes = ", ".join(f"{node:.10g}" for node in res["divided_difference_nodes"])
    lines.append("\n3. VECTOR HỆ SỐ TỶ SAI PHÂN NEWTON:")
    lines.append(f"   * Thứ tự mốc x_i: [{nodes}]")
    lines.append(f"   * Vector [f[x_0], f[x_0,x_1], ..., f[x_0,...,x_n]]:")
    lines.append(f"     [{coeff_vector}]")

    lines.append("\n4. THIẾT LẬP CÔNG THỨC LẶP:")
    lines.append(f"   * Giá trị xấp xỉ ban đầu: t^(0) = {res['t0']:.10f}")
    lines.append(f"   * Công thức hàm lặp:      {res['phi_formula']}")

    lines.append("\n5. BẢNG QUÁ TRÌNH LẶP TÌM NGHIỆM:")
    lines.append(res["df_history"].to_string(index=False))

    status_str = f"HỘI TỤ sau {res['iterations']} bước lặp" if res["converged"] else f"CHƯA HỘI TỤ sau {res['iterations']} bước"
    lines.append(f"\n-> Tình trạng: {status_str}")

    lines.append("\n6. KẾT QUẢ NGHIỆM VÀ ĐÁNH GIÁ ĐỘ CHÍNH XÁC:")
    lines.append(f"   * Giá trị tham số hội tụ: t* = {res['t_final']:.10f}")
    lines.append(f"   * NGHIỆM x* TÌM ĐƯỢC:     x* = {res['x_final']:.10f}")
    lines.append(f"   * Thử lại: P(x*) = {res['check_y']:.10f}")
    lines.append(f"   * Sai số dư |P(x*) - y*| = {res['res_error']:.2e}")
    lines.append("=" * 80)
    return "\n".join(lines)

def format_all_loop_report(solutions, eval_y):
    lines = [
        "=" * 80,
        "BÁO CÁO LẶP NỘI SUY NGƯỢC TRÊN TẤT CẢ CÁC KHOẢNG CÁCH LY",
        f"Giá trị mục tiêu: y* = {eval_y}",
        f"Số nghiệm/khoảng cách ly tìm thấy: {len(solutions)}",
        "\n1. TỔNG HỢP CÁC KHOẢNG KẸP VÀ NGHIỆM:"
    ]

    if not solutions:
        lines.append("Không tìm thấy khoảng kẹp nào cho y* trong dữ liệu.")

    for index, solution in enumerate(solutions, 1):
        (x_left, y_left), (x_right, y_right) = solution["bracket"]
        if solution["kind"] == "exact":
            lines.append(
                f"{index}. Nghiệm tại nút dữ liệu: x* = {solution['root']:.10g}, "
                f"(x_i, y_i) = ({x_left:.10g}, {y_left:.10g})."
            )
            continue

        direction_summaries = []
        for attempt in solution["direction_results"]:
            name = "Newton tiến" if attempt["is_forward"] else "Newton lùi"
            result = attempt["result"]
            if result is None:
                status = f"lỗi: {attempt['error']}"
            elif not math.isfinite(float(result["x_final"])):
                status = "không cho nghiệm hữu hạn"
            elif not result["converged"]:
                status = f"chưa hội tụ, x ≈ {result['x_final']:.10g}"
            elif not (x_left <= result["x_final"] <= x_right):
                status = f"hội tụ ngoài đoạn, x ≈ {result['x_final']:.10g}"
            elif not attempt["valid"]:
                status = (
                    f"sai số dư lớn, x ≈ {result['x_final']:.10g}, "
                    f"|P(x)-y*| = {result['res_error']:.2e}"
                )
            else:
                status = (
                    f"x ≈ {result['x_final']:.10g} "
                    f"({result['iterations']} bước)"
                )
            direction_summaries.append(f"{name}: {status}")
        lines.append(
            f"{index}. Đoạn [{x_left:.10g}, {x_right:.10g}], "
            f"Y = [{y_left:.10g}, {y_right:.10g}]: "
            + "; ".join(direction_summaries)
        )

    for index, solution in enumerate(solutions, 1):
        if solution["kind"] == "exact":
            continue
        for attempt in solution["direction_results"]:
            result = attempt["result"]
            if result is None:
                continue
            method = "Newton Tiến" if attempt["is_forward"] else "Newton Lùi"
            lines.extend([
                "",
                f"{'=' * 20} KHOẢNG {index} - {method} {'=' * 20}",
                format_loop_report(result),
                f"Kết quả có nằm trong khoảng kẹp và đạt tiêu chuẩn không: "
                f"{'Có' if attempt['valid'] else 'Không'}."
            ])

    return "\n".join(lines)

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Phương pháp Lặp nội suy ngược")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-y", "--eval_y", type=float, help="Điểm cần tính ngược y*")
    p.add_argument("-n", "--n_points", help="Số mốc lân cận dùng để lặp (mặc định 4)")
    p.add_argument("--type", choices=["forward", "backward", "auto"], default="auto", help="Mặc định chạy cả Newton tiến và lùi; chọn hướng nếu chỉ muốn chạy một phương pháp")
    p.add_argument("--tol", type=float, default=1e-8, help="Sai số hội tụ (tolerance)")
    p.add_argument("--max_iter", type=int, default=100, help="Số lần lặp tối đa cho mỗi khoảng cách ly")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        p.error("Không đọc được dữ liệu X, Y. Hãy nhập X và Y hoặc x_interval và Y trong input.txt.")

    eval_y = args.eval_y if args.eval_y is not None else params.get("eval_y", -1.43)
    n_pts_cfg = args.n_points or params.get("n_points", 4)
    solutions = solve_all_loop_inverse(
        points, eval_y, n_points=n_pts_cfg, direction=args.type,
        tol=args.tol, max_iter=args.max_iter
    )
    final_output = format_report_output(
        format_all_loop_report(solutions, eval_y), args.input_file
    )

    print(final_output)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_output + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
