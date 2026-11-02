"""
Nội suy Newton sai phân hữu hạn (mốc cách đều)
Tự động in TẤT CẢ các trường hợp:
- PHẦN A: Newton Sai phân hữu hạn TIẾN (Forward Finite Difference - gốc x0, biến t)
- PHẦN B: Newton Sai phân hữu hạn LÙI (Backward Finite Difference - gốc xn, biến tb)
- PHẦN C: TỔNG KẾT & SO SÁNH (Quy đổi về biến x, kiểm chứng tính duy nhất & các cấp đạo hàm)

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

def build_finite_diff_tables(points):
    """
    Xây dựng bảng sai phân hữu hạn mốc cách đều:
    - Bảng tiến Δ (tam giác trên)
    - Bảng lùi ∇ (tam giác dưới)
    """
    xv = [p[0] for p in points]
    yv = [p[1] for p in points]
    n = len(yv)

    tbl_upper = np.full((n, n), np.nan, dtype=float)
    tbl_upper[:, 0] = yv
    for j in range(1, n):
        for i in range(n - j):
            tbl_upper[i, j] = tbl_upper[i + 1, j - 1] - tbl_upper[i, j - 1]

    data_up = {"i": list(range(n)), "x_i": xv, "y_i": yv}
    for j in range(1, n):
        data_up[f"Δ^{j} y"] = tbl_upper[:, j]
    df_upper = pd.DataFrame(data_up)

    tbl_lower = np.full((n, n), np.nan, dtype=float)
    for i in range(n):
        for j in range(i + 1):
            tbl_lower[i, j] = tbl_upper[i - j, j]
    data_low = {"i": list(range(n)), "x_i": xv, "y_i": tbl_lower[:, 0]}
    for j in range(1, n):
        data_low[f"∇^{j} y"] = tbl_lower[:, j]
    df_lower = pd.DataFrame(data_low)

    return df_upper, tbl_upper, df_lower, tbl_lower

def poly_str(coeffs, var="t"):
    terms = []
    for deg, c in enumerate(coeffs):
        if abs(c) < 1e-12: continue
        sign = "+ " if c >= 0 and len(terms) > 0 else ("- " if c < 0 and len(terms) > 0 else "")
        val = abs(c) if len(terms) > 0 else c
        if deg == 0: term = f"{val:.8g}"
        elif deg == 1: term = f"{val:.8g}*{var}"
        else: term = f"{val:.8g}*{var}^{deg}"
        terms.append(f"{sign}{term}")
    return " ".join(terms) if terms else "0"

def convert_pt_to_px(p_t_coeffs, x_base, h):
    """
    Quy đổi đa thức P(t) = sum(a_k * t^k) thành P(x) = sum(c_k * x^k)
    với t = (x - x_base) / h = (1/h)*x - (x_base/h)
    """
    m = len(p_t_coeffs)
    p_x = np.zeros(m, dtype=float)
    t_lin = np.array([-x_base / h, 1.0 / h], dtype=float)
    t_pow = np.array([1.0], dtype=float)
    for k in range(m):
        contrib = p_t_coeffs[k] * t_pow
        p_x[:len(contrib)] += contrib
        t_pow = np.convolve(t_pow, t_lin)
    return p_x

def compute_all_derivatives(coeffs, x_val):
    deriv_list = []
    curr_poly = np.array(coeffs, dtype=float)
    n = len(coeffs) - 1
    for k in range(1, n + 1):
        m = len(curr_poly) - 1
        if m < 1: break
        next_poly = np.array([(j + 1) * curr_poly[j + 1] for j in range(m)], dtype=float)
        val_k = float(sum(c * (x_val ** j) for j, c in enumerate(next_poly))) if x_val is not None else 0.0
        deriv_list.append((k, next_poly, val_k))
        curr_poly = next_poly
    return deriv_list

def run_newton_finite_single(points, tbl_upper, tbl_lower, mode="forward", eval_x=None, section_prefix=""):
    """
    Chạy nội suy Newton sai phân hữu hạn cho 1 chiều:
    - forward: mốc gốc x0 = points[0][0], t = (x - x0)/h
    - backward: mốc gốc xn = points[-1][0], tb = (x - xn)/h
    """
    xv = [p[0] for p in points]
    yv = [p[1] for p in points]
    n_pts = len(points)
    n = n_pts - 1
    h = xv[1] - xv[0]
    is_fwd = (mode == "forward")
    tname = "TIẾN (FORWARD)" if is_fwd else "LÙI (BACKWARD)"
    var_name = "t" if is_fwd else "t_b"

    if is_fwd:
        x_base = xv[0]
        c_coeffs = tbl_upper[0, :].tolist()
        trace_desc = "Lấy hàng đầu tiên của bảng sai phân tiến Δ: y_0, Δ y_0, Δ^2 y_0, ..."
    else:
        x_base = xv[-1]
        c_coeffs = [tbl_upper[n - j, j] for j in range(n + 1)]
        trace_desc = "Lấy hàng cuối cùng của bảng sai phân lùi ∇: y_n, ∇ y_n, ∇^2 y_n, ..."

    # Xây dựng các đa thức cơ sở tích B_i(t)
    # Forward: B_0 = 1, B_1 = t, B_i = B_{i-1} * (t - (i-1))
    # Backward: B_0 = 1, B_1 = t_b, B_i = B_{i-1} * (t_b + (i-1))
    B_list = [np.array([1.0])]
    for i in range(1, n + 1):
        prev = B_list[-1]
        k = i - 1
        factor = [-float(k), 1.0] if is_fwd else [float(k), 1.0]
        Bi = np.convolve(prev, factor)
        B_list.append(Bi)

    # Đa thức P(t) = sum(D_i * B_i(t)) với D_i = C_i / i!
    P_t_total = np.zeros(n + 1, dtype=float)
    step_rows = []
    d_coeffs = []
    for i in range(n + 1):
        Ci = c_coeffs[i]
        fact_i = math.factorial(i)
        Di = Ci / fact_i
        d_coeffs.append(Di)
        Ni = Di * B_list[i]
        P_t_total[:len(Ni)] += Ni

        if i == 0:
            b_expr = "1"
        else:
            factors = [f"({var_name} - {j})" if is_fwd else f"({var_name} + {j})" for j in range(i)]
            # Chuẩn hóa (t - 0) thành t
            factors = [var_name if (f == f"({var_name} - 0)" or f == f"({var_name} + 0)") else f for f in factors]
            b_expr = "*".join(factors)

        step_rows.append({
            "i": i,
            "C_i (Sai phân)": Ci,
            "D_i = C_i / i!": Di,
            f"Biểu thức tích B_i({var_name})": b_expr,
            f"Đa thức số hạng N_i({var_name})": [round(v, 8) for v in Ni.tolist()]
        })
    df_steps = pd.DataFrame(step_rows)

    # Dạng tích chuỗi
    product_terms = []
    for i in range(n + 1):
        Di = d_coeffs[i]
        if abs(Di) < 1e-12: continue
        sign = "+ " if Di >= 0 and len(product_terms) > 0 else ("- " if Di < 0 and len(product_terms) > 0 else "")
        val = abs(Di) if len(product_terms) > 0 else Di
        b_expr = step_rows[i][f"Biểu thức tích B_i({var_name})"]
        if i == 0:
            product_terms.append(f"{val:.8g}")
        else:
            product_terms.append(f"{sign}{val:.8g}*{b_expr}")
    product_form_str = " ".join(product_terms) if product_terms else "0"

    # Quy đổi sang biến x
    P_x_total = convert_pt_to_px(P_t_total, x_base, h)

    # Soạn thảo báo cáo 5 phần
    lines = []
    lines.append("=" * 85)
    lines.append(f"{section_prefix}NỘI SUY NEWTON MỐC CÁCH ĐỀU - {tname}")
    lines.append("=" * 85)
    lines.append(f"Số mốc tính toán: {n_pts}, Bước h: {h:.6f}, Mốc gốc: x_{'0' if is_fwd else 'n'} = {x_base}")
    lines.append(f"Biến chuẩn hóa: {var_name} = (x - {x_base}) / {h:.6f}")
    lines.append(f"Quy tắc chọn hệ số: {trace_desc}\n")

    lines.append(f"1. BẢNG HỆ SỐ NỘI SUY {tname}:")
    lines.append("-" * 85)
    coeff_rows = []
    for i in range(n + 1):
        sym = f"Δ^{i} y_0" if is_fwd else f"∇^{i} y_n"
        coeff_rows.append({
            "Bậc i": i,
            "Ký hiệu sai phân": sym if i > 0 else "y_0" if is_fwd else "y_n",
            "Sai phân C_i": c_coeffs[i],
            "i!": math.factorial(i),
            "Hệ số D_i = C_i/i!": d_coeffs[i],
            f"Biểu thức tích B_i({var_name})": step_rows[i][f"Biểu thức tích B_i({var_name})"]
        })
    lines.append(pd.DataFrame(coeff_rows).to_string(index=False))

    lines.append(f"\n2. BẢNG TÍCH VÀ CÁC BƯỚC DỰNG ĐA THỨC SỐ HẠNG N_i({var_name}) = D_i * B_i({var_name}):")
    lines.append("-" * 85)
    lines.append(df_steps[["i", "C_i (Sai phân)", "D_i = C_i / i!", f"Biểu thức tích B_i({var_name})", f"Đa thức số hạng N_i({var_name})"]].to_string(index=False))

    lines.append(f"\n3. BẢNG HỆ SỐ ĐA THỨC NỘI SUY THEO BIẾN {var_name} VÀ BIẾN x:")
    lines.append("-" * 85)
    df_coeff_comp = pd.DataFrame({
        "Bậc k": np.arange(n + 1),
        f"Hệ số theo {var_name}^k": P_t_total,
        "Hệ số theo x^k": P_x_total
    })
    lines.append(df_coeff_comp.to_string(index=False))

    lines.append(f"\n4. ĐA THỨC NỘI SUY:")
    lines.append("-" * 85)
    lines.append(f"a) Dạng tích Newton (chưa khai triển theo {var_name}):")
    lines.append(f"   P_{mode}({var_name}) = {product_form_str}")
    lines.append(f"\nb) Dạng chính tắc theo biến chuẩn hóa {var_name}:")
    lines.append(f"   P_{mode}({var_name}) = {poly_str(P_t_total, var=var_name)}")
    lines.append(f"\nc) Dạng chính tắc quy đổi theo biến thực x:")
    lines.append(f"   P_{mode}(x) = {poly_str(P_x_total, var='x')}")

    val_final = None
    dict_derivs_x = {}
    lines.append(f"\n5. GIÁ TRỊ NỘI SUY VÀ ĐẠO HÀM TẠI ĐIỂM XÁC ĐỊNH{' (x* = ' + str(eval_x) + ')' if eval_x is not None else ''}:")
    lines.append("-" * 85)
    if eval_x is not None:
        t_eval = (eval_x - x_base) / h
        eval_rows = []
        cum_sum = 0.0
        for i in range(n + 1):
            bi_val = float(sum(c * (t_eval ** deg) for deg, c in enumerate(B_list[i])))
            di = d_coeffs[i]
            ti = di * bi_val
            cum_sum += ti
            eval_rows.append({
                "i": i,
                f"B_i({var_name}*)": bi_val,
                "D_i": di,
                "Số hạng T_i": ti,
                f"Tổng dồn P_i({var_name}*)": cum_sum
            })
        lines.append(f"Giá trị biến chuẩn hóa: {var_name}* = ({eval_x} - {x_base}) / {h:.6f} = {t_eval:.6f}\n")
        lines.append(pd.DataFrame(eval_rows).to_string(index=False))
        val_final = cum_sum
        lines.append(f"\n=> GIÁ TRỊ NỘI SUY:")
        lines.append(f"   P_{mode}({eval_x}) = P_{mode}({var_name}* = {t_eval:.6f}) = {val_final:.10f}")

        # Tính đạo hàm
        derivs_x = compute_all_derivatives(P_x_total, eval_x)
        derivs_t = compute_all_derivatives(P_t_total, t_eval)
        dict_derivs_x = {k: v for k, p, v in derivs_x}
        lines.append(f"\n=> KẾT QUẢ ĐẠO HÀM CÁC CẤP TẠI x* = {eval_x}:")
        d_rows = []
        for (k, p_x, v_x), (_, p_t, v_t) in zip(derivs_x, derivs_t):
            sym = "f'(x*)" if k == 1 else ("f''(x*)" if k == 2 else f"f^({k})(x*)")
            d_rows.append({
                "Cấp (k)": k,
                "Ký hiệu": sym,
                "Giá trị f^(k)(x*)": f"{v_x:.10f}",
                f"P^(k)({var_name}*)": f"{v_t:.10f}",
                "Đa thức đạo hàm P^(k)(x)": poly_str(p_x, var="x")
            })
        lines.append(pd.DataFrame(d_rows).to_string(index=False))
    else:
        lines.append("Không có điểm x* (eval_x) được chỉ định trong input.txt hoặc tham số.")

    return "\n".join(lines), val_final, dict_derivs_x, P_x_total

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Newton mốc cách đều - In toàn bộ trường hợp Tiến & Lùi")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("--type", "-t", choices=["forward", "backward", "tien", "lui", "both", "all"], default="both", help="Loại nội suy (forward/backward/both)")
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức (ví dụ 4)")
    p.add_argument("-n", "--n_points", help="Số mốc nội suy (ví dụ 5 hoặc 'all')")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(1.4, 0.9523), (1.5, 0.9661), (1.6, 0.9753), (1.7, 0.9838), (1.8, 0.9891), (1.9, 0.9928), (2.0, 0.9)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, _, msg = select_window_for_method(points, eval_x, "newton_finite", n_pts_cfg, degree_cfg)

    # Xây dựng bảng sai phân tiến Δ và lùi ∇
    df_upper, tbl_upper, df_lower, tbl_lower = build_finite_diff_tables(sub_points)

    mode_arg = args.type if args.type else params.get("type", "both")
    mode_arg = mode_arg.lower()

    out_sections = []
    if msg:
        out_sections.append(msg + "\n")

    # Bảng sai phân hữu hạn ở đầu
    out_sections.append("=" * 85)
    out_sections.append("BẢNG SAI PHÂN HỮU HẠN MỐC CÁCH ĐỀU")
    out_sections.append("=" * 85)
    out_sections.append("a) BẢNG SAI PHÂN TIẾN Δ (DẠNG TAM GIÁC TRÊN - DÙNG CHO TIẾN):")
    out_sections.append(df_upper.to_string(index=False) + "\n")
    out_sections.append("b) BẢNG SAI PHÂN LÙI ∇ (DẠNG TAM GIÁC DƯỚI - DÙNG CHO LÙI):")
    out_sections.append(df_lower.to_string(index=False) + "\n")

    if mode_arg in ["forward", "tien"]:
        res, _, _, _ = run_newton_finite_single(sub_points, tbl_upper, tbl_lower, mode="forward", eval_x=eval_x)
        out_sections.append(res)
    elif mode_arg in ["backward", "lui"]:
        res, _, _, _ = run_newton_finite_single(sub_points, tbl_upper, tbl_lower, mode="backward", eval_x=eval_x)
        out_sections.append(res)
    else:
        # In ĐỒNG THỜI CẢ TIẾN VÀ LÙI
        res_fwd, val_fwd, d_fwd, p_fwd = run_newton_finite_single(sub_points, tbl_upper, tbl_lower, mode="forward", eval_x=eval_x, section_prefix="PHẦN A: ")
        res_bwd, val_bwd, d_bwd, p_bwd = run_newton_finite_single(sub_points, tbl_upper, tbl_lower, mode="backward", eval_x=eval_x, section_prefix="PHẦN B: ")

        out_sections.append(res_fwd)
        out_sections.append("\n\n" + ("#" * 85) + "\n\n")
        out_sections.append(res_bwd)
        out_sections.append("\n\n" + ("#" * 85) + "\n\n")

        # PHẦN C: TỔNG KẾT & SO SÁNH
        summary_lines = []
        summary_lines.append("=" * 85)
        summary_lines.append("PHẦN C: TỔNG KẾT & SO SÁNH (KIỂM CHỨNG TÍNH ĐỒNG NHẤT GIỮA TIẾN VÀ LÙI)")
        summary_lines.append("=" * 85)

        summary_lines.append("1. So sánh hệ số đa thức quy đổi về biến thực P(x) = sum(c_k * x^k):")
        comp_rows = []
        max_coeff_diff = 0.0
        for deg in range(len(p_fwd)):
            cf = p_fwd[deg]
            cb = p_bwd[deg]
            diff_c = abs(cf - cb)
            max_coeff_diff = max(max_coeff_diff, diff_c)
            comp_rows.append({
                "Bậc (x^k)": deg,
                "Hệ số Tiến (Forward)": cf,
                "Hệ số Lùi (Backward)": cb,
                "Độ lệch tuyệt đối": f"{diff_c:.2e}"
            })
        summary_lines.append(pd.DataFrame(comp_rows).to_string(index=False))
        summary_lines.append(f"=> Độ lệch hệ số cực đại: {max_coeff_diff:.2e} (hoàn toàn triệt tiêu, sai số làm tròn số học)\n")

        if eval_x is not None and val_fwd is not None and val_bwd is not None:
            summary_lines.append(f"2. So sánh tại điểm x* = {eval_x}:")
            summary_lines.append(f"   a) Giá trị hàm f({eval_x}):")
            summary_lines.append(f"      - Newton Tiến (Forward):  P_tien({eval_x}) = {val_fwd:.10f}")
            summary_lines.append(f"      - Newton Lùi  (Backward): P_lui({eval_x})  = {val_bwd:.10f}")
            diff_val = abs(val_fwd - val_bwd)
            summary_lines.append(f"      => Độ lệch: |P_tien - P_lui| = {diff_val:.2e}")

            prime_names = {1: "Đạo hàm cấp 1 f'", 2: "Đạo hàm cấp 2 f''", 3: "Đạo hàm cấp 3 f'''"}
            sub_lbls = ["b", "c", "d"]
            for idx, k in enumerate([1, 2, 3]):
                if k in d_fwd and k in d_bwd:
                    lbl = sub_lbls[idx]
                    t_name = prime_names.get(k, f"Đạo hàm cấp {k} f^({k})")
                    summary_lines.append(f"\n   {lbl}) {t_name}({eval_x}):")
                    summary_lines.append(f"      - Newton Tiến: {d_fwd[k]:.10f}")
                    summary_lines.append(f"      - Newton Lùi:  {d_bwd[k]:.10f}")
                    diff_k = abs(d_fwd[k] - d_bwd[k])
                    summary_lines.append(f"      => Độ lệch: {diff_k:.2e}")

            if diff_val < 1e-9:
                summary_lines.append("\n   => KẾT LUẬN: Đa thức nội suy P(x) và các cấp đạo hàm từ 2 hướng TIẾN và LÙI hoàn toàn trùng khớp!")

        summary_lines.append("\n3. Nhận xét sư phạm về việc lựa chọn công thức khi làm bài thi:")
        summary_lines.append("   - Cả 2 công thức Newton Tiến và Newton Lùi đều cho ra cùng một đa thức nội suy P(x) duy nhất.")
        summary_lines.append("   - Khi nào nên dùng Newton Tiến (Forward)?")
        summary_lines.append("     + Khi điểm x* nằm ở NỬA ĐẦU khoảng mốc (gần x0), tức 0 <= t <= 1.")
        summary_lines.append("     + Lúc này các lũy thừa t^k nhỏ dần, chuỗi Newton hội tụ nhanh, sai số làm tròn nhỏ nhất.")
        summary_lines.append("   - Khi nào nên dùng Newton Lùi (Backward)?")
        summary_lines.append("     + Khi điểm x* nằm ở NỬA CUỐI khoảng mốc (gần xn), tức -1 <= t_b <= 0.")
        summary_lines.append("     + Lúc này các lũy thừa t_b^k nhỏ dần, chuỗi Newton lùi hội tụ nhanh và ít sai số nhất.")
        out_sections.append("\n".join(summary_lines))

    final_content = format_report_output("\n".join(out_sections), args.input_file)
    print("\n" + final_content)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
