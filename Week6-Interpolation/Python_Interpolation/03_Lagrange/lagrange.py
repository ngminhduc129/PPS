"""
Nội suy Lagrange (Lagrange Interpolation)
Tự động in TẤT CẢ các trường hợp:
- PHẦN A: Nội suy Lagrange trên Cửa sổ cục bộ lân cận quanh x* (Local Windowing)
- PHẦN B: Nội suy Lagrange trên Toàn bộ mốc dữ liệu (Full Nodes - nếu N > target)
- PHẦN C: Nội suy Lagrange Ngược (Inverse Lagrange - nếu đề có y*)
- PHẦN D: TỔNG KẾT & SO SÁNH (Đạo hàm các cấp, phân tích hiện tượng Runge & sai số)

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

def lagrange_core(points, eval_x=None, var="x"):
    """
    Tính toán chi tiết đa thức nội suy Lagrange:
    - Bảng các đa thức cơ sở L_i(x)
    - Hệ số đa thức tổng P(x) = sum(y_i * L_i(x))
    - Bảng tính giá trị từng L_i(x*) và tổng dồn P(x*) (nếu có eval_x)
    """
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    n = len(points)
    L_polys = []
    P_coeffs = np.zeros(n, dtype=float)
    step_data = []
    eval_rows = []
    cum_sum = 0.0

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

        # Công thức chuỗi cơ sở
        factors = [f"({var} - {r})" for r in roots]
        l_expr = f"({' * '.join(factors)}) / ({denom:.8g})"

        step_data.append({
            "i": i, f"{var}_i": x[i], "y_i": y[i],
            "Mẫu số w'_i": denom,
            f"Đa thức cơ sở L_i({var})": [round(v, 8) for v in basis.tolist()]
        })

        if eval_x is not None:
            l_val = float(sum(c * (eval_x ** deg) for deg, c in enumerate(basis)))
            term_val = y[i] * l_val
            cum_sum += term_val
            eval_rows.append({
                "i": i,
                f"{var}_i": x[i],
                "y_i": y[i],
                f"L_i({var}*)": l_val,
                f"y_i * L_i({var}*)": term_val,
                f"Tổng dồn P({var}*)": cum_sum
            })

    df_steps = pd.DataFrame(step_data)
    df_coeffs = pd.DataFrame({"Bậc": np.arange(n), "Hệ số": P_coeffs})
    df_eval = pd.DataFrame(eval_rows) if eval_rows else None

    return df_steps, df_coeffs, P_coeffs, df_eval, cum_sum

def poly_str(coeffs, var="x"):
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

def run_lagrange_report(points, eval_x=None, title="", section_prefix="", var="x"):
    df_steps, df_coeffs, P_coeffs, df_eval, val_final = lagrange_core(points, eval_x=eval_x, var=var)
    n = len(points)

    lines = []
    lines.append("=" * 85)
    lines.append(f"{section_prefix}{title}")
    lines.append("=" * 85)
    lines.append(f"Số mốc tính toán: {n}, Bậc đa thức: {n - 1}")
    lines.append(f"Tập mốc {var}: {[p[0] for p in points]}")
    if eval_x is not None:
        lines.append(f"Điểm xác định: {var}* = {eval_x}\n")

    lines.append(f"1. CÁC BƯỚC TÍNH TOÁN VÀ ĐA THỨC CƠ SỞ L_i({var}):")
    lines.append("-" * 85)
    lines.append(df_steps.to_string(index=False))

    lines.append(f"\n2. BẢNG HỆ SỐ ĐA THỨC NỘI SUY P({var}) = sum(a_k * {var}^k):")
    lines.append("-" * 85)
    lines.append(df_coeffs.to_string(index=False))

    lines.append(f"\n3. ĐA THỨC NỘI SUY P({var}):")
    lines.append("-" * 85)
    lines.append(f"P({var}) = {poly_str(P_coeffs, var=var)}")

    dict_derivs = {}
    if eval_x is not None:
        lines.append(f"\n4. BẢNG TÍNH GIÁ TRỊ CƠ SỞ TẠI {var}* = {eval_x}:")
        lines.append("-" * 85)
        lines.append(df_eval.to_string(index=False))
        lines.append(f"\n=> GIÁ TRỊ NỘI SUY: P({eval_x}) = {val_final:.10f}")

        # Đạo hàm
        derivs = compute_all_derivatives(P_coeffs, eval_x)
        dict_derivs = {k: v for k, p, v in derivs}
        lines.append(f"\n5. KẾT QUẢ ĐẠO HÀM CÁC CẤP TẠI {var}* = {eval_x}:")
        lines.append("-" * 85)
        d_rows = []
        for k, p_poly, v_val in derivs:
            sym = f"P'({var}*)" if k == 1 else (f"P''({var}*)" if k == 2 else f"P^({k})({var}*)")
            d_rows.append({
                "Cấp (k)": k, "Ký hiệu": sym, f"Giá trị tại {var}*": f"{v_val:.10f}",
                f"Đa thức đạo hàm P^(k)({var})": poly_str(p_poly, var=var)
            })
        lines.append(pd.DataFrame(d_rows).to_string(index=False))

    return "\n".join(lines), val_final, dict_derivs, P_coeffs

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Lagrange - In toàn bộ trường hợp Cục bộ, Toàn cục, Ngược & So sánh")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("-y", "--eval_y", type=float, help="Điểm cần giải ngược y*")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức")
    p.add_argument("-n", "--n_points", help="Số mốc nội suy (hoặc 'all')")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(0.87, 2.40), (1.24, 2.20), (2.99, 2.00), (3.67, 1.80), (4.23, 1.60)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    eval_y = args.eval_y if args.eval_y is not None else params.get("eval_y", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, _, msg = select_window_for_method(points, eval_x, "lagrange", n_pts_cfg, degree_cfg)

    sections = []
    if msg:
        sections.append(msg + "\n")

    # 1. PHẦN A: Nội suy Lagrange trên Cửa sổ cục bộ
    res_local, val_local, d_local, p_local = run_lagrange_report(
        sub_points, eval_x=eval_x,
        title="NỘI SUY LAGRANGE TRÊN CỬA SỔ CỤC BỘ (LOCAL WINDOWING)",
        section_prefix="PHẦN A: "
    )
    sections.append(res_local)

    # 2. PHẦN B: Nội suy Lagrange trên Toàn bộ mốc (nếu khác cửa sổ cục bộ)
    val_full = None
    if len(points) > len(sub_points):
        res_full, val_full, d_full, p_full = run_lagrange_report(
            points, eval_x=eval_x,
            title="NỘI SUY LAGRANGE TRÊN TOÀN BỘ MỐC DỮ LIỆU (FULL NODES)",
            section_prefix="PHẦN B: "
        )
        sections.append("\n\n" + ("#" * 85) + "\n\n")
        sections.append(res_full)

    # 3. PHẦN C: Nội suy Lagrange Ngược (nếu có eval_y)
    if eval_y is not None:
        inv_points = [(p[1], p[0]) for p in sub_points]
        res_inv, val_inv, _, _ = run_lagrange_report(
            inv_points, eval_x=eval_y,
            title="NỘI SUY LAGRANGE NGƯỢC (INVERSE LAGRANGE: x = P(y))",
            section_prefix="PHẦN C: ",
            var="y"
        )
        sections.append("\n\n" + ("#" * 85) + "\n\n")
        sections.append(res_inv)

    # 4. PHẦN D: TỔNG KẾT & SO SÁNH
    sum_lines = []
    sum_lines.append("\n" + ("#" * 85))
    sum_lines.append("PHẦN D: TỔNG KẾT & NHẬN XÉT LÝ THUYẾT")
    sum_lines.append("#" * 85)

    if eval_x is not None:
        sum_lines.append(f"1. So sánh giá trị nội suy tại x* = {eval_x}:")
        sum_lines.append(f"   - Đa thức cục bộ ({len(sub_points)} mốc): P_local({eval_x}) = {val_local:.10f}")
        if val_full is not None:
            sum_lines.append(f"   - Đa thức toàn cục ({len(points)} mốc): P_full({eval_x})  = {val_full:.10f}")
            diff_lf = abs(val_local - val_full)
            sum_lines.append(f"   => Độ lệch giữa cục bộ và toàn cục: |P_local - P_full| = {diff_lf:.2e}")
        else:
            sum_lines.append("   (Số mốc cục bộ bằng số mốc toàn bộ, đa thức duy nhất).")

    if eval_y is not None:
        sum_lines.append(f"\n2. Bài toán nội suy ngược tại y* = {eval_y}:")
        sum_lines.append(f"   - Nghiệm x* tìm được: x* = P_inv({eval_y}) = {val_inv:.10f}")

    sum_lines.append("\n3. Nhận xét lý thuyết về phương pháp Lagrange:")
    sum_lines.append("   - Ưu điểm: Đa thức cơ sở L_i(x) chỉ phụ thuộc vào các mốc x_i, không phụ thuộc vào giá trị hàm y_i.")
    sum_lines.append("   - Nhược điểm: Khối lượng tính toán lớn khi thêm mốc mới (phải tính lại toàn bộ các L_i(x)).")
    sum_lines.append("   - Hiện tượng Runge: Khi tăng số mốc cách đều lên quá lớn (N > 10), đa thức Lagrange toàn cục dễ dao động mạnh ở 2 biên, do đó trong phòng thi nên ưu tiên cắt cửa sổ cục bộ 4 hoặc 5 mốc quanh x*.")
    sections.append("\n".join(sum_lines))

    final_content = "\n".join(sections)
    print("\n" + final_content)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
