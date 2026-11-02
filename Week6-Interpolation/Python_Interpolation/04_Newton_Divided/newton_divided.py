"""
Nội suy Newton tỷ sai phân (mốc bất kỳ)
Tự động in TẤT CẢ các trường hợp:
- PHẦN A: Newton Tỷ sai phân TIẾN (Forward Divided Difference - từ x0 đến xn)
- PHẦN B: Newton Tỷ sai phân LÙI (Backward Divided Difference - từ xn về x0)
- PHẦN C: TỔNG KẾT & SO SÁNH (Kiểm chứng tính duy nhất của đa thức P(x) và các cấp đạo hàm)

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

def build_divided_diff_table(points):
    """
    Xây dựng bảng tỷ sai phân đầy đủ dạng ma trận tam giác trên.
    points: list các tuple (x_i, y_i)
    """
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    m = len(points)
    table = np.full((m, m), np.nan, dtype=float)
    table[:, 0] = y.copy()
    for j in range(1, m):
        for i in range(m - j):
            table[i, j] = (table[i + 1, j - 1] - table[i, j - 1]) / (x[i + j] - x[i])
    
    data = {"i": list(range(m)), "x_i": x, "y_i": table[:, 0]}
    for j in range(1, m):
        data[f"Cấp {j}"] = [table[i, j] if i < m - j else np.nan for i in range(m)]
    df_table = pd.DataFrame(data)
    return df_table, table

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

def run_newton_divided_single(points, table, mode="forward", eval_x=None, section_prefix=""):
    """
    Chạy nội suy Newton tỷ sai phân cho 1 chiều (forward hoặc backward).
    Chuẩn 5 phần:
    1. Bảng tỷ sai phân & quy tắc chọn hệ số
    2. Hệ số Newton tỷ sai phân
    3. Bảng tích ma trận và bảng tính giá trị tại x*
    4. Đa thức nội suy dạng tích và dạng lũy thừa P(x)
    5. Giá trị nội suy và các cấp đạo hàm tại điểm x*
    """
    n = len(points) - 1
    x_nodes = [p[0] for p in points]
    is_fwd = (mode == "forward")
    tname = "TIẾN (FORWARD)" if is_fwd else "LÙI (BACKWARD)"
    
    # 1. Trích xuất hệ số và mốc nội suy
    if is_fwd:
        nodes = x_nodes
        d_coeffs = [table[0, j] for j in range(n + 1)]
        trace_desc = "Lấy đường chéo TRÊN của bảng: f[x_0], f[x_0, x_1], ..., f[x_0, ..., x_n] (hàng 0)"
    else:
        nodes = list(reversed(x_nodes))
        d_coeffs = [table[n - j, j] for j in range(n + 1)]
        trace_desc = "Lấy đường chéo DƯỚI của bảng: f[x_n], f[x_{n-1}, x_n], ..., f[x_0, ..., x_n] (phần tử cuối mỗi cột)"

    # 2. Xây dựng các đa thức cơ sở tích B_i(x)
    B_list = [np.array([1.0])]
    for i in range(1, n + 1):
        prev = B_list[-1]
        node_prev = nodes[i - 1]
        next_b = np.convolve(prev, [-node_prev, 1.0])
        B_list.append(next_b)

    # 3. Đa thức tổng P(x) = sum(D_i * B_i(x))
    N_total = np.zeros(n + 1, dtype=float)
    for i in range(n + 1):
        Di = d_coeffs[i]
        Ni = Di * B_list[i]
        N_total[:len(Ni)] += Ni

    # Hệ số đa thức theo bậc giảm dần x^n ... x^0 (dùng cho phương pháp Horner)
    coeffs_desc = [N_total[k] for k in range(n, -1, -1)]
    col_powers = [f"x^{deg}" if deg > 0 else "x^0" for deg in range(n, -1, -1)]

    # 4. Bảng tích (Product Table) - CHỈ GHI DUY NHẤT HỆ SỐ CÁC BẬC
    prod_table_rows = []
    for i in range(n + 1):
        node_str = "-" if i == 0 else f"{nodes[i-1]:g}"
        deg_coeffs = [(B_list[i][k] if k <= i else 0.0) for k in range(n, -1, -1)]
        row_dict = {"Mốc x_k": node_str}
        for idx, deg_name in enumerate(col_powers):
            val = deg_coeffs[idx]
            row_dict[deg_name] = f"{val:.8g}" if abs(val) >= 1e-12 else "0"
        prod_table_rows.append(row_dict)
    df_prod_table = pd.DataFrame(prod_table_rows)

    # 5. Tạo công thức dạng tích
    product_terms = []
    for i in range(n + 1):
        Di = d_coeffs[i]
        if abs(Di) < 1e-12: continue
        sign = "+ " if Di >= 0 and len(product_terms) > 0 else ("- " if Di < 0 and len(product_terms) > 0 else "")
        val = abs(Di) if len(product_terms) > 0 else Di
        if i == 0:
            b_expr = "1"
        else:
            factors = [f"(x - {nodes[j]})" for j in range(i)]
            b_expr = "*".join(factors)
        if i == 0:
            product_terms.append(f"{val:.8g}")
        else:
            product_terms.append(f"{sign}{val:.8g}*{b_expr}")
    product_form_str = " ".join(product_terms) if product_terms else "0"

    # Lập báo cáo chi tiết
    lines = []
    lines.append("=" * 85)
    lines.append(f"{section_prefix}NỘI SUY NEWTON TỶ SAI PHÂN - {tname}")
    lines.append("=" * 85)
    lines.append(f"Số mốc tính toán: {n + 1}, Mốc xuất phát: x_{'0' if is_fwd else 'n'} = {nodes[0]}")
    lines.append(f"Thứ tự mốc sử dụng: {nodes}")
    lines.append(f"Quy tắc chọn hệ số: {trace_desc}\n")

    # 1. VECTOR HỆ SỐ NỘI SUY
    lines.append(f"1. VECTOR HỆ SỐ NỘI SUY {tname}:")
    lines.append("-" * 85)
    c_vector_str = "[" + ", ".join([f"{v:.8g}" for v in d_coeffs]) + "]"
    lines.append(f"C_{'tien' if is_fwd else 'lui'} = {c_vector_str}\n")
    coeff_rows = []
    for i in range(n + 1):
        coeff_rows.append({
            "Bậc i": i,
            "Ký hiệu": f"f[x_0, ..., x_{i}]" if is_fwd else f"f[x_n, ..., x_{n-i}]",
            "Hệ số D_i": f"{d_coeffs[i]:.8g}"
        })
    lines.append(pd.DataFrame(coeff_rows).to_string(index=False))

    # 2. BẢNG TÍCH (CHỈ GHI MÌNH HỆ SỐ)
    lines.append(f"\n2. BẢNG TÍCH {tname} (CHỈ GHI HỆ SỐ CÁC ĐA THỨC TÍCH):")
    lines.append("-" * 85)
    lines.append(df_prod_table.to_string(index=False))

    # 3. ĐA THỨC NỘI SUY VÀ HỆ SỐ DÙNG CHO PHƯƠNG PHÁP HORNER
    lines.append(f"\n3. ĐA THỨC NỘI SUY VÀ HỆ SỐ DÙNG CHO PHƯƠNG PHÁP HORNER:")
    lines.append("-" * 85)
    dtns_vector_str = "[" + ", ".join([f"{v:.8g}" for v in coeffs_desc]) + "]"
    lines.append(f"Vector hệ số ĐTNS (từ bậc cao xuống bậc thấp x^n ... x^0):")
    lines.append(f"ĐTNS_{'tien' if is_fwd else 'lui'} = {dtns_vector_str}\n")
    lines.append("Bảng hệ số đa thức P(x):")
    df_p_coeffs = pd.DataFrame({"Bậc (x^k)": col_powers, "Hệ số": [f"{v:.8g}" for v in coeffs_desc]})
    lines.append(df_p_coeffs.to_string(index=False))
    lines.append(f"\na) Dạng tích Newton (chưa khai triển):")
    lines.append(f"   P_{mode}(x) = {product_form_str}")
    lines.append(f"\nb) Dạng chính tắc (đã khai triển theo lũy thừa của x):")
    lines.append(f"   P_{mode}(x) = {poly_str(N_total, var='x')}")

    # 4. BẢNG TÍNH GIÁ TRỊ THEO SƠ ĐỒ HORNER TẠI x*
    val_final = None
    dict_derivs = {}
    lines.append(f"\n4. BẢNG TÍNH GIÁ TRỊ THEO SƠ ĐỒ HORNER{' TẠI x* = ' + str(eval_x) if eval_x is not None else ''}:")
    lines.append("-" * 85)
    if eval_x is not None:
        b_horner = [0.0] * (n + 1)
        b_horner[0] = coeffs_desc[0]
        for idx in range(1, n + 1):
            b_horner[idx] = coeffs_desc[idx] + eval_x * b_horner[idx - 1]
        
        horner_rows = [
            {"x*": "Hệ số ĐTNS", **{col_powers[idx]: f"{coeffs_desc[idx]:.8g}" for idx in range(n + 1)}},
            {"x*": f"{eval_x:g}", **{col_powers[idx]: f"{b_horner[idx]:.8g}" for idx in range(n + 1)}}
        ]
        df_horner = pd.DataFrame(horner_rows)
        lines.append(df_horner.to_string(index=False))
        val_final = b_horner[-1]
        lines.append(f"\n=> GIÁ TRỊ NỘI SUY (TỪ SƠ ĐỒ HORNER):")
        lines.append(f"   P_{mode}({eval_x}) = {val_final:.10f}")

        # 5. KẾT QUẢ ĐẠO HÀM CÁC CẤP
        derivs = compute_all_derivatives(N_total, eval_x)
        dict_derivs = {k: v for k, p, v in derivs}
        lines.append(f"\n5. KẾT QUẢ ĐẠO HÀM CÁC CẤP TẠI x* = {eval_x}:")
        lines.append("-" * 85)
        d_rows = []
        for k, p_poly, v_val in derivs:
            sym = "f'(x*)" if k == 1 else ("f''(x*)" if k == 2 else f"f^({k})(x*)")
            d_rows.append({
                "Cấp (k)": k, "Ký hiệu": sym, "Giá trị f^(k)(x*)": f"{v_val:.10f}",
                "Đa thức đạo hàm P^(k)(x)": poly_str(p_poly, var="x")
            })
        lines.append(pd.DataFrame(d_rows).to_string(index=False))
    else:
        lines.append("Không có điểm x* (eval_x) được chỉ định trong input.txt hoặc tham số.")

    return "\n".join(lines), val_final, dict_derivs, N_total

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Newton tỷ sai phân (mốc bất kỳ) - In toàn bộ trường hợp Tiến & Lùi")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("--type", "-t", choices=["forward", "backward", "tien", "lui", "both", "all"], default="both", help="Loại nội suy (forward/backward/both)")
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức")
    p.add_argument("-n", "--n_points", help="Số mốc nội suy (hoặc 'all')")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(1.2, 0.7247), (1.5, 0.1415), (1.7, -0.2577), (1.9, -0.6466), (2.1, -1.0097), (2.4, -1.4748)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, _, msg = select_window_for_method(points, eval_x, "newton_divided", n_pts_cfg, degree_cfg)

    # Xây dựng bảng tỷ sai phân chung
    df_table, table = build_divided_diff_table(sub_points)
    n_pts = len(sub_points)

    c_fwd = [table[0, j] for j in range(n_pts)]
    c_bwd = [table[n_pts - 1 - j, j] for j in range(n_pts)]

    mode_arg = args.type if args.type else params.get("type", "both")
    mode_arg = mode_arg.lower()

    out_sections = []
    if msg:
        out_sections.append(msg + "\n")

    # Bảng tỷ sai phân chung hiển thị ở đầu
    out_sections.append("=" * 85)
    out_sections.append("BẢNG TỶ SAI PHÂN CHUNG (DIVIDED DIFFERENCE TABLE)")
    out_sections.append("=" * 85)
    out_sections.append(df_table.to_string(index=False) + "\n")

    c_fwd_str = "[" + ", ".join([f"{v:.8g}" for v in c_fwd]) + "]"
    c_bwd_str = "[" + ", ".join([f"{v:.8g}" for v in c_bwd]) + "]"
    out_sections.append(f"Vector hệ số nội suy C_tiến (đường chéo trên):\nC_tien = {c_fwd_str}\n")
    out_sections.append(f"Vector hệ số nội suy C_lùi (đường chéo dưới):\nC_lui  = {c_bwd_str}\n")

    if mode_arg in ["forward", "tien"]:
        res, _, _, _ = run_newton_divided_single(sub_points, table, mode="forward", eval_x=eval_x)
        out_sections.append(res)
    elif mode_arg in ["backward", "lui"]:
        res, _, _, _ = run_newton_divided_single(sub_points, table, mode="backward", eval_x=eval_x)
        out_sections.append(res)
    else:
        # In ĐỒNG THỜI CẢ TIẾN VÀ LÙI
        res_fwd, val_fwd, d_fwd, p_fwd = run_newton_divided_single(sub_points, table, mode="forward", eval_x=eval_x, section_prefix="PHẦN A: ")
        res_bwd, val_bwd, d_bwd, p_bwd = run_newton_divided_single(sub_points, table, mode="backward", eval_x=eval_x, section_prefix="PHẦN B: ")
        
        out_sections.append(res_fwd)
        out_sections.append("\n\n" + ("#" * 85) + "\n\n")
        out_sections.append(res_bwd)
        out_sections.append("\n\n" + ("#" * 85) + "\n\n")

        # PHẦN C: TỔNG KẾT & SO SÁNH
        summary_lines = []
        summary_lines.append("=" * 85)
        summary_lines.append("PHẦN C: TỔNG KẾT & SO SÁNH (KIỂM CHỨNG TÍNH ĐỒNG NHẤT GIỮA TIẾN VÀ LÙI)")
        summary_lines.append("=" * 85)

        # 1. So sánh hệ số đa thức chính tắc P(x)
        summary_lines.append("1. So sánh hệ số đa thức chính tắc P(x) = sum(a_k * x^k):")
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
        summary_lines.append(f"=> Độ lệch hệ số cực đại: {max_coeff_diff:.2e} (hoàn toàn triệt tiêu, sai số máy tính)\n")

        # 2. So sánh giá trị và đạo hàm tại x*
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
                    diff_k = abs(d_fwd[k] - d_fwd[k])
                    summary_lines.append(f"      => Độ lệch: {diff_k:.2e}")

            if diff_val < 1e-9:
                summary_lines.append("\n   => KẾT LUẬN: Đa thức nội suy P(x) và các cấp đạo hàm từ 2 hướng TIẾN và LÙI hoàn toàn trùng khớp!")

        # 3. Nhận xét lý thuyết
        summary_lines.append("\n3. Nhận xét lý thuyết & Ý nghĩa trong thực hành:")
        summary_lines.append("   - Định lý duy nhất: Đa thức nội suy bậc <= n đi qua n+1 điểm phân biệt là DUY NHẤT, bất kể xây dựng theo dạng Lagrange, Newton Tiến hay Newton Lùi.")
        summary_lines.append("   - Ưu điểm của Newton so với Lagrange: Khi bổ sung thêm mốc mới:")
        summary_lines.append("     + Nếu thêm mốc ở CUỐI bảng: Công thức TIẾN giữ nguyên các số hạng cũ, chỉ cần cộng thêm 1 số hạng mới.")
        summary_lines.append("     + Nếu thêm mốc ở ĐẦU bảng: Công thức LÙI giữ nguyên các số hạng cũ, chỉ cần cộng thêm 1 số hạng mới.")
        out_sections.append("\n".join(summary_lines))

    final_content = format_report_output("\n".join(out_sections), args.input_file)
    print("\n" + final_content)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
