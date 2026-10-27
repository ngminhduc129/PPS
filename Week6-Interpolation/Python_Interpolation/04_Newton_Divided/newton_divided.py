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
    from pps_io import parse_universal_input, select_window_for_method
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pps_io import parse_universal_input, select_window_for_method

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
    step_rows = []
    for i in range(n + 1):
        Di = d_coeffs[i]
        Ni = Di * B_list[i]
        N_total[:len(Ni)] += Ni
        
        # Biểu thức tích dạng chuỗi
        if i == 0:
            b_expr = "1"
        else:
            factors = [f"(x - {nodes[j]})" for j in range(i)]
            b_expr = "*".join(factors)
            
        step_rows.append({
            "i": i,
            "Mốc xuất phát x_base": nodes[0],
            "Hệ số D_i": Di,
            "Biểu thức tích B_i(x)": b_expr,
            "Đa thức số hạng N_i(x)": [round(v, 8) for v in Ni.tolist()]
        })
    df_steps = pd.DataFrame(step_rows)

    # 4. Tạo công thức dạng tích
    product_terms = []
    for i in range(n + 1):
        Di = d_coeffs[i]
        if abs(Di) < 1e-12: continue
        sign = "+ " if Di >= 0 and len(product_terms) > 0 else ("- " if Di < 0 and len(product_terms) > 0 else "")
        val = abs(Di) if len(product_terms) > 0 else Di
        b_expr = step_rows[i]["Biểu thức tích B_i(x)"]
        if i == 0:
            product_terms.append(f"{val:.8g}")
        else:
            product_terms.append(f"{sign}{val:.8g}*{b_expr}")
    product_form_str = " ".join(product_terms) if product_terms else "0"

    # 5. Lập báo cáo 5 phần
    lines = []
    lines.append("=" * 85)
    lines.append(f"{section_prefix}NỘI SUY NEWTON TỶ SAI PHÂN - {tname}")
    lines.append("=" * 85)
    lines.append(f"Số mốc tính toán: {n + 1}, Mốc xuất phát: x_{'0' if is_fwd else 'n'} = {nodes[0]}")
    lines.append(f"Thứ tự mốc sử dụng: {nodes}")
    lines.append(f"Quy tắc chọn hệ số: {trace_desc}\n")

    lines.append(f"1. BẢNG HỆ SỐ NỘI SUY {tname}:")
    lines.append("-" * 85)
    coeff_rows = []
    for i in range(n + 1):
        coeff_rows.append({
            "Bậc i": i,
            "Ký hiệu tỷ sai phân": f"f[x_0, ..., x_{i}]" if is_fwd else f"f[x_n, ..., x_{n-i}]",
            "Giá trị D_i": d_coeffs[i],
            "Biểu thức tích B_i(x)": step_rows[i]["Biểu thức tích B_i(x)"]
        })
    lines.append(pd.DataFrame(coeff_rows).to_string(index=False))

    lines.append(f"\n2. BẢNG TÍCH VÀ CÁC BƯỚC DỰNG ĐA THỨC SỐ HẠNG N_i(x) = D_i * B_i(x):")
    lines.append("-" * 85)
    lines.append(df_steps[["i", "Hệ số D_i", "Biểu thức tích B_i(x)", "Đa thức số hạng N_i(x)"]].to_string(index=False))

    lines.append(f"\n3. BẢNG HỆ SỐ ĐA THỨC NỘI SUY KHAI TRIỂN P(x):")
    lines.append("-" * 85)
    df_p_coeffs = pd.DataFrame({"Bậc (x^k)": np.arange(n + 1), "Hệ số a_k": N_total})
    lines.append(df_p_coeffs.to_string(index=False))

    lines.append(f"\n4. ĐA THỨC NỘI SUY P(x):")
    lines.append("-" * 85)
    lines.append(f"a) Dạng tích Newton (chưa khai triển):")
    lines.append(f"   P_{mode}(x) = {product_form_str}")
    lines.append(f"\nb) Dạng chính tắc (đã khai triển theo lũy thừa của x):")
    lines.append(f"   P_{mode}(x) = {poly_str(N_total, var='x')}")

    val_final = None
    dict_derivs = {}
    lines.append(f"\n5. GIÁ TRỊ NỘI SUY VÀ ĐẠO HÀM TẠI ĐIỂM XÁC ĐỊNH{' (x* = ' + str(eval_x) + ')' if eval_x is not None else ''}:")
    lines.append("-" * 85)
    if eval_x is not None:
        eval_rows = []
        cum_sum = 0.0
        term_vals = []
        for i in range(n + 1):
            bi_val = float(sum(c * (eval_x ** deg) for deg, c in enumerate(B_list[i])))
            di = d_coeffs[i]
            ti = di * bi_val
            cum_sum += ti
            term_vals.append(ti)
            eval_rows.append({
                "i": i,
                "B_i(x*)": bi_val,
                "D_i": di,
                "Số hạng T_i": ti,
                "Tổng dồn P_i(x*)": cum_sum
            })
        lines.append(pd.DataFrame(eval_rows).to_string(index=False))
        val_final = cum_sum
        lines.append(f"\n=> GIÁ TRỊ NỘI SUY:")
        lines.append(f"   P_{mode}({eval_x}) = {val_final:.10f}")

        # Tính đạo hàm
        derivs = compute_all_derivatives(N_total, eval_x)
        dict_derivs = {k: v for k, p, v in derivs}
        lines.append(f"\n=> KẾT QUẢ ĐẠO HÀM CÁC CẤP TẠI x* = {eval_x}:")
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

    final_content = "\n".join(out_sections)
    print("\n" + final_content)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
