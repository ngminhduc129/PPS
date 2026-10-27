"""
Nội suy Stirling (Stirling Interpolation)
Xuất kết quả đúng theo 5 phần chuẩn tương tự Newton Finite và Gauss:
1. BẢNG SAI PHÂN HỮU HẠN & QUY TẮC CHỌN SAI PHÂN
2. HỆ SỐ NỘI SUY STIRLING (C CHẴN, C LẺ & CÁC DÃY HỆ SỐ)
3. BẢNG TÍCH STIRLING (MA TRẬN THEO t & BẢNG TÍNH TẠI t*)
4. ĐA THỨC NỘI SUY STIRLING (DẠNG TÍCH, THEO t VÀ QUY ĐỔI THEO x)
5. ĐA THỨC NỘI SUY VÀ ĐẠO HÀM TẠI ĐIỂM XÁC ĐỊNH (VALUE & DERIVATIVES AT EVALUATION POINT)

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

def build_finite_diff_table(points, x0_idx):
    """
    Xây dựng bảng sai phân hữu hạn tiến Δ cho mốc cách đều kèm cột chỉ số t_i = i - x0_idx.
    """
    xv = [p[0] for p in points]
    yv = [p[1] for p in points]
    n_pts = len(yv)
    t_indices = [i - x0_idx for i in range(n_pts)]

    tbl = np.full((n_pts, n_pts), np.nan, dtype=float)
    tbl[:, 0] = yv
    for j in range(1, n_pts):
        for i in range(n_pts - j):
            tbl[i, j] = tbl[i + 1, j - 1] - tbl[i, j - 1]

    data = {"t_i": t_indices, "x_i": xv, "y_i": yv}
    for j in range(1, n_pts):
        col = np.full(n_pts, np.nan)
        col[:(n_pts - j)] = tbl[:(n_pts - j), j]
        data[f"Cấp {j}"] = col
    df_tbl = pd.DataFrame(data)
    return df_tbl, tbl

def extract_stirling_coeffs(tbl, x0_idx, n):
    """
    Trích xuất các hệ số sai phân trung tâm C_i và hệ số rút gọn D_i = C_i / i!
    Phân tách rõ ràng C chẵn và C lẻ:
    - Bậc chẵn (i = 2k): C_{2k} = Δ^{2k} y_{-k} (hàng x0_idx - k, cột Cấp 2k)
    - Bậc lẻ (i = 2k-1, k = (i+1)//2): C_{2k-1} = (Δ^{2k-1} y_{-k} + Δ^{2k-1} y_{-k+1}) / 2
    """
    c_all = []
    d_all = []
    c_even, d_even, details_even = [], [], []
    c_odd, d_odd, details_odd = [], [], []

    for i in range(n + 1):
        fact_i = math.factorial(i)
        if i == 0:
            val = tbl[x0_idx, 0]
            di = val / fact_i
            c_all.append(val); d_all.append(di)
            c_even.append(val); d_even.append(di)
            details_even.append({
                "i": 0, "k": 0, "val": val, "di": di, "fact": fact_i,
                "desc": f"C_0 = y_0 = {val:.8g} => D_0 = {di:.8g}  (hàng {x0_idx}, t = 0)"
            })
        elif i % 2 == 0:
            # Chẵn: i = 2k
            k = i // 2
            row = x0_idx - k
            val = tbl[row, i]
            di = val / fact_i
            c_all.append(val); d_all.append(di)
            c_even.append(val); d_even.append(di)
            details_even.append({
                "i": i, "k": k, "val": val, "di": di, "fact": fact_i,
                "desc": f"C_{i} = Δ^{i} y_{{-{k}}} = {val:.8g} => D_{i} = {val:.8g} / {fact_i}! = {di:.8g}  (hàng {row}, t = -{k})"
            })
        else:
            # Lẻ: i = 2k - 1 (với k = (i + 1) // 2)
            k = (i + 1) // 2
            row1 = x0_idx - k
            row2 = x0_idx - k + 1
            c1 = tbl[row1, i]
            c2 = tbl[row2, i]
            val = (c1 + c2) / 2.0
            di = val / fact_i
            c_all.append(val); d_all.append(di)
            c_odd.append(val); d_odd.append(di)
            k_minus_1 = k - 1
            sub2_str = f"y_{{-{k_minus_1}}}" if k_minus_1 > 0 else "y_0"
            details_odd.append({
                "i": i, "k": k, "val": val, "di": di, "fact": fact_i, "c1": c1, "c2": c2,
                "desc": f"C_{i} = (Δ^{i} y_{{-{k}}} + Δ^{i} {sub2_str}) / 2 = ({c1:.8g} + {c2:.8g}) / 2 = {val:.8g} => D_{i} = {val:.8g} / {fact_i}! = {di:.8g}  (hàng {row1} & {row2})"
            })

    return c_all, d_all, c_even, d_even, c_odd, d_odd, details_even, details_odd

def get_stirling_basis_str(i):
    """
    Trả về biểu thức ký hiệu của đa thức cơ sở tích Stirling B_i(t).
    """
    if i == 0: return "1"
    if i == 1: return "t"
    if i == 2: return "t^2"
    if i % 2 == 1:
        # Lẻ: i = 2k + 1
        k = (i - 1) // 2
        factors = ["t"] + [f"(t^2 - {j}^2)" for j in range(1, k + 1)]
        return "*".join(factors)
    else:
        # Chẵn: i = 2k
        k = i // 2
        factors = ["t^2"] + [f"(t^2 - {j}^2)" for j in range(1, k)]
        return "*".join(factors)

def build_stirling_basis_polynomials(n):
    """
    Xây dựng danh sách các đa thức cơ sở Stirling B_i(t) theo hệ số tăng dần của t.
    B_0(t) = 1
    B_1(t) = t
    B_2(t) = t^2
    B_{2k}(t)   = B_{2k-1}(t) * t
    B_{2k+1}(t) = B_{2k-1}(t) * (t^2 - k^2)
    """
    B_list = [np.array([1.0]), np.array([0.0, 1.0])]
    for i in range(2, n + 1):
        if i % 2 == 0:
            Bi = np.convolve(B_list[i - 1], [0.0, 1.0])
        else:
            k = (i - 1) // 2
            Bi = np.convolve(B_list[i - 2], [-float(k**2), 0.0, 1.0])
        B_list.append(Bi)
    return B_list

def fmt_val(v):
    if abs(v) < 1e-11: return "0"
    if abs(v - round(v)) < 1e-11: return str(int(round(v)))
    formatted = f"{v:.8f}".rstrip("0").rstrip(".")
    return "0" if formatted == "-0" else formatted

def build_product_basis_table_stirling(n, B_list):
    """
    Xây dựng BẢNG TÍCH dạng ma trận theo biến chuẩn hóa t:
    Mỗi hàng i biểu diễn hệ số của đa thức cơ sở tích B_i(t) theo lũy thừa giảm dần:
    i | Biểu thức tích B_i(t) | t^n | t^(n-1) | ... | t^1 | t^0
    """
    matrix_rows = []
    for i in range(n + 1):
        poly_asc = B_list[i]
        poly_desc = poly_asc[::-1]
        row_coeffs = np.zeros(n + 1, dtype=float)
        row_coeffs[(n + 1) - len(poly_desc):] = poly_desc
        b_str = get_stirling_basis_str(i)
        row_display = [str(i), b_str] + [fmt_val(v) for v in row_coeffs]
        matrix_rows.append(row_display)

    cols = ["i", "Biểu thức tích B_i(t)"] + [f"t^{p}" for p in range(n, -1, -1)]
    df_prod = pd.DataFrame(matrix_rows, columns=cols)
    return df_prod

def build_eval_product_table_stirling(n, c_all, d_all, B_list, eval_t):
    """
    Lập bảng tính giá trị các thừa số tích Stirling và tổng dồn tại điểm t*:
    i | Biểu thức B_i(t) | Giá trị B_i(t*) | Hệ số C_i | Hệ số D_i | Số hạng T_i | Tổng dồn P_i(t*)
    """
    rows = []
    cum_sum = 0.0
    term_vals = []
    for i in range(n + 1):
        b_str = get_stirling_basis_str(i)
        bi_poly = B_list[i]
        b_val = float(sum(c * (eval_t ** p) for p, c in enumerate(bi_poly)))
        di = d_all[i]
        ti = di * b_val
        cum_sum += ti
        term_vals.append(ti)
        rows.append({
            "i": i,
            "Biểu thức B_i(t)": b_str,
            "Giá trị B_i(t*)": b_val,
            "Hệ số C_i": c_all[i],
            "Hệ số D_i": di,
            "Số hạng T_i": ti,
            "Tổng dồn P_i(t*)": cum_sum
        })
    df_eval = pd.DataFrame(rows)
    return df_eval, cum_sum, term_vals

def format_product_form_stirling(d_all):
    """
    Tạo chuỗi công thức dạng tích theo biến chuẩn hóa t (chưa khai triển):
    P(t) = D_0 + D_1*t + D_2*t^2 + D_3*t*(t^2-1) + D_4*t^2*(t^2-1) + ...
    """
    terms = []
    for i, c in enumerate(d_all):
        if i == 0:
            terms.append(f"{c:g}")
            continue
        if abs(c) < 1e-12: continue
        sign = "+ " if c >= 0 and len(terms) > 0 else ("- " if c < 0 and len(terms) > 0 else "")
        val_str = f"{abs(c):g}" if len(terms) > 0 else f"{c:g}"
        b_str = get_stirling_basis_str(i)
        terms.append(f"{sign}{val_str}*{b_str}")
    return " ".join(terms) if terms else "0"

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
    """
    Tính đạo hàm tất cả các cấp k = 1, 2, ..., n của đa thức P(x) = sum(coeffs[j] * x^j).
    """
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

def stirling(points, x0_idx, eval_x=None):
    """
    Hàm tính toán chính cho Nội suy Stirling:
    Trả về toàn bộ bảng biểu, các dãy hệ số và kết quả theo 5 phần chuẩn.
    """
    n = len(points) - 1
    if (n + 1) % 2 == 0:
        raise ValueError("Stirling yêu cầu số mốc lẻ.")

    h = points[1][0] - points[0][0]
    x0_val = points[x0_idx][0]

    # 1. Bảng sai phân hữu hạn
    df_tbl, tbl = build_finite_diff_table(points, x0_idx)

    # 2. Hệ số Stirling (C chẵn, C lẻ, D_i)
    c_all, d_all, c_even, d_even, c_odd, d_odd, details_even, details_odd = extract_stirling_coeffs(tbl, x0_idx, n)

    # 3. Đa thức cơ sở B_i(t) và Bảng tích
    B_list = build_stirling_basis_polynomials(n)
    df_prod_basis = build_product_basis_table_stirling(n, B_list)

    # Khai triển rút gọn đa thức P(t) = sum(D_i * B_i(t))
    P_total_t = np.zeros(n + 1, dtype=float)
    steps_data = []
    for i in range(n + 1):
        di = d_all[i]
        ni = di * B_list[i]
        P_total_t[:len(ni)] += ni
        steps_data.append({
            "i": i,
            "C_i": c_all[i],
            "D_i": di,
            "B_i(t) (0->cao)": [round(v, 8) for v in B_list[i].tolist()],
            "N_i(t) (0->cao)": [round(v, 8) for v in ni.tolist()]
        })
    step_df = pd.DataFrame(steps_data)
    coeff_t_df = pd.DataFrame({"Bậc (t^k)": np.arange(n + 1), "Hệ số": P_total_t})

    # Dạng tích & chuỗi đa thức P(t)
    product_form_str = format_product_form_stirling(d_all)
    expanded_t_str = poly_str(P_total_t, var="t")

    # Quy đổi sang đa thức P(x)
    P_total_x = convert_pt_to_px(P_total_t, x0_val, h)
    coeff_x_df = pd.DataFrame({"Bậc (x^k)": np.arange(n + 1), "Hệ số": P_total_x})
    expanded_x_str = poly_str(P_total_x, var="x")

    # Tính toán tại điểm eval_x (nếu có)
    eval_t = (eval_x - x0_val) / h if eval_x is not None else None
    eval_tbl = None
    val_final = None
    term_vals = None
    derivs_x = []
    derivs_t = []
    if eval_x is not None:
        eval_tbl, val_final, term_vals = build_eval_product_table_stirling(n, c_all, d_all, B_list, eval_t)
        derivs_x = compute_all_derivatives(P_total_x, eval_x)
        derivs_t = compute_all_derivatives(P_total_t, eval_t)

    return {
        "n": n, "h": h, "x0_idx": x0_idx, "x0_val": x0_val,
        "df_tbl": df_tbl,
        "c_all": c_all, "d_all": d_all,
        "c_even": c_even, "d_even": d_even,
        "c_odd": c_odd, "d_odd": d_odd,
        "details_even": details_even, "details_odd": details_odd,
        "df_prod_basis": df_prod_basis,
        "step_df": step_df,
        "coeff_t_df": coeff_t_df, "P_total_t": P_total_t,
        "product_form_str": product_form_str, "expanded_t_str": expanded_t_str,
        "coeff_x_df": coeff_x_df, "P_total_x": P_total_x, "expanded_x_str": expanded_x_str,
        "eval_x": eval_x, "eval_t": eval_t,
        "eval_tbl": eval_tbl, "val_final": val_final, "term_vals": term_vals,
        "derivs_x": derivs_x, "derivs_t": derivs_t
    }

def format_stirling_report(res, msg=""):
    """
    Định dạng báo cáo kết quả Stirling 5 phần hoàn chỉnh.
    """
    lines = []
    lines.append("=" * 90)
    lines.append("KẾT QUẢ NỘI SUY STIRLING")
    lines.append("=" * 90)
    if msg: lines.append(msg)
    lines.append(f"Số mốc tính toán: {res['n'] + 1}, Bước h: {res['h']:.6f}, Mốc trung tâm x0 (idx {res['x0_idx']}): {res['x0_val']}")
    lines.append(f"Chuẩn hóa: t = (x - {res['x0_val']}) / {res['h']:.6f}")
    if res['eval_x'] is not None:
        lines.append(f"Điểm xác định: x* = {res['eval_x']} (t* = {res['eval_t']:.6f})")

    # 1. BẢNG SAI PHÂN HỮU HẠN & QUY TẮC CHỌN SAI PHÂN
    lines.append("\n" + "-" * 90)
    lines.append("1. BẢNG SAI PHÂN HỮU HẠN & QUY TẮC CHỌN SAI PHÂN:")
    lines.append("-" * 90)
    lines.append(res['df_tbl'].to_string(index=True))
    lines.append("\n* QUY TẮC CHỌN SAI PHÂN TRUNG TÂM STIRLING (QUANH MỐC x0):")
    lines.append("  - Bậc chẵn (C_{2k}): Lấy trực tiếp sai phân trên hàng ngang đi qua x0 (Δ^{2k} y_{-k}).")
    lines.append("  - Bậc lẻ (C_{2k-1}): Lấy trung bình cộng của 2 sai phân nằm đối xứng 2 bên hàng x0:")
    lines.append("    C_{2k-1} = (Δ^{2k-1} y_{-k} + Δ^{2k-1} y_{-(k-1)}) / 2")

    # 2. HỆ SỐ NỘI SUY STIRLING (C CHẴN, C LẺ & CÁC DÃY HỆ SỐ)
    lines.append("\n" + "-" * 90)
    lines.append("2. HỆ SỐ NỘI SUY STIRLING (C CHẴN, C LẺ & CÁC DÃY HỆ SỐ):")
    lines.append("-" * 90)
    lines.append("Công thức tính hệ số:")
    lines.append("  - Bậc chẵn (i = 2k):   C_{2k} = Δ^{2k} y_{-k}")
    lines.append("  - Bậc lẻ (i = 2k-1): C_{2k-1} = (Δ^{2k-1} y_{-k} + Δ^{2k-1} y_{-(k-1)}) / 2")
    lines.append("  - Hệ số rút gọn:     D_i = C_i / i!")

    lines.append("\na) Chi tiết các hệ số C chẵn (i = 0, 2, 4, ...):")
    for d in res['details_even']:
        lines.append(f"  * {d['desc']}")
    lines.append(f"  => Dãy hệ số C chẵn: [{', '.join([f'{v:.8g}' for v in res['c_even']])}]")
    lines.append(f"  => Dãy hệ số D chẵn (D_{{2k}} = C_{{2k}} / (2k)!): [{', '.join([f'{v:.8g}' for v in res['d_even']])}]")

    lines.append("\nb) Chi tiết các hệ số C lẻ (i = 1, 3, 5, ...):")
    for d in res['details_odd']:
        lines.append(f"  * {d['desc']}")
    lines.append(f"  => Dãy hệ số C lẻ: [{', '.join([f'{v:.8g}' for v in res['c_odd']])}]")
    lines.append(f"  => Dãy hệ số D lẻ (D_{{2k-1}} = C_{{2k-1}} / (2k-1)!): [{', '.join([f'{v:.8g}' for v in res['d_odd']])}]")

    lines.append("\nc) Tổng hợp toàn bộ các dãy hệ số:")
    lines.append(f"  * Dãy sai phân gốc đầy đủ (C_i): [{', '.join([f'{v:.8g}' for v in res['c_all']])}]")
    lines.append(f"  * Dãy hệ số rút gọn đầy đủ (D_i = C_i / i!): [{', '.join([f'{v:.8g}' for v in res['d_all']])}]")
    lines.append(f"  * Dãy hệ số C chẵn: [{', '.join([f'{v:.8g}' for v in res['c_even']])}]")
    lines.append(f"  * Dãy hệ số C lẻ:   [{', '.join([f'{v:.8g}' for v in res['c_odd']])}]")
    lines.append(f"  * Dãy hệ số D chẵn: [{', '.join([f'{v:.8g}' for v in res['d_even']])}]")
    lines.append(f"  * Dãy hệ số D lẻ:   [{', '.join([f'{v:.8g}' for v in res['d_odd']])}]")

    # 3. BẢNG TÍCH STIRLING
    lines.append("\n" + "-" * 90)
    lines.append("3. BẢNG TÍCH (STIRLING):")
    lines.append("-" * 90)
    lines.append("a) Bảng tích Stirling (Ma trận hệ số đa thức cơ sở tích B_i(t) theo lũy thừa giảm dần của t):")
    lines.append(res['df_prod_basis'].to_string(index=False))

    if res['eval_x'] is not None and res['eval_tbl'] is not None:
        lines.append(f"\nb) Bảng tính giá trị các thừa số tích Stirling và tổng dồn tại t* = {res['eval_t']:.6f}:")
        lines.append(res['eval_tbl'].to_string(index=False))
    else:
        lines.append("\nb) Bảng các bước dựng đa thức cơ sở tích B_i(t):")
        lines.append(res['step_df'].to_string(index=False))

    # 4. ĐA THỨC NỘI SUY
    lines.append("\n" + "-" * 90)
    lines.append("4. ĐA THỨC NỘI SUY (STIRLING):")
    lines.append("-" * 90)
    lines.append(f"a) Dạng tích các thừa số theo biến chuẩn hóa t (chưa khai triển):")
    lines.append(f"P(t) = {res['product_form_str']}")

    lines.append(f"\nb) Bảng hệ số đa thức P(t) sau khi khai triển rút gọn theo t:")
    lines.append(res['coeff_t_df'].to_string(index=False))
    lines.append(f"-> Dãy hệ số theo lũy thừa tăng dần của t [c_0, c_1, ..., c_n]: [{', '.join([f'{v:.8g}' for v in res['P_total_t']])}]")
    lines.append(f"-> Dãy hệ số theo lũy thừa giảm dần của t [c_n, ..., c_1, c_0]: [{', '.join([f'{v:.8g}' for v in reversed(res['P_total_t'])])}]")

    lines.append(f"\nc) Dạng đa thức rút gọn theo lũy thừa của t:")
    lines.append(f"P(t) = {res['expanded_t_str']}")

    lines.append(f"\nd) Bảng hệ số đa thức P(x) quy đổi theo lũy thừa của x (với t = (x - {res['x0_val']})/{res['h']:.6f}):")
    lines.append(res['coeff_x_df'].to_string(index=False))
    lines.append(f"-> Dãy hệ số theo lũy thừa tăng dần của x [a_0, a_1, ..., a_n]: [{', '.join([f'{v:.8g}' for v in res['P_total_x']])}]")
    lines.append(f"-> Dãy hệ số theo lũy thừa giảm dần của x [a_n, ..., a_1, a_0] (dùng cho Horner/Casio): [{', '.join([f'{v:.8g}' for v in reversed(res['P_total_x'])])}]")

    lines.append(f"\ne) Dạng đa thức rút gọn theo lũy thừa của x:")
    lines.append(f"P(x) = {res['expanded_x_str']}")

    # 5. ĐA THỨC NỘI SUY VÀ ĐẠO HÀM TẠI ĐIỂM XÁC ĐỊNH
    lines.append("\n" + "-" * 90)
    lines.append(f"5. ĐA THỨC NỘI SUY VÀ ĐẠO HÀM TẠI ĐIỂM XÁC ĐỊNH{' (x* = ' + str(res['eval_x']) + ')' if res['eval_x'] is not None else ''}:")
    lines.append("-" * 90)
    if res['eval_x'] is not None:
        eval_x = res['eval_x']
        eval_t = res['eval_t']
        lines.append(f"Giá trị biến chuẩn hóa:")
        lines.append(f"   t* = (x* - x0) / h = ({eval_x} - {res['x0_val']}) / {res['h']:.6f} = {eval_t:.6f}")

        terms_display = []
        for t_v in res['term_vals']:
            if len(terms_display) == 0:
                terms_display.append(f"{t_v:.10f}")
            else:
                sign = "+ " if t_v >= 0 else "- "
                terms_display.append(f"{sign}{abs(t_v):.10f}")

        sum_str = " ".join(terms_display)
        lines.append(f"\nKhai triển giá trị tại t* = {eval_t:.6f}:")
        lines.append(f"P(t* = {eval_t:.6f}) = {sum_str}")
        lines.append(f"\n=> GIÁ TRỊ NỘI SUY CUỐI CÙNG:")
        lines.append(f"   P({eval_x}) = {res['val_final']:.10f}")

        derivs_x = res['derivs_x']
        derivs_t = res['derivs_t']
        if derivs_x:
            lines.append(f"\n=> KẾT QUẢ ĐẠO HÀM CÁC CẤP TẠI ĐIỂM x* = {eval_x}:")
            deriv_rows = []
            for (k, p_x, v_x), (_, p_t, v_t) in zip(derivs_x, derivs_t):
                sym = "f'(x*)" if k == 1 else ("f''(x*)" if k == 2 else (f"f'''(x*)" if k == 3 else f"f^({k})(x*)"))
                p_x_str = poly_str(p_x, var="x")
                deriv_rows.append({
                    "Cấp (k)": k, "Ký hiệu": sym,
                    "Giá trị f^(k)(x*)": f"{v_x:.10f}",
                    "Giá trị P^(k)(t*)": f"{v_t:.10f}",
                    "Đa thức đạo hàm P^(k)(x)": p_x_str
                })
            df_deriv = pd.DataFrame(deriv_rows)
            lines.append(df_deriv.to_string(index=False))

            lines.append("\nChi tiết các cấp đạo hàm trọng tâm:")
            prime_map = {1: "'", 2: "''", 3: "'''"}
            for k in range(1, min(len(derivs_x) + 1, 4)):
                (_, p_x, v_x) = derivs_x[k - 1]
                (_, p_t, v_t) = derivs_t[k - 1]
                p_prime = prime_map.get(k, f"^({k})")
                p_x_str = poly_str(p_x, var="x")
                p_t_str = poly_str(p_t, var="t")
                h_power_str = "h" if k == 1 else f"h^{k}"
                lines.append(f"   * Đạo hàm cấp {k} (f{p_prime}({eval_x})):")
                lines.append(f"     + Theo biến chuẩn hóa t: P{p_prime}(t* = {eval_t:.6f}) = {v_t:.10f}")
                lines.append(f"       (Đa thức P{p_prime}(t) = {p_t_str})")
                lines.append(f"     + Quy đổi theo biến x:   f{p_prime}({eval_x}) ≈ P{p_prime}({eval_x}) = P{p_prime}(t*) / {h_power_str} = {v_x:.10f}")
                lines.append(f"       (Đa thức P{p_prime}(x) = {p_x_str})")
    else:
        lines.append("Không có điểm x* (eval_x) được chỉ định trong input.txt hoặc tham số.")

    return "\n".join(lines)

def run_gauss_for_stirling_check(points, x0_idx, eval_x):
    """
    Hàm phụ trợ tính nhanh Gauss I và Gauss II trên cùng tập mốc để kiểm chứng lý thuyết:
    P_stirling(x*) = (P_gauss1(x*) + P_gauss2(x*)) / 2
    """
    xv = [p[0] for p in points]; yv = [p[1] for p in points]; m = len(yv); h = xv[1] - xv[0]
    x0_val = xv[x0_idx]
    eval_t = (eval_x - x0_val) / h

    # Bảng sai phân tiến
    tbl = np.full((m, m), np.nan, dtype=float)
    tbl[:, 0] = yv
    for j in range(1, m):
        for i in range(m - j):
            tbl[i, j] = tbl[i + 1, j - 1] - tbl[i, j - 1]

    # Gauss 1
    # B_i(t): B_0=1, B_1=t, B_2=t(t-1), B_3=t(t-1)(t+1), ...
    # Hệ số hàng: r = x0_idx - i//2
    cum_g1 = 0.0
    cum_prod_g1 = 1.0
    for i in range(m):
        r = x0_idx - (i // 2)
        ci = tbl[r, i]
        di = ci / math.factorial(i)
        if i == 0:
            fact_val = 1.0
        elif i == 1:
            fact_val = eval_t
            cum_prod_g1 *= fact_val
        else:
            node = i // 2 if i % 2 == 0 else -((i - 1) // 2)
            fact_val = eval_t - node
            cum_prod_g1 *= fact_val
        cum_g1 += di * cum_prod_g1

    # Gauss 2
    # B_i(t): B_0=1, B_1=t, B_2=t(t+1), B_3=t(t+1)(t-1), ...
    # Hệ số trên bảng tiến: r = x0_idx - (i+1)//2 nếu i > 0
    cum_g2 = 0.0
    cum_prod_g2 = 1.0
    for i in range(m):
        r = x0_idx if i == 0 else x0_idx - ((i + 1) // 2)
        ci = tbl[r, i]
        di = ci / math.factorial(i)
        if i == 0:
            fact_val = 1.0
        elif i == 1:
            fact_val = eval_t
            cum_prod_g2 *= fact_val
        else:
            node = -((i) // 2) if i % 2 == 0 else ((i - 1) // 2)
            fact_val = eval_t - node
            cum_prod_g2 *= fact_val
        cum_g2 += di * cum_prod_g2

    return cum_g1, cum_g2

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Stirling - Chuẩn 5 phần, in toàn bộ trường hợp & đối chiếu Gauss")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("--x0_index", type=int, help="Chỉ số mốc trung tâm x0")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức (ví dụ 4, 6)")
    p.add_argument("-n", "--n_points", help="Số mốc nội suy (số lẻ: 5, 7 hoặc 'all')")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(9.2, 9.4341319), (9.3, 9.4307764), (9.4, 9.4261142), (9.5, 9.4211191),
                  (9.6, 9.4170553), (9.7, 9.4147476), (9.8, 9.4152900), (9.9, 9.4196762), (10.0, 9.4288617)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    n_pts_cfg = "all" if args.all else (args.n_points or params.get("n_points", None))
    degree_cfg = args.degree or params.get("degree", None)

    sub_points, x0_idx, msg = select_window_for_method(points, eval_x, "stirling", n_pts_cfg, degree_cfg)
    if args.x0_index is not None:
        x0_idx = args.x0_index
    elif "x0_index" in params and len(points) == len(sub_points):
        x0_idx = params["x0_index"]

    if len(sub_points) % 2 == 0:
        sub_points = sub_points[:-1]
        x0_idx = min(x0_idx, len(sub_points) - 1)

    if args.x0_index is None and "x0_index" not in params:
        x0_idx = (len(sub_points) - 1) // 2

    # 1. Trường hợp 1: Stirling với mốc tâm x0 tối ưu
    res = stirling(sub_points, x0_idx, eval_x)
    report_opt = format_stirling_report(res, msg)

    sections = [report_opt]

    # 2. Trường hợp 2: Lựa chọn mốc tâm liền kề (nếu x* nằm giữa 2 mốc)
    h = sub_points[1][0] - sub_points[0][0]
    n_pts = len(sub_points)
    max_half = (n_pts - 1) // 2
    if eval_x is not None:
        # Kiểm tra mốc liền kề
        alt_idx = None
        if eval_x > sub_points[x0_idx][0] and x0_idx + 1 < n_pts - 1:
            alt_idx = x0_idx + 1
        elif eval_x < sub_points[x0_idx][0] and x0_idx - 1 >= 1:
            alt_idx = x0_idx - 1

        # Nếu có thể cắt cửa sổ đối xứng cho mốc liền kề từ points gốc
        if alt_idx is not None:
            # Tìm chỉ số trong points gốc
            alt_val = sub_points[alt_idx][0]
            orig_indices = [i for i, p in enumerate(points) if abs(p[0] - alt_val) < 1e-9]
            if orig_indices:
                orig_k = orig_indices[0]
                m_rad = (len(sub_points) - 1) // 2
                if orig_k - m_rad >= 0 and orig_k + m_rad < len(points):
                    alt_sub = points[orig_k - m_rad : orig_k + m_rad + 1]
                    res_alt = stirling(alt_sub, m_rad, eval_x)
                    alt_t = (eval_x - alt_val) / h
                    alt_lines = []
                    alt_lines.append("\n" + ("#" * 90))
                    alt_lines.append(f"PHẦN B: NỘI SUY STIRLING VỚI MỐC TÂM LIỀN KỀ x0 = {alt_val} (idx {alt_idx})")
                    alt_lines.append("#" * 90)
                    alt_lines.append(f"Điểm x* = {eval_x} nằm giữa {min(sub_points[x0_idx][0], alt_val)} và {max(sub_points[x0_idx][0], alt_val)}.")
                    alt_lines.append(f"Tại tâm x0 = {alt_val}: t* = ({eval_x} - {alt_val}) / {h:.6f} = {alt_t:.6f}")
                    alt_lines.append(f"Đa thức: P_alt(x) = {res_alt['expanded_x_str']}")
                    alt_lines.append(f"=> Giá trị nội suy: P_alt({eval_x}) = {res_alt['val_final']:.10f}")
                    diff_centers = abs(res['val_final'] - res_alt['val_final'])
                    alt_lines.append(f"=> Độ lệch giữa 2 lựa chọn mốc tâm: |P_opt - P_alt| = {diff_centers:.2e}")
                    sections.append("\n".join(alt_lines))

        # 3. Phần C: Đối chiếu kiểm chứng quan hệ giữa Stirling và Gauss
        val_g1, val_g2 = run_gauss_for_stirling_check(sub_points, x0_idx, eval_x)
        avg_gauss = (val_g1 + val_g2) / 2.0
        val_stirl = res['val_final']
        diff_theory = abs(val_stirl - avg_gauss)

        verif_lines = []
        verif_lines.append("\n" + ("#" * 90))
        verif_lines.append("PHẦN C: ĐỐI CHIẾU KIỂM CHỨNG QUAN HỆ GIỮA STIRLING VÀ GAUSS (GAUSS I & GAUSS II)")
        verif_lines.append("#" * 90)
        verif_lines.append(f"Lý thuyết toán học: Công thức nội suy Stirling chính là TRUNG BÌNH CỘNG SỐ HỌC giữa Gauss I (Thuận) và Gauss II (Nghịch):")
        verif_lines.append("   P_stirling(x*) = [ P_gauss1(x*) + P_gauss2(x*) ] / 2")
        verif_lines.append(f"\nKiểm chứng tại điểm x* = {eval_x}:")
        verif_lines.append(f"   - Giá trị theo Gauss I  (Thuận):  P_gauss1({eval_x})  = {val_g1:.10f}")
        verif_lines.append(f"   - Giá trị theo Gauss II (Nghịch): P_gauss2({eval_x})  = {val_g2:.10f}")
        verif_lines.append(f"   - Trung bình cộng (Gauss 1 + Gauss 2) / 2 = {avg_gauss:.10f}")
        verif_lines.append(f"   - Giá trị tính theo Stirling:     P_stirling({eval_x}) = {val_stirl:.10f}")
        verif_lines.append(f"   => Độ lệch kiểm chứng: |P_stirling - (P_gauss1 + P_gauss2)/2| = {diff_theory:.2e}")
        if diff_theory < 1e-9:
            verif_lines.append("   => KẾT LUẬN: Đẳng thức toán học hoàn toàn chính xác! Stirling là trung bình cộng của 2 đa thức Gauss.")
        sections.append("\n".join(verif_lines))

    final_output = "\n".join(sections)
    print("\n" + final_output)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_output + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
