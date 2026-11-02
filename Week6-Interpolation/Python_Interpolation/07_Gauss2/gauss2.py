"""
Nội suy Gauss II (Gauss nghịch) & Gauss I (Gauss thuận)
Xuất kết quả đúng theo 5 phần chuẩn tương tự Newton Divided / Newton Finite:
1. BẢNG SAI PHÂN HỮU HẠN (FINITE DIFFERENCE TABLE - CẢ TIẾN VÀ LÙI)
2. HỆ SỐ GAUSS (GAUSS COEFFICIENTS: C_i & D_i = C_i / i!)
3. BẢNG TÍCH GAUSS (PRODUCT BASIS TABLE: MA TRẬN THEO t & BẢNG TÍNH TẠI t*)
4. ĐA THỨC NỘI SUY GAUSS (GAUSS POLYNOMIAL THEO BIẾN t VÀ BIẾN x)
5. ĐA THỨC NỘI SUY TẠI ĐIỂM XÁC ĐỊNH (VALUE AT EVALUATION POINT & DERIVATIVES)

Hỗ trợ chạy riêng Gauss II (-t gauss2) hoặc Gauss I (-t gauss1).
Mặc định: tính đồng thời cả GAUSS II (NGHỊCH) và GAUSS I (THUẬN), sau đó TỔNG KẾT & SO SÁNH kiểm chứng tính đồng nhất.
Nhập: input.txt hoặc file đề thi (.txt, .csv, .xlsx)
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

def build_finite_diff_tables(points, x0_idx=None):
    """
    Xây dựng bảng sai phân hữu hạn dạng ma trận:
    - Bảng tam giác trên (dạng sai phân tiến Δ): dùng cho Gauss I (Thuận).
    - Bảng tam giác dưới (dạng sai phân lùi ∇): dùng cho Gauss II (Nghịch).
    Cột t_i: Chỉ số tương đối so với mốc trung tâm x0 (..., -2, -1, 0, 1, 2, ...)
    """
    xv = [p[0] for p in points]
    yv = [p[1] for p in points]
    n = len(yv)
    if x0_idx is None:
        x0_idx = (n - 1) // 2

    t_indices = [i - x0_idx for i in range(n)]

    tbl_upper = np.full((n, n), np.nan, dtype=float)
    tbl_upper[:, 0] = yv
    for j in range(1, n):
        for i in range(n - j):
            tbl_upper[i, j] = tbl_upper[i + 1, j - 1] - tbl_upper[i, j - 1]
    data_up = {"t_i": t_indices, "x_i": xv, "y_i": yv}
    for j in range(1, n):
        data_up[f"Cấp {j}"] = tbl_upper[:, j]
    df_upper = pd.DataFrame(data_up)

    tbl_lower = np.full((n, n), np.nan, dtype=float)
    for i in range(n):
        for j in range(i + 1):
            tbl_lower[i, j] = tbl_upper[i - j, j]
    data_low = {"t_i": t_indices, "x_i": xv, "y_i": tbl_lower[:, 0]}
    for j in range(1, n):
        data_low[f"Cấp {j}"] = tbl_lower[:, j]
    df_lower = pd.DataFrame(data_low)

    return df_upper, tbl_upper, df_lower, tbl_lower

def extract_gauss_coeffs(tbl_upper, tbl_lower, x0_idx, m, mode="gauss2"):
    """
    Trích xuất các hệ số sai phân C_i và hệ số Newton/Gauss D_i = C_i / i!
    Đồng thời trả về thông tin đối chiếu trên cả 2 bảng Tiến (Δ) và Lùi (∇):
    - Gauss I (Thuận):
        + Trên bảng tiến Δ: hàng r_delta = x0_idx - i // 2
        + Trên bảng lùi ∇: hàng r_nabla = x0_idx + (i + 1) // 2
    - Gauss II (Nghịch):
        + Trên bảng lùi ∇: hàng r_nabla = x0_idx + i // 2
        + Trên bảng tiến Δ: hàng r_delta = x0_idx - (i + 1) // 2
    """
    c_coeffs = []
    d_coeffs = []
    details = []
    for i in range(m):
        fact_i = math.factorial(i)
        if mode == "gauss1":
            r_delta = x0_idx - (i // 2)
            r_nabla = x0_idx + ((i + 1) // 2) if i > 0 else x0_idx
            val_delta = tbl_upper[r_delta, i]
            val = val_delta
            if i == 0:
                sym_delta = "y_0"
                sym_nabla = "y_0"
            else:
                sub_delta = -(i // 2) if i % 2 == 0 else -((i - 1) // 2)
                sub_delta_str = f"y_{{{sub_delta}}}" if sub_delta != 0 else "y_0"
                sym_delta = f"Δ^{i} {sub_delta_str}" if i > 1 else f"Δ {sub_delta_str}"

                sub_nabla = (i + 1) // 2 if i % 2 == 1 else i // 2
                sym_nabla = f"∇^{i} y_{{{sub_nabla}}}" if i > 1 else f"∇ y_{{{sub_nabla}}}"
        else:
            r_nabla = x0_idx + (i // 2)
            r_delta = x0_idx - ((i + 1) // 2) if i > 0 else x0_idx
            val_nabla = tbl_lower[r_nabla, i]
            val = val_nabla
            if i == 0:
                sym_nabla = "y_0"
                sym_delta = "y_0"
            else:
                sub_nabla = i // 2 if i % 2 == 0 else (i - 1) // 2
                sub_nabla_str = f"y_{{{sub_nabla}}}" if sub_nabla != 0 else "y_0"
                sym_nabla = f"∇^{i} {sub_nabla_str}" if i > 1 else f"∇ {sub_nabla_str}"

                sub_delta = -((i + 1) // 2) if i % 2 == 1 else -(i // 2)
                sym_delta = f"Δ^{i} y_{{{sub_delta}}}" if i > 1 else f"Δ y_{{{sub_delta}}}"

        c_coeffs.append(val)
        d_coeffs.append(val / fact_i)
        details.append({
            "i": i,
            "val": val,
            "fact": fact_i,
            "r_delta": r_delta,
            "t_delta": r_delta - x0_idx,
            "sym_delta": sym_delta,
            "r_nabla": r_nabla,
            "t_nabla": r_nabla - x0_idx,
            "sym_nabla": sym_nabla,
        })
    return c_coeffs, d_coeffs, details

def build_gauss_polynomial(m, d_coeffs, mode="gauss2"):
    """
    Xây dựng các đa thức cơ sở B_i(t) và đa thức rút gọn P(t) theo biến chuẩn hóa t.
    """
    sdata = []
    P_total_t = np.zeros(m, dtype=float)
    B_poly = np.array([1.0], dtype=float)

    for i in range(m):
        Di = d_coeffs[i]
        if i == 0:
            added_factor = "1"
            Bc = np.array([1.0], dtype=float)
        elif i == 1:
            added_factor = "t"
            Bc = np.array([0.0, 1.0], dtype=float)
        else:
            if mode == "gauss1":
                k = i // 2 if i % 2 == 0 else (i - 1) // 2
                conv_term = np.array([-float(k), 1.0]) if i % 2 == 0 else np.array([float(k), 1.0])
                added_factor = f"(t - {k})" if i % 2 == 0 else f"(t + {k})"
            else:
                k = (i - 1) // 2 if i % 2 == 1 else i // 2
                conv_term = np.array([float(k), 1.0]) if i % 2 == 0 else np.array([-float(k), 1.0])
                added_factor = f"(t + {k})" if i % 2 == 0 else f"(t - {k})"
            Bc = np.convolve(B_poly, conv_term)

        B_poly = Bc.copy()
        Ni = Di * B_poly
        P_total_t[:len(Ni)] += Ni

        sdata.append({
            "i": i,
            "D_i": Di,
            "Nhân thêm": added_factor,
            "Tích B_i(t) [0->cao]": [round(v, 8) for v in B_poly.tolist()],
            "Số hạng N_i(t) [0->cao]": [round(v, 8) for v in Ni.tolist()]
        })

    step_df = pd.DataFrame(sdata)
    coeff_t_df = pd.DataFrame({"Bậc (t^k)": np.arange(m), "Hệ số": P_total_t})
    return step_df, coeff_t_df, P_total_t

def convert_pt_to_px(p_t_coeffs, x_base, h):
    m = len(p_t_coeffs)
    p_x = np.zeros(m, dtype=float)
    t_lin = np.array([-x_base / h, 1.0 / h], dtype=float)
    t_pow = np.array([1.0], dtype=float)
    for k in range(m):
        contrib = p_t_coeffs[k] * t_pow
        p_x[:len(contrib)] += contrib
        t_pow = np.convolve(t_pow, t_lin)
    return p_x

def build_product_basis_table_gauss(m, mode="gauss2"):
    rows = []
    poly = np.array([1.0], dtype=float)

    r0 = np.zeros(m, dtype=float)
    r0[-1] = 1.0
    rows.append(r0)

    for i in range(1, m):
        if mode == "gauss1":
            node_t = 0 if i == 1 else (i // 2 if i % 2 == 0 else -((i - 1) // 2))
        else:
            node_t = 0 if i == 1 else (-((i) // 2) if i % 2 == 0 else ((i - 1) // 2))
        poly = np.convolve(poly, [1.0, -float(node_t)])
        r = np.zeros(m, dtype=float)
        r[m - len(poly):] = poly
        rows.append(r)

    def fmt_val(v):
        if abs(v) < 1e-11: return "0"
        if abs(v - round(v)) < 1e-11: return str(int(round(v)))
        formatted = f"{v:.8f}".rstrip("0").rstrip(".")
        return "0" if formatted == "-0" else formatted

    formatted_rows = []
    for i in range(m):
        if i == 0:
            moc_str = ""
        elif i == 1:
            moc_str = "0"
        else:
            if mode == "gauss1":
                moc_str = f"{i // 2}" if i % 2 == 0 else f"{-((i - 1) // 2)}"
            else:
                moc_str = f"{-((i) // 2)}" if i % 2 == 0 else f"{(i - 1) // 2}"
        row_strs = [moc_str] + [fmt_val(v) for v in rows[i]]
        formatted_rows.append(row_strs)

    cols = ["Mốc"] + [f"t^{p}" for p in range(m - 1, -1, -1)]
    df_prod = pd.DataFrame(formatted_rows, columns=cols)
    return df_prod, np.array(rows)

def build_eval_product_table_gauss(m, d_coeffs, eval_t, mode="gauss2"):
    rows = []
    cum_prod = 1.0
    cum_sum = 0.0

    name_prod = "Tích Pi_i(t*)"
    col_node = "Mốc nhân thêm"
    col_factor = "Thừa số nhân thêm"
    col_coeff = "Hệ số D_i" if mode == "gauss1" else "Hệ số D'_i"
    term_values = []

    for i in range(m):
        Di = d_coeffs[i]
        if i == 0:
            node_str = "-"
            factor_str = "1"
            val_factor = 1.0
        elif i == 1:
            node_str = "0"
            val_factor = eval_t
            factor_str = f"({eval_t:.6g} - 0) = {val_factor:.6g}"
            cum_prod *= val_factor
        else:
            if mode == "gauss1":
                node_t = i // 2 if i % 2 == 0 else -((i - 1) // 2)
            else:
                node_t = -((i) // 2) if i % 2 == 0 else ((i - 1) // 2)

            node_str = f"{node_t}"
            val_factor = eval_t - node_t
            sign_part = f"- {node_t}" if node_t >= 0 else f"+ {abs(node_t)}"
            factor_str = f"({eval_t:.6g} {sign_part}) = {val_factor:.6g}"
            cum_prod *= val_factor

        term_val = Di * cum_prod
        cum_sum += term_val
        term_values.append(term_val)

        rows.append({
            "i": i,
            col_node: node_str,
            col_factor: factor_str,
            name_prod: cum_prod,
            col_coeff: Di,
            "Số hạng T_i": term_val,
            "Tổng dồn P_i(t*)": cum_sum
        })

    return pd.DataFrame(rows), cum_sum, term_values

def format_product_form_gauss(d_coeffs, mode="gauss2"):
    m = len(d_coeffs)
    terms = []
    factors = []

    for i in range(m):
        c = d_coeffs[i]
        if i == 0:
            terms.append(f"{c:g}")
        elif i == 1:
            factors.append("t")
            prod_str = "t"
            sign = "+ " if c >= 0 and len(terms) > 0 else ("- " if c < 0 and len(terms) > 0 else "")
            val_str = f"{abs(c):g}" if len(terms) > 0 else f"{c:g}"
            terms.append(f"{sign}{val_str}*{prod_str}")
        else:
            if mode == "gauss1":
                node = i // 2 if i % 2 == 0 else -((i - 1) // 2)
                f_str = f"(t - {node})" if node > 0 else f"(t + {abs(node)})"
            else:
                node = -((i) // 2) if i % 2 == 0 else ((i - 1) // 2)
                f_str = f"(t + {abs(node)})" if node < 0 else f"(t - {node})"
            factors.append(f_str)
            prod_str = "*".join(factors)
            sign = "+ " if c >= 0 and len(terms) > 0 else ("- " if c < 0 and len(terms) > 0 else "")
            val_str = f"{abs(c):g}" if len(terms) > 0 else f"{c:g}"
            terms.append(f"{sign}{val_str}*{prod_str}")

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

def run_gauss(sub_points, x0_idx, eval_x=None, mode="gauss2", section_prefix=""):
    m = len(sub_points)
    h = sub_points[1][0] - sub_points[0][0]
    x0_val = sub_points[x0_idx][0]
    mode_name = "GAUSS II (GAUSS NGHỊCH)" if mode == "gauss2" else "GAUSS I (GAUSS THUẬN)"
    m_name_short = "Gauss II" if mode == "gauss2" else "Gauss I"

    df_upper, tbl_upper, df_lower, tbl_lower = build_finite_diff_tables(sub_points, x0_idx=x0_idx)
    c_coeffs, d_coeffs, details = extract_gauss_coeffs(tbl_upper, tbl_lower, x0_idx, m, mode=mode)

    step_df, coeff_t_df, P_total_t = build_gauss_polynomial(m, d_coeffs, mode=mode)
    product_form_str = format_product_form_gauss(d_coeffs, mode=mode)
    expanded_t_str = poly_str(P_total_t, var="t")

    P_total_x = convert_pt_to_px(P_total_t, x0_val, h)
    coeff_x_df = pd.DataFrame({"Bậc (x^k)": np.arange(m), "Hệ số": P_total_x})
    expanded_x_str = poly_str(P_total_x, var="x")

    eval_t = (eval_x - x0_val) / h if eval_x is not None else None

    lines = []
    lines.append("=" * 90)
    lines.append(f"KẾT QUẢ NỘI SUY {section_prefix}{mode_name}")
    lines.append("=" * 90)
    lines.append(f"Số mốc tính toán: {m}, Bước h: {h:.6f}, Mốc trung tâm x0 (idx {x0_idx}): {x0_val}")
    lines.append(f"Chuẩn hóa: t = (x - {x0_val}) / {h:.6f}")
    if eval_x is not None:
        lines.append(f"Điểm xác định: x* = {eval_x} (t* = {eval_t:.6f})")

    lines.append("\n" + "-" * 90)
    lines.append(f"1. BẢNG SAI PHÂN HỮU HẠN & ĐƯỜNG ĐI ZIG-ZAG (CẢ TIẾN VÀ LÙI):")
    lines.append("-" * 90)
    lines.append("a) BẢNG SAI PHÂN TIẾN Δ (DẠNG TAM GIÁC TRÊN):")
    lines.append(df_upper.to_string(index=True))

    lines.append("\n\nb) BẢNG SAI PHÂN LÙI ∇ (DẠNG TAM GIÁC DƯỚI):")
    lines.append(df_lower.to_string(index=True))

    lines.append("\n\nc) QUY TẮC CHỌN SAI PHÂN & ĐƯỜNG ĐI ZIG-ZAG:")
    if mode == "gauss2":
        lines.append("   * Trên bảng LÙI ∇ (Tam giác dưới) - Bước đầu đi LÊN theo chiều âm (t < 0):")
        lines.append(f"     y_0 -> ∇ y_0 -> ∇^2 y_1 -> ∇^3 y_1 -> ∇^4 y_2 -> ∇^5 y_2 ... (hàng {x0_idx} + i // 2)")
        lines.append("   * Tương đương trên bảng TIẾN Δ (Tam giác trên):")
        lines.append(f"     y_0 -> Δ y_{{-1}} -> Δ^2 y_{{-1}} -> Δ^3 y_{{-2}} -> Δ^4 y_{{-2}} -> Δ^5 y_{{-3}} ... (hàng {x0_idx} - (i + 1) // 2)")
        lines.append("   => Hai cách tra cứu trên bảng LÙI và TIẾN cho ra các hệ số C'_i HOÀN TOÀN TRÙNG KHỚP NHAU!")
    else:
        lines.append("   * Trên bảng TIẾN Δ (Tam giác trên) - Bước đầu đi XUỐNG theo chiều dương (t > 0):")
        lines.append(f"     y_0 -> Δ y_0 -> Δ^2 y_{{-1}} -> Δ^3 y_{{-1}} -> Δ^4 y_{{-2}} -> Δ^5 y_{{-2}} ... (hàng {x0_idx} - i // 2)")
        lines.append("   * Tương đương trên bảng LÙI ∇ (Tam giác dưới):")
        lines.append(f"     y_0 -> ∇ y_1 -> ∇^2 y_1 -> ∇^3 y_2 -> ∇^4 y_2 -> ∇^5 y_3 ... (hàng {x0_idx} + (i + 1) // 2)")
        lines.append("   => Hai cách tra cứu trên bảng TIẾN và LÙI cho ra các hệ số C_i HOÀN TOÀN TRÙNG KHỚP NHAU!")

    lines.append("\n" + "-" * 90)
    coeff_prefix = "D'" if mode == "gauss2" else "D"
    label_coeffs = f"{coeff_prefix}_0, {coeff_prefix}_1, ..., {coeff_prefix}_n"
    lines.append(f"2. HỆ SỐ NỘI SUY ({label_coeffs}) - ĐỐI CHIẾU CẢ SAI PHÂN TIẾN VÀ LÙI:")
    lines.append("-" * 90)
    lines.append(f"Công thức tính hệ số: {coeff_prefix}_i = C_i / i!")
    coeff_lines = []
    for d in details:
        j = d["i"]
        c_val = d["val"]
        fact_j = d["fact"]
        d_val = c_val / fact_j
        s_del = d["sym_delta"]
        r_del = d["r_delta"]
        t_del = d["t_delta"]
        s_nab = d["sym_nabla"]
        r_nab = d["r_nabla"]
        t_nab = d["t_nabla"]
        sign_tdel = f"{t_del:+d}"
        sign_tnab = f"{t_nab:+d}"
        if mode == "gauss2":
            coeff_lines.append(f"  {coeff_prefix}_{j} = {s_nab} / {j}! = {s_del} / {j}! = {c_val:.10g} / {fact_j} = {d_val:.10g}  (Lùi: t = {sign_tnab} [hàng {r_nab}] | Tiến: t = {sign_tdel} [hàng {r_del}])")
        else:
            coeff_lines.append(f"  {coeff_prefix}_{j} = {s_del} / {j}! = {s_nab} / {j}! = {c_val:.10g} / {fact_j} = {d_val:.10g}  (Tiến: t = {sign_tdel} [hàng {r_del}] | Lùi: t = {sign_tnab} [hàng {r_nab}])")
    lines.extend(coeff_lines)
    lines.append(f"\n-> Dãy sai phân gốc (C_i): [{', '.join([f'{v:.8g}' for v in c_coeffs])}]")
    lines.append(f"-> Dãy hệ số rút gọn ({coeff_prefix}_i): [{', '.join([f'{v:.8g}' for v in d_coeffs])}]")

    lines.append("\n" + "-" * 90)
    lines.append(f"3. BẢNG TÍCH ({m_name_short.upper()}):")
    lines.append("-" * 90)
    df_prod_basis, _ = build_product_basis_table_gauss(m, mode=mode)
    lines.append(f"Bảng tích {m_name_short} (Ma trận hệ số đa thức tích theo lũy thừa giảm dần của t):")
    lines.append(df_prod_basis.to_string(index=False))

    val_final = None
    dict_derivs_x = {}
    if eval_x is not None:
        eval_tbl, val_final, term_vals = build_eval_product_table_gauss(m, d_coeffs, eval_t, mode=mode)
        lines.append(f"\n(Bảng tính giá trị các thừa số tích {m_name_short} và tổng dồn tại t* = {eval_t:.6f}):")
        lines.append(eval_tbl.to_string(index=False))
    else:
        lines.append(f"\n(Bảng các đa thức cơ sở tích {m_name_short} B_i(t)):")
        lines.append(step_df[["i", "D_i", "Nhân thêm", "Tích B_i(t) [0->cao]"]].to_string(index=False))

    lines.append("\n" + "-" * 90)
    lines.append(f"4. ĐA THỨC NỘI SUY ({m_name_short.upper()}):")
    lines.append("-" * 90)
    lines.append(f"a) Dạng tích các thừa số theo biến chuẩn hóa t (chưa khai triển):")
    lines.append(f"P_{mode}(t) = {product_form_str}")

    lines.append(f"\nb) Bảng hệ số đa thức P_{mode}(t) sau khi khai triển rút gọn theo t:")
    lines.append(coeff_t_df.to_string(index=False))
    lines.append(f"-> Dãy hệ số theo lũy thừa tăng dần của t [c_0, c_1, ..., c_n]: [{', '.join([f'{v:.8g}' for v in P_total_t])}]")

    lines.append(f"\nc) Dạng đa thức rút gọn theo lũy thừa của t:")
    lines.append(f"P_{mode}(t) = {expanded_t_str}")

    lines.append(f"\nd) Bảng hệ số đa thức P_{mode}(x) quy đổi theo lũy thừa của x (với t = (x - {x0_val})/{h:.6f}):")
    lines.append(coeff_x_df.to_string(index=False))
    lines.append(f"-> Dãy hệ số theo lũy thừa tăng dần [a_0, a_1, ..., a_n]: [{', '.join([f'{v:.8g}' for v in P_total_x])}]")
    lines.append(f"-> Dãy hệ số theo lũy thừa giảm dần [a_n, ..., a_1, a_0] (dùng cho Horner/Casio): [{', '.join([f'{v:.8g}' for v in reversed(P_total_x)])}]")

    lines.append(f"\ne) Dạng đa thức rút gọn theo lũy thừa của x:")
    lines.append(f"P_{mode}(x) = {expanded_x_str}")

    lines.append("\n" + "-" * 90)
    lines.append(f"5. ĐA THỨC NỘI SUY TẠI ĐIỂM XÁC ĐỊNH{' (x* = ' + str(eval_x) + ')' if eval_x is not None else ''}:")
    lines.append("-" * 90)
    if eval_x is not None:
        eval_tbl, val_eval, term_vals = build_eval_product_table_gauss(m, d_coeffs, eval_t, mode=mode)
        lines.append(f"Giá trị biến chuẩn hóa:")
        lines.append(f"   t* = (x* - x0) / h = ({eval_x} - {x0_val}) / {h:.6f} = {eval_t:.6f}")

        terms_display = []
        for t_v in term_vals:
            if len(terms_display) == 0:
                terms_display.append(f"{t_v:.10f}")
            else:
                sign = "+ " if t_v >= 0 else "- "
                terms_display.append(f"{sign}{abs(t_v):.10f}")

        sum_str = " ".join(terms_display)
        lines.append(f"\nKhai triển giá trị tại t* = {eval_t:.6f}:")
        lines.append(f"P_{mode}(t* = {eval_t:.6f}) = {sum_str}")
        lines.append(f"\n=> GIÁ TRỊ NỘI SUY CUỐI CÙNG:")
        lines.append(f"   P_{mode}({eval_x}) = {val_eval:.10f}")

        derivs_x = compute_all_derivatives(P_total_x, eval_x)
        derivs_t = compute_all_derivatives(P_total_t, eval_t)
        dict_derivs_x = {k: v for k, p, v in derivs_x}

        lines.append(f"\n=> KẾT QUẢ ĐẠO HÀM CÁC CẤP TẠI ĐIỂM x* = {eval_x}:")
        deriv_table_rows = []
        for (k, p_x, v_x), (_, p_t, v_t) in zip(derivs_x, derivs_t):
            symbol = "f'(x*)" if k == 1 else ("f''(x*)" if k == 2 else (f"f'''(x*)" if k == 3 else f"f^({k})(x*)"))
            p_x_str = poly_str(p_x, var="x")
            deriv_table_rows.append({
                "Cấp (k)": k,
                "Ký hiệu": symbol,
                "Giá trị f^(k)(x*)": f"{v_x:.10f}",
                "Giá trị P^(k)(t*)": f"{v_t:.10f}",
                "Đa thức đạo hàm P^(k)(x)": p_x_str
            })
        df_deriv_summary = pd.DataFrame(deriv_table_rows)
        lines.append(df_deriv_summary.to_string(index=False))

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
            lines.append(f"     + Theo biến chuẩn hóa t: P{p_prime}_{mode}(t* = {eval_t:.6f}) = {v_t:.10f}")
            lines.append(f"       (Đa thức P{p_prime}(t) = {p_t_str})")
            lines.append(f"     + Quy đổi theo biến x:   f{p_prime}({eval_x}) ≈ P{p_prime}_{mode}({eval_x}) = P{p_prime}(t*) / {h_power_str} = {v_x:.10f}")
            lines.append(f"       (Đa thức P{p_prime}(x) = {p_x_str})")
    else:
        lines.append("Không có điểm x* (eval_x) được chỉ định trong input.txt hoặc tham số.")

    return "\n".join(lines), val_final, dict_derivs_x

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Nội suy Gauss II (Nghịch) & Gauss I (Thuận) - Chuẩn 5 phần")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-t", "--type", choices=["gauss1", "gauss2", "both", "1", "2", "thuan", "nghich"], help="Loại nội suy (gauss2/gauss1/both)")
    p.add_argument("-x", "--eval_x", type=float, help="Điểm cần tính x*")
    p.add_argument("--x0_index", type=int, help="Chỉ số mốc trung tâm x0")
    p.add_argument("-d", "--degree", type=int, help="Bậc đa thức mong muốn (số chẵn: 4, 6)")
    p.add_argument("-n", "--n_points", help="Số mốc nội suy (số lẻ: 5, 7 hoặc 'all')")
    p.add_argument("--all", action="store_true", help="Dùng toàn bộ mốc dữ liệu")
    args = p.parse_args()

    points, params = parse_universal_input(args.input_file)
    if not points:
        points = [(9.2, 9.4341319), (9.3, 9.4307764), (9.4, 9.4261142), (9.5, 9.4211191),
                  (9.6, 9.4170553), (9.7, 9.4147476), (9.8, 9.4152900), (9.9, 9.4196762), (10.0, 9.4288617)]

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)
    cfg_n = "all" if args.all else (args.n_points or params.get("n_points", None))
    cfg_deg = args.degree or params.get("degree", None)

    sub_points, x0_idx, msg = select_window_for_method(points, eval_x, "gauss2", cfg_n, cfg_deg)
    if args.x0_index is not None:
        x0_idx = args.x0_index
    elif "x0_index" in params and len(points) == len(sub_points):
        x0_idx = params["x0_index"]

    if len(sub_points) % 2 == 0:
        sub_points = sub_points[:-1]
        x0_idx = min(x0_idx, len(sub_points) - 1)

    if args.x0_index is None and "x0_index" not in params:
        x0_idx = (len(sub_points) - 1) // 2

    mode_arg = args.type or params.get("type", "both")
    mode_arg = mode_arg.lower() if mode_arg else "both"

    if ("2" in mode_arg or "nghich" in mode_arg) and "both" not in mode_arg and "1" not in mode_arg and "thuan" not in mode_arg:
        final_output, _, _ = run_gauss(sub_points, x0_idx, eval_x, mode="gauss2")
        if msg:
            final_output = msg + "\n\n" + final_output
    elif ("1" in mode_arg or "thuan" in mode_arg) and "both" not in mode_arg and "2" not in mode_arg and "nghich" not in mode_arg:
        final_output, _, _ = run_gauss(sub_points, x0_idx, eval_x, mode="gauss1")
        if msg:
            final_output = msg + "\n\n" + final_output
    else:
        # In cả GAUSS II (NGHỊCH) và GAUSS I (THUẬN)
        res_g2, val_g2, d_g2 = run_gauss(sub_points, x0_idx, eval_x, mode="gauss2", section_prefix="PHẦN A: ")
        res_g1, val_g1, d_g1 = run_gauss(sub_points, x0_idx, eval_x, mode="gauss1", section_prefix="PHẦN B: ")

        summary_lines = []
        summary_lines.append("=" * 90)
        summary_lines.append("PHẦN C: TỔNG KẾT & SO SÁNH (KIỂM CHỨNG TÍNH ĐỒNG NHẤT GIỮA GAUSS II VÀ GAUSS I)")
        summary_lines.append("=" * 90)
        if eval_x is not None and val_g1 is not None and val_g2 is not None:
            summary_lines.append(f"1. So sánh tại điểm x* = {eval_x}:")
            summary_lines.append(f"   a) Giá trị hàm f({eval_x}):")
            summary_lines.append(f"      - Gauss II (Nghịch): P_gauss2({eval_x}) = {val_g2:.10f}")
            summary_lines.append(f"      - Gauss I  (Thuận):  P_gauss1({eval_x}) = {val_g1:.10f}")
            diff_val = abs(val_g1 - val_g2)
            summary_lines.append(f"      => Độ lệch tuyệt đối: |P_gauss2 - P_gauss1| = {diff_val:.2e}")

            prime_map = {1: "'", 2: "''", 3: "'''"}
            prime_names = {1: "Đạo hàm cấp 1 f'", 2: "Đạo hàm cấp 2 f''", 3: "Đạo hàm cấp 3 f'''"}
            sub_lbls = ["b", "c", "d", "e", "f"]
            for idx, k in enumerate([1, 2, 3]):
                if k in d_g1 and k in d_g2:
                    lbl = sub_lbls[idx]
                    t_name = prime_names.get(k, f"Đạo hàm cấp {k} f^({k})")
                    s_sym = prime_map.get(k, f"^({k})")
                    summary_lines.append(f"\n   {lbl}) {t_name}({eval_x}):")
                    summary_lines.append(f"      - Gauss II (Nghịch): P{s_sym}_gauss2({eval_x}) = {d_g2[k]:.10f}")
                    summary_lines.append(f"      - Gauss I  (Thuận):  P{s_sym}_gauss1({eval_x}) = {d_g1[k]:.10f}")
                    diff_k = abs(d_g1[k] - d_g2[k])
                    summary_lines.append(f"      => Độ lệch tuyệt đối: |P{s_sym}_gauss2 - P{s_sym}_gauss1| = {diff_k:.2e}")

            if diff_val < 1e-9:
                summary_lines.append("\n   => Kết luận: Cả giá trị nội suy và các cấp đạo hàm từ 2 hướng NGHỊCH và THUẬN hoàn toàn trùng khớp!")
        summary_lines.append("\n2. Nhận xét lý thuyết:")
        summary_lines.append("   - Đa thức nội suy P(x) là duy nhất dù xây dựng theo Gauss II hay Gauss I trên cùng tập mốc đối xứng quanh x0.")
        summary_lines.append("   - Khác biệt quy trình chọn sai phân:")
        summary_lines.append("     + GAUSS II (Nghịch): Đường đi sai phân zig-zag đi lên trước qua mốc x0 (y0 -> Δy_{-1} -> Δ^2 y_{-1} -> Δ^3 y_{-2} ...), nhân thừa số t(t+1)(t-1)(t+2)... (tối ưu khi -0.5 < t < 0).")
        summary_lines.append("     + GAUSS I (Thuận): Đường đi sai phân zig-zag đi xuống trước qua mốc x0 (y0 -> Δy0 -> Δ^2 y_{-1} -> Δ^3 y_{-1} ...), nhân thừa số t(t-1)(t+1)(t-2)... (tối ưu khi 0 < t < 0.5).")
        summary_lines.append("     + Cả hai phương pháp đều có thể tra cứu tương đương trên cả Bảng sai phân tiến Δ (tam giác trên) và Bảng sai phân lùi ∇ (tam giác dưới).")

        sep = "\n\n" + ("#" * 90) + "\n\n"
        final_output = res_g2 + sep + res_g1 + "\n\n" + "\n".join(summary_lines)
        if msg:
            final_output = msg + "\n\n" + final_output

    final_output = format_report_output(final_output, args.input_file)
    print("\n" + final_output)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(final_output + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
