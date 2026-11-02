"""
Lược đồ Horner (Horner Division Table & Repeated Horner for Derivatives)
Nhập: input.txt hoặc tham số dòng lệnh
Xuất: output.txt và in ra màn hình
"""
import os, sys, math, argparse
import numpy as np
import pandas as pd

try:
    from pps_io import format_report_output
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pps_io import format_report_output

def poly_str_from_coeffs(coeffs_high_to_low, var="x"):
    """
    Tạo chuỗi đa thức đẹp từ dãy hệ số bậc giảm dần [a_n, a_(n-1), ..., a_0]
    """
    n = len(coeffs_high_to_low) - 1
    terms = []
    for idx, c in enumerate(coeffs_high_to_low):
        deg = n - idx
        if abs(c) < 1e-12:
            continue
        sign = "+ " if c >= 0 and len(terms) > 0 else ("- " if c < 0 and len(terms) > 0 else "")
        val = abs(c) if len(terms) > 0 else c
        if deg == 0:
            term = f"{val:.8g}"
        elif deg == 1:
            term = f"{val:.8g}*{var}"
        else:
            term = f"{val:.8g}*{var}^{deg}"
        terms.append(f"{sign}{term}")
    return " ".join(terms) if terms else "0"

def synthetic_division(a_coeffs_desc, c: float):
    """
    Thực hiện phép chia Horner 1 lần cho đa thức P(x) tại điểm c.
    a_coeffs_desc: danh sách hệ số từ bậc cao xuống bậc thấp [a_n, a_(n-1), ..., a_0]
    c: điểm tính giá trị / nghiệm thử
    """
    n = len(a_coeffs_desc) - 1
    # a_by_deg[i] là hệ số của x^i
    a_by_deg = list(reversed(a_coeffs_desc))
    b = [0.0] * (n + 1)
    added_terms = ["-"] * (n + 1)

    b[n] = a_by_deg[n]
    for i in range(n - 1, -1, -1):
        add_val = c * b[i + 1]
        added_terms[i] = add_val
        b[i] = a_by_deg[i] + add_val

    # 1. Bảng chi tiết dạng dọc
    rows_vertical = []
    for i in range(n, -1, -1):
        add_str = "-" if i == n else f"{added_terms[i]:.8g}"
        rows_vertical.append({
            "Bậc (i)": i,
            "Hệ số a_i": a_by_deg[i],
            f"c * b_(i+1) (với c={c:g})": add_str,
            "b_i = a_i + c*b_(i+1)": b[i]
        })
    df_vertical = pd.DataFrame(rows_vertical)

    # 2. Bảng ngang truyền thống (Horizontal Horner Table)
    deg_labels = [f"x^{i}" if i > 0 else "Hằng số" for i in range(n, -1, -1)]
    row_a = [f"{a_by_deg[i]:.8g}" for i in range(n, -1, -1)]
    row_add = ["-"] + [f"{added_terms[i]:.8g}" for i in range(n - 1, -1, -1)]
    row_b = [f"{b[i]:.8g}" for i in range(n, -1, -1)]
    
    df_horizontal = pd.DataFrame({
        "Thành phần": ["Hệ số ban đầu a_i", f"Nhân thêm (c={c:g})", "Kết quả b_i"],
        **{deg_labels[idx]: [row_a[idx], row_add[idx], row_b[idx]] for idx in range(n + 1)}
    })

    p_c = b[0]
    q_coeffs_desc = [b[i] for i in range(n, 0, -1)]
    return df_horizontal, df_vertical, p_c, q_coeffs_desc, b

def repeated_horner_derivatives(a_coeffs_desc, c: float):
    """
    Lược đồ Horner lặp để tính tất cả các cấp đạo hàm P^(k)(c) (k = 0, 1, ..., n).
    Tại mỗi lượt k, synthetic_division cho ra số dư b_0.
    Theo định lý Taylor: b_0 = P^(k)(c) / k! => P^(k)(c) = b_0 * k!
    """
    n = len(a_coeffs_desc) - 1
    curr_a = a_coeffs_desc.copy()
    deriv_results = []
    horner_pyramid_rows = []

    for k in range(n + 1):
        m = len(curr_a) - 1
        b = [0.0] * (m + 1)
        b[0] = curr_a[0] # curr_a[0] là hệ số bậc cao nhất
        for j in range(1, m + 1):
            b[j] = curr_a[j] + c * b[j - 1]
        
        rem = b[-1] # số dư của lượt này
        fact_k = math.factorial(k)
        deriv_val = rem * fact_k
        
        symbol = f"P({c:g})" if k == 0 else (f"P'({c:g})" if k == 1 else (f"P''({c:g})" if k == 2 else (f"P'''({c:g})" if k == 3 else f"P^({k})({c:g})")))
        deriv_results.append({
            "Cấp (k)": k,
            "Ký hiệu": symbol,
            "Số dư Horner b_0^(k)": rem,
            "k!": fact_k,
            f"Đạo hàm P^(k)({c:g})": deriv_val
        })
        
        horner_pyramid_rows.append(b.copy())
        curr_a = b[:-1] # lấy thương làm đa thức cho lượt tiếp theo

    df_derivs = pd.DataFrame(deriv_results)
    return df_derivs, horner_pyramid_rows

def read_input_txt(filepath: str):
    """
    Đọc file input.txt linh hoạt:
    Hỗ trợ cả dấu chấm thập phân, dấu phẩy danh sách hoặc phân cách phẩy.
    """
    if not os.path.exists(filepath):
        return None, None, "desc", None
    coeffs, c, derivative_order = None, None, None
    order = "desc" # Mặc định bậc giảm dần a_n ... a_0

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip().lower()
                v = v.strip()
                if k in ["c", "x", "eval_x", "point"]:
                    c = float(v.replace(",", "."))
                elif k in ["k", "derivative_order", "derivative", "cap_dao_ham"]:
                    derivative_order = int(v)
                elif k in ["order", "sort", "chieu"]:
                    order = v.lower()
                elif k in ["coeffs", "a", "p", "he_so", "poly"]:
                    clean_v = v.replace(";", " ").replace(",", " ")
                    coeffs = [float(x) for x in clean_v.split() if x]
            else:
                parts = line.replace(";", " ").replace(",", " ").split()
                if len(parts) == 1 and c is None:
                    try:
                        c = float(parts[0])
                    except ValueError:
                        pass
                elif len(parts) > 1 and coeffs is None:
                    try:
                        coeffs = [float(x) for x in parts]
                    except ValueError:
                        pass

    return coeffs, c, order, derivative_order

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description="Lược đồ Horner (Synthetic Division & Repeated Horner for Derivatives)")
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    p.add_argument("-c", "--point", type=float, help="Điểm cần tính c")
    p.add_argument("-k", "--derivative_order", type=int, help="Cấp đạo hàm cần lấy P^(k)(c), 0 <= k <= bậc đa thức")
    p.add_argument("--coeffs", nargs="+", type=float, help="Dãy hệ số của đa thức")
    p.add_argument("--order", choices=["asc", "desc"], default=None, help="Thứ tự hệ số (desc: a_n..a_0, asc: a_0..a_n)")
    args = p.parse_args()

    file_coeffs, file_c, file_order, file_derivative_order = read_input_txt(args.input_file)
    coeffs = args.coeffs if args.coeffs is not None else file_coeffs
    c = args.point if args.point is not None else file_c
    derivative_order = args.derivative_order if args.derivative_order is not None else file_derivative_order
    order = args.order or file_order or "desc"

    if coeffs is None or c is None:
        coeffs = [-0.011833333, 0.62711667, -11.426082, 87.839112, -235.13741]
        c = 9.68

    if "asc" in order or "tang" in order:
        coeffs_desc = list(reversed(coeffs))
    else:
        coeffs_desc = list(coeffs)

    n = len(coeffs_desc) - 1
    if derivative_order is not None and not 0 <= derivative_order <= n:
        p.error(f"Cấp đạo hàm k phải nằm trong khoảng 0..{n}.")

    poly_display = poly_str_from_coeffs(coeffs_desc, var="x")
    coeffs_asc = list(reversed(coeffs_desc))

    # Thực hiện chia Horner 1 lần
    df_horiz, df_vert, p_c, q_coeffs_desc, b_full = synthetic_division(coeffs_desc, c)
    q_poly_display = poly_str_from_coeffs(q_coeffs_desc, var="x")

    # Thực hiện chia Horner lặp tính tất cả đạo hàm
    df_derivs, pyramid = repeated_horner_derivatives(coeffs_desc, c)

    lines = []
    lines.append("=" * 90)
    lines.append("LƯỢC ĐỒ HORNER (SYNTHETIC DIVISION & TÍNH ĐẠO HÀM LẶP)")
    lines.append("=" * 90)
    lines.append(f"1. THÔNG TIN ĐA THỨC VÀ ĐIỂM TÍNH:")
    lines.append(f"   - Đa thức ban đầu: P(x) = {poly_display}")
    lines.append(f"   - Bậc của đa thức: n = {n}")
    lines.append(f"   - Điểm đánh giá:   c = {c:g}")
    if derivative_order is not None:
        lines.append(f"   - Cấp đạo hàm yêu cầu: k = {derivative_order}")
    lines.append(f"   - Dãy hệ số lũy thừa giảm dần [a_n, ..., a_0]: [{', '.join([f'{v:.8g}' for v in coeffs_desc])}]")
    lines.append(f"   - Dãy hệ số lũy thừa tăng dần [a_0, ..., a_n]: [{', '.join([f'{v:.8g}' for v in coeffs_asc])}]")

    lines.append("\n" + "-" * 90)
    lines.append(f"2. BẢNG CHIA HORNER TRUYỀN THỐNG (DẠNG NGANG):")
    lines.append("-" * 90)
    lines.append(df_horiz.to_string(index=False))

    lines.append("\n" + "-" * 90)
    lines.append(f"3. BẢNG CHI TIẾT TỪNG BƯỚC (DẠNG DỌC):")
    lines.append("-" * 90)
    lines.append(df_vert.to_string(index=False))

    lines.append("\n" + "-" * 90)
    lines.append(f"4. KẾT QUẢ PHÉP CHIA VÀ GIÁ TRỊ TẠI c = {c:g}:")
    lines.append("-" * 90)
    lines.append(f"   a) Giá trị đa thức tại c: P({c:g}) = {p_c:.10f}")
    lines.append(f"   b) Đa thức thương:        q(x) = {q_poly_display}")
    lines.append(f"      (Hệ số q(x) giảm dần: [{', '.join([f'{v:.8g}' for v in q_coeffs_desc])}])")
    lines.append(f"   c) Đẳng thức chia:        P(x) = (x - {c:g}) * q(x) + {p_c:.10f}")

    lines.append("\n" + "-" * 90)
    lines.append(f"5. BẢNG HORNER LẶP - TÍNH TẤT CẢ CÁC CẤP ĐẠO HÀM TẠI c = {c:g}:")
    lines.append("-" * 90)
    lines.append(df_derivs.to_string(index=False))

    if derivative_order is not None:
        derivative_row = df_derivs.iloc[derivative_order]
        derivative_value = derivative_row[f"Đạo hàm P^(k)({c:g})"]
        lines.append(
            f"\n   >>> ĐẠO HÀM ĐƯỢC YÊU CẦU: "
            f"P^({derivative_order})({c:g}) = {derivative_value:.10f} <<<"
        )

    lines.append("\nChi tiết các cấp đạo hàm trọng tâm:")
    prime_map = {1: "'", 2: "''", 3: "'''"}
    for k in range(1, min(n + 1, 4)):
        row_k = df_derivs.iloc[k]
        val_k = row_k[f"Đạo hàm P^(k)({c:g})"]
        b0_k = row_k["Số dư Horner b_0^(k)"]
        p_sym = prime_map.get(k, f"^({k})")
        lines.append(f"   * Đạo hàm cấp {k}: P{p_sym}({c:g}) = {k}! * {b0_k:.10g} = {val_k:.10f}")

    # 6. Khai triển Taylor quanh c
    taylor_terms = []
    for k in range(n + 1):
        b0_k = df_derivs.iloc[k]["Số dư Horner b_0^(k)"]
        if abs(b0_k) < 1e-12:
            continue
        sign = "+ " if b0_k >= 0 and len(taylor_terms) > 0 else ("- " if b0_k < 0 and len(taylor_terms) > 0 else "")
        val = abs(b0_k) if len(taylor_terms) > 0 else b0_k
        if k == 0:
            term = f"{val:.8g}"
        elif k == 1:
            term = f"{val:.8g}*(x - {c:g})"
        else:
            term = f"{val:.8g}*(x - {c:g})^{k}"
        taylor_terms.append(f"{sign}{term}")
    taylor_str = " ".join(taylor_terms) if taylor_terms else "0"

    lines.append("\n" + "-" * 90)
    lines.append(f"6. KHAI TRIỂN TAYLOR CỦA P(x) TẠI LÂN CẬN c = {c:g}:")
    lines.append("-" * 90)
    lines.append(f"P(x) = {taylor_str}")

    content = format_report_output("\n".join(lines), args.input_file)
    print("\n" + content)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
