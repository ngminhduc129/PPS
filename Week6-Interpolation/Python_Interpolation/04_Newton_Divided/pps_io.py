"""
pps_io.py - Universal Smart Input/Output Parser for Numerical Analysis (PPS)
Hỗ trợ tự động:
1. Copy-paste nguyên xi file đề thi Midterm (.txt, .csv, .xlsx) của HUST (20241, 20251,...).
2. Tự động nhận diện và trích xuất điểm cần tính x* (eval_x) và y* (eval_y).
3. Hỗ trợ nhập tự do: hàng ngang x/y, rút gọn x0+h+y, hoặc 2 cột dọc.
4. Auto-Windowing: Tự động trích xuất cửa sổ điểm lân cận quanh x* cho các bài toán dữ liệu lớn (100+ mốc).
5. Làm tròn số hiển thị mặc định 7 chữ số; tùy chỉnh bằng round_digits = N trong input.txt.
"""
import os, sys, re
import numpy as np
import pandas as pd

_MINUS_TRANSLATION = str.maketrans({
    "\u2212": "-",  # MINUS SIGN
    "\u2010": "-",  # HYPHEN
    "\u2011": "-",  # NON-BREAKING HYPHEN
    "\u2012": "-",  # FIGURE DASH
    "\u2013": "-",  # EN DASH
    "\u2014": "-",  # EM DASH
    "\ufe63": "-",  # SMALL HYPHEN-MINUS
    "\uff0d": "-",  # FULLWIDTH HYPHEN-MINUS
})

def _parse_float(value):
    text = str(value).strip().translate(_MINUS_TRANSLATION)
    text = re.sub(r"^([+-])\s+(?=\d|\.)", r"\1", text)
    return float(text.replace(",", "."))

def format_report_output(report, input_file, default_digits=7):
    """Round displayed decimal values using round_digits from the input file."""
    digits = default_digits
    try:
        with open(input_file, "r", encoding="utf-8", errors="ignore") as source:
            for line in source:
                match = re.match(r"\s*round_digits\s*[:=]\s*(\d+)", line, re.IGNORECASE)
                if match:
                    digits = int(match.group(1))
                    break
    except FileNotFoundError:
        pass

    if not 0 <= digits <= 15:
        raise ValueError("round_digits phải nằm trong khoảng 0..15.")

    number_pattern = re.compile(r"(?<![\w.])[+-]?(?:\d+\.\d+|\.\d+|\d+)(?:[eE][+-]?\d+)?(?![\w])")

    def round_number(match):
        value = match.group(0)
        line_prefix = report[:match.start()].rsplit("\n", 1)[-1]
        line = line_prefix + report[match.start():].split("\n", 1)[0]
        if not line_prefix.strip() and re.match(rf"\s*{re.escape(value)}\.\s", line):
            return value
        numeric_value = float(value)
        if "e" in value.lower() or (numeric_value != 0 and round(numeric_value, digits) == 0):
            return f"{numeric_value:.{digits}e}"
        if "." not in value:
            return value
        return f"{numeric_value:.{digits}f}"

    return number_pattern.sub(round_number, report)

def _resolve_filepath(path_str):
    if not isinstance(path_str, str):
        return None
    if os.path.exists(path_str) and os.path.isfile(path_str):
        return os.path.abspath(path_str)
    
    cur = os.getcwd()
    for prefix in ["", "../../", "../", "../../../", "e:/2026.1/PPS"]:
        candidate = os.path.normpath(os.path.join(cur, prefix, path_str))
        if os.path.exists(candidate) and os.path.isfile(candidate):
            return candidate
    return None

def parse_universal_input(filepath_or_text):
    """
    Đọc dữ liệu từ đường dẫn file hoặc chuỗi text trực tiếp.
    Trả về:
        points: list các tuple (x_i, y_i)
        params: dict chứa các tham số (eval_x, eval_y, degree, n_points, type, x0_index,...)
    """
    params = {}
    lines = []
    
    # 1. Kiểm tra xem là đường dẫn file hay là text
    resolved_file = _resolve_filepath(filepath_or_text) if isinstance(filepath_or_text, str) else None
    if resolved_file:
        ext = os.path.splitext(resolved_file)[1].lower()
        if ext in [".xlsx", ".xls"]:
            try:
                df = pd.read_excel(resolved_file)
                return _parse_from_dataframe(df, params)
            except Exception as e:
                print(f"[Cảnh báo] Không thể đọc Excel qua pandas: {e}")
        
        with open(resolved_file, "r", encoding="utf-8", errors="ignore") as f:
            lines = [l.strip() for l in f.readlines()]
    else:
        if isinstance(filepath_or_text, str):
            # Nếu người dùng truyền một tên file nhưng không tìm thấy
            if any(filepath_or_text.endswith(ext) for ext in [".csv", ".txt", ".xlsx", ".xls"]):
                print(f"[Cảnh báo] Không tìm thấy file: '{filepath_or_text}'! Kiểm tra lại đường dẫn.")
            lines = [l.strip() for l in filepath_or_text.splitlines()]
        else:
            lines = list(filepath_or_text)

    # 2. Phân tích các dòng text
    explicit_x, explicit_y = None, None
    raw_rows = []
    
    for line in lines:
        if not line:
            continue
        clean = line.strip().translate(_MINUS_TRANSLATION)
        clean = re.sub(r"([+-])\s+(?=\d|\.)", r"\1", clean)
        
        # Bắt các cấu hình tham số dạng key = value hoặc key: value
        if ("=" in clean or ":" in clean) and not clean.startswith("#"):
            # Kiểm tra xem có phải dòng gán tham số không
            delim = "=" if "=" in clean else ":"
            parts = clean.split(delim, 1)
            k = parts[0].strip().lower()
            v = parts[1].strip()
            
            if k in ["eval_x", "x_eval", "x*", "x_bar", "xbar", "x_calc"]:
                try: params["eval_x"] = _parse_float(v)
                except: pass
                continue
            elif k in ["eval_y", "y_eval", "y*", "y_bar", "ybar", "y_calc"]:
                try: params["eval_y"] = _parse_float(v)
                except: pass
                continue
            elif k in ["n_points", "num_points", "so_moc", "points_count"]:
                if v.lower() == "all":
                    params["n_points"] = "all"
                else:
                    try: params["n_points"] = int(v)
                    except: pass
                continue
            elif k in ["degree", "order", "bac", "cap"]:
                try: params["degree"] = int(v)
                except: pass
                continue
            elif k == "type":
                params["type"] = "backward" if ("back" in v.lower() or "lui" in v.lower()) else "forward"
                continue
            elif k in ["x0_index", "x0_idx"]:
                try: params["x0_index"] = int(v)
                except: pass
                continue
            elif k == "x0":
                try: params["x0"] = _parse_float(v)
                except: pass
                continue
            elif k == "h":
                try: params["h"] = _parse_float(v)
                except: pass
                continue
            elif k in ["x", "x_i", "xs", "x_values"]:
                toks = [_parse_float(t) for t in re.split(r"[\s,;\t]+", v) if t]
                if len(toks) == 1 and "eval_x" not in params:
                    params["eval_x"] = toks[0]
                elif len(toks) > 1:
                    explicit_x = toks
                continue
            elif k in ["y", "y_i", "ys", "y_values"]:
                toks = [_parse_float(t) for t in re.split(r"[\s,;\t]+", v) if t]
                if len(toks) == 1 and "eval_y" not in params:
                    params["eval_y"] = toks[0]
                elif len(toks) > 1:
                    explicit_y = toks
                continue

        # Kiểm tra comment
        if clean.startswith("#"):
            # Kiểm tra xem comment có chứa tham số không (vd: # n_points: 5)
            comm_body = clean.lstrip("#").strip()
            if ("=" in comm_body or ":" in comm_body):
                delim = "=" if "=" in comm_body else ":"
                k, v = [p.strip() for p in comm_body.split(delim, 1)]
                k_low = k.lower()
                if k_low in ["n_points", "num_points", "so_moc"]:
                    params["n_points"] = "all" if v.lower() == "all" else int(v)
                elif k_low in ["degree", "order", "bac"]:
                    params["degree"] = int(v)
                elif k_low in ["eval_x", "x*", "x_bar"]:
                    params["eval_x"] = _parse_float(v)
                elif k_low in ["eval_y", "y*", "y_bar"]:
                    params["eval_y"] = _parse_float(v)
                elif k_low == "type":
                    params["type"] = "backward" if ("back" in v.lower() or "lui" in v.lower()) else "forward"
            continue
            
        raw_rows.append(clean)

    # 3. Phân tích định dạng bảng số liệu (Table Parsing)
    if not raw_rows:
        if explicit_x is not None and explicit_y is not None:
            points = [(explicit_x[i], explicit_y[i]) for i in range(min(len(explicit_x), len(explicit_y)))]
            return points, params
        if "x0" in params and "h" in params and explicit_y is not None:
            points = [(params["x0"] + i * params["h"], explicit_y[i]) for i in range(len(explicit_y))]
            return points, params
        return [], params

    # Kiểm tra dòng đầu tiên xem có phải Header của đề Midterm hay không
    header_row = raw_rows[0]
    # Phân tách tokens theo dấu phẩy, tab, chấm phẩy hoặc khoảng trắng
    header_toks = [t.strip() for t in re.split(r"[\t,;]+", header_row) if t.strip()]
    
    # TRƯỜNG HỢP A: Đề thi Midterm 20241 dạng "x_i, y_i, x, 4, y, 5" hoặc tab
    # Ví dụ: ['x_i', 'y_i', 'x', '4', 'y', '5']
    if any(h.lower() in ["x_i", "xi"] for h in header_toks):
        # Duyệt trong header_toks tìm x và y truy vấn
        for idx, token in enumerate(header_toks):
            t_low = token.lower()
            if t_low == "x" and idx + 1 < len(header_toks):
                try:
                    val = _parse_float(header_toks[idx + 1])
                    if "eval_x" not in params: params["eval_x"] = val
                except ValueError: pass
            elif t_low == "y" and idx + 1 < len(header_toks):
                try:
                    val = _parse_float(header_toks[idx + 1])
                    if "eval_y" not in params: params["eval_y"] = val
                except ValueError: pass
        
        # TRƯỜNG HỢP B: Đề thi Midterm 20251 dạng header có "x_bar", "y_bar"
        x_bar_col = None
        y_bar_col = None
        for c_idx, h in enumerate(header_toks):
            if "x_bar" in h.lower() or "xbar" in h.lower():
                x_bar_col = c_idx
            elif "y_bar" in h.lower() or "ybar" in h.lower():
                y_bar_col = c_idx
                
        # Dòng 2 có thể chứa giá trị x_bar và y_bar
        start_row_idx = 1
        if (x_bar_col is not None or y_bar_col is not None) and len(raw_rows) > 1:
            row2_toks = [t.strip() for t in re.split(r"[\t,;]+", raw_rows[1]) if t.strip()]
            # Nếu dòng 2 có giá trị ở các cột x_bar, y_bar
            if len(row2_toks) >= 3:
                # Kiểm tra 2 số đầu tiên là x_i, y_i
                try:
                    if x_bar_col is not None and "eval_x" not in params and len(row2_toks) >= 3:
                        params["eval_x"] = _parse_float(row2_toks[2])
                    if y_bar_col is not None and "eval_y" not in params and len(row2_toks) >= 4:
                        params["eval_y"] = _parse_float(row2_toks[3])
                except ValueError: pass
        
        # Lấy dữ liệu 2 cột x_i, y_i từ start_row_idx trở đi
        points = []
        for r in raw_rows[start_row_idx:]:
            # Phân tách theo dấu tab, phẩy, chấm phẩy hoặc khoảng trắng
            cols = [t.strip() for t in re.split(r"[\t,;]+", r) if t.strip()]
            if len(cols) >= 2:
                try:
                    xi = _parse_float(cols[0])
                    yi = _parse_float(cols[1])
                    points.append((xi, yi))
                except ValueError:
                    pass
        return points, params

    # TRƯỜNG HỢP C & E: Kiểm tra ma trận số trong raw_rows
    def _split_line_cols(line_str):
        if "\t" in line_str:
            return [t.strip() for t in line_str.split("\t") if t.strip()]
        elif ";" in line_str:
            return [t.strip() for t in line_str.split(";") if t.strip()]
        elif "," in line_str and "." in line_str:
            return [t.strip() for t in line_str.split(",") if t.strip()]
        else:
            return [t.strip() for t in re.split(r"\s+", line_str) if t.strip()]

    lines_of_nums = []
    for r in raw_rows:
        toks = _split_line_cols(r)
        try:
            nums = [_parse_float(t) for t in toks]
            if nums: lines_of_nums.append(nums)
        except ValueError:
            pass

    table_points = []
    if len(lines_of_nums) == 2 and len(lines_of_nums[0]) > 2 and len(lines_of_nums[0]) == len(lines_of_nums[1]):
        table_points = [(lines_of_nums[0][i], lines_of_nums[1][i]) for i in range(len(lines_of_nums[0]))]
    elif len(lines_of_nums) >= 2 and all(len(r) == 2 for r in lines_of_nums):
        table_points = [(r[0], r[1]) for r in lines_of_nums]
    elif len(lines_of_nums) >= 2 and len(lines_of_nums[0]) >= 2:
        for r in lines_of_nums:
            if len(r) >= 2: table_points.append((r[0], r[1]))

    # Nếu có bảng 2 cột số nhiều điểm (ví dụ dán đề thi phía dưới), ưu tiên lấy bảng số đó
    if len(table_points) > len(explicit_x or []):
        return table_points, params

    # Nếu không có bảng dọc hoặc explicit x, y có số điểm lớn hơn
    if explicit_x is not None and explicit_y is not None:
        points = [(explicit_x[i], explicit_y[i]) for i in range(min(len(explicit_x), len(explicit_y)))]
        return points, params

    # Nhập rút gọn x0 + h + y
    if "x0" in params and "h" in params and explicit_y is not None:
        points = [(params["x0"] + i * params["h"], explicit_y[i]) for i in range(len(explicit_y))]
        return points, params

    if table_points:
        return table_points, params

    return [], params

def _parse_from_dataframe(df, params):
    """Trích xuất mốc nội suy và x_bar/y_bar từ pandas DataFrame (Excel/CSV)"""
    col_names = [str(c).strip() for c in df.columns]
    
    # Tìm eval_x, eval_y trong tên cột (vd 20241: ['x_i', 'y_i', 'x', 4, 'y', 5])
    for idx, c in enumerate(col_names):
        c_low = c.lower()
        if c_low == "x" and idx + 1 < len(col_names):
            try: params["eval_x"] = _parse_float(col_names[idx + 1])
            except: pass
        elif c_low == "y" and idx + 1 < len(col_names):
            try: params["eval_y"] = _parse_float(col_names[idx + 1])
            except: pass
        elif "x_bar" in c_low or "xbar" in c_low:
            val = df[df.columns[idx]].dropna().iloc[0] if not df[df.columns[idx]].dropna().empty else None
            if val is not None:
                try: params["eval_x"] = _parse_float(val)
                except: pass
        elif "y_bar" in c_low or "ybar" in c_low:
            val = df[df.columns[idx]].dropna().iloc[0] if not df[df.columns[idx]].dropna().empty else None
            if val is not None:
                try: params["eval_y"] = _parse_float(val)
                except: pass

    # Lấy 2 cột đầu làm x_i, y_i
    points = []
    x_col = df.columns[0]
    y_col = df.columns[1]
    sub = df[[x_col, y_col]].dropna()
    for _, row in sub.iterrows():
        try:
            points.append((_parse_float(row[0]), _parse_float(row[1])))
        except:
            pass
    return points, params

def select_window_for_method(points, eval_x=None, method="gauss1", n_points=None, degree=None, eval_y=None):
    """
    Tự động chọn mốc trung tâm x0 và cắt cửa sổ điểm tối ưu quanh x* (hoặc y*) cho phương pháp tương ứng.
    
    Args:
        points: list (x_i, y_i)
        eval_x: float điểm cần nội suy
        method: 'gauss1', 'gauss2', 'stirling', 'bessel', 'newton_finite', 'newton_divided', 'lagrange',
                'inverse_function', 'loop_inverse', 'spline'
        n_points: số điểm mong muốn (int hoặc 'all')
        degree: bậc mong muốn (n_points = degree + 1)
        eval_y: float điểm cần nội suy ngược (nếu có)
        
    Returns:
        subset_points: danh sách điểm được chọn
        x0_idx_in_subset: chỉ số của mốc x0 trong subset (quan trọng cho Gauss/Stirling/Bessel/Lặp)
        msg: thông báo tóm tắt quá trình chọn cửa sổ
    """
    total_pts = len(points)
    if total_pts == 0:
        return [], 0, "Không có dữ liệu điểm."

    # Nếu người dùng yêu cầu dùng toàn bộ
    if n_points == "all":
        # Tìm x0 gần x* nhất trong toàn bộ
        if eval_x is not None:
            diffs = [abs(p[0] - eval_x) for p in points]
            best_idx = int(np.argmin(diffs))
        elif eval_y is not None:
            diffs = [abs(p[1] - eval_y) for p in points]
            best_idx = int(np.argmin(diffs))
        else:
            best_idx = total_pts // 2
        return points, best_idx, f"Sử dụng toàn bộ {total_pts} mốc điểm (x0_idx = {best_idx})."

    # Nếu số điểm sẵn có nhỏ (<= 11 điểm), không cần cắt trừ khi người dùng chỉ định rõ
    if total_pts <= 11 and n_points is None and degree is None:
        if eval_x is not None:
            diffs = [abs(p[0] - eval_x) for p in points]
            best_idx = int(np.argmin(diffs))
        elif eval_y is not None:
            diffs = [abs(p[1] - eval_y) for p in points]
            best_idx = int(np.argmin(diffs))
        else:
            best_idx = total_pts // 2
        return points, best_idx, f"Số mốc nhỏ ({total_pts} điểm), giữ nguyên toàn bộ."

    # Xác định số mốc mục tiêu
    method_low = method.lower()
    
    if degree is not None:
        target_pts = degree + 1
    elif n_points is not None:
        target_pts = int(n_points)
    else:
        if method_low in ["gauss1", "gauss2", "stirling"]:
            target_pts = 5  # Bậc 4
        elif method_low == "bessel":
            target_pts = 6  # 6 mốc (mỗi bên 3 mốc)
        elif method_low in ["loop_inverse", "inverse_function"]:
            target_pts = 4  # Chuẩn đề thi 4 mốc bậc 3
        elif method_low == "spline":
            target_pts = min(total_pts, 8)
        else:
            target_pts = 5  # Bậc 4

    # Đảm bảo điều kiện chẵn/lẻ của từng phương pháp
    if method_low in ["gauss1", "gauss2", "stirling"]:
        if target_pts % 2 == 0:
            target_pts += 1  # Tăng lên số lẻ gần nhất
    elif method_low == "bessel":
        if target_pts % 2 != 0:
            target_pts += 1  # Tăng lên số chẵn gần nhất

    target_pts = min(target_pts, total_pts)

    # Xử lý trường hợp phương pháp ngược (Inverse / Loop) dựa trên eval_y
    if method_low in ["inverse_function", "loop_inverse"] and eval_y is not None:
        bracket_idx = None
        for i in range(total_pts - 1):
            if (points[i][1] - eval_y) * (points[i + 1][1] - eval_y) <= 0:
                bracket_idx = i
                break
        if bracket_idx is not None:
            start = max(0, min(bracket_idx, total_pts - target_pts))
            k_y = bracket_idx
        else:
            diffs_y = [abs(p[1] - eval_y) for p in points]
            k_y = int(np.argmin(diffs_y))
            start = max(0, min(k_y - target_pts // 2, total_pts - target_pts))

        end = start + target_pts
        sub = points[start:end]
        x0_in_sub = k_y - start
        msg = (f"[Auto-Windowing] Đã phát hiện y* = {eval_y}. Chọn khoảng kẹp [{points[k_y][1]:.6f}, {points[min(k_y+1, total_pts-1)][1]:.6f}] "
               f"(x = {points[k_y][0]}..{points[min(k_y+1, total_pts-1)][0]}), "
               f"trích xuất cửa sổ {len(sub)} mốc lân cận (từ x = {sub[0][0]} đến {sub[-1][0]}).")
        return sub, x0_in_sub, msg

    if eval_x is None:
        if eval_y is not None:
            diffs_y = [abs(p[1] - eval_y) for p in points]
            k_y = int(np.argmin(diffs_y))
            start = max(0, min(k_y - target_pts // 2, total_pts - target_pts))
            end = start + target_pts
            sub = points[start:end]
            return sub, k_y - start, f"Cắt {len(sub)} mốc quanh y* = {eval_y}."
        start = max(0, (total_pts - target_pts) // 2)
        end = start + target_pts
        sub = points[start:end]
        return sub, len(sub) // 2, f"Không có eval_x/eval_y, lấy {target_pts} mốc chính giữa."

    # Tìm mốc gần eval_x nhất
    diffs = [abs(p[0] - eval_x) for p in points]
    k = int(np.argmin(diffs))
    
    # 1. Với GAUSS I, GAUSS II, STIRLING (cần số mốc lẻ, đối xứng quanh mốc k)
    if method_low in ["gauss1", "gauss2", "stirling"]:
        m = (target_pts - 1) // 2  # Số mốc mỗi bên
        start = k - m
        end = k + m + 1
        
        # Xử lý biên nếu chạm đầu hoặc cuối mốc
        if start < 0:
            start = 0
            end = min(total_pts, target_pts)
        elif end > total_pts:
            end = total_pts
            start = max(0, total_pts - target_pts)
            
        sub = points[start:end]
        x0_in_sub = k - start
        msg = (f"[Auto-Windowing] Phát hiện dữ liệu lớn ({total_pts} mốc). "
               f"Đã chọn mốc trung tâm x0 = {points[k][0]} (gần x* = {eval_x} nhất), "
               f"trích xuất cửa sổ {len(sub)} mốc lân cận đối xứng (từ x = {sub[0][0]} đến {sub[-1][0]}).")
        return sub, x0_in_sub, msg

    # 2. Với BESSEL (cần số mốc chẵn, đối xứng quanh khoảng [x_k, x_{k+1}] chứa x*)
    elif method_low == "bessel":
        # Xác định đoạn [x_k, x_{k+1}]
        if points[k][0] > eval_x and k > 0:
            k = k - 1
        # Giờ x_k <= eval_x (hoặc k là mốc dưới)
        m = target_pts // 2
        start = k - m + 1
        end = k + m + 1
        
        if start < 0:
            start = 0
            end = min(total_pts, target_pts)
        elif end > total_pts:
            end = total_pts
            start = max(0, total_pts - target_pts)
            
        sub = points[start:end]
        x0_in_sub = (len(sub) // 2) - 1
        msg = (f"[Auto-Windowing] Phát hiện dữ liệu lớn ({total_pts} mốc). "
               f"Bessel: Chọn khoảng trung tâm [{sub[x0_in_sub][0]}, {sub[x0_in_sub+1][0]}], "
               f"trích xuất cửa sổ {len(sub)} mốc đối xứng.")
        return sub, x0_in_sub, msg

    # 3. Với NEWTON TIẾN / LÙI (Finite)
    elif method_low == "newton_finite":
        # Tìm mốc k bên trái ngay trước eval_x (points[k][0] <= eval_x)
        # để đảm bảo mốc gốc x0 <= eval_x (t >= 0 cho Newton tiến, tránh ngoại suy ở đầu)
        left_indices = [i for i, p in enumerate(points) if p[0] <= eval_x]
        if left_indices:
            k_left = left_indices[-1]
        else:
            k_left = 0

        start = max(0, min(k_left, total_pts - target_pts))
        end = start + target_pts
        sub = points[start:end]
        x0_in_sub = 0
        msg = (f"[Auto-Windowing] Trích xuất cửa sổ {len(sub)} mốc lân cận cho Newton Finite "
               f"(từ x = {sub[0][0]} đến {sub[-1][0]}, mốc gốc x0 = {sub[0][0]} <= x* = {eval_x}).")
        return sub, x0_in_sub, msg

    # 4. Với LAGRANGE / NEWTON DIVIDED
    else:
        # Lấy target_pts điểm gần eval_x nhất
        indexed_diffs = sorted(enumerate(diffs), key=lambda x: x[1])
        selected_indices = sorted([idx for idx, _ in indexed_diffs[:target_pts]])
        sub = [points[i] for i in selected_indices]
        msg = (f"[Auto-Windowing] Trích xuất {len(sub)} mốc gần x* = {eval_x} nhất.")
        return sub, 0, msg
