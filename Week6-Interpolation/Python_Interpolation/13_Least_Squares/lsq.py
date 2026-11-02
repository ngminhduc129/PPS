"""
lsq.py - Phương pháp Bình Phương Bé Nhất (Least Squares Fitting)
Gộp toàn bộ nội dung từ pp6_LSQ.ipynb và mở rộng cho bài thi Giữa kỳ (Midterm).

Bao gồm:
1. Universal File/Data Parser: Đọc .csv, .txt, .xlsx, hỗ trợ phân cách dấu phẩy, chấm phẩy, tab, dấu cách, và số thực dùng dấu phẩy ','.
2. Thuật toán Bình phương bé nhất tổng quát (Normal Equations):
   - Ma trận thiết kế Phi (n x m)
   - Ma trận Gram M = Phi^T * Phi (m x m)
   - Véctơ b' = Phi^T * y (m x 1)
   - Giải hệ Ma = b'
   - Sai số SSE, MSE, RMSE, MAE, MaxError, MAPE, R^2
3. Các mô hình hồi quy / hàm thực nghiệm:
   - Linear Basis: y = a*x + b*sin(x) (Notebook Cell 20)
   - Trigonometric / Fourier: y = a0 + a1*cos(x) + a2*cos(2x) + b1*sin(x) + b2*sin(2x) (Notebook Cell 7)
   - Exponential: y = a * exp(b*x) (Notebook Cell 25, hỗ trợ cả y > 0 và y < 0)
   - Power Law: y = a * x^b (Notebook Cell 31)
   - Polynomial: y = a0 + a1*x + ... + am*x^m
   - Hyperbolic / Rational: y = a + b/x, y = 1/(a + bx), y = x/(a + bx)
   - Logarithmic: y = a + b*ln(x)
   - Midterm K70 Model: y = a * exp(b*x^2 + c/x) (20251GKtest1)
4. Chế độ so sánh tự động (Compare All Models) để tìm hàm khớp tốt nhất theo R^2 và RMSE.
5. Dự báo giá trị y_hat(x*) khi có mốc x*.
"""

import os
import sys
import re
import argparse
import numpy as np
import pandas as pd
from typing import List, Tuple, Callable, Dict, Any, Optional

# Use the shared interpolation parser when this module is run from its folder.
try:
    from pps_io import format_report_output, parse_universal_input
except ImportError:
    interpolation_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..")
    )
    if interpolation_dir not in sys.path:
        sys.path.insert(0, interpolation_dir)
    from pps_io import format_report_output, parse_universal_input

HAS_PPS_IO = True


# ==============================================================================
# 1. PARSER DỮ LIỆU ĐA NĂNG (UNIVERSAL DATA PARSER)
# ==============================================================================
def parse_xy_data(filepath: str, delimiter: Optional[str] = None) -> List[Tuple[float, float]]:
    """
    Đọc file dữ liệu (CSV/TXT), tự động nhận diện dấu phân cách (; , tab, space)
    và hỗ trợ dấu phẩy thập phân ',' theo đúng định dạng trong pp6_LSQ.ipynb.
    """
    data_points = []
    detected_delimiter = delimiter

    if not os.path.exists(filepath):
        print(f"Lỗi: Không tìm thấy file '{filepath}'")
        return []

    # Delimiter Sniffing
    if detected_delimiter is None:
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                first_line = ""
                for line in f:
                    stripped = line.strip().lstrip('\ufeff')
                    if stripped and not stripped.startswith('#'):
                        first_line = stripped
                        break

                if ';' in first_line:
                    detected_delimiter = ';'
                elif '\t' in first_line:
                    detected_delimiter = '\t'
                elif ' ' in first_line:
                    parts = re.split(r'\s+', first_line)
                    if len(parts) >= 2:
                        try:
                            float(parts[0].replace(',', '.'))
                            float(parts[1].replace(',', '.'))
                            detected_delimiter = ' '
                        except ValueError:
                            pass

                if detected_delimiter is None and ',' in first_line:
                    detected_delimiter = ','

                if detected_delimiter is None:
                    detected_delimiter = ','
        except Exception as e:
            print(f"Lỗi khi đọc file để phát hiện dấu phân cách: {e}")
            return []

    # File Parsing
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line_number, raw_line in enumerate(f, 1):
                line = raw_line.strip().lstrip('\ufeff')
                if not line or line.startswith('#'):
                    continue

                if detected_delimiter == ' ':
                    parts = re.split(r'\s+', line)
                else:
                    parts = line.split(detected_delimiter)

                if len(parts) >= 2:
                    x_str, y_str = parts[0].strip(), parts[1].strip()
                    try:
                        x_val = float(x_str.replace(',', '.'))
                        y_val = float(y_str.replace(',', '.'))
                        data_points.append((x_val, y_val))
                    except ValueError:
                        # Bỏ qua dòng tiêu đề như 'x_i,y_i' hoặc ',2.95'
                        continue
    except Exception as e:
        print(f"Lỗi khi phân tích dữ liệu: {e}")
        return []

    return data_points


def load_dataset(input_source: str) -> Tuple[List[Tuple[float, float]], Dict[str, Any]]:
    """
    Đọc dữ liệu từ file hoặc đường dẫn, trả về danh sách điểm và các tham số (eval_x, eval_y).
    """
    params = {}
    points = []

    # Nếu có pps_io, dùng pps_io để trích xuất cả eval_x, eval_y (đề thi Midterm)
    if HAS_PPS_IO:
        try:
            pts, p = parse_universal_input(input_source)
            if pts and len(pts) > 0:
                return pts, p
        except Exception:
            pass

    # Nếu không hoặc pps_io không đọc được, dùng parser dự phòng
    if os.path.exists(input_source):
        points = parse_xy_data(input_source)
    else:
        # Thử đọc trực tiếp chuỗi text
        lines = input_source.strip().splitlines()
        for l in lines:
            parts = re.split(r'[,\t; ]+', l.strip().lstrip('\ufeff'))
            if len(parts) >= 2:
                try:
                    pts = (float(parts[0].replace(',', '.')), float(parts[1].replace(',', '.')))
                    points.append(pts)
                except ValueError:
                    continue

    return points, params


# ==============================================================================
# 2. CHUYỂN ĐỔI BIẾN (DATA TRANSFORMATION)
# ==============================================================================
def modify_value(
    original_points: List[Tuple[float, float]],
    transform: Callable[[float, float], Tuple[float, float]] = lambda x, y: (x, y)
) -> List[Tuple[float, float]]:
    """
    Áp dụng phép biến đổi (X, Y) = transform(x, y) cho từng điểm dữ liệu.
    """
    transformed_points = []
    for x, y in original_points:
        try:
            X, Y = transform(x, y)
            transformed_points.append((float(X), float(Y)))
        except Exception as e:
            raise ValueError(f"Lỗi khi áp dụng transform tại điểm (x={x}, y={y}): {e}")
    return transformed_points


# ==============================================================================
# 3. TÍNH TOÁN SAI SỐ (ERROR METRICS)
# ==============================================================================
def calculate_least_squares_error(
    y_vals: np.ndarray,
    b_prime: np.ndarray,
    a_coeffs: np.ndarray,
    y_pred: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """
    Tính toán toàn diện các chỉ số sai số:
    - SSE (Sum of Squared Errors / S_min): y^T y - a^T b' hoặc sum((y - y_hat)^2)
    - MSE (Mean Squared Error): SSE / n
    - RMSE (Root Mean Square Error): sqrt(MSE)
    - MAE (Mean Absolute Error): mean(|y - y_hat|)
    - MaxError: max(|y - y_hat|)
    - MAPE (Mean Absolute Percentage Error): mean(|(y - y_hat)/y|) * 100%
    - R2 (Coefficient of Determination): 1 - SSE / SS_tot
    """
    n = len(y_vals)
    y_flat = y_vals.flatten()

    # Tính SSE từ công thức đại số tuyến tính: SSE = y^T*y - a^T*b'
    # và đối chiếu với sum((y - y_hat)^2) để đảm bảo ổn định số học
    sse_matrix = float(y_vals.T.dot(y_vals)[0, 0] - a_coeffs.T.dot(b_prime)[0, 0])
    sse_matrix = max(0.0, sse_matrix)

    if y_pred is not None:
        y_pred_flat = y_pred.flatten()
        diff = y_flat - y_pred_flat
        sse = float(np.sum(diff**2))
        mae = float(np.mean(np.abs(diff)))
        max_err = float(np.max(np.abs(diff)))

        # MAPE (tránh chia 0)
        non_zero = np.abs(y_flat) > 1e-12
        if np.any(non_zero):
            mape = float(np.mean(np.abs(diff[non_zero] / y_flat[non_zero])) * 100.0)
        else:
            mape = 0.0
    else:
        sse = sse_matrix
        mae = 0.0
        max_err = 0.0
        mape = 0.0

    mse = sse / n if n > 0 else 0.0
    rmse = float(np.sqrt(mse))

    # R^2
    ss_tot = float(np.sum((y_flat - np.mean(y_flat))**2))
    r2 = 1.0 - (sse / ss_tot) if ss_tot > 1e-14 else 1.0

    return {
        "sse": sse,
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "max_err": max_err,
        "mape": mape,
        "r2": r2,
        "sse_matrix": sse_matrix
    }


# ==============================================================================
# 4. HỆ PHƯƠNG TRÌNH CHUẨN (CORE NORMAL EQUATIONS SOLVER)
# ==============================================================================
def least_squares_fit(
    points: List[Tuple[float, float]],
    basis_funcs: List[Callable[[np.ndarray], np.ndarray]],
    basis_names: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray, Dict[str, float]]:
    """
    Tìm các hệ số tối ưu [a1, a2, ..., am] cho hệ hàm cơ sở [phi_1, ..., phi_m]
    bằng phương pháp Phương Trình Chuẩn (Normal Equations).

    Returns:
        df_phi: Ma trận thiết kế Phi (n x m)
        df_M: Ma trận Gram M = Phi^T * Phi (m x m)
        df_b_prime: Véctơ vế phải b' = Phi^T * y (m x 1)
        df_a: DataFrame hệ số a_j
        a_coeffs: Mảng numpy (m x 1)
        metrics: Dict chứa SSE, MSE, RMSE, R^2,...
    """
    n = len(points)
    m = len(basis_funcs)

    if n < m:
        print(f"Cảnh báo: Số điểm dữ liệu (n = {n}) nhỏ hơn số hàm cơ sở (m = {m}). Hệ thiếu định!")

    x_vals = np.array([p[0] for p in points], dtype=np.float64)
    y_vals = np.array([p[1] for p in points], dtype=np.float64).reshape(-1, 1)

    # 1. Xây dựng ma trận thiết kế Phi (n x m)
    phi_matrix = np.zeros((n, m), dtype=np.float64)
    for j, func in enumerate(basis_funcs):
        phi_matrix[:, j] = func(x_vals)

    col_names = basis_names if basis_names and len(basis_names) == m else [f"\u03C6_{j+1}(x)" for j in range(m)]
    df_phi = pd.DataFrame(phi_matrix, columns=col_names)

    # 2. Ma trận Gram M = Phi^T * Phi (m x m)
    phi_t = phi_matrix.T
    M = phi_t.dot(phi_matrix)
    df_M = pd.DataFrame(M, index=[f"Row_{j+1}" for j in range(m)], columns=[f"Col_{j+1}" for j in range(m)])

    # 3. Véctơ vế phải b' = Phi^T * y (m x 1)
    b_prime = phi_t.dot(y_vals)
    df_b_prime = pd.DataFrame(b_prime, index=[f"Row_{j+1}" for j in range(m)], columns=["b'_j"])

    # 4. Giải hệ phương trình chuẩn Ma = b'
    try:
        a_coeffs = np.linalg.solve(M, b_prime)
    except np.linalg.LinAlgError:
        print("Lỗi: Ma trận M bị suy biến (singular). Các hàm cơ sở phụ thuộc tuyến tính!")
        a_coeffs = np.linalg.lstsq(M, b_prime, rcond=None)[0]

    df_a = pd.DataFrame(a_coeffs, columns=["a_j"], index=[f"a_{j+1}" for j in range(m)])

    # 5. Giá trị dự báo trên không gian tuyến tính
    y_pred_linear = phi_matrix.dot(a_coeffs)
    metrics = calculate_least_squares_error(y_vals, b_prime, a_coeffs, y_pred=y_pred_linear)

    return df_phi, df_M, df_b_prime, df_a, a_coeffs, metrics


# ==============================================================================
# 5. CÁC MÔ HÌNH HỒI QUY CHUYÊN BIỆT (SPECIALIZED MODELS)
# ==============================================================================
class ModelResult:
    def __init__(
        self,
        name: str,
        formula_str: str,
        predict_fn: Callable[[float], float],
        coeffs: Dict[str, float],
        df_phi: pd.DataFrame,
        df_M: pd.DataFrame,
        df_b_prime: pd.DataFrame,
        df_a: pd.DataFrame,
        linear_coeffs: np.ndarray,
        metrics: Dict[str, float],
        transformed_points: Optional[List[Tuple[float, float]]] = None
    ):
        self.name = name
        self.formula_str = formula_str
        self.predict_fn = predict_fn
        self.coeffs = coeffs
        self.df_phi = df_phi
        self.df_M = df_M
        self.df_b_prime = df_b_prime
        self.df_a = df_a
        self.linear_coeffs = linear_coeffs
        self.metrics = metrics
        self.transformed_points = transformed_points

    def print_report(self, eval_x: Optional[float] = None, show_phi_rows: int = 10):
        print("\n" + "=" * 80)
        print(f"MÔ HÌNH: {self.name.upper()}")
        print("=" * 80)
        print(f"-> Phương trình hàm khớp: {self.formula_str}\n")

        print("1. Hệ số tuyến tính hóa giải từ hệ chuẩn Ma = b':")
        print(self.df_a.to_string())

        print("\n2. Ma trận Gram M = \u03A6^T * \u03A6 (kích thước {}x{}):".format(self.df_M.shape[0], self.df_M.shape[1]))
        print(self.df_M.to_string())

        print("\n3. Véctơ vế phải b' = \u03A6^T * y:")
        print(self.df_b_prime.to_string())

        if show_phi_rows > 0:
            total_n = len(self.df_phi)
            print(f"\n4. Ma trận thiết kế \u03A6 ({total_n} hàng x {self.df_phi.shape[1]} cột) [Hiển thị tối đa {show_phi_rows} hàng đầu]:")
            print(self.df_phi.head(show_phi_rows).to_string())

        print("\n5. Các hệ số của mô hình:")
        for k, v in self.coeffs.items():
            print(f"   {k} = {v:.10f}")

        print("\n6. Đánh giá sai số (Error Metrics):")
        print(f"   - Tổng bình phương sai số (SSE / S_min) : {self.metrics['sse']:.10f}")
        print(f"   - Sai số trung bình bình phương (MSE)   : {self.metrics['mse']:.10f}")
        print(f"   - Căn sai số trung bình (RMSE)          : {self.metrics['rmse']:.10f}")
        print(f"   - Sai số tuyệt đối trung bình (MAE)      : {self.metrics['mae']:.10f}")
        print(f"   - Sai số tuyệt đối lớn nhất (Max Error)  : {self.metrics['max_err']:.10f}")
        print(f"   - Sai số tương đối TB (MAPE)            : {self.metrics['mape']:.4f}%")
        print(f"   - Hệ số xác định (R^2)                  : {self.metrics['r2']:.10f}")

        if eval_x is not None:
            pred_y = self.predict_fn(eval_x)
            print(f"\n7. Dự báo tại x* = {eval_x}:")
            print(f"   y_hat({eval_x}) = {pred_y:.10f}")


def fit_linear_combination(
    points: List[Tuple[float, float]],
    basis_funcs: List[Callable[[np.ndarray], np.ndarray]],
    basis_names: List[str],
    formula_template: str,
    model_name: str = "Tuyến tính kết hợp (Linear Basis)"
) -> ModelResult:
    """
    Khớp mô hình tổ hợp tuyến tính dạng y = a1*phi_1(x) + a2*phi_2(x) + ...
    Ví dụ: y = a*x + b*sin(x) trong Cell 20 của notebook.
    """
    df_phi, df_M, df_b_prime, df_a, a_coeffs, metrics = least_squares_fit(points, basis_funcs, basis_names)

    coeffs = {basis_names[j]: float(a_coeffs[j, 0]) for j in range(len(basis_funcs))}

    # Tạo hàm dự đoán
    def predict_fn(x_val: float) -> float:
        val = 0.0
        x_arr = np.array([x_val], dtype=np.float64)
        for j, fn in enumerate(basis_funcs):
            val += a_coeffs[j, 0] * fn(x_arr)[0]
        return float(val)

    # Tính sai số thực tế trên dữ liệu gốc
    x_orig = np.array([p[0] for p in points])
    y_orig = np.array([p[1] for p in points]).reshape(-1, 1)
    y_preds = np.array([predict_fn(x) for x in x_orig]).reshape(-1, 1)
    metrics = calculate_least_squares_error(y_orig, df_b_prime.values, a_coeffs, y_pred=y_preds)

    # Tạo chuỗi công thức
    formula_str = formula_template
    for name, val in coeffs.items():
        formula_str = formula_str.replace(f"{{{name}}}", f"{val:+.6f}")

    return ModelResult(
        name=model_name,
        formula_str=formula_str,
        predict_fn=predict_fn,
        coeffs=coeffs,
        df_phi=df_phi,
        df_M=df_M,
        df_b_prime=df_b_prime,
        df_a=df_a,
        linear_coeffs=a_coeffs,
        metrics=metrics
    )


def fit_notebook_linear_x_sinx(points: List[Tuple[float, float]]) -> ModelResult:
    """
    Khớp mô hình: y = a*x + b*sin(x) theo đúng Cell 20 của pp6_LSQ.ipynb.
    """
    basis_funcs = [
        lambda x: x,
        lambda x: np.sin(x)
    ]
    basis_names = ["x", "sin(x)"]
    return fit_linear_combination(
        points,
        basis_funcs,
        basis_names,
        formula_template="y = {x} * x + {sin(x)} * sin(x)",
        model_name="Mô hình y = a*x + b*sin(x) (Notebook Cell 20)"
    )


def fit_trigonometric(points: List[Tuple[float, float]]) -> ModelResult:
    """
    Khớp mô hình chuỗi lượng giác (Fourier-like) theo Cell 7 của pp6_LSQ.ipynb:
    y = a1 + a2*cos(x) + a3*cos(2x) + a4*sin(x) + a5*sin(2x)
    """
    basis_funcs = [
        lambda x: np.ones_like(x),
        lambda x: np.cos(x),
        lambda x: np.cos(2 * x),
        lambda x: np.sin(x),
        lambda x: np.sin(2 * x)
    ]
    basis_names = ["1", "cos(x)", "cos(2x)", "sin(x)", "sin(2x)"]
    formula_tmpl = "y = {1} + {cos(x)}*cos(x) + {cos(2x)}*cos(2x) + {sin(x)}*sin(x) + {sin(2x)}*sin(2x)"
    return fit_linear_combination(points, basis_funcs, basis_names, formula_tmpl, "Mô hình lượng giác (Notebook Cell 7)")


def fit_polynomial(points: List[Tuple[float, float]], degree: int = 2) -> ModelResult:
    """
    Khớp đa thức bậc m: y = a0 + a1*x + a2*x^2 + ... + am*x^m
    """
    basis_funcs = [lambda x, p=p: x**p for p in range(degree + 1)]
    basis_names = [f"x^{p}" if p > 1 else ("x" if p == 1 else "1") for p in range(degree + 1)]

    terms = []
    for p in range(degree + 1):
        name = basis_names[p]
        terms.append(f"{{{name}}}*{name}" if name != "1" else f"{{{name}}}")
    formula_tmpl = "y = " + " + ".join(terms)

    return fit_linear_combination(points, basis_funcs, basis_names, formula_tmpl, f"Đa thức bậc {degree} (Polynomial)")


def fit_exponential(points: List[Tuple[float, float]]) -> ModelResult:
    """
    Khớp mô hình hàm mũ: y = a * exp(b*x)
    Tuyến tính hóa:
    - Nếu y > 0: ln(y) = ln(a) + b*x   => Y = A + B*X với A = ln(a) => a = exp(A), b = B.
    - Nếu y < 0: ln(-y) = ln(-a) + b*x => A = ln(-a) => a = -exp(A), b = B (như trong đề thi 20251GKtest1!).
    """
    # Tuyến tính hóa qua ln|y|
    transformed_points = []
    for x, y in points:
        val = abs(y) if abs(y) > 1e-12 else 1e-12
        transformed_points.append((x, np.log(val)))

    basis_funcs = [
        lambda x: np.ones_like(x),
        lambda x: x
    ]
    basis_names = ["1", "x"]

    df_phi, df_M, df_b_prime, df_a, a_coeffs, _ = least_squares_fit(
        transformed_points, basis_funcs, basis_names
    )

    A = float(a_coeffs[0, 0])
    B = float(a_coeffs[1, 0])

    # Tự động chọn dấu của a (+exp(A) hoặc -exp(A)) để SSE nhỏ nhất
    x_orig = np.array([p[0] for p in points])
    y_orig = np.array([p[1] for p in points]).reshape(-1, 1)

    a_pos = np.exp(A)
    a_neg = -np.exp(A)
    b = B

    sse_pos = np.sum((y_orig.flatten() - a_pos * np.exp(b * x_orig))**2)
    sse_neg = np.sum((y_orig.flatten() - a_neg * np.exp(b * x_orig))**2)
    a = a_pos if sse_pos <= sse_neg else a_neg

    def predict_fn(x_val: float) -> float:
        return float(a * np.exp(b * x_val))

    y_preds = np.array([predict_fn(x) for x in x_orig]).reshape(-1, 1)
    metrics = calculate_least_squares_error(y_orig, df_b_prime.values, a_coeffs, y_pred=y_preds)

    formula_str = f"y = ({a:+.8f}) * exp({b:+.8f} * x)"

    return ModelResult(
        name="Hàm mũ y = a * exp(b*x) (Notebook Cell 25)",
        formula_str=formula_str,
        predict_fn=predict_fn,
        coeffs={"A (ln|a|)": A, "B (b)": B, "a": a, "b": b},
        df_phi=df_phi,
        df_M=df_M,
        df_b_prime=df_b_prime,
        df_a=df_a,
        linear_coeffs=a_coeffs,
        metrics=metrics,
        transformed_points=transformed_points
    )


def fit_power_law(points: List[Tuple[float, float]]) -> ModelResult:
    """
    Khớp mô hình hàm lũy thừa: y = a * x^b (Cell 31)
    Tuyến tính hóa: ln|y| = ln|a| + b * ln|x| => Y = A + B*X với X = ln|x|, Y = ln|y|.
    """
    # Tuyến tính hóa qua ln|y| = ln|a| + b * ln|x|
    transformed_points = []
    for x, y in points:
        val_x = abs(x) if abs(x) > 1e-12 else 1e-12
        val_y = abs(y) if abs(y) > 1e-12 else 1e-12
        transformed_points.append((np.log(val_x), np.log(val_y)))

    basis_funcs = [
        lambda x: np.ones_like(x),
        lambda x: x
    ]
    basis_names = ["1", "ln(x)"]

    df_phi, df_M, df_b_prime, df_a, a_coeffs, _ = least_squares_fit(
        transformed_points, basis_funcs, basis_names
    )

    A = float(a_coeffs[0, 0])
    B = float(a_coeffs[1, 0])

    x_orig = np.array([p[0] for p in points])
    y_orig = np.array([p[1] for p in points]).reshape(-1, 1)

    a_pos = np.exp(A)
    a_neg = -np.exp(A)
    b = B

    sse_pos = np.sum((y_orig.flatten() - a_pos * (np.abs(x_orig) ** b))**2)
    sse_neg = np.sum((y_orig.flatten() - a_neg * (np.abs(x_orig) ** b))**2)
    a = a_pos if sse_pos <= sse_neg else a_neg

    def predict_fn(x_val: float) -> float:
        val_x = abs(x_val) if abs(x_val) > 1e-12 else 1e-12
        return float(a * (val_x ** b))

    y_preds = np.array([predict_fn(x) for x in x_orig]).reshape(-1, 1)
    metrics = calculate_least_squares_error(y_orig, df_b_prime.values, a_coeffs, y_pred=y_preds)

    formula_str = f"y = ({a:+.8f}) * x^({b:+.8f})"

    return ModelResult(
        name="Hàm lũy thừa y = a * x^b (Notebook Cell 31)",
        formula_str=formula_str,
        predict_fn=predict_fn,
        coeffs={"A (ln|a|)": A, "B (b)": B, "a": a, "b": b},
        df_phi=df_phi,
        df_M=df_M,
        df_b_prime=df_b_prime,
        df_a=df_a,
        linear_coeffs=a_coeffs,
        metrics=metrics,
        transformed_points=transformed_points
    )


def fit_midterm_k70(points: List[Tuple[float, float]]) -> ModelResult:
    """
    Khớp mô hình đề thi Giữa kỳ K70 (20251GKtest1):
    y = a * exp(b * x^2 + c / x)
    Tuyến tính hóa:
    ln|y| = ln|a| + b * x^2 + c * (1/x)
    => Y = A0 + b * X1 + c * X2 với X1 = x^2, X2 = 1/x, A0 = ln|a|.
    """
    transformed_points = []
    x_vals = []
    Y_vals = []
    for x, y in points:
        if abs(x) < 1e-12:
            continue
        val_y = abs(y) if abs(y) > 1e-12 else 1e-12
        x_vals.append(x)
        Y_vals.append(np.log(val_y))

    N = len(x_vals)
    x_arr = np.array(x_vals, dtype=np.float64)
    Y_arr = np.array(Y_vals, dtype=np.float64).reshape(-1, 1)

    phi_matrix = np.column_stack([np.ones(N), x_arr**2, 1.0 / x_arr])
    basis_names = ["1", "x^2", "1/x"]
    df_phi = pd.DataFrame(phi_matrix, columns=basis_names)

    phi_t = phi_matrix.T
    M = phi_t.dot(phi_matrix)
    df_M = pd.DataFrame(M, index=["Row_1", "Row_2", "Row_3"], columns=["Col_1", "Col_2", "Col_3"])

    b_prime = phi_t.dot(Y_arr)
    df_b_prime = pd.DataFrame(b_prime, index=["Row_1", "Row_2", "Row_3"], columns=["b'_j"])

    a_coeffs = np.linalg.solve(M, b_prime)
    df_a = pd.DataFrame(a_coeffs, columns=["a_j"], index=["A0 (ln|a|)", "b", "c"])

    A0 = float(a_coeffs[0, 0])
    b = float(a_coeffs[1, 0])
    c = float(a_coeffs[2, 0])

    # Tự động chọn dấu tối ưu của a (+exp(A0) hoặc -exp(A0))
    x_orig = np.array([p[0] for p in points])
    y_orig = np.array([p[1] for p in points]).reshape(-1, 1)

    a_pos = np.exp(A0)
    a_neg = -np.exp(A0)

    safe_exp = np.exp(b * (x_orig**2) + c / np.where(np.abs(x_orig) < 1e-12, 1e-12, x_orig))
    sse_pos = np.sum((y_orig.flatten() - a_pos * safe_exp)**2)
    sse_neg = np.sum((y_orig.flatten() - a_neg * safe_exp)**2)
    a = a_pos if sse_pos <= sse_neg else a_neg

    def predict_fn(x_val: float) -> float:
        if abs(x_val) < 1e-12:
            return 0.0
        return float(a * np.exp(b * (x_val**2) + c / x_val))

    x_orig = np.array([p[0] for p in points])
    y_orig = np.array([p[1] for p in points]).reshape(-1, 1)
    y_preds = np.array([predict_fn(x) for x in x_orig]).reshape(-1, 1)
    metrics = calculate_least_squares_error(y_orig, b_prime, a_coeffs, y_pred=y_preds)

    formula_str = f"y = ({a:+.8f}) * exp({b:+.8f}*x^2 + ({c:+.8f})/x)"

    return ModelResult(
        name="Hàm thực nghiệm Midterm 20251 y = a * exp(b*x^2 + c/x)",
        formula_str=formula_str,
        predict_fn=predict_fn,
        coeffs={"A0 (ln|a|)": A0, "b": b, "c": c, "a": a},
        df_phi=df_phi,
        df_M=df_M,
        df_b_prime=df_b_prime,
        df_a=df_a,
        linear_coeffs=a_coeffs,
        metrics=metrics
    )


def fit_hyperbolic(points: List[Tuple[float, float]], model_type: str = "a+b/x") -> ModelResult:
    """
    Khớp các mô hình hypebolic/hữu tỉ:
    - Type 'a+b/x': y = a + b/x => basis [1, 1/x]
    - Type '1/(a+bx)': 1/y = a + bx => Y = 1/y, basis [1, x]
    - Type 'x/(a+bx)': y = x/(a+bx) => x/y = a + bx => Y = x/y, basis [1, x]
    """
    if model_type == "a+b/x":
        basis_funcs = [lambda x: np.ones_like(x), lambda x: 1.0 / np.where(np.abs(x) < 1e-12, 1e-12, x)]
        basis_names = ["1", "1/x"]
        return fit_linear_combination(points, basis_funcs, basis_names, "y = {1} + {1/x} / x", "Hypebolic y = a + b/x")

    elif model_type == "1/(a+bx)":
        trans_pts = []
        for x, y in points:
            val_y = y if abs(y) > 1e-12 else 1e-12
            trans_pts.append((x, 1.0 / val_y))

        basis_funcs = [lambda x: np.ones_like(x), lambda x: x]
        basis_names = ["1", "x"]
        df_phi, df_M, df_b_prime, df_a, a_coeffs, _ = least_squares_fit(trans_pts, basis_funcs, basis_names)
        a = float(a_coeffs[0, 0])
        b = float(a_coeffs[1, 0])

        def predict_fn(x_val: float) -> float:
            denom = a + b * x_val
            return float(1.0 / denom) if abs(denom) > 1e-12 else 0.0

        x_orig = np.array([p[0] for p in points])
        y_orig = np.array([p[1] for p in points]).reshape(-1, 1)
        y_preds = np.array([predict_fn(x) for x in x_orig]).reshape(-1, 1)
        metrics = calculate_least_squares_error(y_orig, df_b_prime.values, a_coeffs, y_pred=y_preds)

        return ModelResult(
            name="Hypebolic y = 1 / (a + bx)",
            formula_str=f"y = 1 / ({a:+.6f} + {b:+.6f}*x)",
            predict_fn=predict_fn,
            coeffs={"a": a, "b": b},
            df_phi=df_phi,
            df_M=df_M,
            df_b_prime=df_b_prime,
            df_a=df_a,
            linear_coeffs=a_coeffs,
            metrics=metrics
        )

    else:  # x/(a+bx)
        trans_pts = []
        for x, y in points:
            val_y = y if abs(y) > 1e-12 else 1e-12
            trans_pts.append((x, x / val_y))
        basis_funcs = [lambda x: np.ones_like(x), lambda x: x]
        basis_names = ["1", "x"]
        df_phi, df_M, df_b_prime, df_a, a_coeffs, _ = least_squares_fit(trans_pts, basis_funcs, basis_names)
        a = float(a_coeffs[0, 0])
        b = float(a_coeffs[1, 0])

        def predict_fn(x_val: float) -> float:
            denom = a + b * x_val
            return float(x_val / denom) if abs(denom) > 1e-12 else 0.0

        x_orig = np.array([p[0] for p in points])
        y_orig = np.array([p[1] for p in points]).reshape(-1, 1)
        y_preds = np.array([predict_fn(x) for x in x_orig]).reshape(-1, 1)
        metrics = calculate_least_squares_error(y_orig, df_b_prime.values, a_coeffs, y_pred=y_preds)

        return ModelResult(
            name="Hypebolic y = x / (a + bx)",
            formula_str=f"y = x / ({a:+.6f} + {b:+.6f}*x)",
            predict_fn=predict_fn,
            coeffs={"a": a, "b": b},
            df_phi=df_phi,
            df_M=df_M,
            df_b_prime=df_b_prime,
            df_a=df_a,
            linear_coeffs=a_coeffs,
            metrics=metrics
        )


def fit_logarithmic(points: List[Tuple[float, float]]) -> ModelResult:
    """
    Khớp mô hình logarit: y = a + b * ln(x)
    """
    trans_pts = [(np.log(abs(x) if abs(x) > 1e-12 else 1e-12), y) for x, y in points]
    basis_funcs = [lambda x: np.ones_like(x), lambda x: x]
    basis_names = ["1", "ln(x)"]
    return fit_linear_combination(trans_pts, basis_funcs, basis_names, "y = {1} + {ln(x)}*ln(x)", "Logarit y = a + b*ln(x)")


# ==============================================================================
# 6. SO SÁNH TẤT CẢ MÔ HÌNH (COMPARE ALL MODELS)
# ==============================================================================
def compare_all_models(points: List[Tuple[float, float]], eval_x: Optional[float] = None) -> pd.DataFrame:
    """
    Chạy tất cả các mô hình có thể áp dụng trên tập điểm và lập bảng so sánh theo R^2 và RMSE.
    """
    results: List[ModelResult] = []

    # Danh sách các hàm thử nghiệm
    model_runners = [
        ("Đa thức bậc 1 (y = a0 + a1*x)", lambda: fit_polynomial(points, degree=1)),
        ("Đa thức bậc 2 (Parabol)", lambda: fit_polynomial(points, degree=2)),
        ("Đa thức bậc 3 (Cubic)", lambda: fit_polynomial(points, degree=3)),
        ("Hàm y = a*x + b*sin(x)", lambda: fit_notebook_linear_x_sinx(points)),
        ("Hàm mũ y = a*exp(b*x)", lambda: fit_exponential(points)),
        ("Hàm lũy thừa y = a*x^b", lambda: fit_power_law(points)),
        ("Hypebolic y = a + b/x", lambda: fit_hyperbolic(points, "a+b/x")),
        ("Hypebolic y = 1/(a + bx)", lambda: fit_hyperbolic(points, "1/(a+bx)")),
        ("Hàm Midterm y = a*exp(b*x^2 + c/x)", lambda: fit_midterm_k70(points)),
    ]

    for name, runner in model_runners:
        try:
            res = runner()
            results.append(res)
        except Exception as e:
            # Bỏ qua mô hình nếu dữ liệu không thỏa điều kiện (ví dụ x=0 cho 1/x hoặc ln)
            pass

    # Lập bảng so sánh
    rows = []
    for r in results:
        pred_val = r.predict_fn(eval_x) if eval_x is not None else np.nan
        rows.append({
            "Mô hình": r.name,
            "R^2 (Độ phù hợp)": r.metrics["r2"],
            "RMSE": r.metrics["rmse"],
            "SSE (S_min)": r.metrics["sse"],
            "MAE": r.metrics["mae"],
            f"y_hat(x*={eval_x})": pred_val if eval_x is not None else None,
            "Phương trình khớp": r.formula_str
        })

    df_comp = pd.DataFrame(rows)
    # Sắp xếp theo R^2 giảm dần (càng gần 1 càng tốt)
    df_comp = df_comp.sort_values(by="R^2 (Độ phù hợp)", ascending=False).reset_index(drop=True)

    print("\n" + "=" * 90)
    print("BẢNG TỔNG HỢP SO SÁNH CÁC MÔ HÌNH HỒI QUY (SẮP XẾP THEO R^2 GIẢM DẦN)")
    print("=" * 90)
    cols_to_show = ["Mô hình", "R^2 (Độ phù hợp)", "RMSE", "SSE (S_min)"]
    if eval_x is not None:
        cols_to_show.append(f"y_hat(x*={eval_x})")
    cols_to_show.append("Phương trình khớp")

    print(df_comp[cols_to_show].to_string(index=True))
    best_model = df_comp.iloc[0]
    print(f"\n-> MÔ HÌNH TỐT NHẤT: '{best_model['Mô hình']}' với R^2 = {best_model['R^2 (Độ phù hợp)']:.8f}, RMSE = {best_model['RMSE']:.8f}")

    return df_comp


# ==============================================================================
# 7. CHƯƠNG TRÌNH CHÍNH & CLI (MAIN ENTRY POINT)
# ==============================================================================
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description="Chương trình Giải Thuật toán Bình Phương Bé Nhất (Least Squares Fit) - Numerical Analysis PPS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ví dụ sử dụng:
  # 1. Chạy trên file dữ liệu mặc định 0.csv với mô hình tuyến tính x + sin(x):
  python lsq.py 0.csv --model x_sinx

  # 2. Khớp hàm mũ y = a*exp(b*x) và dự báo tại x = 2.5:
  python lsq.py 0.csv --model exp --eval_x 2.5

  # 3. So sánh TẤT CẢ các mô hình và chọn mô hình tối ưu nhất:
  python lsq.py 0.csv --compare

  # 4. Khớp đa thức bậc 3:
  python lsq.py 0.csv --model poly --degree 3

  # 5. Khớp mô hình đề thi Giữa kỳ 20251GKtest1:
  python lsq.py ../../Midterm/20251GKtest1/20251GKtest1.csv --model midterm
"""
    )
    parser.add_argument("file", nargs="?", default=None, help="Đường dẫn tới file dữ liệu (.csv, .txt, .xlsx). Mặc định dùng input.txt cạnh script")
    parser.add_argument("--model", choices=["all", "x_sinx", "exp", "power", "poly", "trig", "hyper1", "hyper2", "hyper3", "midterm", "log"], default="all", help="Mô hình hồi quy cần tính toán (mặc định: all để so sánh tất cả)")
    parser.add_argument("--degree", type=int, default=2, help="Bậc của đa thức nếu chọn --model poly (mặc định: 2)")
    parser.add_argument("--eval_x", type=float, default=None, help="Mốc x* cần dự báo giá trị y_hat(x*)")
    parser.add_argument("--compare", action="store_true", help="So sánh toàn bộ mô hình và tìm mô hình có R^2 cao nhất")
    parser.add_argument("--output", type=str, default=os.path.join(script_dir, "output.txt"), help="Đường dẫn file lưu kết quả báo cáo (mặc định output.txt cạnh script)")

    args = parser.parse_args()

    # Xác định file đầu vào
    input_file = args.file
    if input_file is None:
        candidates = ["input.txt", "0.csv", "data.csv"]
        for c in candidates:
            candidate_path = os.path.join(script_dir, c)
            if os.path.exists(candidate_path):
                input_file = candidate_path
                break

    if input_file is None:
        print("Không tìm thấy file 0.csv hoặc input.txt!")
        print("Đang chạy chế độ kiểm thử mẫu từ notebook pp6_LSQ.ipynb...")
        # Lấy dữ liệu mẫu Cell 26
        points = [
            (1.0, 7.1),
            (2.0, 27.8),
            (3.0, 63.1),
            (4.0, 110.0),
            (5.0, 161.0)
        ]
        params = {"eval_x": 3.5}
        print(f"Dữ liệu mẫu ({len(points)} điểm): {points}")
    else:
        print(f"-> Đang tải dữ liệu từ: {input_file}")
        points, params = load_dataset(input_file)
        if not points:
            print(f"Không thể đọc được điểm dữ liệu nào từ '{input_file}'. Vui lòng kiểm tra lại!")
            return

    eval_x = args.eval_x if args.eval_x is not None else params.get("eval_x", None)

    print(f"-> Tổng số mốc dữ liệu n = {len(points)}")
    print(f"-> Phạm vi x: [{min(p[0] for p in points):.4f}, {max(p[0] for p in points):.4f}]")
    print(f"-> Phạm vi y: [{min(p[1] for p in points):.4f}, {max(p[1] for p in points):.4f}]")
    if eval_x is not None:
        print(f"-> Mốc dự báo eval_x (x*) = {eval_x}")

    import io
    from contextlib import redirect_stdout

    report_buffer = io.StringIO()
    with redirect_stdout(report_buffer):
        if args.compare or args.model == "all":
            compare_all_models(points, eval_x=eval_x)
            print("\n" + "=" * 80)
            print("CHI TIẾT MÔ HÌNH GỐC TRONG NOTEBOOK pp6_LSQ.ipynb:")
            print("=" * 80)
            fit_notebook_linear_x_sinx(points).print_report(eval_x=eval_x)
        else:
            res = None
            if args.model == "x_sinx":
                res = fit_notebook_linear_x_sinx(points)
            elif args.model == "exp":
                res = fit_exponential(points)
            elif args.model == "power":
                res = fit_power_law(points)
            elif args.model == "poly":
                res = fit_polynomial(points, degree=args.degree)
            elif args.model == "trig":
                res = fit_trigonometric(points)
            elif args.model == "hyper1":
                res = fit_hyperbolic(points, "a+b/x")
            elif args.model == "hyper2":
                res = fit_hyperbolic(points, "1/(a+bx)")
            elif args.model == "hyper3":
                res = fit_hyperbolic(points, "x/(a+bx)")
            elif args.model == "midterm":
                res = fit_midterm_k70(points)
            elif args.model == "log":
                res = fit_logarithmic(points)

            if res:
                res.print_report(eval_x=eval_x)

    rounding_input = input_file or os.path.join(script_dir, "input.txt")
    report = format_report_output(report_buffer.getvalue(), rounding_input)
    print(report)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"-> Đã lưu thành công vào '{args.output}'!")


if __name__ == "__main__":
    main()
