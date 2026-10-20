"""
Script giải chi tiết đề thi Giữa kỳ 20251 Đề 1 (sử dụng 20251GKtest1.csv)
Câu a: Hồi quy phi tuyến thực nghiệm y = a * exp(b*x^2 + c/x)
Câu b: Phương pháp lặp giải f(x) = y_bar
Câu c: Nội suy Bessel 8 mốc tính f(x_bar) và f'(x_bar)
"""
import os, math
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, "e:/2026.1/PPS/Week6-Interpolation/Python_Interpolation")
from pps_io import parse_universal_input

csv_path = "e:/2026.1/PPS/Midterm/20251GKtest1/20251GKtest1.csv"
points, params = parse_universal_input(csv_path)

x_bar = params["eval_x"]
y_bar = params["eval_y"]
x_pts = np.array([p[0] for p in points])
y_pts = np.array([p[1] for p in points])

print(f"Tổng số mốc dữ liệu: N = {len(x_pts)}")
print(f"x_bar = {x_bar}")
print(f"y_bar = {y_bar}")
print(f"Khoảng của x: [{x_pts[0]}, {x_pts[-1]}]")
print(f"Khoảng của y: [{y_pts.min()}, {y_pts.max()}]")

# ==============================================================================
# CÂU A: TÌM HÀM THỰC NGHIỆM y = a * exp(b*x^2 + c/x) VÀ XÁC ĐỊNH SAI SỐ
# ==============================================================================
print("\n" + "="*80)
print("CÂU A: TÌM HÀM THỰC NGHIỆM y = a * exp(b*x^2 + c/x)")
print("="*80)
# Nhận xét: Các giá trị y_i trong 20251GKtest1 đều ÂM (y_i < 0)!
# Hàm mũ exp(...) luôn dương, do đó hệ số a phải ÂM (a < 0).
# Ta đặt: -y = (-a) * exp(b*x^2 + c/x) với -a > 0.
# Lấy ln hai vế:
# ln(-y) = ln(-a) + b * x^2 + c * (1/x)
# Đặt:
# Y = ln(-y)
# A0 = ln(-a) => a = -exp(A0)
# X1 = x^2
# X2 = 1/x
# Ta có mô hình hồi quy tuyến tính bội: Y = A0 + b * X1 + c * X2
# Phương pháp bình phương bé nhất (OLS): (Phi^T * Phi) * [A0, b, c]^T = Phi^T * Y

N = len(x_pts)
all_negative = np.all(y_pts < 0)
print(f"Kiểm tra dấu y_i: Toàn bộ y_i < 0 ? -> {all_negative}")

Y = np.log(-y_pts) if all_negative else np.log(np.abs(y_pts))
Phi = np.column_stack([np.ones(N), x_pts**2, 1.0 / x_pts])

# Giải hệ phương trình chuẩn
theta, residuals, rank, s = np.linalg.lstsq(Phi, Y, rcond=None)
A0, b, c = theta
a = -np.exp(A0) if all_negative else np.exp(A0)

print(f"Hệ số tuyến tính hóa:")
print(f"  A0 = ln(|a|) = {A0:.8f}")
print(f"  b  = {b:.8f}")
print(f"  c  = {c:.8f}")
print(f"Hệ số của hàm thực nghiệm y = a * exp(b*x^2 + c/x):")
print(f"  a = {a:.8f}")
print(f"  b = {b:.8f}")
print(f"  c = {c:.8f}")
print(f"Phương trình hàm thực nghiệm:")
print(f"  y = ({a:.8f}) * exp({b:.8f}*x^2 + ({c:.8f})/x)")

# Tính giá trị dự báo y_pred và các loại sai số
y_pred = a * np.exp(b * (x_pts**2) + c / x_pts)

# Sai số tuyệt đối tại từng điểm e_i = |y_i - y_pred|
abs_errors = np.abs(y_pts - y_pred)
rel_errors = np.abs((y_pts - y_pred) / y_pts) * 100 # phần trăm

# Các chỉ số sai số tổng hợp:
sse = np.sum((y_pts - y_pred)**2) # Tổng bình phương sai số (Sum of Squared Errors)
mse = np.mean((y_pts - y_pred)**2) # Trung bình bình phương sai số
rmse = np.sqrt(mse) # Sai số căn quân phương (Root Mean Square Error)
mae = np.mean(abs_errors) # Sai số tuyệt đối trung bình (Mean Absolute Error)
max_err = np.max(abs_errors) # Sai số tuyệt đối cực đại
mape = np.mean(rel_errors) # Sai số tương đối phần trăm trung bình

# Hệ số xác định R^2
ss_tot = np.sum((y_pts - np.mean(y_pts))**2)
r2 = 1 - (sse / ss_tot)

print("\n--- Đánh giá sai số của hàm thực nghiệm ---")
print(f"  1. Tổng bình phương sai số (SSE)      : {sse:.8f}")
print(f"  2. Sai số tuyệt đối trung bình (MAE)   : {mae:.8f}")
print(f"  3. Sai số căn quân phương (RMSE)       : {rmse:.8f}")
print(f"  4. Sai số tuyệt đối lớn nhất (Max Err) : {max_err:.8f}")
print(f"  5. Sai số tương đối trung bình (MAPE) : {mape:.4f}%")
print(f"  6. Hệ số xác định R^2                 : {r2:.8f}")

# ==============================================================================
# CÂU C: NỘI SUY BESSEL VỚI 8 MỐC ĐỂ TÍNH f(x_bar) VÀ f'(x_bar)
# ==============================================================================
print("\n" + "="*80)
print(f"CÂU C: ĐA THỨC NỘI SUY BESSEL 8 MỐC TẠI x_bar = {x_bar}")
print("="*80)

# 1. Chọn 8 mốc nội suy:
# Ta cần 8 mốc chẵn đối xứng quanh x_bar = 4.15.
# Tìm mốc k sao cho x_k <= x_bar < x_{k+1}
h = x_pts[1] - x_pts[0]
diffs = np.abs(x_pts - x_bar)
k_nearest = np.argmin(diffs)
if x_pts[k_nearest] > x_bar and k_nearest > 0:
    k_base = k_nearest - 1
else:
    k_base = k_nearest

# Bessel với 8 mốc (bậc 7): m = 8 // 2 = 4 mốc mỗi bên.
# Cửa sổ từ k_base - 3 đến k_base + 4 (tổng 8 mốc: x_{-3}, x_{-2}, x_{-1}, x_0, x_1, x_2, x_3, x_4)
idx_start = k_base - 3
idx_end = idx_start + 8
bessel_x = x_pts[idx_start:idx_end]
bessel_y = y_pts[idx_start:idx_end]
x0 = bessel_x[3] # mốc x_0
x1 = bessel_x[4] # mốc x_1

print(f"Bước nhảy h = {h:.6f}")
print(f"Khoảng chứa x_bar = {x_bar}: [{x0:.4f}, {x1:.4f}]")
print(f"8 mốc được chọn (từ chỉ số {idx_start} đến {idx_end-1} trong file):")
for i, (bx, by) in enumerate(zip(bessel_x, bessel_y)):
    idx_rel = i - 3
    print(f"  x_{{{idx_rel}}} = {bx:.4f}, y_{{{idx_rel}}} = {by:.6f}")

# Lập bảng sai phân hữu hạn cho 8 mốc
n_pts = len(bessel_x)
diff_table = np.full((n_pts, n_pts), np.nan)
diff_table[:, 0] = bessel_y
for j in range(1, n_pts):
    for i in range(n_pts - j):
        diff_table[i, j] = diff_table[i + 1, j - 1] - diff_table[i, j - 1]

diff_cols = ["x_i", "y_i"] + [f"Delta^{j}" for j in range(1, n_pts)]
bessel_df_data = {"x_i": bessel_x, "y_i": bessel_y}
for j in range(1, n_pts):
    bessel_df_data[f"Delta^{j}"] = [diff_table[i, j] if i < n_pts - j else np.nan for i in range(n_pts)]
df_bessel_tbl = pd.DataFrame(bessel_df_data)
print("\nBảng sai phân hữu hạn 8 mốc:")
print(df_bessel_tbl.to_string(index=True))

# 2. Xây dựng đa thức Bessel:
# Chuẩn hóa biến: t = (x - x0) / h
# Biến Bessel: p = t - 0.5 (hoặc viết theo t)
# Công thức Bessel theo t:
# P(t) = y_{1/2} + (t - 1/2) Delta y_0 + (t(t-1)/2) Delta^2 y_{1/2} + ((t-1/2)t(t-1)/6) Delta^3 y_0
#       + ((t+1)t(t-1)(t-2)/24) Delta^4 y_{1/2} + ((t-1/2)(t+1)t(t-1)(t-2)/120) Delta^5 y_0
#       + ((t+2)(t+1)t(t-1)(t-2)(t-3)/720) Delta^6 y_{1/2} + ((t-1/2)(t+2)(t+1)t(t-1)(t-2)(t-3)/5040) Delta^7 y_0
# Với:
# y_{1/2} = (y_0 + y_1) / 2
# Delta y_0 = y_1 - y_0
# Delta^2 y_{1/2} = (Delta^2 y_{-1} + Delta^2 y_0) / 2
# Delta^3 y_0 = Delta^3 y_{-1}
# Delta^4 y_{1/2} = (Delta^4 y_{-2} + Delta^4 y_{-1}) / 2
# Delta^5 y_0 = Delta^5 y_{-2}
# Delta^6 y_{1/2} = (Delta^6 y_{-3} + Delta^6 y_{-2}) / 2
# Delta^7 y_0 = Delta^7 y_{-3}

# Ở đây chỉ số hàng trong diff_table:
# x_{-3} là hàng 0
# x_{-2} là hàng 1
# x_{-1} là hàng 2
# x_0    là hàng 3
# x_1    là hàng 4
# x_2    là hàng 5
# x_3    là hàng 6
# x_4    là hàng 7

row_x0 = 3
row_x1 = 4

# Thu thập các hệ số sai phân C_j:
C = []
C.append(0.5 * (diff_table[row_x0, 0] + diff_table[row_x1, 0]))       # j=0
C.append(diff_table[row_x0, 1])                                       # j=1
C.append(0.5 * (diff_table[row_x0-1, 2] + diff_table[row_x0, 2]))     # j=2
C.append(diff_table[row_x0-1, 3])                                     # j=3
C.append(0.5 * (diff_table[row_x0-2, 4] + diff_table[row_x0-1, 4]))   # j=4
C.append(diff_table[row_x0-2, 5])                                     # j=5
C.append(0.5 * (diff_table[row_x0-3, 6] + diff_table[row_x0-2, 6]))   # j=6
C.append(diff_table[row_x0-3, 7])                                     # j=7

print("\nCác hệ số sai phân trung tâm Bessel C_j:")
for j, c_val in enumerate(C):
    print(f"  C_{j} = {c_val:.8f}")

# Xây dựng các đa thức cơ sở B_j(t) theo biến t:
# B_0(t) = 1
# B_1(t) = t - 1/2
# B_2(t) = t(t-1) / 2 = 1/2 t^2 - 1/2 t
# B_3(t) = (t - 1/2) t(t-1) / 6 = (t - 1/2) * B_2(t) / 3
# B_{2k}(t)   = B_{2k-2}(t) * (t + k - 1)(t - k) / (2k * (2k-1))
# B_{2k+1}(t) = (t - 1/2) / (2k+1) * B_{2k}(t)

# Sử dụng numpy polynomial để biểu diễn chính xác P(t):
B_poly = [np.array([1.0]), np.array([-0.5, 1.0])]
for i in range(2, 8):
    if i % 2 == 0:
        k = i // 2
        fact = np.convolve(np.array([float(k - 1), 1.0]), np.array([-float(k), 1.0]))
        Bi = np.convolve(B_poly[i - 2], fact)
    else:
        Bi = np.convolve(B_poly[i - 1], [-0.5, 1.0])
    B_poly.append(Bi)

# Đa thức P(t) = sum_{j=0}^7 (C_j / j!) * B_j(t)
P_t_coeffs = np.zeros(8, dtype=float)
bessel_steps = []
for j in range(8):
    Dj = C[j] / math.factorial(j)
    term_poly = Dj * B_poly[j]
    P_t_coeffs[:len(term_poly)] += term_poly
    bessel_steps.append({
        "j": j,
        "C_j": C[j],
        "D_j = C_j/j!": Dj,
        "B_j(t) (0->cao)": [round(v, 6) for v in B_poly[j].tolist()],
        "Term(t)": [round(v, 6) for v in term_poly.tolist()]
    })

print("\nBảng chi tiết các số hạng Bessel:")
df_bessel_steps = pd.DataFrame(bessel_steps)
print(df_bessel_steps[["j", "C_j", "D_j = C_j/j!"]].to_string(index=False))

print("\nHệ số đa thức P(t) = a_0 + a_1*t + a_2*t^2 + ... + a_7*t^7:")
for p_deg, p_coeff in enumerate(P_t_coeffs):
    print(f"  t^{p_deg}: {p_coeff:.10f}")

# Tính giá trị tại x_bar = 4.15:
t_bar = (x_bar - x0) / h
f_xbar = sum(P_t_coeffs[deg] * (t_bar ** deg) for deg in range(8))

# Tính đạo hàm P'(t) và f'(x_bar):
# Vì x = x0 + t*h nên df/dx = (df/dt) * (dt/dx) = (1/h) * P'(t)
P_prime_coeffs = [deg * P_t_coeffs[deg] for deg in range(1, 8)]
P_prime_t_bar = sum(P_prime_coeffs[deg-1] * (t_bar ** (deg - 1)) for deg in range(1, 8))
df_xbar = P_prime_t_bar / h

print(f"\n--- Kết quả tính gần đúng tại x_bar = {x_bar} ---")
print(f"  t_bar = (x_bar - x0)/h = ({x_bar} - {x0:.4f}) / {h:.6f} = {t_bar:.8f}")
print(f"  P(t_bar) = f(x_bar)   = {f_xbar:.10f}")
print(f"  P'(t_bar)              = {P_prime_t_bar:.10f}")
print(f"  f'(x_bar) = P'/h      = {df_xbar:.10f}")

# ==============================================================================
# CÂU B: PHƯƠNG PHÁP LẶP GIẢI PHƯƠNG TRÌNH f(x) = y_bar
# ==============================================================================
print("\n" + "="*80)
print(f"CÂU B: PHƯƠNG PHÁP LẶP GIẢI PHƯƠNG TRÌNH f(x) = y_bar ({y_bar})")
print("="*80)

# 1. Tìm khoảng phân ly nghiệm [a, b]:
# Dò trong bảng dữ liệu:
idx_iso = None
for i in range(len(y_pts) - 1):
    if (y_pts[i] - y_bar) * (y_pts[i+1] - y_bar) <= 0:
        idx_iso = i
        break

if idx_iso is not None:
    xa, xb = x_pts[idx_iso], x_pts[idx_iso+1]
    ya, yb = y_pts[idx_iso], y_pts[idx_iso+1]
    print(f"Khoảng phân ly nghiệm trên bảng dữ liệu: [{xa}, {xb}]")
    print(f"  y({xa}) = {ya} <= y_bar = {y_bar} <= y({xb}) = {yb}")

# Cách tiếp cận chuẩn 1: Phương pháp lặp nội suy ngược (Inverse Interpolation Loop)
# Sử dụng đa thức nội suy Newton tiến bậc 2 hoặc 3 trên 3-4 mốc lân cận [xa, xb]
# P(x) = y_0 + t * Delta y_0 + t(t-1)/2 Delta^2 y_0 + ... = y_bar
# Chuyển về dạng lặp: t_{k+1} = phi(t_k)
# t = (y_bar - y_0 - t(t-1)/2 Delta^2 y_0 - ... ) / Delta y_0
print("\n--- Cách 1: Phương pháp lặp nội suy ngược (chuẩn theo giáo trình PPS) ---")
# Lấy 4 mốc quanh xa:
idx_start_iter = max(0, idx_iso - 1)
iter_x = x_pts[idx_start_iter:idx_start_iter+4]
iter_y = y_pts[idx_start_iter:idx_start_iter+4]
print("4 mốc lân cận dùng để lặp nội suy ngược:")
for ix, iy in zip(iter_x, iter_y):
    print(f"  x = {ix:.4f}, y = {iy:.6f}")

# Sai phân hữu hạn cho 4 mốc:
dy0 = iter_y[1] - iter_y[0]
d2y0 = (iter_y[2] - iter_y[1]) - (iter_y[1] - iter_y[0])
d3y0 = ((iter_y[3] - iter_y[2]) - (iter_y[2] - iter_y[1])) - d2y0

print(f"  y_0 = {iter_y[0]:.6f}")
print(f"  Delta y_0   = {dy0:.6f}")
print(f"  Delta^2 y_0 = {d2y0:.6f}")
print(f"  Delta^3 y_0 = {d3y0:.6f}")

# Hàm lặp phi(t):
# t = (y_bar - y_0) / Delta y_0 - (Delta^2 y_0 / (2 * Delta y_0)) * t(t-1) - (Delta^3 y_0 / (6 * Delta y_0)) * t(t-1)(t-2)
t0_init = (y_bar - iter_y[0]) / dy0
print(f"Giá trị lặp ban đầu t^(0) = (y_bar - y_0) / Delta y_0 = {t0_init:.8f}")

t_curr = t0_init
eps = 1e-8
iter_history = []
iter_history.append({"k": 0, "t_k": t_curr, "x_k": iter_x[0] + t_curr * h, "diff": np.nan})

for it in range(1, 20):
    t_prev = t_curr
    # Công thức lặp phi(t):
    correction = (d2y0 / (2 * dy0)) * t_prev * (t_prev - 1) + (d3y0 / (6 * dy0)) * t_prev * (t_prev - 1) * (t_prev - 2)
    t_curr = t0_init - correction
    diff = abs(t_curr - t_prev)
    x_k = iter_x[0] + t_curr * h
    iter_history.append({"k": it, "t_k": t_curr, "x_k": x_k, "diff": diff})
    if diff < eps:
        break

df_iter = pd.DataFrame(iter_history)
print("\nBảng các bước lặp nội suy ngược:")
print(df_iter.to_string(index=False))
x_root_inv = iter_x[0] + t_curr * h
print(f"-> Nghiệm x tìm được từ lặp nội suy ngược: x* = {x_root_inv:.8f}")

# Cách tiếp cận chuẩn 2: Phương pháp lặp đơn trên hàm thực nghiệm (ở câu a)
# Ta có hàm thực nghiệm: y(x) = a * exp(b*x^2 + c/x) = y_bar
# ln(-y_bar) = ln(-a) + b*x^2 + c/x
# => b*x^2 + c/x + (ln(-a) - ln(-y_bar)) = 0
# Chia hai vế cho x hoặc rút x = phi(x) để giải phương trình
print("\n--- Cách 2: Giải f(x) = y_bar trên hàm thực nghiệm câu a ---")
# ln(-y_bar) = A0 + b*x^2 + c/x
const_val = A0 - np.log(-y_bar)
print(f"Phương trình tương đương: b*x^2 + c/x + {const_val:.8f} = 0")
# Đặt g(x) = b*x^2 + c/x + const_val = 0
# Dùng phương pháp lặp đơn: x = phi(x) = c / (-const_val - b*x^2)
# hoặc phương pháp Newton-Raphson: x_{k+1} = x_k - g(x_k)/g'(x_k)
g = lambda x: b * (x**2) + c / x + const_val
g_prime = lambda x: 2 * b * x - c / (x**2)

x_newton = x_root_inv
newton_steps = []
newton_steps.append({"k": 0, "x_k": x_newton, "g(x_k)": g(x_newton)})
for it in range(1, 15):
    x_prev = x_newton
    x_newton = x_prev - g(x_prev) / g_prime(x_prev)
    newton_steps.append({"k": it, "x_k": x_newton, "g(x_k)": g(x_newton)})
    if abs(x_newton - x_prev) < 1e-8:
        break
df_newton = pd.DataFrame(newton_steps)
print("Bảng lặp Newton trên hàm thực nghiệm:")
print(df_newton.to_string(index=False))
print(f"-> Nghiệm x trên hàm thực nghiệm: x* = {x_newton:.8f}")
