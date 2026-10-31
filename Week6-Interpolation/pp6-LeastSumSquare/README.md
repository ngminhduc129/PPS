# Phương Pháp Bình Phương Bé Nhất (Least Squares Fitting)

Mã nguồn Python độc lập được chuyển đổi và nâng cấp toàn diện từ file notebook [`pp6_LSQ.ipynb`](file:///e:/2026.1/PPS/Week6-Interpolation/pp6-LeastSumSquare/pp6_LSQ.ipynb).

---

## 1. File Mã Nguồn
- **[`lsq.py`](file:///e:/2026.1/PPS/Week6-Interpolation/pp6-LeastSumSquare/lsq.py)**: Script chính chứa toàn bộ thuật toán, bộ đọc dữ liệu đa năng, các mô hình thực nghiệm và giao diện dòng lệnh (CLI).
- **[`pp6_lsq.py`](file:///e:/2026.1/PPS/Week6-Interpolation/pp6-LeastSumSquare/pp6_lsq.py)**: Điểm gọi phụ (alias entry point).
- **[`0.csv`](file:///e:/2026.1/PPS/Week6-Interpolation/pp6-LeastSumSquare/0.csv)**: Dữ liệu thực nghiệm mẫu từ notebook gốc (141 điểm).

---

## 2. Các Tính Năng Đã Tích Hợp

### 2.1. Đọc dữ liệu vạn năng (Universal Input Parser)
- Hỗ trợ định dạng `.csv`, `.txt`, `.xlsx` từ các đề thi Giữa kỳ HUST (20241, 20251,...).
- Tự động nhận diện dấu phân cách (`;`, `,`, tab `\t`, dấu cách).
- Tự động xử lý số thực có phần thập phân dùng dấu phẩy `,` (ví dụ `2,95`).
- Tự động trích xuất điểm dự báo $x^*$ (`eval_x`) và $y^*$ (`eval_y`).

### 2.2. Thuật toán Phương Trình Chuẩn (Normal Equations)
- Ma trận thiết kế $\Phi$ kích thước $n \times m$.
- Ma trận Gram đối xứng $M = \Phi^T \Phi$ kích thước $m \times m$.
- Véctơ vế phải $b' = \Phi^T y$ kích thước $m \times 1$.
- Nghiệm $a = M^{-1} b'$ giải bằng phép phân rã đại số tuyến tính chính xác cao.
- Toàn bộ các chỉ số đánh giá sai số:
  - $SSE = \sum (y_i - \hat{y}_i)^2 = y^Ty - a^T b'$ (Tổng bình phương sai số / $S_{\min}$)
  - $MSE = \frac{SSE}{n}$ (Trung bình bình phương sai số)
  - $RMSE = \sqrt{MSE}$ (Căn sai số trung bình)
  - $MAE = \frac{1}{n} \sum |y_i - \hat{y}_i|$ (Sai số tuyệt đối trung bình)
  - $\text{Max Error} = \max |y_i - \hat{y}_i|$ (Sai số cực đại)
  - $MAPE = \frac{1}{n} \sum \left|\frac{y_i - \hat{y}_i}{y_i}\right| \times 100\%$ (Sai số tương đối phần trăm)
  - $R^2 = 1 - \frac{SSE}{SS_{\text{tot}}}$ (Hệ số xác định độ phù hợp)

### 2.3. Danh sách các mô hình hồi quy (Models Supported)
| Mô hình | Phương trình toán học | Tham số dòng lệnh | Gốc từ notebook |
| :--- | :--- | :--- | :--- |
| **Linear Basis** | $y = a \cdot x + b \cdot \sin(x)$ | `--model x_sinx` | Cell 20 |
| **Lượng giác Fourier** | $y = a_1 + a_2 \cos x + a_3 \cos 2x + a_4 \sin x + a_5 \sin 2x$ | `--model trig` | Cell 7 |
| **Hàm mũ (Exponential)** | $y = a \cdot e^{bx}$ (tự động tối ưu dấu $a > 0$ hoặc $a < 0$) | `--model exp` | Cell 25 |
| **Hàm lũy thừa (Power)** | $y = a \cdot x^b$ | `--model power` | Cell 31 |
| **Đa thức (Polynomial)** | $y = a_0 + a_1 x + a_2 x^2 + \dots + a_k x^k$ | `--model poly --degree k` | Tổng quát |
| **Hypebolic loại 1** | $y = a + \frac{b}{x}$ | `--model hyper1` | Thực nghiệm |
| **Hypebolic loại 2** | $y = \frac{1}{a + bx}$ | `--model hyper2` | Thực nghiệm |
| **Hypebolic loại 3** | $y = \frac{x}{a + bx}$ | `--model hyper3` | Thực nghiệm |
| **Logarit** | $y = a + b \ln(x)$ | `--model log` | Thực nghiệm |
| **Đề thi Midterm K70** | $y = a \cdot \exp(b x^2 + c/x)$ | `--model midterm` | 20251GKtest1 |
| **So sánh toàn bộ** | Chạy toàn bộ mô hình và xếp hạng theo $R^2$ | `--model all` hoặc `--compare` | Toàn diện |

---

## 3. Hướng Dẫn Sử Dụng

### Chạy so sánh toàn bộ các mô hình trên file mặc định 0.csv:
```bash
python lsq.py
# hoặc
python lsq.py 0.csv --compare
```

### Chạy mô hình cụ thể và dự báo giá trị tại $x^* = 2.5$:
```bash
python lsq.py 0.csv --model exp --eval_x 2.5
```

### Chạy khớp đa thức bậc 3:
```bash
python lsq.py 0.csv --model poly --degree 3
```

### Chạy trực tiếp file đề thi Giữa kỳ:
```bash
python lsq.py ../../Midterm/20251GKtest1/20251GKtest1.csv --model midterm
```

### Xuất báo cáo chi tiết ra file text:
```bash
python lsq.py 0.csv --model all --output report.txt
```
