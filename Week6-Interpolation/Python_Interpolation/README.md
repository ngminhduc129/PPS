# Các phương pháp Nội suy (Interpolation Methods)

Tất cả các phương pháp nội suy được tổ chức theo cấu trúc module độc lập, tích hợp **Universal Smart Parser (`pps_io.py`)** cho phép:
1. **Dán thẳng file đề thi Midterm (`.txt`, `.csv`, `.xlsx`)** mà không cần sửa file hay xóa header.
2. **Tự động nhận diện điểm cần tính $x^*$ (`eval_x`) và $y^*$ (`eval_y`)** trực tiếp từ đề bài.
3. **Auto-Windowing thông minh**: Với các đề thi có hơn 100 điểm mốc, code tự động chọn mốc trung tâm $x_0 \approx x^*$ và cắt cửa sổ 5 hoặc 7 điểm đối xứng quanh $x^*$ (chuẩn theo yêu cầu tính đa thức bậc 4, 6 trong phòng thi).

---

## 📁 Danh sách các thư mục phương pháp

| Thư mục | Phương pháp | Script thực thi | File Input | File Output |
| :--- | :--- | :--- | :--- | :--- |
| [`01_Chebyshev/`](01_Chebyshev/) | Mốc Chebyshev | `chebyshev.py` | `input.txt` | `output.txt` |
| [`02_Horner/`](02_Horner/) | Lược đồ Horner | `horner.py` | `input.txt` | `output.txt` |
| [`03_Lagrange/`](03_Lagrange/) | Nội suy Lagrange | `lagrange.py` | `input.txt` | `output.txt` |
| [`04_Newton_Divided/`](04_Newton_Divided/) | Newton tỷ sai phân (mốc bất kỳ) | `newton_divided.py` | `input.txt` | `output.txt` |
| [`05_Newton_Finite/`](05_Newton_Finite/) | Newton sai phân hữu hạn (mốc cách đều) | `newton_finite.py` | `input.txt` | `output.txt` |
| [`06_Gauss1/`](06_Gauss1/) | Gauss I (Gauss thuận) | `gauss1.py` | `input.txt` | `output.txt` |
| [`07_Gauss2/`](07_Gauss2/) | Gauss II (Gauss nghịch) | `gauss2.py` | `input.txt` | `output.txt` |
| [`08_Stirling/`](08_Stirling/) | Stirling | `stirling.py` | `input.txt` | `output.txt` |
| [`09_Bessel/`](09_Bessel/) | Bessel | `bessel.py` | `input.txt` | `output.txt` |
| [`10_Inverse_Function/`](10_Inverse_Function/) | Nội suy Hàm Ngược (Đổi biến $x=P(y)$) | `inverse_function.py` | `input.txt` | `output.txt` |
| [`11_Loop_Inverse/`](11_Loop_Inverse/) | Lặp Nội suy Ngược ($t=\varphi(t)$) | `loop_inverse.py` | `input.txt` | `output.txt` |
| [`12_Spline/`](12_Spline/) | Ghép trơn Spline (Linear, Quad, Cubic) | `spline.py` | `input.txt` | `output.txt` |
| [`13_Least_Squares/`](13_Least_Squares/) | Bình phương tối thiểu (Least Squares) | `lsq.py` | `input.txt` | `output.txt` |

---

## ⚡ Các cách nhập dữ liệu tiện nhất khi đi thi

### Cách 1: Copy-Paste nguyên xi file đề thi vào `input.txt` (Tiện nhất & Khuyên dùng)
Người dùng chỉ cần mở file đề thi (`20241GKtest1.txt` / `.csv` hoặc `20251GKtest1`), nhấn `Ctrl+A` -> `Ctrl+C`, dán thẳng vào `input.txt`:
```text
x_i	y_i	x	4	y	5
1	-0.2067
1.115	0.1728
...
```
- Không cần xóa header `x_i y_i x 4 y 5` hay `x_bar y_bar`.
- Code tự trích xuất `eval_x = 4.0` (và `eval_y = 5.0`).
- Nếu muốn chỉnh số mốc hoặc bậc đa thức, có thể thêm vào đầu file:
  ```text
  # n_points: 5   (hoặc degree: 4)
  ```

### Cách 2: Chạy trực tiếp chỉ định đường dẫn file đề thi
Không cần copy file vào thư mục, chỉ cần truyền đường dẫn:
```bash
python stirling.py ../../../Midterm/20251GKtest1/20251GKtest1.csv
```

### Cách 3: Nhập dữ liệu tự do truyền thống (cho bài tập ngắn)
- **Hàng ngang:**
  ```text
  eval_x = 1.65
  x = 1.40  1.50  1.60  1.70  1.80  1.90  2.00
  y = 0.9523 0.9661 0.9753 0.9838 0.9891 0.9928 0.9000
  ```
- **Rút gọn mốc cách đều:**
  ```text
  eval_x = 1.65
  x0 = 1.40
  h = 0.10
  y = 0.9523 0.9661 0.9753 0.9838 0.9891 0.9928 0.9000
  ```
- **2 cột số dọc:**
  ```text
  eval_x = 1.65
  1.40  0.9523
  1.50  0.9661
  1.60  0.9753
  ```

---

## 🎯 Tham số dòng lệnh nâng cao (CLI Options)

Tất cả các script phương pháp (`<method>.py`) đều hỗ trợ các tham số:
- `-x <val>` hoặc `--eval_x <val>`: Điểm cần tính giá trị nội suy.
- `-d <val>` hoặc `--degree <val>`: Bậc đa thức cần dựng (ví dụ 4 hoặc 6).
- `-n <val>` hoặc `--n_points <val>`: Số mốc dữ liệu cần dùng.
- `--all`: Sử dụng toàn bộ mốc dữ liệu mà không cắt cửa sổ.
- `-t <type>` hoặc `--type <type>`: Chỉ định trường hợp cụ thể (ví dụ `tien`, `lui`, `both`), mặc định luôn là `both` (tự động in tất cả các trường hợp).

### Bình phương tối thiểu (Least Squares)
Chạy từ thư mục `13_Least_Squares/`; mặc định đọc `input.txt`, so sánh các mô hình và ghi báo cáo vào `output.txt`:
```bash
python lsq.py
```
Chọn mô hình hoặc tập dữ liệu khác:
```bash
python lsq.py input.txt --model poly --degree 3 --eval_x 2.5
python lsq.py path/to/data.csv --model exp --eval_x 2.5 --output report.txt
```
Các lựa chọn `--model`: `all`, `x_sinx`, `exp`, `power`, `poly`, `trig`, `hyper1`, `hyper2`, `hyper3`, `midterm`, `log`.

---

## 🌟 Tính năng: Tự động in TẤT CẢ các trường hợp (All Cases In One Run)

Đối với các phương pháp từ **Lagrange** đến **Bessel**, khi chạy mặc định (không cần truyền thêm cờ), chương trình tự động phân tích và in ra toàn bộ mọi khía cạnh:
- **`03_Lagrange`**: In đồng thời **Cửa sổ cục bộ lân cận quanh $x^*$** + **Toàn bộ mốc** (đánh giá hiện tượng Runge) + **Nội suy ngược $x=P(y^*)$** (nếu có $y^*$) + Đạo hàm các cấp.
- **`04_Newton_Divided`**: In đồng thời **Newton Tỷ sai phân TIẾN** (từ $x_0$) + **Newton Tỷ sai phân LÙI** (từ $x_n$) + Bảng so sánh kiểm chứng tính duy nhất của $P(x)$.
- **`05_Newton_Finite`**: In đồng thời **Newton Sai phân hữu hạn TIẾN** ($\Delta^k y_0, t$) + **Newton Sai phân hữu hạn LÙI** ($\nabla^k y_n, t_b$) + Quy đổi sang $P(x)$ + Đạo hàm + So sánh vùng hội tụ tối ưu.
- **`06_Gauss1` & `07_Gauss2`**: In đồng thời **Gauss I (Thuận)** + **Gauss II (Nghịch)** + So sánh tính đồng nhất tuyệt đối.
- **`08_Stirling`**: In kết quả với **Mốc tâm $x_0$ tối ưu** + **Mốc tâm liền kề** + Đối chiếu kiểm chứng công thức lý thuyết $P_{stirling}(x^*) \equiv \frac{P_{gauss1}(x^*) + P_{gauss2}(x^*)}{2}$.
- **`09_Bessel`**: In đồng thời dạng theo **Biến chuẩn hóa $t$** ($0 \le t \le 1$) và dạng theo **Biến đối xứng $u = t - 1/2$** ($-0.5 \le u \le 0.5$) kinh điển trong phòng thi + So sánh kiểm chứng tính đồng nhất 100%.
