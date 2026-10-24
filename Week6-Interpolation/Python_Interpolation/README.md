# Các phương pháp Nội suy (Interpolation Methods)

Tất cả các phương pháp nội suy được tổ chức theo cấu trúc module độc lập, tích hợp **Universal Smart Parser (`pps_io.py`)** cho phép:
1. **Dán thẳng file đề thi Midterm (`.txt`, `.csv`, `.xlsx`)** mà không cần sửa file hay xóa header.
2. **Tự động nhận diện điểm cần tính $x^*$ (`eval_x`) và $y^*$ (`eval_y`)** trực tiếp từ đề bài.
3. **Auto-Windowing thông minh**: Với các đề thi có hơn 100 điểm mốc, code tự động chọn mốc trung tâm $x_0 \approx x^*$ và cắt cửa sổ 5 hoặc 7 điểm đối xứng quanh $x^*$ (chuẩn theo yêu cầu tính đa thức bậc 4, 6 trong phòng thi).
4. **Công cụ giải nhanh đề thi `solve_midterm.py`**: Chạy và so sánh tất cả các phương pháp chỉ với một dòng lệnh duy nhất.

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

## 🏆 Công cụ giải nhanh đề thi tổng hợp `solve_midterm.py`

Nằm tại thư mục gốc `Python_Interpolation/`, cho phép giải bất kỳ phương pháp nào hoặc so sánh tất cả các phương pháp cùng lúc:

```bash
# 1. Giải phương pháp Stirling trên đề thi 20251:
python solve_midterm.py --file ../../Midterm/20251GKtest1/20251GKtest1.csv --method stirling

# 2. Giải phương pháp Gauss I trên đề thi 20241 với đa thức bậc 4:
python solve_midterm.py --file ../../Midterm/20241GKtest1/20241GKtest1.txt --method gauss1 --degree 4

# 3. So sánh tất cả các phương pháp trên cùng 1 đề thi:
python solve_midterm.py --file ../../Midterm/20251GKtest1/20251GKtest1.csv --compare
```

---

## 🎯 Tham số dòng lệnh nâng cao (CLI Options)

Tất cả các script phương pháp (`<method>.py`) đều hỗ trợ các tham số:
- `-x <val>` hoặc `--eval_x <val>`: Điểm cần tính giá trị nội suy.
- `-d <val>` hoặc `--degree <val>`: Bậc đa thức cần dựng (ví dụ 4 hoặc 6).
- `-n <val>` hoặc `--n_points <val>`: Số mốc dữ liệu cần dùng.
- `--all`: Sử dụng toàn bộ mốc dữ liệu mà không cắt cửa sổ.
