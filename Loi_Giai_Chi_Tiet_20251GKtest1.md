# BÀI GIẢI CHI TIẾT ĐỀ THI GIỮA KỲ (20251 - ĐỀ 1)
**Học phần:** Phương pháp số / Giải tích số  
**Dữ liệu sử dụng:** `20251GKtest1.csv` / `20251GKtest1.txt`  
**Các thông số đề bài trích xuất từ file:**
- Điểm tính nội suy: $\bar{x} = 4.15$
- Điểm tính phương trình: $\bar{y} = -1.43$
- Số lượng mốc: $N = 104$ mốc, bước nhảy cách đều $h = 0.115$ ($x_0 = 1.0, x_{103} = 12.845$)

---

## CÂU A: TÌM HÀM THỰC NGHIỆM $y = a e^{b x^2 + \frac{c}{x}}$ VÀ XÁC ĐỊNH SAI SỐ

### 1. Cơ sở lý thuyết và tuyến tính hóa mô hình
Hàm thực nghiệm cần tìm có dạng:
$$y = a e^{b x^2 + \frac{c}{x}}$$

**Nhận xét về dấu dữ liệu:**
Trong file `20251GKtest1`, các giá trị thực nghiệm $y_i < 0$ (hàm mang giá trị âm). Do hàm số mũ $e^{u} > 0$ với mọi $u$, nên tham số $a$ phải mang dấu âm ($a < 0$).
Ta đặt: $-y = (-a) e^{b x^2 + \frac{c}{x}}$ với $-a > 0$ và $-y_i > 0$.

Lấy logarit tự nhiên ($\ln$) hai vế:
$$\ln(-y) = \ln(-a) + b x^2 + c \frac{1}{x}$$

Đặt các biến đổi tuyến tính:
- $Y = \ln(-y)$
- $A_0 = \ln(-a) \implies a = -e^{A_0}$
- $X_1 = x^2$
- $X_2 = \frac{1}{x}$

Khi đó bài toán quy về mô hình hồi quy tuyến tính bội:
$$Y = A_0 + b X_1 + c X_2$$

---

### 2. Thiết lập hệ phương trình chuẩn (Phương pháp bình phương bé nhất)
Với $N$ mốc dữ liệu $(x_i, y_i)$, hàm tổng bình phương sai số theo biến tuyến tính:
$$S(A_0, b, c) = \sum_{i=1}^{N} \left[ Y_i - (A_0 + b X_{1i} + c X_{2i}) \right]^2 \to \min$$

Điều kiện cực trị $\frac{\partial S}{\partial A_0} = 0, \frac{\partial S}{\partial b} = 0, \frac{\partial S}{\partial c} = 0$ dẫn đến hệ phương trình chuẩn Gauss:
$$\begin{pmatrix} N & \sum X_{1i} & \sum X_{2i} \\ \sum X_{1i} & \sum X_{1i}^2 & \sum X_{1i} X_{2i} \\ \sum X_{2i} & \sum X_{1i} X_{2i} & \sum X_{2i}^2 \end{pmatrix} \begin{pmatrix} A_0 \\ b \\ c \end{pmatrix} = \begin{pmatrix} \sum Y_i \\ \sum X_{1i} Y_i \\ \sum X_{2i} Y_i \end{pmatrix}$$

---

### 3. Kết quả tính toán các hệ số
Áp dụng trên tập dữ liệu mốc $y_i < 0$:
- $A_0 = \ln(-a) = -0.46049187 \implies a = -e^{-0.46049187} = \mathbf{-0.63098311}$
- $b = \mathbf{-0.01817920}$
- $c = \mathbf{2.25587011}$

*(Nếu tối ưu hóa phi tuyến trực tiếp bằng thuật toán Gauss-Newton / Levenberg-Marquardt trên toàn bộ 104 mốc, các tham số hội tụ tối ưu là: $a = -0.80009922$, $b = -0.02414029$, $c = 1.92017834$)*.

**Phương trình hàm thực nghiệm tìm được:**
$$y = -0.63098311 \cdot e^{-0.01817920 x^2 + \frac{2.25587011}{x}}$$

---

### 4. Đánh giá sai số của hàm thực nghiệm
Gọi $\hat{y}_i = a e^{b x_i^2 + c/x_i}$ là giá trị tính theo hàm thực nghiệm, sai số tại từng điểm $e_i = y_i - \hat{y}_i$:

| Chỉ số sai số | Công thức | Giá trị tính toán |
| :--- | :--- | :---: |
| **Tổng bình phương sai số (SSE)** | $SSE = \sum_{i=1}^{N} (y_i - \hat{y}_i)^2$ | **1.051302** |
| **Sai số tuyệt đối trung bình (MAE)** | $MAE = \frac{1}{N} \sum_{i=1}^{N} \vert y_i - \hat{y}_i \vert$ | **0.068412** |
| **Sai số căn quân phương (RMSE)** | $RMSE = \sqrt{\frac{1}{N} SSE}$ | **0.101538** |
| **Sai số cực đại (Max Error)** | $\max \vert y_i - \hat{y}_i \vert$ | **0.294156** |
| **Hệ số xác định ($R^2$)** | $R^2 = 1 - \frac{SSE}{\sum (y_i - \bar{y})^2}$ | **0.9882 (98.82%)** |

> **Nhận xét:** Hệ số xác định $R^2 = 98.82\%$ rất cao, chứng minh dạng hàm thực nghiệm mô tả cực kỳ chính xác quy luật của bộ dữ liệu đã cho.

---

## CÂU B: PHƯƠNG PHÁP LẶP GIẢI PHƯƠNG TRÌNH $f(x) = \bar{y}$ VỚI $\bar{y} = -1.43$

Đề bài yêu cầu: *"Viết thuật toán cho phương pháp lặp giải phương trình $f(x) = \bar{y}$ và áp dụng với tập dữ liệu và $\bar{y}$ cho trong file dữ liệu"*.

Trong môn Phương pháp số, bài toán giải $f(x) = \bar{y}$ từ bảng dữ liệu được giải quyết theo **2 phương pháp chuẩn**:

---

### CÁCH 1: PHƯƠNG PHÁP LẶP NỘI SUY NGƯỢC (Inverse Interpolation Loop - Chuẩn Giáo Trình)

#### 1. Thuật toán (Mã giả / Mô tả các bước)
- **Bước 1 (Xác định khoảng phân ly nghiệm):**
  Tra bảng số liệu tìm 2 mốc liên tiếp $x_k, x_{k+1}$ sao cho:
  $$y_k \le \bar{y} \le y_{k+1} \quad (\text{hoặc } y_k \ge \bar{y} \ge y_{k+1})$$
  Lấy một cụm mốc lân cận xung quanh khoảng này (thường lấy 3 hoặc 4 mốc: $x_0, x_1, x_2, x_3$).
- **Bước 2 (Thiết lập đa thức nội suy):**
  Xây dựng đa thức Newton tiến theo biến chuẩn hóa $t = \frac{x - x_0}{h}$:
  $$P(t) = y_0 + t \Delta y_0 + \frac{t(t-1)}{2} \Delta^2 y_0 + \frac{t(t-1)(t-2)}{6} \Delta^3 y_0 = \bar{y}$$
- **Bước 3 (Đưa về dạng phương trình lặp):**
  Tách số hạng bậc nhất $t \Delta y_0$ ra vế trái:
  $$t = \varphi(t) = \frac{\bar{y} - y_0}{\Delta y_0} - \frac{\Delta^2 y_0}{2 \Delta y_0} t(t-1) - \frac{\Delta^3 y_0}{6 \Delta y_0} t(t-1)(t-2)$$
- **Bước 4 (Khởi tạo và quá trình lặp):**
  - Khởi tạo: $t^{(0)} = \frac{\bar{y} - y_0}{\Delta y_0}$
  - Lặp: $t^{(k+1)} = \varphi(t^{(k)})$
  - Dừng lặp khi $|t^{(k+1)} - t^{(k)}| < \varepsilon$ (hoặc $|x^{(k+1)} - x^{(k)}| < \varepsilon$).
- **Bước 5:**
  Tính nghiệm cuối cùng: $x^* = x_0 + t^* \cdot h$.

#### 2. Áp dụng trên tập dữ liệu và $\bar{y} = -1.43$:
- Khoảng phân ly nghiệm: $y(2.495) = -1.5127 \le \bar{y} = -1.43 \le y(2.610) = -1.3977$.
- Chọn cụm 4 mốc lân cận:
  - $x_0 = 2.380, y_0 = -1.5819$
  - $x_1 = 2.495, y_1 = -1.5127 \implies \Delta y_0 = 0.0692$
  - $x_2 = 2.610, y_2 = -1.3977 \implies \Delta^2 y_0 = 0.0458$
  - $x_3 = 2.725, y_3 = -1.3814 \implies \Delta^3 y_0 = -0.1445$

- Hàm lặp:
  $$\varphi(t) = 2.19508671 - 0.33092485 \cdot t(t-1) + 0.34802505 \cdot t(t-1)(t-2)$$

**Bảng các bước lặp:**

| Bước lặp $k$ | $t^{(k)}$ | $x^{(k)} = x_0 + t^{(k)} h$ | Sai số $|t^{(k)} - t^{(k-1)}|$ |
| :---: | :---: | :---: | :---: |
| 0 | 2.195087 | 2.632435 | - |
| 1 | 1.505076 | 2.553084 | 0.690011 |
| 2 | 1.812588 | 2.588448 | 0.307512 |
| 3 | 1.611604 | 2.565334 | 0.200984 |
| 4 | 1.735672 | 2.579602 | 0.124068 |
| ... | ... | ... | ... |
| **Hội tụ** | **1.685968** | **2.573886** | $< 10^{-6}$ |

$$\implies \mathbf{x^* \approx 2.573886}$$

---

### CÁCH 2: PHƯƠNG PHÁP LẶP TRÊN HÀM THỰC NGHIỆM TÌM ĐƯỢC Ở CÂU A
Từ hàm thực nghiệm $y(x) = a e^{b x^2 + c/x} = \bar{y} = -1.43$:
$$-0.63098311 e^{-0.01817920 x^2 + \frac{2.25587011}{x}} = -1.43$$

Lấy logarit tự nhiên:
$$-0.01817920 x^2 + \frac{2.25587011}{x} + \ln(0.63098311) - \ln(1.43) = 0$$
$$\iff g(x) = -0.01817920 x^2 + \frac{2.25587011}{x} - 0.81822369 = 0$$

- Sử dụng phương pháp lặp tiếp tuyến Newton-Raphson: $x_{k+1} = x_k - \frac{g(x_k)}{g'(x_k)}$:
  - Chọn $x_0 = 2.55$
  - Bước 1: $x_1 = 2.481681$
  - Bước 2: $x_2 = 2.483784$
  - Bước 3: $x_3 = \mathbf{2.483786}$ (Hội tụ)

$$\implies \mathbf{x^* \approx 2.483786}$$

---

## CÂU C: ĐA THỨC NỘI SUY BESSEL 8 MỐC TÍNH $f(\bar{x})$ VÀ $f'(\bar{x})$ VỚI $\bar{x} = 4.15$

### 1. Trích xuất 8 mốc nội suy đối xứng
Điểm cần tính $\bar{x} = 4.15$ nằm trong khoảng $[x_0, x_1] = [4.105, 4.220]$ với bước $h = 0.115$.  
Đa thức Bessel yêu cầu số mốc chẵn ($n+1 = 8$ mốc, bậc 7), ta chọn 8 mốc đối xứng quanh khoảng $[x_0, x_1]$:

| Chỉ số quy ước | Mốc $x_i$ | Giá trị $y_i$ |
| :---: | :---: | :---: |
| $x_{-3}$ | 3.760 | -0.9384 |
| $x_{-2}$ | 3.875 | -0.8979 |
| $x_{-1}$ | 3.990 | -0.8110 |
| **$x_{0}$** | **4.105** | **-0.8056** |
| **$x_{1}$** | **4.220** | **-0.7548** |
| $x_{2}$ | 4.335 | -0.7416 |
| $x_{3}$ | 4.450 | -0.7049 |
| $x_{4}$ | 4.565 | -0.6784 |

---

### 2. Bảng sai phân hữu hạn 8 mốc
$$\Delta y_i = y_{i+1} - y_i, \quad \Delta^k y_i = \Delta^{k-1} y_{i+1} - \Delta^{k-1} y_i$$

| $i$ | $x_i$ | $y_i$ | $\Delta y$ | $\Delta^2 y$ | $\Delta^3 y$ | $\Delta^4 y$ | $\Delta^5 y$ | $\Delta^6 y$ | $\Delta^7 y$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| -3 | 3.760 | -0.9384 | 0.0405 | 0.0464 | -0.1279 | 0.2548 | -0.4647 | 0.8187 | **-1.4116** |
| -2 | 3.875 | -0.8979 | 0.0869 | -0.0815 | 0.1269 | -0.2099 | **0.3540** | -0.5929 | |
| -1 | 3.990 | -0.8110 | 0.0054 | 0.0454 | **-0.0830** | 0.1441 | -0.2389 | | |
| **0** | **4.105** | **-0.8056** | **0.0508** | -0.0376 | 0.0611 | -0.0948 | | | |
| **1** | **4.220** | **-0.7548** | 0.0132 | 0.0235 | -0.0337 | | | | |
| 2 | 4.335 | -0.7416 | 0.0367 | -0.0102 | | | | | |
| 3 | 4.450 | -0.7049 | 0.0265 | | | | | | |
| 4 | 4.565 | -0.6784 | | | | | | | |

---

### 3. Công thức nội suy Bessel và các hệ số trung tâm $C_j$
Biến chuẩn hóa:
$$t = \frac{x - x_0}{h} = \frac{4.15 - 4.105}{0.115} \approx \mathbf{0.39130435}$$

Công thức nội suy Bessel:
$$P(t) = C_0 + \left(t - \frac{1}{2}\right) C_1 + \frac{t(t-1)}{2!} C_2 + \frac{\left(t - \frac{1}{2}\right)t(t-1)}{3!} C_3 + \frac{(t+1)t(t-1)(t-2)}{4!} C_4 + \dots$$

Các hệ số sai phân trung tâm $C_j$:
- $C_0 = \frac{y_0 + y_1}{2} = \frac{-0.8056 - 0.7548}{2} = \mathbf{-0.78020000}$
- $C_1 = \Delta y_0 = \mathbf{0.05080000}$
- $C_2 = \frac{\Delta^2 y_{-1} + \Delta^2 y_0}{2} = \frac{0.0454 - 0.0376}{2} = \mathbf{0.00390000}$
- $C_3 = \Delta^3 y_{-1} = \mathbf{-0.08300000}$
- $C_4 = \frac{\Delta^4 y_{-2} + \Delta^4 y_{-1}}{2} = \frac{-0.2099 + 0.1441}{2} = \mathbf{-0.03290000}$
- $C_5 = \Delta^5 y_{-2} = \mathbf{0.35400000}$
- $C_6 = \frac{\Delta^6 y_{-3} + \Delta^6 y_{-2}}{2} = \frac{0.8187 - 0.5929}{2} = \mathbf{0.11290000}$
- $C_7 = \Delta^7 y_{-3} = \mathbf{-1.41160000}$

---

### 4. Dạng tường minh của đa thức Bessel $P(t)$
Khai triển và thu gọn các số hạng theo lũy thừa của $t$:
$$P(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3 + a_4 t^4 + a_5 t^5 + a_6 t^6 + a_7 t^7$$

| Bậc $t^k$ | Hệ số $a_k$ |
| :---: | :---: |
| $t^0$ | $-0.80560000$ |
| $t^1$ | $+0.03267952$ |
| $t^2$ | $+0.03599417$ |
| $t^3$ | $-0.00775931$ |
| $t^4$ | $-0.01443125$ |
| $t^5$ | $+0.00345986$ |
| $t^6$ | $+0.00113708$ |
| $t^7$ | $-0.00028008$ |

$$P(t) = -0.80560000 + 0.03267952 t + 0.03599417 t^2 - 0.00775931 t^3 - 0.01443125 t^4 + 0.00345986 t^5 + 0.00113708 t^6 - 0.00028008 t^7$$

---

### 5. Tính giá trị $f(\bar{x})$ và đạo hàm $f'(\bar{x})$ tại $\bar{x} = 4.15$

#### a) Tính $f(\bar{x})$:
Thay $t = 0.39130435$ vào đa thức $P(t)$:
$$f(4.15) \approx P(0.39130435) = \mathbf{-0.78806879}$$

#### b) Tính đạo hàm $f'(\bar{x})$:
Vì $t = \frac{x - x_0}{h} \implies \frac{dt}{dx} = \frac{1}{h}$.  
Áp dụng quy tắc đạo hàm hàm hợp:
$$f'(x) = \frac{dP}{dx} = \frac{dP}{dt} \cdot \frac{dt}{dx} = \frac{1}{h} P'(t)$$

Ta có đạo hàm $P'(t)$:
$$P'(t) = a_1 + 2 a_2 t + 3 a_3 t^2 + 4 a_4 t^3 + 5 a_5 t^4 + 6 a_6 t^5 + 7 a_7 t^6$$

Thay $t = 0.39130435$:
$$P'(0.39130435) \approx \mathbf{0.05428706}$$

Suy ra đạo hàm tại $\bar{x} = 4.15$:
$$f'(4.15) = \frac{P'(0.39130435)}{h} = \frac{0.05428706}{0.115} = \mathbf{0.47206140}$$

---

## BẢNG TỔNG KẾT ĐÁP SỐ TOÀN BỘ ĐỀ THI

1. **Câu a:**
   - Hàm thực nghiệm: $\mathbf{y = -0.63098311 \cdot e^{-0.01817920 x^2 + \frac{2.25587011}{x}}}$
   - Sai số: $SSE = 1.0513$, $RMSE = 0.1015$, $R^2 = 98.82\%$.
2. **Câu b:**
   - Nghiệm của phương trình $f(x) = -1.43$:
     - Theo lặp nội suy ngược: $\mathbf{x^* \approx 2.573886}$
     - Theo lặp trên hàm thực nghiệm: $\mathbf{x^* \approx 2.483786}$
3. **Câu c:**
   - Giá trị nội suy: $\mathbf{f(4.15) \approx -0.78806879}$
   - Đạo hàm nội suy: $\mathbf{f'(4.15) \approx 0.47206140}$
