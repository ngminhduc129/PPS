# TỔNG HỢP THUẬT TOÁN PHƯƠNG PHÁP SỐ (TỪ HORNER ĐẾN LEAST SQUARES)

---

## 1. Sơ đồ Horner (Tính giá trị đa thức và đạo hàm cấp $k$ tại $x = c$)
**Bài toán:** Cho đa thức $P(x) = a_0 + a_1 x + \dots + a_N x^N$ và điểm $x = c$. Tính $P(c)$ và đạo hàm cấp $k$ là $P^{(k)}(c)$.

**Các bước thực hiện:**
1. **Bước 1 (Khởi tạo mảng hệ số):** Gán mảng $b$ chứa các hệ số của đa thức theo thứ tự từ bậc cao nhất đến bậc tự do: $b_i = a_{N-i}$ với $i = 0, 1, \dots, N$ (tức là $b = [a_N, a_{N-1}, \dots, a_0]$).
2. **Bước 2 (Lặp chia tổng hợp $k+1$ lần):** Cho biến đếm lần chia $m$ chạy từ $0$ đến $k$:
   * Giữ nguyên phần tử đầu tiên: $b_0 = b_0$.
   * Cho $i$ chạy từ $1$ đến $N - m$: Cập nhật $b_i = b_{i-1} \cdot c + b_i$.
   * Lưu kết quả phần dư cuối mỗi lần chia: Hệ số dư tại vị trí $b_{N-m}$ chính là giá trị $\frac{P^{(m)}(c)}{m!}$.
3. **Bước 3 (Tính đạo hàm cấp $k$):** 
   * Nếu $m = 0$: Giá trị đa thức là $P(c) = b_N$.
   * Nếu cần đạo hàm cấp $k$ ($k \le N$): Tính $P^{(k)}(c) = k! \cdot b_{N-k}$. (Nếu $k > N$ thì $P^{(k)}(c) = 0$).

---

## 2. Đa thức nội suy Lagrange (Mốc bất kỳ)
**Bài toán:** Cho bộ $(n+1)$ điểm $(x_i, y_i)$ với $i = 0, \dots, n$ ($x_i$ đôi một khác nhau). Tìm đa thức nội suy $L_n(x)$ bậc không quá $n$.

**Các bước thực hiện:**
1. **Bước 1 (Nhập dữ liệu):** Nhập mảng mốc $x = [x_0, x_1, \dots, x_n]$ và mảng giá trị $y = [y_0, y_1, \dots, y_n]$.
2. **Bước 2 (Xây dựng tích cơ sở hoặc đa thức Lagrange cơ bản $l_i(x)$):** Với mỗi $i = 0, \dots, n$:
   * Tính mẫu số: $D_i = \prod_{j=0, j \neq i}^{n} (x_i - x_j)$.
   * Xây dựng tử số là đa thức tích: $w_i(x) = \prod_{j=0, j \neq i}^{n} (x - x_j)$ (sử dụng phép nhân đa thức / bảng tính tích qua sơ đồ Horner ngược).
   * Đa thức cơ bản thứ $i$: $l_i(x) = \frac{w_i(x)}{D_i}$.
3. **Bước 3 (Tổng hợp đa thức nội suy):** 
   * Tính tổ hợp tuyến tính các hệ số tương ứng: $L_n(x) = \sum_{i=0}^{n} y_i \cdot l_i(x)$.
   * Thu gọn hệ số theo bậc từ cao đến thấp (hoặc từ thấp đến cao).

---

## 3. Đa thức nội suy Newton với mốc bất kỳ
**Bài toán:** Cho bộ $(n+1)$ điểm $(x_i, y_i)$ với $i = 0, \dots, n$. Xây dựng đa thức nội suy Newton tiến (hoặc lùi).

**Các bước thực hiện:**
1. **Bước 1 (Nhập dữ liệu):** Nhập mảng $x_i, y_i$ ($i = 0, \dots, n$). Khởi tạo ma trận bảng tỷ sai phân $F$ kích thước $(n+1) \times (n+1)$, gán cột đầu tiên $F_{i,0} = y_i$.
2. **Bước 2 (Lập bảng tỷ sai phân):** 
   * Cho cột cấp tỷ sai phân $j$ chạy từ $1$ đến $n$:
   * Cho hàng $i$ chạy từ $j$ đến $n$: Tính $F_{i,j} = \frac{F_{i, j-1} - F_{i-1, j-1}}{x_i - x_{i-j}}$.
3. **Bước 3 (Trích xuất hệ số và lập bảng tính tích):**
   * Với công thức **Newton tiến** xuất phát từ $x_0$: Lấy đường chéo trên $c_j = F_{j,j}$ (tức $f[x_0, \dots, x_j]$) với $j = 0, \dots, n$.
   * Khởi tạo đa thức $P(x) = c_0$ và đa thức tích cơ sở $\omega_0(x) = 1$.
4. **Bước 4 (Cộng dồn đa thức):** Cho $j$ chạy từ $1$ đến $n$:
   * Nhân thêm nhân tử mới vào bảng tính tích: $\omega_j(x) = \omega_{j-1}(x) \cdot (x - x_{j-1})$.
   * Cập nhật đa thức: $P(x) = P(x) + c_j \cdot \omega_j(x)$.

---

## 4. Đa thức nội suy Newton với mốc cách đều (Tiến & Lùi)
**Bài toán:** Cho bộ điểm $(x_k, y_k)$ cách đều với bước nhảy $h = x_{k+1} - x_k = \text{const}$ ($k = 0, \dots, n$). Tính giá trị nội suy tại $\bar{x}$.

**Các bước thực hiện:**
1. **Bước 1 (Kiểm tra điều kiện và lập bảng sai phân):**
   * Xác định bước nhảy $h = x_1 - x_0$.
   * Khởi tạo cột 0 của bảng sai phân $\Delta^0 y_i = y_i$ ($i = 0, \dots, n$).
   * Cho cấp sai phân $j$ chạy từ $1$ đến $n$, tính $\Delta^j y_i = \Delta^{j-1} y_{i+1} - \Delta^{j-1} y_i$ (lưu ý: $\Delta y_k = \nabla y_{k+1}$).
2. **Bước 2 (Đổi biến sang biến phụ $t$):**
   * **Nếu $\bar{x}$ gần đầu bảng (dùng Newton tiến từ $x_0$):** Đặt $t = \frac{\bar{x} - x_0}{h}$. Các sai phân sử dụng là $\Delta^k y_0$.
   * **Nếu $\bar{x}$ gần cuối bảng (dùng Newton lùi từ $x_n$):** Đặt $t = \frac{\bar{x} - x_n}{h}$. Các sai phân sử dụng là $\nabla^k y_n$ (tức $\Delta^k y_{n-k}$).
3. **Bước 3 (Xây dựng đa thức theo $t$):**
   * **Newton tiến:** $P_n(t) = y_0 + \frac{\Delta y_0}{1!} t + \frac{\Delta^2 y_0}{2!} t(t-1) + \dots + \frac{\Delta^n y_0}{n!} t(t-1)\dots(t-n+1)$.
   * **Newton lùi:** $P_n(t) = y_n + \frac{\nabla y_n}{1!} t + \frac{\nabla^2 y_n}{2!} t(t+1) + \dots + \frac{\nabla^n y_n}{n!} t(t+1)\dots(t+n-1)$.
4. **Bước 4 (Tính giá trị):** Rút gọn đa thức theo lũy thừa của $t$ qua bảng tính tích, sau đó dùng sơ đồ Horner tính $P_n(t)$ tại giá trị $t$ vừa tìm.

---

## 5. Đa thức nội suy trung tâm: Gauss I, Gauss II, Stirling và Bessel
**Bài toán:** Cho bảng mốc cách đều bước $h$, cần tính gần đúng $f(\bar{x})$ và $f'(\bar{x})$ khi $\bar{x}$ nằm ở khoảng giữa bảng dữ liệu.

### 5.1. Công thức Gauss I, Gauss II và Stirling (Sử dụng $2N+1$ mốc nội suy lẻ)
1. **Bước 1 (Chọn mốc trung tâm và trích xuất điểm):**
   * Tìm mốc $x_0$ trong bảng gần $\bar{x}$ nhất.
   * Trích xuất $2N+1$ mốc đối xứng quanh $x_0$: $\{x_{-N}, \dots, x_{-1}, x_0, x_1, \dots, x_N\}$ với điều kiện bảng dữ liệu đủ điểm hai phía.
   * Đổi biến: $t = \frac{\bar{x} - x_0}{h}$ (thường $|t| \le 0.5$, chọn Stirling tốt nhất khi $|t| \le 0.25$).
2. **Bước 2 (Lập bảng sai phân trung tâm):** Tính các cấp sai phân đến cấp $2N$ cho bộ $2N+1$ điểm vừa trích xuất.
3. **Bước 3 (Xác định hệ số và xây dựng đa thức):**
   * **Với Gauss I:** Đi theo đường ziczac xuống trước: $y_0 \to \Delta y_0 \to \Delta^2 y_{-1} \to \Delta^3 y_{-1} \to \Delta^4 y_{-2} \dots$
     $$P(t) = y_0 + t\Delta y_0 + \frac{\Delta^2 y_{-1}}{2!}t(t-1) + \frac{\Delta^3 y_{-1}}{3!}(t+1)t(t-1) + \frac{\Delta^4 y_{-2}}{4!}(t+1)t(t-1)(t-2) + \dots$$
   * **Với Gauss II:** Đi theo đường ziczac lên trước: $y_0 \to \Delta y_{-1} \to \Delta^2 y_{-1} \to \Delta^3 y_{-2} \to \Delta^4 y_{-2} \dots$
   * **Với Stirling (Trung bình cộng của Gauss I và Gauss II):** Tách riêng cụm bậc chẵn và bậc lẻ:
     * Hệ số chẵn: $a_{2i} = \frac{\Delta^{2i} y_{-i}}{(2i)!}$ đi kèm tích $t^2(t^2 - 1^2)\dots(t^2 - (i-1)^2)$.
     * Hệ số lẻ: $a_{2i-1} = \frac{\Delta^{2i-1} y_{-i} + \Delta^{2i-1} y_{-(i-1)}}{2 \cdot (2i-1)!}$ đi kèm tích $t(t^2 - 1^2)\dots(t^2 - (i-1)^2)$.
4. **Bước 4 (Tính $f(\bar{x})$ và $f'(\bar{x})$):**
   * Thu gọn đa thức $P(t)$ bằng bảng tính tích theo biến $t^2$.
   * Tính $f(\bar{x}) \approx P(t)$ và $f'(\bar{x}) \approx \frac{1}{h} P'(t)$ bằng sơ đồ Horner.

### 5.2. Công thức nội suy Bessel (Sử dụng $2N$ mốc nội suy chẵn)
1. **Bước 1 (Trích xuất mốc):** Tìm khoảng $[x_0, x_1]$ chứa $\bar{x}$ (thích hợp nhất khi $t = \frac{\bar{x}-x_0}{h} \in [0.25, 0.75]$). Trích xuất $2N$ mốc đối xứng qua khoảng này: $\{x_{-(N-1)}, \dots, x_0, x_1, \dots, x_N\}$.
2. **Bước 2 (Đổi biến phụ $u$):** Đặt $u = \frac{\bar{x} - x_0}{h} - \frac{1}{2} = t - \frac{1}{2}$.
3. **Bước 3 (Lập bảng sai phân và tách chẵn, lẻ):**
   * Tính bảng sai phân đến cấp $2N-1$.
   * Hệ số chẵn: Trung bình cộng hai sai phân cấp chẵn đối xứng $\frac{\Delta^{2i} y_{-i} + \Delta^{2i} y_{-(i-1)}}{2 \cdot (2i)!}$ nhân với tích $(u^2 - \frac{1}{4})(u^2 - \frac{9}{4})\dots(u^2 - \frac{(2i-1)^2}{4})$.
   * Hệ số lẻ: Sai phân cấp lẻ ở giữa $\frac{\Delta^{2i+1} y_{-i}}{(2i+1)!}$ nhân với $u(u^2 - \frac{1}{4})\dots(u^2 - \frac{(2i-1)^2}{4})$.
4. **Bước 4 (Tính giá trị và đạo hàm):** Lập đa thức theo $u$, dùng Horner tính $f(\bar{x}) \approx P(u)$ và $f'(\bar{x}) \approx \frac{1}{h} P'(u)$.

---

## 6. Nội suy ngược (Phương pháp lặp giải phương trình $f(x) = \bar{y}$)
**Bài toán:** Cho bảng dữ liệu $(x_i, y_i)$ có mốc cách đều bước $h$ và giá trị $\bar{y}$. Tìm $x$ sao cho $f(x) = \bar{y}$.

**Các bước thực hiện:**
1. **Bước 1 (Tìm khoảng cách ly nghiệm):** Quét bảng dữ liệu để tìm chỉ số $k$ sao cho $(y_k - \bar{y})(y_{k+1} - \bar{y}) < 0$. Khoảng cách ly nghiệm là $(x_k, x_{k+1})$ và $\Delta y_k = y_{k+1} - y_k \neq 0$. Trích xuất bộ điểm phù hợp quanh khoảng này và lập bảng sai phân.
2. **Bước 2 (Xây dựng hàm lặp theo vị trí khoảng cách ly):**
   * **Giải pháp 1 (Dùng Newton tiến xuất phát từ $x_k$ - khi khoảng cách ly ở đầu/giữa bảng):**
     * Từ $\bar{y} = y_k + \Delta y_k t + \frac{\Delta^2 y_k}{2!} t(t-1) + \dots$, rút $t$ từ số hạng bậc nhất:
       $$t = \varphi(t) = \frac{\bar{y} - y_k}{\Delta y_k} - \frac{1}{\Delta y_k} \left[ \frac{\Delta^2 y_k}{2!} t(t-1) + \frac{\Delta^3 y_k}{3!} t(t-1)(t-2) + \dots \right]$$
     * Khởi tạo xấp xỉ đầu: $t_0 = \frac{\bar{y} - y_k}{\Delta y_k}$.
   * **Giải pháp 2 (Dùng Newton lùi xuất phát từ $x_{k+1}$ - khi khoảng cách ly ở cuối bảng):**
     * Từ đa thức Newton lùi tại $x_{k+1}$, rút $t$:
       $$t = \psi(t) = \frac{\bar{y} - y_{k+1}}{\nabla y_{k+1}} - \frac{1}{\nabla y_{k+1}} \left[ \frac{\nabla^2 y_{k+1}}{2!} t(t+1) + \dots \right]$$
     * Khởi tạo xấp xỉ đầu: $t_0 = \frac{\bar{y} - y_{k+1}}{\nabla y_{k+1}}$.
3. **Bước 3 (Thực hiện vòng lặp đơn):**
   * Tính $t_{j+1} = \varphi(t_j)$ (hoặc $\psi(t_j)$) với $j = 0, 1, 2, \dots$
   * Kiểm tra điều kiện dừng: $|t_{j+1} - t_j| < \varepsilon$ (hoặc hết số lần lặp tối đa).
4. **Bước 4 (Quy đổi về nghiệm $\bar{x}$):**
   * Nếu dùng Newton tiến: $\bar{x} = x_k + t_{j+1} \cdot h$.
   * Nếu dùng Newton lùi: $\bar{x} = x_{k+1} + t_{j+1} \cdot h$.

---

## 7. Hàm ghép trơn (Spline tuyến tính, Spline cấp 2, Spline cấp 3)
**Bài toán:** Cho bộ $n$ điểm $(x_i, y_i)_{i=1}^n$ sắp xếp tăng dần $x_1 < x_2 < \dots < x_n$. Đặt $h_k = x_{k+1} - x_k$ với $k = 1, \dots, n-1$. Tìm hàm ghép $S(x) = S_k(x)$ trên từng đoạn $[x_k, x_{k+1}]$ có bậc $p$ và trơn đến cấp $p-1$.

### 7.1. Spline tuyến tính (Cấp 1)
1. **Bước 1:** Tính độ dài các đoạn $h_k = x_{k+1} - x_k$ với $k = 1, \dots, n-1$.
2. **Bước 2:** Trên mỗi đoạn $[x_k, x_{k+1}]$, tính hệ số của $S_k(x) = a_k x + b_k$:
   * $a_k = \frac{y_{k+1} - y_k}{h_k}$
   * $b_k = \frac{y_k x_{k+1} - y_{k+1} x_k}{h_k}$

### 7.2. Spline cấp 2 (Quadratic Spline)
1. **Bước 1:** Đặt ẩn phụ đạo hàm cấp 1 tại các nút $m_k = S'_k(x_k)$. Tính $h_k = x_{k+1} - x_k$ và vế phải $\gamma_k = 2 \frac{y_{k+1} - y_k}{h_k}$ với $k = 1, \dots, n-1$.
2. **Bước 2:** Bổ sung 1 điều kiện biên (ví dụ cho trước $m_1 = y'_1$ hoặc $m_1 = 0$) để giải hệ $(n-1)$ phương trình tuyến tính tìm $m_2, \dots, m_n$:
   $$m_k + m_{k+1} = \gamma_k \quad (k = 1, \dots, n-1)$$
3. **Bước 3:** Tính hệ số đa thức $S_k(x) = a_k x^2 + b_k x + c_k$ trên từng đoạn $[x_k, x_{k+1}]$:
   * $a_k = \frac{m_{k+1} - m_k}{2h_k}$
   * $b_k = \frac{m_k x_{k+1} - m_{k+1} x_k}{h_k}$
   * $c_k = \frac{-m_k x_{k+1}^2 + m_{k+1} x_k^2}{2h_k} + y_k + \frac{m_k h_k}{2}$

### 7.3. Spline cấp 3 (Cubic Spline)
1. **Bước 1 (Tính bước nhảy và vế phải):**
   * Tính $h_k = x_{k+1} - x_k$ với $k = 1, \dots, n-1$.
   * Tính $\gamma_k = \frac{y_{k+2} - y_{k+1}}{h_{k+1}} - \frac{y_{k+1} - y_k}{h_k}$ với $k = 1, \dots, n-2$.
2. **Bước 2 (Lập và giải hệ phương trình 3 đường chéo tìm $\alpha_k = S''(x_k)$):**
   * Hệ $(n-2)$ phương trình với $n$ ẩn $(\alpha_1, \dots, \alpha_n)$:
     $$\frac{h_k}{6}\alpha_k + \frac{h_k + h_{k+1}}{3}\alpha_{k+1} + \frac{h_{k+1}}{6}\alpha_{k+2} = \gamma_k \quad (k = 1, \dots, n-2)$$
   * Áp dụng 2 điều kiện biên (ví dụ Spline tự nhiên: $\alpha_1 = 0, \alpha_n = 0$) để đưa về hệ $n$ phương trình $n$ ẩn và giải tìm toàn bộ vector $\alpha$.
3. **Bước 3 (Tính hệ số $S_k(x) = a_k x^3 + b_k x^2 + c_k x + d_k$ với $k = 1, \dots, n-1$):**
   * $a_k = \frac{\alpha_{k+1} - \alpha_k}{6h_k}$
   * $b_k = \frac{3\alpha_k x_{k+1} - 3\alpha_{k+1} x_k}{6h_k}$
   * $c_k = \frac{-3\alpha_k x_{k+1}^2 + 3\alpha_{k+1} x_k^2}{6h_k} + \frac{y_{k+1} - y_k}{h_k} - \frac{(\alpha_{k+1} - \alpha_k)h_k}{6}$
   * $d_k = \frac{\alpha_k x_{k+1}^3 - \alpha_{k+1} x_k^3}{6h_k} + \frac{y_k x_{k+1} - y_{k+1} x_k}{h_k} + \frac{(\alpha_{k+1} x_k - \alpha_k x_{k+1})h_k}{6}$

---

## 8. Phương pháp Bình phương tối thiểu (Least Squares)
**Bài toán:** Cho bộ dữ liệu thực nghiệm $(x_i, y_i)_{i=1}^n$. Tìm hàm xấp xỉ $g(x)$ sao cho sai số toàn phương trung bình $\sigma = \sqrt{\frac{1}{n}\sum_{i=1}^n (y_i - g(x_i))^2}$ (hoặc tổng bình phương sai lệch $S = \sum_{i=1}^n (y_i - g(x_i))^2$) đạt cực tiểu.

### 8.1. Mô hình tuyến tính theo tham số: $g(x) = c_1 \varphi_1(x) + c_2 \varphi_2(x) + \dots + c_m \varphi_m(x)$
1. **Bước 1 (Nhập dữ liệu và chọn hệ hàm cơ sở):**
   * Nhập mảng điểm $(x_i, y_i)$ với $i = 1, \dots, n$ (lập vector cột $y = [y_1, y_2, \dots, y_n]^T$).
   * Xác định các hàm cơ sở $\varphi_1(x), \varphi_2(x), \dots, \varphi_m(x)$ (ví dụ: hàm đa thức có $\varphi_1(x)=1, \varphi_2(x)=x, \dots$; hoặc $g(x) = ax + b\sin x$ có $\varphi_1(x)=x, \varphi_2(x)=\sin x$).
2. **Bước 2 (Xây dựng ma trận thiết kế $\Phi$ và hệ phương trình chuẩn):**
   * Tính ma trận $\Phi$ kích thước $n \times m$, trong đó phần tử hàng $i$ cột $j$ là $\Phi_{i,j} = \varphi_j(x_i)$, tức $\Phi = [\varphi_1 \quad \varphi_2 \quad \dots \quad \varphi_m]$.
   * Tính ma trận vuông đối xứng $M = \Phi^T \Phi$ (kích thước $m \times m$).
   * Tính vector vế phải $b = \Phi^T y$ (kích thước $m \times 1$).
3. **Bước 3 (Giải hệ phương trình tuyến tính):**
   * Giải hệ phương trình chuẩn $M c = b$ để tìm vector tham số $c = [c_1, c_2, \dots, c_m]^T$.
4. **Bước 4 (Kết luận và đánh giá sai số):**
   * Viết hàm thực nghiệm $g(x) = \sum_{j=1}^m c_j \varphi_j(x)$.
   * Tính sai số trung bình phương: $\sigma_n = \sqrt{\frac{1}{n}\sum_{i=1}^n (y_i - g(x_i))^2}$ (hoặc tổng bình phương phần dư $SSE = \sum_{i=1}^n (y_i - g(x_i))^2$).

### 8.2. Mô hình phi tuyến đưa được về dạng tuyến tính (Ví dụ: $y = a e^{b_1 \varphi_1(x) + \dots + b_m \varphi_m(x)} + C_0$)
1. **Bước 1 (Tuyến tính hóa mô hình):**
   * Chuyển hằng số đã biết (nếu có, ví dụ $C_0 = -2$ trong $y = ae^{bx} - 2$) sang vế trái: $\tilde{y}_i = y_i - C_0$.
   * Kiểm tra dấu dữ liệu và lấy logarit tự nhiên hai vế:
     * Nếu $\tilde{y}_i > 0, \forall i$: Đặt $Y_i = \ln(\tilde{y}_i)$ và $A = \ln(a) \Rightarrow Y = A + b_1 \varphi_1(x) + \dots + b_m \varphi_m(x)$.
     * Nếu $\tilde{y}_i < 0, \forall i$: Đặt $Y_i = \ln(-\tilde{y}_i)$ và $A = \ln(-a) \Rightarrow Y = A + b_1 \varphi_1(x) + \dots + b_m \varphi_m(x)$.
2. **Bước 2 (Biến đổi bảng dữ liệu):** Tính mảng giá trị mới $Y = [Y_1, \dots, Y_n]^T$ và đặt thêm hàm cơ sở hằng số $\varphi_0(x) = 1$ đi kèm tham số $A$.
3. **Bước 3 (Áp dụng thuật toán tuyến tính):**
   * Lập ma trận $\Phi = [1 \quad \varphi_1 \quad \dots \quad \varphi_m]$, tính $M = \Phi^T \Phi$ và $B = \Phi^T Y$.
   * Giải hệ $M X = B$ để thu được nghiệm $X = [A, b_1, \dots, b_m]^T$.
4. **Bước 4 (Khôi phục tham số gốc và tính sai số):**
   * Nếu $\tilde{y}_i > 0$: Tính $a = e^A$. Nếu $\tilde{y}_i < 0$: Tính $a = -e^A$.
   * Thay $a, b_1, \dots, b_m$ trở lại hàm gốc $g(x)$ và tính sai số trực tiếp trên bộ dữ liệu gốc $(x_i, y_i)$: $\sigma_n = \sqrt{\frac{1}{n}\sum_{i=1}^n (y_i - g(x_i))^2}$.