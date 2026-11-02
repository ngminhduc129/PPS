# PPS Simulator — Web mô phỏng thuật toán phương pháp số

Web tĩnh mô phỏng **10 nhóm thuật toán & bảng so sánh** theo tài liệu `Algorithm.md`: từ sơ đồ Horner, mốc Chebyshev, nội suy Lagrange/Newton, nội suy trung tâm, nội suy ngược, spline đến bình phương tối thiểu.

- Giao diện **song ngữ Việt / Anh** (nút `VI | EN` trên header).
- **Dark / Light theme**, responsive cho điện thoại.
- **Không cần cài đặt, không cần build, không phụ thuộc thư viện** — chạy offline.
- Hiển thị **bảng trung gian của từng bước**, tô sáng ô/hệ số đang dùng và **đồ thị SVG**.
- Có **animation** cho bảng chia tổng hợp Horner, Newton, và vòng lặp nội suy ngược.

## Cách chạy

Mở trực tiếp:

```
PPS-Simulator/index.html
```

Nhấp đúp vào file là chạy trên Chrome / Edge / Firefox. Không cần `npm install` hay server.

> Nếu muốn chạy qua HTTP (tùy chọn):
> ```powershell
> cd PPS-Simulator
> npx serve .
> ```

## Các thuật toán (dựa trên Algorithm.md)

| # | Thuật toán | Nội dung mô phỏng |
|---|---|---|
| 1 | Sơ đồ Horner | Bảng chia tổng hợp lặp, `P(c)`, đạo hàm các cấp, khai triển Taylor, đa thức thương `Q(x)`, nhân ngược Horner `P(x) = (x−c)Q(x) + r` |
| 2 | Mốc Chebyshev | Sinh mốc tối ưu Chebyshev trên `[a, b]`, đánh giá chặn sai số lý thuyết, mô phỏng triệt tiêu hiện tượng Runge so với mốc cách đều |
| 3 | Nội suy Lagrange | Mẫu số `Dᵢ`, tử số `wᵢ(x)`, đa thức cơ sở `lᵢ(x)`, hệ số đa thức chính tắc, kiểm tra tại mốc |
| 4 | Newton tỷ sai phân | Bảng tỷ sai phân đầy đủ, vector tiến `C_tiến` (đường chéo), vector lùi `C_lùi` (hàng cuối), bảng tích lũy thừa `Bᵢ(x)`, nhân lồng Horner, đạo hàm cấp `r` |
| 5 | Newton mốc cách đều | Bảng sai phân `Δʲyᵢ`, tự chọn tiến/lùi theo `x̄`, hệ số theo `t`, giá trị và đạo hàm |
| 6 | Nội suy trung tâm | Gauss I, Gauss II, Stirling, Bessel; tự gợi ý công thức theo `t`; bảng sai phân trung tâm tô ô dùng; `f(x̄)`, `f′(x̄)` |
| 7 | Nội suy ngược | **Cách 1**: Đổi biến hàm ngược `x = P(y)` qua tỷ sai phân; **Cách 2**: Lặp điểm bất động `t = φ(t)` với bảng trace chi tiết, hệ số co `q`, residual và chặn sai số hậu nghiệm |
| 8 | Spline | Tuyến tính (bậc 1), Bậc 2 (`z₀ = S′(x₀)`), Bậc 3: Tự nhiên (`M₀ = Mₙ = 0`), Moment biên (`M₀, Mₙ`), Biên kẹp Clamped (`S′(x₀), S′(xₙ)`) giải bằng Thomas TDMA |
| 9 | Bình phương tối thiểu | Mô hình tuyến tính theo tham số: `Φ`, `ΦᵀΦ`, `Φᵀy`, nghiệm; Mô hình mũ `y = a·e^(bφ₁(x)+cφ₂(x)) + C₀` qua logarit; Mô hình lũy thừa `y = a·xᵇ`; Tính đầy đủ bộ chỉ số SSE, MAE, RMSE, Max Error, R², SST |
| 10 | Bảng so sánh thuật toán | Bảng tổng hợp đối chiếu toàn diện 9 nhóm thuật toán theo Mục 10 tài liệu `Algorithm.md` |

### Định dạng nhập dữ liệu

- Các dãy số cách nhau bởi **dấu phẩy, dấu cách, dấu chấm phẩy hoặc xuống dòng**.
- Ví dụ: `0, 1, 2, 3, 4` hoặc `0 1 2 3 4`.
- Hỗ trợ số âm và ký hiệu khoa học: `-1.5`, `1e-3`.
- Nút **“Ví dụ mẫu”** nạp sẵn bộ dữ liệu phù hợp cho từng thuật toán.

## Cấu trúc dự án

```
PPS-Simulator/
├─ index.html                  # khung ứng dụng
├─ assets/styles.css           # giao diện, theme, bảng, biểu đồ, responsive
├─ src/
│  ├─ i18n.js                  # từ điển VI/EN
│  ├─ format.js                # định dạng số, đa thức, parse dữ liệu
│  ├─ linalg.js                # khử Gauss có pivot, Thomas, ma trận
│  ├─ polynomial.js            # cộng/nhân/đánh giá đa thức, Horner
│  ├─ horner.js
│  ├─ lagrange.js
│  ├─ newton.js                # tỷ sai phân + sai phân cách đều
│  ├─ central.js               # Gauss I/II, Stirling, Bessel
│  ├─ inverse.js               # nội suy ngược
│  ├─ spline.js                # spline cấp 1/2/3
│  ├─ leastSquares.js          # hệ phương trình chuẩn + logarit hóa
│  ├─ chart.js                 # vẽ SVG
│  ├─ ui.js                    # helper DOM
│  ├─ views/                   # giao diện từng thuật toán
│  └─ app.js                   # router, sidebar, theme, i18n
├─ tests/
│  ├─ run-tests.mjs            # unit test engine (Node)
│  └─ browser-smoke.mjs        # smoke test giao diện bằng Edge headless
└─ README.md
```

## Kiểm thử

```powershell
cd PPS-Simulator

# Unit test engine: 28 test
node tests/run-tests.mjs

# Smoke test giao diện: mở Edge headless, đi qua 8 thuật toán,
# bấm ví dụ mẫu, đổi tab, đổi ngôn ngữ và theme
node tests/browser-smoke.mjs
```

Bộ test kiểm tra các trường hợp tiêu biểu:

- Horner với đa thức bậc 3 và đạo hàm cấp 1–3.
- Newton tỷ sai phân trùng kết quả Lagrange tại điểm bất kỳ.
- Newton tiến/lùi tái tạo chính xác đa thức bậc 3.
- Gauss I, Gauss II, Stirling, Bessel tái tạo chính xác đa thức bậc 3.
- Nội suy ngược giải `x² = 6.25`.
- Spline qua đúng các mốc; spline cấp 2 liên tục đạo hàm; spline cấp 3 tự nhiên có `α₁ = αₙ = 0`.
- Least squares khôi phục chính xác dữ liệu tuyến tính/bậc hai và mô hình `y = 2e^{0.5x} + 1`.

## Ghi chú kỹ thuật

- Hệ số đa thức trong code lưu theo **bậc tăng dần**, riêng ô nhập Horner nhận **bậc cao → thấp** đúng như tài liệu.
- Các công thức trung tâm dùng chỉ số tương đối `y₋ᵢ … yᵢ`; bảng hiển thị nhãn tương ứng và tô đúng các ô được dùng (kể cả hai ô trung bình của Stirling/Bessel).
- Spline cấp 2 giải lần lượt `mₖ₊₁ = γₖ − mₖ` từ điều kiện biên `m₁`; spline cấp 3 giải hệ ba đường chéo bằng **Thomas algorithm** với biên tự nhiên.
- Least squares giải hệ chuẩn `ΦᵀΦc = Φᵀy` bằng **khử Gauss có chọn pivot**; mô hình logarit kiểm tra tất cả `y − C₀` cùng dấu trước khi lấy `ln`.
- Đồ thị là SVG thuần, tự co giãn theo khung; không dùng thư viện ngoài.
- Dữ liệu và lựa chọn ngôn ngữ/theme được lưu cục bộ (`localStorage`).

## Nguồn

Nội dung thuật toán dựa trên file `Midterm - PPS.md` trong cùng thư mục vault. Giao diện và phần mô phỏng được viết mới cho mục đích học tập.
