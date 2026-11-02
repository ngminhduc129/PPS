"""
apply_universal_io.py - Script cập nhật đồng bộ tất cả các thư mục phương pháp với pps_io và auto-windowing
"""
import os, shutil

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PPS_IO_SRC = os.path.join(BASE_DIR, "pps_io.py")

folders = [
    "01_Chebyshev",
    "02_Horner",
    "03_Lagrange",
    "04_Newton_Divided",
    "05_Newton_Finite",
    "06_Gauss1",
    "07_Gauss2",
    "08_Stirling",
    "09_Bessel",
    "10_Inverse_Function",
    "11_Loop_Inverse",
    "12_Spline",
    "13_Least_Squares"
]

# Copy pps_io.py vào từng folder để mỗi phương pháp tự chứa (self-contained)
for f in folders:
    dest = os.path.join(BASE_DIR, f, "pps_io.py")
    shutil.copy2(PPS_IO_SRC, dest)
    print(f"Đã copy pps_io.py vào {f}")

print(f"Hoàn tất copy pps_io.py vào {len(folders)} thư mục phương pháp!")
