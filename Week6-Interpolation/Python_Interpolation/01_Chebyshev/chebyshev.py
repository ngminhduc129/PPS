"""
Mốc Chebyshev (Chebyshev Nodes)
Nhập: input.txt
Xuất: output.txt
"""
import os, sys, argparse
import numpy as np
import pandas as pd

def choose_chebyshev_nodes(a: float, b: float, n: int):
    condition = 1 if (a == -1.0 and b == 1.0) else 0
    i_values = np.arange(0, n + 1)
    t_values = np.cos(np.pi * (2 * i_values + 1) / (2 * n + 2))
    x_values = t_values if condition == 1 else 0.5 * ((b - a) * t_values + (b + a))
    df = pd.DataFrame({"i": i_values, "t_i (cos)": t_values, "x_i (mốc chọn)": x_values})
    return df, x_values

def read_input_txt(filepath: str):
    if not os.path.exists(filepath):
        return None, None, None
    a, b, n = None, None, None
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"): continue
            if "=" in line:
                k, v = line.split("=", 1)
                k, v = k.strip().lower(), v.strip().replace(",", ".")
                if k == "a": a = float(v)
                elif k == "b": b = float(v)
                elif k == "n": n = int(float(v))
            else:
                parts = line.replace(";", " ").replace(",", " ").split()
                if len(parts) >= 3:
                    a, b, n = float(parts[0]), float(parts[1]), int(parts[2])
                    break
    return a, b, n

def main():
    s_dir = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser()
    p.add_argument("input_file", nargs="?", default=os.path.join(s_dir, "input.txt"))
    p.add_argument("output_file", nargs="?", default=os.path.join(s_dir, "output.txt"))
    args = p.parse_args()

    print(f"-> Đọc file input: {args.input_file}")
    a, b, n = read_input_txt(args.input_file)
    if a is None or b is None or n is None: a, b, n = -1.0, 1.0, 4

    df, x_values = choose_chebyshev_nodes(a, b, n)
    out = ["="*60, "KẾT QUẢ TÍNH MỐC CHEBYSHEV", "="*60,
           f"Khoảng [a, b]: [{a}, {b}]", f"Bậc n: {n} (Tổng số mốc: {n + 1})", "-"*60,
           df.to_string(index=False), "-"*60, f"Danh sách mốc x_i: {x_values.tolist()}"]
    content = "\n".join(out)
    print("\n" + content)
    with open(args.output_file, "w", encoding="utf-8") as f:
        f.write(content + "\n")
    print(f"\n-> Đã lưu kết quả ra: '{args.output_file}'")

if __name__ == "__main__":
    main()
