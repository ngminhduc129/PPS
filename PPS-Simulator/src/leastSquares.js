/* ============================================================
   leastSquares.js — linear-in-parameters least squares,
                     log-linearized exponential, and power-law fitting
                     Updated based on Section 9 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});

  const BASIS = [
    { key: 'one', labelKey: 'basis.one', fn: function () { return 1; } },
    { key: 'x', labelKey: 'basis.x', fn: function (x) { return x; } },
    { key: 'x2', labelKey: 'basis.x2', fn: function (x) { return x * x; } },
    { key: 'x3', labelKey: 'basis.x3', fn: function (x) { return x * x * x; } },
    { key: 'sin', labelKey: 'basis.sin', fn: function (x) { return Math.sin(x); } },
    { key: 'cos', labelKey: 'basis.cos', fn: function (x) { return Math.cos(x); } },
    { key: 'exp', labelKey: 'basis.exp', fn: function (x) { return Math.exp(x); } },
    { key: 'ln', labelKey: 'basis.ln', fn: function (x) { if (!(x > 0)) throw new Error('domain'); return Math.log(x); } },
    { key: 'inv', labelKey: 'basis.inv', fn: function (x) { if (Math.abs(x) < 1e-15) throw new Error('domain'); return 1 / x; } },
    { key: 'sqrt', labelKey: 'basis.sqrt', fn: function (x) { if (x < 0) throw new Error('domain'); return Math.sqrt(x); } }
  ];

  function getBasis(keys) {
    if (!keys || !keys.length) throw new Error('basis');
    return keys.map(function (key) {
      const found = BASIS.find(function (b) { return b.key === key; });
      if (!found) throw new Error('basis');
      return found;
    });
  }

  function validatePoints(xs, ys) {
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 2) throw new Error('tooFew');
    for (let i = 0; i < xs.length; i++) {
      if (!isFinite(xs[i]) || !isFinite(ys[i])) throw new Error('numeric');
    }
  }

  function buildPhi(xs, basis) {
    return xs.map(function (x) {
      return basis.map(function (b) {
        const value = b.fn(x);
        if (!isFinite(value)) throw new Error('domain');
        return value;
      });
    });
  }

  /**
   * Đánh giá các chỉ số sai số theo Mục 9.2 Bước 6 - Algorithm.md:
   * SSE = sum(e_i^2)
   * MAE = (1 / N) * sum(|e_i|)
   * RMSE = sqrt(SSE / N)
   * Max_Error = max(|e_i|)
   * SST = sum((y_i - mean(y))^2)
   * R_squared = 1 - SSE / SST
   */
  function metrics(ys, fitted, numParams) {
    const N = ys.length;
    let sse = 0;
    let sumAbs = 0;
    let maxError = 0;

    for (let i = 0; i < N; i++) {
      const e = ys[i] - fitted[i];
      const absE = Math.abs(e);
      sse += e * e;
      sumAbs += absE;
      if (absE > maxError) maxError = absE;
    }

    const mean = ys.reduce(function (a, b) { return a + b; }, 0) / N;
    let sst = 0;
    for (let i = 0; i < N; i++) {
      sst += (ys[i] - mean) * (ys[i] - mean);
    }

    const mae = sumAbs / N;
    const rmse = Math.sqrt(sse / N);
    const r2 = sst < 1e-15 ? (sse < 1e-15 ? 1 : 0) : 1 - sse / sst;
    const p = (typeof numParams === 'number' && numParams > 0) ? numParams : 1;
    const r2adj = (N > p && sst > 1e-15) ? 1 - ((1 - r2) * (N - 1)) / (N - p) : r2;

    return {
      sse: sse,
      mae: mae,
      rmse: rmse,
      sigma: rmse,
      maxError: maxError,
      sst: sst,
      r2: r2,
      r2adj: r2adj
    };
  }

  /**
   * Linear least squares g(x) = sum c_j phi_j(x).
   */
  function linear(xs, ys, keys) {
    validatePoints(xs, ys);
    const basis = getBasis(keys);
    if (xs.length < basis.length) throw new Error('leastSquares');

    const Phi = buildPhi(xs, basis);
    const MT = PPS.linalg.transpose(Phi);
    const M = PPS.linalg.matmul(MT, Phi);
    const rhs = PPS.linalg.matvec(MT, ys);
    const { solution: coef } = PPS.linalg.solve(M, rhs);

    let invM = null;
    let pinv = null;
    try {
      invM = PPS.linalg.invert(M);
      pinv = PPS.linalg.matmul(invM, MT);
    } catch (e) {
      invM = null;
      pinv = null;
    }

    let invPhi = null;
    if (xs.length === basis.length) {
      try {
        invPhi = PPS.linalg.invert(Phi);
      } catch (e) {
        invPhi = null;
      }
    }

    const evaluate = function (x) {
      let sum = 0;
      for (let j = 0; j < basis.length; j++) sum += coef[j] * basis[j].fn(x);
      return sum;
    };

    const fitted = xs.map(evaluate);
    const stats = metrics(ys, fitted, basis.length);
    const residuals = xs.map(function (x, i) {
      return { x: x, y: ys[i], fitted: fitted[i], residual: ys[i] - fitted[i] };
    });

    return {
      mode: 'linear',
      basis: basis,
      phi: Phi,
      M: M,
      invM: invM,
      invPhi: invPhi,
      pinv: pinv,
      rhs: rhs,
      coef: coef,
      evaluate: evaluate,
      fitted: fitted,
      residuals: residuals,
      sse: stats.sse,
      mae: stats.mae,
      rmse: stats.rmse,
      sigma: stats.sigma,
      maxError: stats.maxError,
      sst: stats.sst,
      r2: stats.r2,
      r2adj: stats.r2adj,
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; })
    };
  }

  /**
   * Fit phi tuyến mũ y = a * exp(sum b_j phi_j(x)) + C0 (Mục 9.2 - Algorithm.md)
   * Tuyến tính hóa: Y_lin = ln|y - C0| = A_0 + sum b_j phi_j(x)
   */
  function logLinear(xs, ys, keys, c0) {
    validatePoints(xs, ys);
    c0 = isFinite(c0) ? Number(c0) : 0;
    const basis = getBasis(keys);

    const yt = ys.map(function (y) { return y - c0; });
    const positive = yt.filter(function (v) { return v > 1e-14; }).length;
    const negative = yt.filter(function (v) { return v < -1e-14; }).length;
    if (positive + negative !== yt.length) throw new Error('sign');
    if (positive !== yt.length && negative !== yt.length) throw new Error('sameSign');
    const sign = positive === yt.length ? 1 : -1;
    const Y = yt.map(function (v) { return Math.log(Math.abs(v)); });

    const extended = [{ key: 'one', labelKey: 'basis.one', fn: function () { return 1; } }].concat(basis);
    if (xs.length < extended.length) throw new Error('leastSquares');

    const Phi = buildPhi(xs, extended);
    const MT = PPS.linalg.transpose(Phi);
    const M = PPS.linalg.matmul(MT, Phi);
    const rhs = PPS.linalg.matvec(MT, Y);
    const { solution: coef } = PPS.linalg.solve(M, rhs);

    let invM = null;
    let pinv = null;
    try {
      invM = PPS.linalg.invert(M);
      pinv = PPS.linalg.matmul(invM, MT);
    } catch (e) {
      invM = null;
      pinv = null;
    }

    let invPhi = null;
    if (xs.length === extended.length) {
      try {
        invPhi = PPS.linalg.invert(Phi);
      } catch (e) {
        invPhi = null;
      }
    }

    const A = coef[0];
    const b = coef.slice(1);
    const a = sign * Math.exp(A);

    const evaluate = function (x) {
      let exponent = 0;
      for (let j = 0; j < basis.length; j++) exponent += b[j] * basis[j].fn(x);
      return a * Math.exp(exponent) + c0;
    };

    const fitted = xs.map(evaluate);
    const stats = metrics(ys, fitted, extended.length);
    const residuals = xs.map(function (x, i) {
      return { x: x, y: ys[i], fitted: fitted[i], residual: ys[i] - fitted[i] };
    });

    return {
      mode: 'log',
      basis: basis,
      extended: extended,
      phi: Phi,
      M: M,
      invM: invM,
      invPhi: invPhi,
      pinv: pinv,
      rhs: rhs,
      coef: coef,
      A: A,
      a: a,
      b: b,
      c0: c0,
      sign: sign,
      evaluate: evaluate,
      fitted: fitted,
      residuals: residuals,
      sse: stats.sse,
      mae: stats.mae,
      rmse: stats.rmse,
      sigma: stats.sigma,
      maxError: stats.maxError,
      sst: stats.sst,
      r2: stats.r2,
      r2adj: stats.r2adj,
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; })
    };
  }

  /**
   * Fit mô hình lũy thừa y = a * x^b (Mục 9.2 - Nhận dạng nhanh)
   * Tuyến tính hóa: ln(y) = ln(a) + b * ln(x) với x > 0, y > 0
   */
  function fitPower(xs, ys) {
    validatePoints(xs, ys);
    for (let i = 0; i < xs.length; i++) {
      if (xs[i] <= 0 || ys[i] <= 0) throw new Error('positiveOnly');
    }

    const lnX = xs.map(Math.log);
    const lnY = ys.map(Math.log);

    // Hồi quy ln(y) theo (1, ln(x))
    const linResult = linear(lnX, lnY, ['one', 'x']);
    const A = linResult.coef[0];
    const b = linResult.coef[1];
    const a = Math.exp(A);

    const evaluate = function (x) {
      if (x <= 0) return NaN;
      return a * Math.pow(x, b);
    };

    const fitted = xs.map(evaluate);
    const stats = metrics(ys, fitted, 2);
    const residuals = xs.map(function (x, i) {
      return { x: x, y: ys[i], fitted: fitted[i], residual: ys[i] - fitted[i] };
    });

    return {
      mode: 'power',
      a: a,
      b: b,
      A: A,
      basis: linResult.basis,
      phi: linResult.phi,
      M: linResult.M,
      invM: linResult.invM,
      invPhi: linResult.invPhi,
      pinv: linResult.pinv,
      rhs: linResult.rhs,
      evaluate: evaluate,
      fitted: fitted,
      residuals: residuals,
      sse: stats.sse,
      mae: stats.mae,
      rmse: stats.rmse,
      sigma: stats.sigma,
      maxError: stats.maxError,
      sst: stats.sst,
      r2: stats.r2,
      r2adj: stats.r2adj,
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; })
    };
  }

  /**
   * Sinh cấu trúc bảng tổng tích lũy (Summation table) phục vụ trình bày bài thi:
   * Cột: i, x_i, y_i, các cột phi_j, phi_j * y_i, và hàng tổng sum ở cuối.
   */
  function buildSumTable(xs, ys, basis) {
    const N = xs.length;
    const cols = [
      { key: 'i', label: 'i' },
      { key: 'x', label: 'xᵢ' },
      { key: 'y', label: 'yᵢ' }
    ];

    basis.forEach(function (b) {
      if (b.key !== 'one' && b.key !== 'x') {
        cols.push({ key: 'phi_' + b.key, label: 'φ(' + b.key + ')' });
      }
    });

    basis.forEach(function (b) {
      cols.push({ key: 'prod_' + b.key, label: 'φ(' + b.key + ')·yᵢ' });
    });

    const rows = [];
    const sums = {};
    cols.forEach(function (c) { sums[c.key] = 0; });

    for (let i = 0; i < N; i++) {
      const row = { i: i, x: xs[i], y: ys[i] };
      sums.x += xs[i];
      sums.y += ys[i];

      basis.forEach(function (b) {
        const val = b.fn(xs[i]);
        if (b.key !== 'one' && b.key !== 'x') {
          row['phi_' + b.key] = val;
          sums['phi_' + b.key] += val;
        }
        const prod = val * ys[i];
        row['prod_' + b.key] = prod;
        sums['prod_' + b.key] += prod;
      });
      rows.push(row);
    }

    return { cols: cols, rows: rows, sums: sums };
  }

  PPS.leastSquares = {
    BASIS: BASIS,
    getBasis: getBasis,
    buildPhi: buildPhi,
    metrics: metrics,
    linear: linear,
    logLinear: logLinear,
    fitPower: fitPower,
    buildSumTable: buildSumTable
  };
})();
