/* ============================================================
   spline.js — Linear, Quadratic, and Cubic Splines
               Cubic Spline supports Natural, Given Moments, and Clamped
               Updated based on Section 8 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const poly = PPS.poly;

  function validate(xs, ys, minPoints) {
    minPoints = minPoints || 2;
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < minPoints) throw new Error('tooFew');
    for (let i = 0; i < xs.length; i++) {
      if (!isFinite(xs[i]) || !isFinite(ys[i])) throw new Error('numeric');
      if (i > 0 && xs[i] <= xs[i - 1]) throw new Error('order');
    }
  }

  function steps(xs, ys) {
    const n = xs.length;
    const h = [];
    const gammaQuad = [];
    const gammaCubic = [];
    for (let i = 0; i < n - 1; i++) {
      h.push(xs[i + 1] - xs[i]);
    }
    for (let i = 0; i < n - 1; i++) {
      gammaQuad.push((ys[i + 1] - ys[i]) / h[i]);
    }
    for (let i = 0; i < n - 2; i++) {
      const g1 = (ys[i + 2] - ys[i + 1]) / h[i + 1];
      const g0 = (ys[i + 1] - ys[i]) / h[i];
      gammaCubic.push(6 * (g1 - g0));
    }
    return { h: h, gammaQuad: gammaQuad, gammaCubic: gammaCubic };
  }

  function buildSegments(xs, ys, coefsList, degree, localDefs) {
    const segs = [];
    for (let i = 0; i < coefsList.length; i++) {
      const c = coefsList[i];
      const local = localDefs ? localDefs[i] : null;
      const highest = c.length >= 1 ? c[c.length - 1] : 0;
      const second = c.length >= 2 ? c[c.length - 2] : 0;
      segs.push({
        index: i,
        x0: xs[i],
        x1: xs[i + 1],
        coef: c,
        a: highest,
        b: second,
        local: local,
        degree: degree,
        evaluate: function (x) {
          if (local) {
            const dx = x - local.xi;
            return local.A + local.B * dx + local.C * dx * dx + (local.D ? local.D * Math.pow(dx, 3) : 0);
          }
          return poly.evaluate(c, x);
        },
        evalDerivative: function (x) {
          if (local) {
            const dx = x - local.xi;
            return local.B + 2 * local.C * dx + (local.D ? 3 * local.D * dx * dx : 0);
          }
          return poly.evalWithDerivative(c, x).derivative;
        }
      });
    }
    return segs;
  }

  /**
   * Linear Spline (bậc 1)
   */
  function linear(xs, ys) {
    validate(xs, ys, 2);
    const n = xs.length;
    const { h } = steps(xs, ys);
    const coefsList = [];
    const localDefs = [];

    for (let i = 0; i < n - 1; i++) {
      const a1 = (ys[i + 1] - ys[i]) / h[i];
      const a0 = ys[i] - a1 * xs[i];
      coefsList.push([a0, a1]);
      localDefs.push({ xi: xs[i], A: ys[i], B: a1, C: 0 });
    }

    return {
      degree: 1,
      h: h,
      segments: buildSegments(xs, ys, coefsList, 1, localDefs),
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; })
    };
  }

  /**
   * Quadratic Spline (bậc 2): z_0 = S'(x_0) cho trước (mặc định z_0 = 0)
   */
  function quadratic(xs, ys, z0) {
    validate(xs, ys, 2);
    z0 = isFinite(z0) ? Number(z0) : 0;
    const n = xs.length;
    const { h, gammaQuad } = steps(xs, ys);
    const z = [z0];

    for (let i = 0; i < n - 1; i++) {
      const zNext = 2 * gammaQuad[i] - z[i];
      z.push(zNext);
    }

    const coefsList = [];
    const localDefs = [];

    for (let i = 0; i < n - 1; i++) {
      const hi = h[i];
      const a2 = (z[i + 1] - z[i]) / (2 * hi);
      const a1 = z[i] - 2 * a2 * xs[i];
      const a0 = ys[i] - a2 * xs[i] * xs[i] - a1 * xs[i];
      coefsList.push([a0, a1, a2]);
      localDefs.push({ xi: xs[i], A: ys[i], B: z[i], C: a2 });
    }

    return {
      degree: 2,
      h: h,
      z: z,
      z0: z0,
      gamma: gammaQuad,
      segments: buildSegments(xs, ys, coefsList, 2, localDefs),
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; })
    };
  }

  /**
   * Cubic Spline (bậc 3) & Thuật toán Thomas (TDMA)
   * opts:
   *   bcType: 'natural' | 'moments' | 'clamped'
   *   M0, Mn: cho bcType = 'moments'
   *   s0, sn: cho bcType = 'clamped' (s0 = S'(x0), sn = S'(xn))
   */
  function cubic(xs, ys, opts) {
    validate(xs, ys, 2);
    opts = opts || {};
    const bcType = opts.bcType || 'natural';
    const n = xs.length;
    const { h, gammaCubic } = steps(xs, ys);
    let M = new Array(n).fill(0);
    const system = [];

    if (n === 2) {
      // 2 mốc: nếu natural hoặc moments, M = [M0, M1]
      const M0 = bcType === 'moments' ? Number(opts.M0 || 0) : 0;
      const M1 = bcType === 'moments' ? Number(opts.Mn || 0) : 0;
      M = [M0, M1];
    } else if (bcType === 'clamped') {
      // Clamped boundary: S'(x0) = s0, S'(xn) = sn
      const s0 = isFinite(opts.s0) ? Number(opts.s0) : 0;
      const sn = isFinite(opts.sn) ? Number(opts.sn) : 0;

      // Hệ (n) ẩn M_0..M_{n-1} gồm n phương trình 3 đường chéo
      const mSize = n;
      const lower = [];
      const diag = [];
      const upper = [];
      const rhs = [];

      // Hàng 0: 2 h_0 M_0 + h_0 M_1 = 6 * ((y_1 - y_0)/h_0 - s0)
      diag.push(2 * h[0]);
      upper.push(h[0]);
      rhs.push(6 * ((ys[1] - ys[0]) / h[0] - s0));
      system.push({ i: 0, left: 0, mid: 2 * h[0], right: h[0], rhs: rhs[0] });

      // Hàng nội điểm 1..n-2
      for (let k = 0; k < n - 2; k++) {
        lower.push(h[k]);
        diag.push(2 * (h[k] + h[k + 1]));
        upper.push(h[k + 1]);
        rhs.push(gammaCubic[k]);
        system.push({
          i: k + 1,
          left: h[k],
          mid: 2 * (h[k] + h[k + 1]),
          right: h[k + 1],
          rhs: gammaCubic[k]
        });
      }

      // Hàng cuối n-1: h_{n-2} M_{n-2} + 2 h_{n-2} M_{n-1} = 6 * (sn - (y_{n-1} - y_{n-2})/h_{n-2})
      const lastH = h[n - 2];
      lower.push(lastH);
      diag.push(2 * lastH);
      const rhsLast = 6 * (sn - (ys[n - 1] - ys[n - 2]) / lastH);
      rhs.push(rhsLast);
      system.push({ i: n - 1, left: lastH, mid: 2 * lastH, right: 0, rhs: rhsLast });

      M = PPS.linalg.thomas(lower, diag, upper, rhs);
    } else {
      // Natural hoặc Given Moments
      const M0 = bcType === 'moments' ? Number(opts.M0 || 0) : 0;
      const Mn = bcType === 'moments' ? Number(opts.Mn || 0) : 0;
      M[0] = M0;
      M[n - 1] = Mn;

      const unknowns = n - 2;
      const lower = [];
      const diag = [];
      const upper = [];
      const rhs = [];

      for (let k = 0; k < unknowns; k++) {
        diag.push(2 * (h[k] + h[k + 1]));
        let r = gammaCubic[k];
        if (k === 0) r -= h[0] * M0;
        if (k === unknowns - 1) r -= h[k + 1] * Mn;
        rhs.push(r);

        if (k > 0) lower.push(h[k]);
        if (k < unknowns - 1) upper.push(h[k + 1]);

        system.push({
          i: k + 1,
          left: k === 0 ? 0 : h[k],
          mid: 2 * (h[k] + h[k + 1]),
          right: k === unknowns - 1 ? 0 : h[k + 1],
          rhs: r
        });
      }

      const solved = PPS.linalg.thomas(lower, diag, upper, rhs);
      for (let i = 0; i < solved.length; i++) {
        M[i + 1] = solved[i];
      }
    }

    // Xây dựng đa thức từng đoạn theo biến cục bộ (x - x_i) (Mục 8.4 C - Algorithm.md)
    // S_i(x) = A_i + B_i(x - x_i) + C_i(x - x_i)^2 + D_i(x - x_i)^3
    const localDefs = [];
    const coefsList = [];

    for (let i = 0; i < n - 1; i++) {
      const hi = h[i];
      const Ai = ys[i];
      const Bi = (ys[i + 1] - ys[i]) / hi - (hi / 6) * (2 * M[i] + M[i + 1]);
      const Ci = M[i] / 2;
      const Di = (M[i + 1] - M[i]) / (6 * hi);

      localDefs.push({ xi: xs[i], A: Ai, B: Bi, C: Ci, D: Di });

      // Chuyển sang dạng chính tắc a0 + a1*x + a2*x^2 + a3*x^3
      const xi = xs[i];
      const a3 = Di;
      const a2 = Ci - 3 * Di * xi;
      const a1 = Bi - 2 * Ci * xi + 3 * Di * xi * xi;
      const a0 = Ai - Bi * xi + Ci * xi * xi - Di * Math.pow(xi, 3);
      coefsList.push([a0, a1, a2, a3]);
    }

    return {
      degree: 3,
      bcType: bcType,
      h: h,
      gamma: gammaCubic,
      M: M,
      alpha: M,
      system: system,
      segments: buildSegments(xs, ys, coefsList, 3, localDefs),
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; })
    };
  }

  function evaluate(spline, x) {
    const segs = spline.segments;
    if (!segs.length) return NaN;
    if (x <= segs[0].x0) return segs[0].evaluate(x);
    if (x >= segs[segs.length - 1].x1) return segs[segs.length - 1].evaluate(x);
    for (const seg of segs) {
      if (x >= seg.x0 - 1e-12 && x <= seg.x1 + 1e-12) return seg.evaluate(x);
    }
    return segs[segs.length - 1].evaluate(x);
  }

  function verify(xs, ys, spline) {
    return xs.map(function (x, i) {
      const value = evaluate(spline, x);
      return { index: i, x: x, y: ys[i], computed: value, error: value - ys[i] };
    });
  }

  PPS.spline = {
    linear: linear,
    quadratic: quadratic,
    cubic: cubic,
    evaluate: evaluate,
    verify: verify,
    steps: steps
  };
})();
