/* ============================================================
   horner.js — synthetic division, P(c), derivatives, Taylor expansion,
               and reverse Horner multiplication
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const { factorial } = PPS.format;

  /**
   * Evaluate P(c) and its k-th derivative using k+1 synthetic divisions.
   * coefsDesc: coefficients from highest degree to constant term.
   */
  function evaluate(coefsDesc, c, k) {
    if (!coefsDesc.length) throw new Error('empty');
    c = Number(c);
    k = Math.max(0, Math.floor(Number(k) || 0));
    if (!isFinite(c)) throw new Error('numeric');

    const N = coefsDesc.length - 1;
    const maxPass = Math.min(k, N);
    const b = coefsDesc.slice();
    const passes = [];
    const derivatives = [];

    for (let m = 0; m <= maxPass; m++) {
      const before = b.slice();
      for (let i = 1; i <= N - m; i++) {
        b[i] = b[i - 1] * c + b[i];
      }
      const remainder = b[N - m];
      const derivativeValue = factorial(m) * remainder;
      derivatives[m] = derivativeValue;
      passes.push({
        m: m,
        before: before,
        after: b.slice(),
        remainder: remainder,
        remainderIndex: N - m,
        derivativeValue: derivativeValue,
        factorial: factorial(m)
      });
    }

    const value = passes[0] ? passes[0].remainder : 0;
    const derivative = k <= N ? derivatives[k] : 0;
    // Đa thức thương Q(x) khi chia P(x) cho (x - c)
    const quotientDesc = passes[0] ? passes[0].after.slice(0, N) : [];

    // Khai triển Taylor: P(x) = sum (rem_m * (x - c)^m) với rem_m = P^(m)(c)/m!
    const taylorTerms = passes.map(function (p) {
      return {
        m: p.m,
        rem: p.remainder,
        derivative: p.derivativeValue
      };
    });

    return {
      N: N,
      c: c,
      k: k,
      value: value,
      derivative: derivative,
      derivatives: derivatives,
      passes: passes,
      quotientDesc: quotientDesc,
      taylorTerms: taylorTerms,
      kTooBig: k > N
    };
  }

  /**
   * Thuật toán Nhân ngược Horner (Mục 1.3 - Algorithm.md):
   * Tính P(x) = (x - c)*Q(x) + r
   * qAsc[0..N-1]: hệ số của Q(x) theo bậc tăng.
   * Trả về mảng a[0..N] theo bậc tăng.
   */
  function reverseMultiplyAsc(qAsc, c, r) {
    c = Number(c);
    r = Number(r || 0);
    const N = qAsc.length;
    if (N === 0) return [r];
    const a = new Array(N + 1);
    a[0] = r - c * qAsc[0];
    for (let k = 1; k < N; k++) {
      a[k] = qAsc[k - 1] - c * qAsc[k];
    }
    a[N] = qAsc[N - 1];
    return a;
  }

  /**
   * Nhân ngược Horner với mảng hệ số Q theo bậc giảm (Desc).
   * qDesc: [q_{N-1}, ..., q_0]
   * Trả về aDesc theo bậc giảm.
   */
  function reverseMultiplyDesc(qDesc, c, r) {
    const qAsc = qDesc.slice().reverse();
    const aAsc = reverseMultiplyAsc(qAsc, c, r);
    return aAsc.slice().reverse();
  }

  PPS.horner = {
    evaluate: evaluate,
    reverseMultiplyAsc: reverseMultiplyAsc,
    reverseMultiplyDesc: reverseMultiplyDesc
  };
})();
