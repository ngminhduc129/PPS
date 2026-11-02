/* ============================================================
   polynomial.js — polynomial arithmetic and Horner evaluation
   Coefficients are stored in ascending order:
   coefs[0] + coefs[1] x + coefs[2] x^2 + ...
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});

  function trim(coefs, eps) {
    eps = eps == null ? 1e-14 : eps;
    let n = coefs.length;
    while (n > 1 && Math.abs(coefs[n - 1]) < eps) n--;
    return coefs.slice(0, n);
  }

  function add(a, b) {
    const n = Math.max(a.length, b.length);
    const out = new Array(n).fill(0);
    for (let i = 0; i < n; i++) out[i] = (a[i] || 0) + (b[i] || 0);
    return out;
  }

  function scale(a, s) {
    return a.map(function (x) { return x * s; });
  }

  function mul(a, b) {
    if (!a.length || !b.length) return [0];
    const out = new Array(a.length + b.length - 1).fill(0);
    for (let i = 0; i < a.length; i++) {
      if (Math.abs(a[i]) < 1e-18) continue;
      for (let j = 0; j < b.length; j++) {
        out[i + j] += a[i] * b[j];
      }
    }
    return out;
  }

  // Multiply polynomial by (x - root)
  function mulLinear(coefs, root) {
    return mul(coefs, [-root, 1]);
  }

  // Multiply polynomial by (x^2 - r2)
  function mulQuadratic(coefs, r2) {
    return mul(coefs, [-r2, 0, 1]);
  }

  function evaluate(coefs, x) {
    let value = 0;
    for (let i = coefs.length - 1; i >= 0; i--) {
      value = value * x + coefs[i];
    }
    return value;
  }

  // Horner evaluation of P(x) and P'(x)
  function evalWithDerivative(coefs, x) {
    const n = coefs.length;
    if (n === 0) return { value: 0, derivative: 0 };
    let value = coefs[n - 1];
    let derivative = 0;
    for (let i = n - 2; i >= 0; i--) {
      derivative = derivative * x + value;
      value = value * x + coefs[i];
    }
    return { value: value, derivative: derivative };
  }

  function derivativeCoefs(coefs) {
    if (coefs.length <= 1) return [0];
    return coefs.slice(1).map(function (c, i) { return c * (i + 1); });
  }

  function fromRoots(roots) {
    let poly = [1];
    for (const root of roots) poly = mulLinear(poly, root);
    return poly;
  }

  function sample(coefs, xmin, xmax, count, varName) {
    // For polynomials expressed in t, caller maps t back manually.
    const points = [];
    const n = count || 240;
    for (let i = 0; i <= n; i++) {
      const x = xmin + (xmax - xmin) * i / n;
      points.push({ x: x, y: evaluate(coefs, x) });
    }
    return points;
  }

  PPS.poly = {
    trim: trim,
    add: add,
    scale: scale,
    mul: mul,
    mulLinear: mulLinear,
    mulQuadratic: mulQuadratic,
    evaluate: evaluate,
    evalWithDerivative: evalWithDerivative,
    derivativeCoefs: derivativeCoefs,
    fromRoots: fromRoots,
    sample: sample
  };
})();
