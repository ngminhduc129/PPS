/* ============================================================
   lagrange.js — Lagrange interpolation at arbitrary nodes
   Updated with explicit product table omega(x) and Horner division table
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const poly = PPS.poly;
  const fmt = (PPS.format && PPS.format.fmt) ? PPS.format.fmt : function (v) { return String(v); };

  function padPolyDesc(desc, totalDegree) {
    const res = new Array(totalDegree).fill(0);
    const start = totalDegree - desc.length;
    for (let i = 0; i < desc.length; i++) {
      res[start + i] = desc[i];
    }
    return res;
  }

  function interpolate(xs, ys) {
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 2) throw new Error('tooFew');

    const n = xs.length;
    for (let i = 0; i < n; i++) {
      if (!isFinite(xs[i]) || !isFinite(ys[i])) throw new Error('numeric');
      for (let j = i + 1; j < n; j++) {
        if (Math.abs(xs[i] - xs[j]) < 1e-12) throw new Error('dupX');
      }
    }

    // 1. LẬP BẢNG TÍCH (Product Table): Tính omega(x) = prod_{j=0}^{n-1} (x - x_j)
    const productTable = [];
    let omegaAsc = [1];
    productTable.push({
      step: 0,
      factor: '1',
      node: null,
      polyAsc: [1],
      polyDesc: padPolyDesc([1], n + 1)
    });

    for (let j = 0; j < n; j++) {
      omegaAsc = poly.mulLinear(omegaAsc, xs[j]);
      const signStr = xs[j] >= 0 ? ' − ' + fmt(xs[j], 4) : ' + ' + fmt(-xs[j], 4);
      productTable.push({
        step: j + 1,
        factor: '(x' + signStr + ')',
        node: xs[j],
        polyAsc: omegaAsc.slice(),
        polyDesc: padPolyDesc(omegaAsc.slice().reverse(), n + 1)
      });
    }

    const omegaDesc = omegaAsc.slice().reverse();

    // 2. LẬP BẢNG CHIA (Division Table): Chia Horner omega(x) cho (x - x_i)
    // w_i(x) = omega(x) / (x - x_i)
    const divisionTable = [];
    const basis = [];
    let coef = new Array(n).fill(0);

    for (let i = 0; i < n; i++) {
      const xi = xs[i];
      // Phép chia Horner: omegaDesc chia cho (x - xi)
      const bDesc = [omegaDesc[0]]; // b_0 = a_n = 1
      for (let k = 1; k < n; k++) {
        bDesc.push(omegaDesc[k] + xi * bDesc[k - 1]);
      }
      const remainder = omegaDesc[n] + xi * bDesc[n - 1]; // luon bang 0 vi omega(xi) = 0
      const wiAsc = bDesc.slice().reverse();

      // Mẫu số D_i = prod_{j != i} (xi - xj) = w_i(xi)
      let Di = 1;
      for (let j = 0; j < n; j++) {
        if (j === i) continue;
        Di *= (xi - xs[j]);
      }
      if (Math.abs(Di) < 1e-13) throw new Error('dupX');

      const ci = ys[i] / Di;
      const li = poly.scale(wiAsc, 1 / Di);
      const termPoly = poly.scale(wiAsc, ci);
      coef = poly.add(coef, termPoly);

      divisionTable.push({
        i: i,
        xi: xi,
        yi: ys[i],
        Di: Di,
        ci: ci,
        wiDesc: bDesc,
        wiAsc: wiAsc,
        remainder: remainder,
        termPoly: termPoly
      });

      basis.push({
        i: i,
        xi: xi,
        yi: ys[i],
        Di: Di,
        ci: ci,
        wi: wiAsc,
        li: li,
        termPoly: termPoly
      });
    }

    coef = poly.trim(coef);
    const evaluate = function (x) { return poly.evaluate(coef, x); };

    return {
      coef: coef,
      basis: basis,
      degree: n - 1,
      evaluate: evaluate,
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; }),
      omega: omegaAsc,
      omegaDesc: omegaDesc,
      productTable: productTable,
      divisionTable: divisionTable
    };
  }

  PPS.lagrange = { interpolate: interpolate };
})();
