/* ============================================================
   newton.js — divided differences (arbitrary nodes) and
               finite differences (equally spaced nodes)
               Updated based on Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const poly = PPS.poly;
  const { factorial, fmt } = PPS.format;

  function validate(xs, ys) {
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 2) throw new Error('tooFew');
    for (let i = 0; i < xs.length; i++) {
      if (!isFinite(xs[i]) || !isFinite(ys[i])) throw new Error('numeric');
      for (let j = i + 1; j < xs.length; j++) {
        if (Math.abs(xs[i] - xs[j]) < 1e-12) throw new Error('dupX');
      }
    }
  }

  function dividedTable(xs, ys) {
    const n = xs.length;
    const F = Array.from({ length: n }, function () { return new Array(n).fill(null); });
    for (let i = 0; i < n; i++) F[i][0] = ys[i];
    for (let j = 1; j < n; j++) {
      for (let i = j; i < n; i++) {
        F[i][j] = (F[i][j - 1] - F[i - 1][j - 1]) / (xs[i] - xs[i - j]);
      }
    }
    return F;
  }

  /**
   * Tính đạo hàm cấp r (r = 0..K) theo biến x của dạng Newton tại eval_x (Mục 4.2 C - Algorithm.md)
   * c: vector hệ số Newton
   * Z: vector mốc tương ứng (Tiến: Z = X; Lùi: Z = X đảo ngược)
   */
  function evalDerivatives(c, Z, eval_x, K) {
    const n = c.length - 1;
    K = Math.min(n, Math.max(0, Math.floor(K || 0)));
    const R = new Array(K + 1).fill(0);
    R[0] = c[n];

    for (let k = n - 1; k >= 0; k--) {
      const diff = eval_x - Z[k];
      for (let r = K; r >= 1; r--) {
        R[r] = R[r] * diff + r * R[r - 1];
      }
      R[0] = R[0] * diff + c[k];
    }
    return R;
  }

  /**
   * Nhân lồng Horner tính P(eval_x) dạng Newton (Mục 4.2 C)
   */
  function nestedEval(c, Z, eval_x) {
    const n = c.length - 1;
    let val = c[n];
    for (let k = n - 1; k >= 0; k--) {
      val = val * (eval_x - Z[k]) + c[k];
    }
    return val;
  }

  function interpolate(xs, ys) {
    validate(xs, ys);
    const n = xs.length;
    const F = dividedTable(xs, ys);

    // Vector hệ số tiến (đường chéo chính: F[j, j])
    const cForward = [];
    for (let j = 0; j < n; j++) cForward.push(F[j][j]);

    // Vector hệ số lùi (hàng cuối: F[n - 1, j])
    const cBackward = [];
    for (let j = 0; j < n; j++) cBackward.push(F[n - 1][j]);

    let P = [cForward[0]];
    let omega = [1];
    const omegas = [[1]];
    const terms = [{
      j: 0,
      coef: cForward[0],
      omega: [1],
      roots: [],
      factor: '1',
      sourceIndex: 0
    }];

    // Bảng tích (Product Table - Mục 4.3 B)
    // Các hàng Bi(x) theo lũy thừa giảm dần
    const productTable = [];
    // Hàng 0: B_0(x) = 1
    productTable.push({
      node: null,
      j: 0,
      polyAsc: [1],
      polyDesc: padPolyDesc([1], n)
    });

    for (let j = 1; j < n; j++) {
      omega = poly.mulLinear(omega, xs[j - 1]);
      const c = cForward[j];
      P = poly.add(P, poly.scale(omega, c));
      omegas.push(omega);
      terms.push({
        j: j,
        coef: c,
        omega: omega.slice(),
        roots: xs.slice(0, j),
        factor: factorString(xs.slice(0, j)),
        sourceIndex: j
      });
      productTable.push({
        node: xs[j - 1],
        j: j,
        polyAsc: omega.slice(),
        polyDesc: padPolyDesc(omega.slice().reverse(), n)
      });
    }

    return {
      table: F,
      coefs: cForward,
      cForward: cForward,
      cBackward: cBackward,
      omegas: omegas,
      terms: terms,
      productTable: productTable,
      coef: poly.trim(P),
      degree: n - 1,
      nodes: xs.map(function (x, i) { return { x: x, y: ys[i] }; }),
      evaluate: function (x) { return nestedEval(cForward, xs, x); },
      evalDerivatives: function (x, K) { return evalDerivatives(cForward, xs, x, K); }
    };
  }

  function padPolyDesc(desc, totalDegree) {
    const res = new Array(totalDegree).fill(0);
    const start = totalDegree - desc.length;
    for (let i = 0; i < desc.length; i++) {
      res[start + i] = desc[i];
    }
    return res;
  }

  function factorString(roots, varName) {
    varName = varName || 't';
    if (!roots.length) return '1';
    return roots.map(function (r) {
      if (Math.abs(r) < 1e-12) return varName;
      if (r > 0) return '(' + varName + ' − ' + fmt(r, 4) + ')';
      return '(' + varName + ' + ' + fmt(-r, 4) + ')';
    }).join('');
  }

  function differenceTable(ys) {
    const n = ys.length;
    const D = Array.from({ length: n }, function () { return new Array(n).fill(null); });
    for (let i = 0; i < n; i++) D[0][i] = ys[i];
    for (let j = 1; j < n; j++) {
      for (let i = 0; i + j < n; i++) {
        D[j][i] = D[j - 1][i + 1] - D[j - 1][i];
      }
    }
    return D;
  }

  function checkEqualSpacing(xs) {
    if (xs.length < 2) throw new Error('tooFew');
    const h = xs[1] - xs[0];
    if (!(h > 0)) throw new Error('spacing');
    let maxDev = 0;
    for (let i = 1; i < xs.length; i++) {
      const diff = (xs[i] - xs[i - 1]) - h;
      maxDev = Math.max(maxDev, Math.abs(diff));
    }
    const tolerance = 1e-8 * Math.max(1, Math.abs(h));
    if (maxDev > tolerance) throw new Error('spacing');
    return { h: h, maxDev: maxDev };
  }

  function autoMode(xs, xbar) {
    const mid = (xs[0] + xs[xs.length - 1]) / 2;
    return xbar <= mid ? 'forward' : 'backward';
  }

  function composeLinear(polyAsc, a, b) {
    if (!polyAsc.length) return [0];
    let acc = [polyAsc[polyAsc.length - 1]];
    const linear = [b, a]; // b + a*x
    for (let k = polyAsc.length - 2; k >= 0; k--) {
      acc = poly.add(poly.mul(acc, linear), [polyAsc[k]]);
    }
    return poly.trim(acc);
  }

  function computeDividedCoeffs(xs, ys) {
    const n = xs.length;
    const diffs = ys.slice();
    const coefs = [diffs[0]];
    for (let order = 1; order < n; order++) {
      for (let i = 0; i < n - order; i++) {
        diffs[i] = (diffs[i + 1] - diffs[i]) / (xs[i + order] - xs[i]);
      }
      coefs.push(diffs[0]);
    }
    return coefs;
  }

  function equalSpacing(xs, ys, xbar, mode) {
    validate(xs, ys);
    const { h } = checkEqualSpacing(xs);
    if (!isFinite(xbar)) throw new Error('xbar');
    const n = xs.length;
    const D = differenceTable(ys);

    // Vector sai phân hữu hạn tiến & lùi
    const cForward = [];
    const cBackward = [];
    for (let j = 0; j < n; j++) {
      cForward.push(D[j][0]);
      cBackward.push(D[j][n - 1 - j]);
    }

    let chosen = mode;
    if (!chosen || chosen === 'auto') chosen = autoMode(xs, xbar);

    if (chosen === 'forward') {
      const t = (xbar - xs[0]) / h;
      let P = [ys[0]];
      let omega = [1];
      const coefT = [ys[0]];
      const terms = [{
        j: 0,
        diff: ys[0],
        coef: ys[0],
        factor: '1',
        poly: [1],
        sourceIndex: 0,
        sourceLabel: 'y₀'
      }];
      const productTable = [{
        j: 0,
        factor: '1',
        factorStr: '1',
        polyAsc: [1],
        polyDesc: padPolyDesc([1], n),
        diff: ys[0],
        coef: ys[0],
        termPoly: [ys[0]]
      }];

      for (let j = 1; j < n; j++) {
        omega = poly.mulLinear(omega, j - 1);
        const diff = D[j][0];
        const c = diff / factorial(j);
        const termPoly = poly.scale(omega, c);
        P = poly.add(P, termPoly);
        coefT.push(c);
        terms.push({
          j: j,
          diff: diff,
          coef: c,
          factor: factorString(range(j)),
          poly: omega.slice(),
          sourceIndex: j,
          sourceLabel: 'Δ' + PPS.format.superscript(j) + 'y₀'
        });
        const factorTerms = ['t'];
        for (let k = 1; k < j; k++) factorTerms.push('(t − ' + k + ')');
        productTable.push({
          j: j,
          factor: '(t − ' + (j - 1) + ')',
          factorStr: factorTerms.join(''),
          polyAsc: omega.slice(),
          polyDesc: padPolyDesc(omega.slice().reverse(), n),
          diff: diff,
          coef: c,
          termPoly: termPoly
        });
      }
      const ev = poly.evalWithDerivative(P, t);
      const polyTrimmedT = poly.trim(P);
      const polyX = composeLinear(polyTrimmedT, 1 / h, -xs[0] / h);
      const cDivided = computeDividedCoeffs(xs, ys);

      return {
        mode: 'forward',
        h: h,
        t: t,
        table: D,
        coefT: coefT,
        polyT: polyTrimmedT,
        polyX: polyX,
        terms: terms,
        cForward: cForward,
        cBackward: cBackward,
        cNormalized: coefT,
        cDivided: cDivided,
        productTable: productTable,
        value: ev.value,
        derivative: ev.derivative / h,
        anchorIndex: 0,
        anchorX: xs[0],
        usedIndices: range(n),
        note: 'forward'
      };
    }

    // backward: anchor at x_n
    const anchor = n - 1;
    const t = (xbar - xs[anchor]) / h;
    let P = [ys[anchor]];
    let omega = [1];
    const coefT = [ys[anchor]];
    const terms = [{
      j: 0,
      diff: ys[anchor],
      coef: ys[anchor],
      factor: '1',
      poly: [1],
      sourceIndex: anchor,
      sourceLabel: 'yₙ'
    }];
    const productTable = [{
      j: 0,
      factor: '1',
      factorStr: '1',
      polyAsc: [1],
      polyDesc: padPolyDesc([1], n),
      diff: ys[anchor],
      coef: ys[anchor],
      termPoly: [ys[anchor]]
    }];

    for (let j = 1; j < n; j++) {
      omega = poly.mulLinear(omega, -(j - 1));
      const diff = D[j][anchor - j];
      const c = diff / factorial(j);
      const termPoly = poly.scale(omega, c);
      P = poly.add(P, termPoly);
      coefT.push(c);
      terms.push({
        j: j,
        diff: diff,
        coef: c,
        factor: factorString(range(j).map(function (s) { return -s; })),
        poly: omega.slice(),
        sourceIndex: anchor - j,
        sourceLabel: '∇' + PPS.format.superscript(j) + 'yₙ'
      });
      const factorTerms = ['t'];
      for (let k = 1; k < j; k++) factorTerms.push('(t + ' + k + ')');
      productTable.push({
        j: j,
        factor: '(t + ' + (j - 1) + ')',
        factorStr: factorTerms.join(''),
        polyAsc: omega.slice(),
        polyDesc: padPolyDesc(omega.slice().reverse(), n),
        diff: diff,
        coef: c,
        termPoly: termPoly
      });
    }
    const ev = poly.evalWithDerivative(P, t);
    const polyTrimmedT = poly.trim(P);
    const polyX = composeLinear(polyTrimmedT, 1 / h, -xs[anchor] / h);
    const cDivided = computeDividedCoeffs(xs.slice().reverse(), ys.slice().reverse());

    return {
      mode: 'backward',
      h: h,
      t: t,
      table: D,
      coefT: coefT,
      polyT: polyTrimmedT,
      polyX: polyX,
      terms: terms,
      cForward: cForward,
      cBackward: cBackward,
      cNormalized: coefT,
      cDivided: cDivided,
      productTable: productTable,
      value: ev.value,
      derivative: ev.derivative / h,
      anchorIndex: anchor,
      anchorX: xs[anchor],
      usedIndices: range(n).reverse(),
      note: 'backward'
    };
  }

  function range(n) {
    const out = [];
    for (let i = 0; i < n; i++) out.push(i);
    return out;
  }

  PPS.newton = {
    validate: validate,
    dividedTable: dividedTable,
    interpolate: interpolate,
    differenceTable: differenceTable,
    checkEqualSpacing: checkEqualSpacing,
    equalSpacing: equalSpacing,
    factorString: factorString,
    autoMode: autoMode,
    nestedEval: nestedEval,
    evalDerivatives: evalDerivatives
  };
})();
