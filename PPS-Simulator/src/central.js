/* ============================================================
   central.js — Gauss I, Gauss II, Stirling, Bessel
   central interpolation and derivatives
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const poly = PPS.poly;
  const { factorial, fmt } = PPS.format;

  const SUB = {
    '-': '₋', '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
    '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'
  };

  function sub(n) {
    return String(n).split('').map(function (ch) { return SUB[ch] || ch; }).join('');
  }

  function padPolyDesc(desc, totalDegree) {
    const res = new Array(totalDegree).fill(0);
    const start = totalDegree - desc.length;
    for (let i = 0; i < desc.length; i++) {
      res[start + i] = desc[i];
    }
    return res;
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

  function diffTable(values) {
    const size = values.length;
    const D = Array.from({ length: size }, function () { return new Array(size).fill(null); });
    for (let i = 0; i < size; i++) D[0][i] = values[i];
    for (let j = 1; j < size; j++) {
      for (let i = 0; i + j < size; i++) {
        D[j][i] = D[j - 1][i + 1] - D[j - 1][i];
      }
    }
    return D;
  }

  function range(from, to) {
    if (to == null) { to = from; from = 0; }
    const out = [];
    for (let i = from; i < to; i++) out.push(i);
    return out;
  }

  function productRoots(roots) {
    let p = [1];
    for (const r of roots) p = poly.mulLinear(p, r);
    return p;
  }

  function productSquares(ks) {
    let p = [1];
    for (const k of ks) p = poly.mulQuadratic(p, k * k);
    return p;
  }

  function productQuarters(i) {
    let p = [1];
    for (let k = 1; k <= i; k++) {
      const r = (2 * k - 1);
      p = poly.mulQuadratic(p, (r * r) / 4);
    }
    return p;
  }

  // Factor string for Bessel: (u² − 1/4)(u² − 9/4)...
  function squareFactorString(i, varName) {
    if (!i) return '1';
    const parts = [];
    for (let k = 1; k <= i; k++) {
      const r = (2 * k - 1);
      parts.push(varName + '² − ' + fmt((r * r) / 4, 4));
    }
    return '(' + parts.join(')(') + ')';
  }

  // Factor string for Stirling: (t² − 1)(t² − 4)...
  function stirlingSquareFactorString(i, varName) {
    if (!i) return '1';
    const parts = [];
    for (let k = 1; k <= i; k++) {
      parts.push(varName + '² − ' + (k * k));
    }
    return '(' + parts.join(')(') + ')';
  }

  function intFactorString(roots, varName) {
    if (!roots.length) return '1';
    return roots.map(function (r) {
      if (Math.abs(r) < 1e-12) return varName;
      return r > 0 ? '(' + varName + ' − ' + fmt(r, 4) + ')' : '(' + varName + ' + ' + fmt(-r, 4) + ')';
    }).join('');
  }

  function nearestIndex(xs, xbar) {
    let best = 0;
    let bestDist = Infinity;
    for (let i = 0; i < xs.length; i++) {
      const d = Math.abs(xs[i] - xbar);
      if (d < bestDist) { bestDist = d; best = i; }
    }
    return best;
  }

  function containingInterval(xs, xbar) {
    const n = xs.length;
    if (xbar <= xs[0]) return 0;
    if (xbar >= xs[n - 1]) return n - 2;
    for (let i = 0; i < n - 1; i++) {
      if (xs[i] <= xbar && xbar <= xs[i + 1]) return i;
    }
    return n - 2;
  }

  function suggestFormula(xs, xbar, h) {
    const center = nearestIndex(xs, xbar);
    const tCenter = (xbar - xs[center]) / h;
    if (Math.abs(tCenter) <= 0.25) return { formula: 'stirling', reason: 'stirling', tCenter: tCenter };
    if (xbar >= xs[0] && xbar <= xs[xs.length - 1]) {
      return { formula: 'bessel', reason: 'bessel', tCenter: tCenter };
    }
    return { formula: tCenter >= 0 ? 'gauss1' : 'gauss2', reason: 'gauss', tCenter: tCenter };
  }

  function computeBlock(xs, ys, start, size) {
    const values = ys.slice(start, start + size);
    return { start: start, size: size, values: values, table: diffTable(values) };
  }

  function makeColumns(xs, start, size, offset) {
    const cols = [];
    for (let local = 0; local < size; local++) {
      cols.push({
        local: local,
        rel: start + local - offset,
        x: xs[start + local]
      });
    }
    return cols;
  }

  function gaussTerms(block, N, formula, centerOffset) {
    const D = block.table;
    const terms = [];
    const cells = [];
    const maxOrder = 2 * N;
    for (let j = 0; j <= maxOrder; j++) {
      let local;
      let r;
      if (j === 0) {
        local = N;
        r = 0;
      } else if (formula === 'gauss1') {
        r = Math.floor(j / 2);
        local = N - r;
      } else {
        r = Math.ceil(j / 2);
        local = N - r;
      }
      if (local < 0 || local > 2 * N) continue;
      const value = D[j][local];
      if (value == null) continue;

      let roots = [];
      if (j > 0) {
        const low = -Math.floor((j - 1) / 2);
        const high = Math.ceil((j - 1) / 2);
        for (let s = low; s <= high; s++) roots.push(s);
      }
      const factorPoly = productRoots(roots);
      const coef = value / factorial(j);
      terms.push({
        j: j,
        coef: coef,
        factorPoly: factorPoly,
        factorString: intFactorString(roots, 't'),
        sources: [{ j: j, local: local, rel: local - N, label: 'Δ' + PPS.format.superscript(j) + 'y' + sub(local - N) }],
        poly: poly.scale(factorPoly, coef)
      });
      cells.push({ j: j, local: local });
    }
    return { terms: terms, cells: cells };
  }

  function stirlingTerms(block, N) {
    const D = block.table;
    const terms = [];
    const cells = [];
    for (let i = 0; i <= N; i++) {
      const j = 2 * i;
      const local = N - i;
      const value = D[j][local];
      if (value == null) continue;
      const factorPoly = productSquares(range(i));
      const coef = value / factorial(j);
      terms.push({
        j: j,
        coef: coef,
        factorPoly: factorPoly,
        factorString: i === 0 ? '1' : 't²' + (i > 1 ? stirlingSquareFactorString(i - 1, 't').replace(/^\(/, '·(') : ''),
        sources: [{ j: j, local: local, rel: -i, label: 'Δ' + PPS.format.superscript(j) + 'y' + sub(-i) }],
        poly: poly.scale(factorPoly, coef)
      });
      cells.push({ j: j, local: local });
    }
    for (let i = 1; i <= N; i++) {
      const j = 2 * i - 1;
      const a = N - i;
      const b = N - i + 1;
      const v1 = D[j][a];
      const v2 = D[j][b];
      if (v1 == null || v2 == null) continue;
      let factorPoly = productSquares(range(1, i));
      factorPoly = poly.mulLinear(factorPoly, 0); // multiply by t
      const coef = (v1 + v2) / (2 * factorial(j));
      terms.push({
        j: j,
        coef: coef,
        factorPoly: factorPoly,
        factorString: i === 1 ? 't' : 't' + (i > 1 ? stirlingSquareFactorString(i - 1, 't').replace(/^\(/, '·(') : ''),
        sources: [
          { j: j, local: a, rel: -i, label: 'Δ' + PPS.format.superscript(j) + 'y' + sub(-i) },
          { j: j, local: b, rel: -(i - 1), label: 'Δ' + PPS.format.superscript(j) + 'y' + sub(-(i - 1)) }
        ],
        averaged: true,
        poly: poly.scale(factorPoly, coef)
      });
      cells.push({ j: j, local: a });
      cells.push({ j: j, local: b });
    }
    return { terms: terms.sort(function (a, b) { return a.j - b.j; }), cells: cells };
  }

  function besselTerms(block, N) {
    const D = block.table;
    const terms = [];
    const cells = [];
    const zeroLocal = N - 1;
    for (let i = 0; i <= N - 1; i++) {
      const j = 2 * i;
      const a = zeroLocal - i;
      const b = zeroLocal - i + 1;
      const v1 = D[j][a];
      const v2 = D[j][b];
      if (v1 == null || v2 == null) continue;
      const factorPoly = productQuarters(i);
      const coef = (v1 + v2) / (2 * factorial(j));
      terms.push({
        j: j,
        coef: coef,
        factorPoly: factorPoly,
        factorString: squareFactorString(i, 'u'),
        sources: [
          { j: j, local: a, rel: a - zeroLocal, label: 'Δ' + PPS.format.superscript(j) + 'y' + sub(a - zeroLocal) },
          { j: j, local: b, rel: b - zeroLocal, label: 'Δ' + PPS.format.superscript(j) + 'y' + sub(b - zeroLocal) }
        ],
        averaged: true,
        poly: poly.scale(factorPoly, coef)
      });
      cells.push({ j: j, local: a });
      cells.push({ j: j, local: b });
    }
    for (let i = 0; i <= N - 1; i++) {
      const j = 2 * i + 1;
      const a = zeroLocal - i;
      const value = D[j][a];
      if (value == null) continue;
      let factorPoly = productQuarters(i);
      factorPoly = poly.mulLinear(factorPoly, 0);
      const coef = value / factorial(j);
      terms.push({
        j: j,
        coef: coef,
        factorPoly: factorPoly,
        factorString: i === 0 ? 'u' : 'u' + squareFactorString(i, 'u').replace(/^\(/, '·('),
        sources: [{ j: j, local: a, rel: -i, label: 'Δ' + PPS.format.superscript(j) + 'y' + sub(-i) }],
        poly: poly.scale(factorPoly, coef)
      });
      cells.push({ j: j, local: a });
    }
    return { terms: terms.sort(function (a, b) { return a.j - b.j; }), cells: cells };
  }

  /**
   * Main entry.
   * opts: { formula: 'auto'|'gauss1'|'gauss2'|'stirling'|'bessel', N: number|'auto' }
   */
  function interpolate(xs, ys, xbar, opts) {
    opts = opts || {};
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 3) throw new Error('central');
    if (!isFinite(xbar)) throw new Error('xbar');
    const { h } = PPS.newton.checkEqualSpacing(xs);
    for (let i = 0; i < xs.length; i++) {
      if (!isFinite(ys[i])) throw new Error('numeric');
    }

    const requested = opts.formula || 'auto';
    const suggestion = suggestFormula(xs, xbar, h);
    let formula = requested === 'auto' ? suggestion.formula : requested;
    if (['gauss1', 'gauss2', 'stirling', 'bessel'].indexOf(formula) === -1) formula = suggestion.formula;

    const n = xs.length;
    let centerIndex, centerX, t, u, block, N, offset, maxN, clamped = false;

    if (formula === 'bessel') {
      const k = containingInterval(xs, xbar);
      centerIndex = k;
      centerX = xs[k];
      t = (xbar - centerX) / h;
      u = t - 0.5;
      maxN = Math.min(k + 1, n - 1 - k);
      if (maxN < 1) throw new Error('central');
      const req = opts.N === 'auto' || opts.N == null || opts.N === '' ? Math.min(maxN, 3) : Math.floor(Number(opts.N));
      N = Math.max(1, Math.min(maxN, isFinite(req) ? req : maxN));
      if (req > maxN) clamped = true;
      const start = k - (N - 1);
      block = computeBlock(xs, ys, start, 2 * N);
      offset = centerIndex;
      const built = besselTerms(block, N);
      const P = built.terms.reduce(function (acc, term) { return poly.add(acc, term.poly); }, [0]);
      const ev = poly.evalWithDerivative(P, u);
      const P_trimmed = poly.trim(P);
      const polyX = composeLinear(P_trimmed, 1 / h, -(centerX / h + 0.5));
      const nodesZ = xs.slice(start, start + 2 * N);
      const cVector = built.terms.map(function (t) { return t.coef * factorial(t.j); });
      const dVector = built.terms.map(function (t) { return t.coef; });
      const productTable = built.terms.map(function (term) {
        return {
          j: term.j,
          factorStr: term.factorString,
          polyAsc: term.factorPoly.slice(),
          polyDesc: padPolyDesc(term.factorPoly.slice().reverse(), 2 * N),
          diff: term.coef * factorial(term.j),
          coef: term.coef,
          termPoly: term.poly
        };
      });

      return {
        formula: formula,
        requested: requested,
        autoReason: requested === 'auto' ? suggestion.reason : null,
        centerIndex: centerIndex,
        centerX: centerX,
        h: h,
        t: t,
        u: u,
        N: N,
        maxN: maxN,
        clamped: clamped,
        blockStart: start,
        blockSize: 2 * N,
        blockValues: block.values,
        table: block.table,
        columns: makeColumns(xs, start, 2 * N, offset),
        terms: built.terms,
        highlight: built.cells,
        polyT: P_trimmed,
        polyX: polyX,
        nodesZ: nodesZ,
        cVector: cVector,
        dVector: dVector,
        productTable: productTable,
        variable: 'u',
        value: ev.value,
        derivative: ev.derivative / h,
        usedIndices: range(2 * N).map(function (i) { return start + i; }),
        evalAt: function (x) {
          return poly.evaluate(P, (x - centerX) / h - 0.5);
        }
      };
    }

    // Gauss I / Gauss II / Stirling: 2N+1 nodes around nearest node
    let center = nearestIndex(xs, xbar);
    if (center < 1) center = 1;
    if (center > n - 2) center = n - 2;
    centerIndex = center;
    centerX = xs[center];
    t = (xbar - centerX) / h;
    maxN = Math.min(center, n - 1 - center);
    if (maxN < 1) throw new Error('central');
    const req = opts.N === 'auto' || opts.N == null || opts.N === '' ? Math.min(maxN, 3) : Math.floor(Number(opts.N));
    N = Math.max(0, Math.min(maxN, isFinite(req) ? req : maxN));
    if (req > maxN) clamped = true;
    const start = center - N;
    block = computeBlock(xs, ys, start, 2 * N + 1);
    offset = center;

    let built;
    if (formula === 'stirling') built = stirlingTerms(block, N);
    else built = gaussTerms(block, N, formula, center);

    const P = built.terms.reduce(function (acc, term) { return poly.add(acc, term.poly); }, [0]);
    const ev = poly.evalWithDerivative(P, t);
    const P_trimmed = poly.trim(P);
    const polyX = composeLinear(P_trimmed, 1 / h, -centerX / h);
    const nodesZ = xs.slice(start, start + 2 * N + 1);
    const cVector = built.terms.map(function (t) { return t.coef * factorial(t.j); });
    const dVector = built.terms.map(function (t) { return t.coef; });
    const productTable = built.terms.map(function (term) {
      return {
        j: term.j,
        factorStr: term.factorString,
        polyAsc: term.factorPoly.slice(),
        polyDesc: padPolyDesc(term.factorPoly.slice().reverse(), 2 * N + 1),
        diff: term.coef * factorial(term.j),
        coef: term.coef,
        termPoly: term.poly
      };
    });

    return {
      formula: formula,
      requested: requested,
      autoReason: requested === 'auto' ? suggestion.reason : null,
      centerIndex: centerIndex,
      centerX: centerX,
      h: h,
      t: t,
      u: null,
      N: N,
      maxN: maxN,
      clamped: clamped,
      blockStart: start,
      blockSize: 2 * N + 1,
      blockValues: block.values,
      table: block.table,
      columns: makeColumns(xs, start, 2 * N + 1, offset),
      terms: built.terms,
      highlight: built.cells,
      polyT: P_trimmed,
      polyX: polyX,
      nodesZ: nodesZ,
      cVector: cVector,
      dVector: dVector,
      productTable: productTable,
      variable: 't',
      value: ev.value,
      derivative: ev.derivative / h,
      usedIndices: range(2 * N + 1).map(function (i) { return start + i; }),
      evalAt: function (x) {
        return poly.evaluate(P, (x - centerX) / h);
      }
    };
  }

  PPS.central = {
    interpolate: interpolate,
    suggestFormula: suggestFormula,
    diffTable: diffTable,
    sub: sub,
    squareFactorString: squareFactorString,
    stirlingSquareFactorString: stirlingSquareFactorString
  };
})();
