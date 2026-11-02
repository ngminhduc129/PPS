/* ============================================================
   inverse.js — inverse interpolation:
   Method 1: Inverse Function Interpolation x = P(y)
   Method 2: Fixed-point iteration loop t = phi(t)
   Updated strictly based on Section 7 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const poly = PPS.poly;
  const { factorial } = PPS.format;

  function padPolyDesc(desc, totalDegree) {
    const res = new Array(totalDegree).fill(0);
    const start = totalDegree - desc.length;
    for (let i = 0; i < desc.length; i++) {
      res[start + i] = desc[i];
    }
    return res;
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

  /**
   * Quét các khoảng đơn điệu từ dữ liệu (Mục 7.2.2 Bước 1 - Algorithm.md)
   */
  function scanMonotonicIntervals(xs, ys) {
    const n = xs.length;
    if (n < 2) return [];
    const intervals = [];
    let currDir = 0;
    let startIdx = 0;

    for (let i = 0; i < n - 1; i++) {
      const diff = ys[i + 1] - ys[i];
      const d = Math.abs(diff) < 1e-14 ? 0 : (diff > 0 ? 1 : -1);

      if (d === 0) {
        if (currDir !== 0) intervals.push({ start: startIdx, end: i, dir: currDir });
        currDir = 0;
        startIdx = i + 1;
      } else if (currDir === 0) {
        startIdx = i;
        currDir = d;
      } else if (d !== currDir) {
        intervals.push({ start: startIdx, end: i, dir: currDir });
        startIdx = i;
        currDir = d;
      }
    }
    if (currDir !== 0) {
      intervals.push({ start: startIdx, end: n - 1, dir: currDir });
    }
    return intervals;
  }

  /**
   * Tìm tất cả các khoảng cách ly nghiệm [x_k, x_{k+1}] thỏa mãn (y_k - ybar)*(y_{k+1} - ybar) <= 0
   */
  function findBrackets(xs, ys, ybar) {
    const n = xs.length;
    const brackets = [];
    let yScale = Math.max(1, Math.abs(ybar));
    for (let i = 0; i < n; i++) {
      yScale = Math.max(yScale, Math.abs(ys[i]));
    }
    const yTol = 1e-12 * yScale;
    for (let i = 0; i < n; i++) {
      if (Math.abs(ys[i] - ybar) <= yTol) {
        brackets.push({ exact: true, index: i, x: xs[i], y: ys[i] });
      }
    }
    for (let k = 0; k < n - 1; k++) {
      const a = ys[k] - ybar;
      const b = ys[k + 1] - ybar;
      if (a * b < 0) {
        brackets.push({
          exact: false,
          k: k,
          x0: xs[k],
          x1: xs[k + 1],
          y0: ys[k],
          y1: ys[k + 1]
        });
      }
    }
    return brackets;
  }

  function findBracket(xs, ys, ybar) {
    const list = findBrackets(xs, ys, ybar);
    return list.length ? list[0] : null;
  }

  /**
   * Lập bảng tỷ sai phân cho hàm ngược x = P(y).
   * Sắp xếp các mốc theo y tăng dần để giảm sai số (giống Python 10_Inverse_Function).
   */
  function buildInverseTable(rawX, rawY) {
    const order = rawX.map(function (_, i) { return i; })
      .sort(function (a, b) { return rawY[a] - rawY[b]; });
    const subX = order.map(function (i) { return rawX[i]; });
    const subY = order.map(function (i) { return rawY[i]; });
    const m = subX.length;

    for (let i = 0; i < m; i++) {
      for (let j = i + 1; j < m; j++) {
        if (Math.abs(subY[i] - subY[j]) < 1e-14) throw new Error('dupY');
      }
    }

    const F = Array.from({ length: m }, function () { return new Array(m).fill(null); });
    for (let i = 0; i < m; i++) F[i][0] = subX[i];
    for (let j = 1; j < m; j++) {
      for (let i = 0; i < m - j; i++) {
        const denom = subY[i + j] - subY[i];
        if (Math.abs(denom) < 1e-14) throw new Error('dupY');
        F[i][j] = (F[i + 1][j - 1] - F[i][j - 1]) / denom;
      }
    }

    const coefs = [];
    for (let j = 0; j < m; j++) coefs.push(F[0][j]);

    const cBackward = [];
    for (let j = 0; j < m; j++) cBackward.push(F[m - 1 - j][j]);

    // Bảng tích lũy thừa B_i(y) (Product Table) theo biến y
    const productTable = [];
    let omegaY = [1];
    productTable.push({
      j: 0,
      node: null,
      factor: '1',
      polyAsc: [1],
      polyDesc: padPolyDesc([1], m)
    });

    let P_y = [coefs[0]];
    for (let j = 1; j < m; j++) {
      omegaY = poly.mulLinear(omegaY, subY[j - 1]);
      productTable.push({
        j: j,
        node: subY[j - 1],
        factor: '(y − ' + subY[j - 1] + ')',
        polyAsc: omegaY.slice(),
        polyDesc: padPolyDesc(omegaY.slice().reverse(), m)
      });
      P_y = poly.add(P_y, poly.scale(omegaY, coefs[j]));
    }

    return {
      subX: subX,
      subY: subY,
      table: F,
      coefs: coefs,
      cForward: coefs,
      cBackward: cBackward,
      productTable: productTable,
      polyCanonical: poly.trim(P_y)
    };
  }

  function hornerInverse(coefs, subY, ybar) {
    let value = coefs[coefs.length - 1];
    for (let k = coefs.length - 2; k >= 0; k--) {
      value = value * (ybar - subY[k]) + coefs[k];
    }
    return value;
  }

  function linearInverse(bracket, ybar) {
    const dy = bracket.y1 - bracket.y0;
    if (Math.abs(dy) < 1e-14) return null;
    return bracket.x0 + (ybar - bracket.y0) * (bracket.x1 - bracket.x0) / dy;
  }

  /**
   * Nội suy ngược trên một cửa sổ mốc.
   * Nếu đa thức P(y) cho nghiệm nằm ngoài khoảng cách ly (hoặc không hữu hạn),
   * dùng nội suy tuyến tính nghịch đảo trên hai đầu khoảng (fallback).
   */
  function evaluateInverseWindow(rawX, rawY, bracket, ybar) {
    let built = null;
    try {
      built = buildInverseTable(rawX, rawY);
    } catch (e) {
      built = null;
    }

    const slope = (bracket.x1 - bracket.x0) / (bracket.y1 - bracket.y0);

    if (!built) {
      const xlin = linearInverse(bracket, ybar);
      return {
        xbar: xlin,
        polyX: null,
        fallback: true,
        inBracket: xlin != null && xlin >= bracket.x0 && xlin <= bracket.x1,
        derivative: slope,
        subX: rawX,
        subY: rawY,
        tableY: null,
        coefsD: null,
        cForward: null,
        cBackward: null,
        productTable: null,
        polyCanonical: null
      };
    }

    const polyX = hornerInverse(built.coefs, built.subY, ybar);
    let xbar = polyX;
    let fallback = false;

    if (!isFinite(polyX) || polyX < bracket.x0 || polyX > bracket.x1) {
      const xlin = linearInverse(bracket, ybar);
      if (xlin != null) {
        xbar = xlin;
        fallback = true;
      }
    }

    const derivative = fallback
      ? slope
      : PPS.newton.evalDerivatives(built.coefs, built.subY, ybar, 1)[1];

    return {
      xbar: xbar,
      polyX: polyX,
      fallback: fallback,
      inBracket: xbar != null && xbar >= bracket.x0 && xbar <= bracket.x1,
      derivative: derivative,
      subX: built.subX,
      subY: built.subY,
      tableY: built.table,
      coefsD: built.coefs,
      cForward: built.cForward,
      cBackward: built.cBackward,
      productTable: built.productTable,
      polyCanonical: built.polyCanonical
    };
  }

  function inverseWindow(xs, ys, bracket, nPoints) {
    const run = monotoneRun(xs, ys, bracket.k);
    const stencil = selectStencil(run.start, run.end, bracket.k, nPoints);
    return {
      run: run,
      stencil: stencil,
      rawX: xs.slice(stencil.start, stencil.end + 1),
      rawY: ys.slice(stencil.start, stencil.end + 1)
    };
  }

  function forwardResidual(xs, ys, xbar, ybar) {
    try {
      const fwdInterp = PPS.newton.interpolate(xs, ys);
      return Math.abs(fwdInterp.evaluate(xbar) - ybar);
    } catch (e) {
      return null;
    }
  }

  /**
   * CÁCH 1: Phương pháp Nội suy Hàm Ngược x = P(y) (Mục 7.2 - Algorithm.md)
   * Đổi vai trò: Y làm biến độc lập, X làm giá trị hàm.
   * Cửa sổ mốc được chọn giống Cách 2 (nhánh đơn điệu + n_points),
   * nên ô n_points điều khiển được cả hai cách.
   */
  function solveInverseFunction(xs, ys, ybar, opts) {
    opts = opts || {};
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 2) throw new Error('tooFew');
    if (!isFinite(ybar)) throw new Error('numeric');

    const intervals = scanMonotonicIntervals(xs, ys);
    const brackets = findBrackets(xs, ys, ybar);
    if (!brackets.length) throw new Error('noBracket');

    const bracket = brackets[0];
    const nPoints = normalizePoints(opts.nPoints);

    if (bracket.exact) {
      return {
        mode: 'function',
        exact: true,
        xbar: bracket.x,
        ybar: ybar,
        residual: 0,
        derivative: null,
        intervals: intervals,
        bracket: bracket,
        nPoints: nPoints,
        stencil: null,
        fallback: false,
        inBracket: true
      };
    }

    const window = inverseWindow(xs, ys, bracket, nPoints);
    const evaluated = evaluateInverseWindow(window.rawX, window.rawY, bracket, ybar);
    const residual = forwardResidual(xs, ys, evaluated.xbar, ybar);

    return Object.assign({
      mode: 'function',
      exact: false,
      ybar: ybar,
      nPoints: nPoints,
      run: window.run,
      stencil: {
        start: window.stencil.start,
        end: window.stencil.end,
        count: window.stencil.count,
        xs: window.rawX,
        ys: window.rawY
      },
      intervals: intervals,
      bracket: bracket
    }, evaluated, { residual: residual });
  }

  // ============================================================
  // CÁCH 2: LẶP NỘI SUY NGƯỢC t = phi(t)
  // Port theo Python_Interpolation/11_Loop_Inverse:
  //   1) Cắt bracket về nhánh đơn điệu lớn nhất chứa nó.
  //   2) Chọn tối đa nPoints mốc liên tiếp quanh bracket trong nhánh.
  //   3) Lập lại bảng sai phân trên đúng cửa sổ đó.
  //   4) Chạy cả Newton tiến và Newton lùi khi method = 'auto',
  //      rồi kiểm chứng nghiệm nằm trong bracket và residual đạt ngưỡng.
  // ============================================================

  function normalizePoints(nPoints) {
    if (typeof nPoints === 'string' && nPoints.trim().toLowerCase() === 'all') return 'all';
    if (nPoints == null || nPoints === '') return 4;
    const value = Number(nPoints);
    if (!isFinite(value)) return 4;
    return Math.max(2, Math.floor(value));
  }

  /**
   * Mở rộng bracket [k, k+1] thành nhánh đơn điệu lớn nhất chứa nó.
   * Tương ứng _bracketing_run của Python 11_Loop_Inverse.
   */
  function monotoneRun(xs, ys, k) {
    const dir = Math.sign(ys[k + 1] - ys[k]);
    let start = k;
    let end = k + 1;
    while (start > 0 && Math.sign(ys[start] - ys[start - 1]) === dir) start--;
    while (end < xs.length - 1 && Math.sign(ys[end + 1] - ys[end]) === dir) end++;
    return { start: start, end: end, dir: dir };
  }

  /**
   * Chọn tối đa nPoints mốc liên tiếp quanh bracket trong nhánh đơn điệu.
   * nPoints === 'all' -> lấy cả nhánh.
   * Tương ứng _select_bracket_stencil của Python.
   */
  function selectStencil(runStart, runEnd, bracketK, nPoints) {
    const len = runEnd - runStart + 1;
    const count = nPoints === 'all' ? len : Math.min(nPoints, len);
    const localK = bracketK - runStart;
    const offset = Math.max(0, Math.min(
      localK - Math.floor((count - 2) / 2),
      len - count
    ));
    return {
      start: runStart + offset,
      end: runStart + offset + count - 1,
      count: count
    };
  }

  function buildAttemptOptions(opts, ybar, h, yScale) {
    const tol = opts.tolT || opts.eps || 1e-8;
    return {
      ybar: ybar,
      h: h,
      yScale: yScale,
      tolT: tol,
      tolRes: opts.tolRes || 1e-8,
      tolX: opts.tolX || 1e-8,
      maxIter: opts.maxIter ? Math.max(1, Math.floor(opts.maxIter)) : 100,
      qCoeff: opts.qCoeff != null && opts.qCoeff > 0 && opts.qCoeff < 1 ? Number(opts.qCoeff) : null
    };
  }

  function failedAttempt(isForward, reason, attemptOpts, bracket) {
    return {
      method: isForward ? 'forward' : 'backward',
      valid: false,
      reason: reason || 'numeric',
      status: 'DIVERGED',
      converged: false,
      xbar: null,
      t: null,
      t0: null,
      h: attemptOpts.h,
      history: [],
      bracket: bracket,
      anchorIndex: null,
      anchorX: null,
      residual: null,
      checkY: null,
      iterations: 0,
      stencil: null,
      diffTable: null,
      c1: null,
      baseY: null,
      cDiff: null,
      cNormalized: null,
      dividedNodes: null,
      cDivided: null,
      productTable: null,
      polyT: null
    };
  }

  /**
   * Chạy một hướng lặp (Newton tiến hoặc lùi) trên cửa sổ mốc đã chọn.
   */
  function runIterationAttempt(xs, ys, bracket, stencil, isForward, attemptOpts) {
    const method = isForward ? 'forward' : 'backward';
    const subX = xs.slice(stencil.start, stencil.end + 1);
    const subY = ys.slice(stencil.start, stencil.end + 1);
    const m = subX.length;

    if (m < 2) return failedAttempt(isForward, 'tooFew', attemptOpts, bracket);

    const D = PPS.newton.differenceTable(subY);
    const c1 = isForward ? D[1][0] : D[1][m - 2];
    const baseY = isForward ? subY[0] : subY[m - 1];
    const xAnchor = isForward ? subX[0] : subX[m - 1];
    const anchorIndex = isForward ? stencil.start : stencil.end;
    const yScale = attemptOpts.yScale;

    if (!isFinite(c1) || Math.abs(c1) < 1e-14) {
      const failed = failedAttempt(isForward, 'ZERO_DIFF', attemptOpts, bracket);
      failed.stencil = { start: stencil.start, end: stencil.end, xs: subX, ys: subY };
      failed.diffTable = D;
      failed.c1 = c1;
      failed.baseY = baseY;
      return failed;
    }

    const illConditioned = Math.abs(c1) < 1e-10 * Math.max(1, yScale);
    const t0 = (attemptOpts.ybar - baseY) / c1;

    // Vector hệ số sai phân hữu hạn & chuẩn hóa D_i = C_i / i!
    const cDiff = [];
    const cNormalized = [];
    for (let j = 0; j < m; j++) {
      const val = isForward ? D[j][0] : D[j][m - 1 - j];
      cDiff.push(val);
      cNormalized.push(val / factorial(j));
    }

    // Vector hệ số tỷ sai phân Newton (Mục 3 Python loop_inverse)
    const dividedNodes = isForward ? subX.slice() : subX.slice().reverse();
    const dividedYs = isForward ? subY.slice() : subY.slice().reverse();
    const cDivided = computeDividedCoeffs(dividedNodes, dividedYs);

    // Bảng tích lũy thừa B_i(t) theo biến chuẩn hóa t
    const productTable = [];
    let omegaT = [1];
    productTable.push({
      j: 0,
      node: null,
      factor: '1',
      factorStr: '1',
      polyAsc: [1],
      polyDesc: padPolyDesc([1], m),
      diff: cDiff[0],
      coef: cNormalized[0],
      termPoly: [cNormalized[0]]
    });

    let P_t = [cNormalized[0]];
    for (let j = 1; j < m; j++) {
      const root = isForward ? (j - 1) : -(j - 1);
      omegaT = poly.mulLinear(omegaT, root);
      const signStr = isForward ? ' − ' : ' + ';
      const factor = '(t' + signStr + (j - 1) + ')';
      const factorTerms = ['t'];
      for (let k = 1; k < j; k++) {
        factorTerms.push('(t' + (isForward ? ' − ' : ' + ') + k + ')');
      }
      const factorStr = factorTerms.join('');
      const termPoly = poly.scale(omegaT, cNormalized[j]);
      P_t = poly.add(P_t, termPoly);

      productTable.push({
        j: j,
        node: root,
        factor: factor,
        factorStr: factorStr,
        polyAsc: omegaT.slice(),
        polyDesc: padPolyDesc(omegaT.slice().reverse(), m),
        diff: cDiff[j],
        coef: cNormalized[j],
        termPoly: termPoly
      });
    }
    const polyT = poly.trim(P_t);

    function phi(tVal) {
      let correction = 0;
      for (let r = 2; r < m; r++) {
        const delta = isForward ? D[r][0] : D[r][m - 1 - r];
        let prod = 1;
        for (let j = 0; j < r; j++) {
          prod *= isForward ? (tVal - j) : (tVal + j);
        }
        correction += (delta / (factorial(r) * c1)) * prod;
      }
      return t0 - correction;
    }

    function evalPoly(tVal) {
      let value = baseY;
      for (let r = 1; r < m; r++) {
        const delta = isForward ? D[r][0] : D[r][m - 1 - r];
        let prod = 1;
        for (let j = 0; j < r; j++) {
          prod *= isForward ? (tVal - j) : (tVal + j);
        }
        value += (delta / factorial(r)) * prod;
      }
      return value;
    }

    const xStart = xAnchor + t0 * attemptOpts.h;
    const history = [{
      step: 0,
      t: t0,
      x: xStart,
      diffT: null,
      diffX: null,
      residual: Math.abs(evalPoly(t0) - attemptOpts.ybar),
      inBracket: xStart >= bracket.x0 && xStart <= bracket.x1,
      conclusion: 'init'
    }];

    let tCurr = t0;
    let status = 'MAX_ITER';
    let boundX = null;
    let iterations = 0;

    for (let step = 1; step <= attemptOpts.maxIter; step++) {
      const tPrev = tCurr;
      const tNext = phi(tPrev);

      if (!isFinite(tNext)) {
        history.push({
          step: step,
          t: tNext,
          x: null,
          diffT: null,
          diffX: null,
          residual: null,
          inBracket: false,
          conclusion: 'diverged'
        });
        tCurr = tNext;
        status = 'DIVERGED';
        iterations = step;
        break;
      }

      tCurr = tNext;
      iterations = step;
      const diffT = Math.abs(tCurr - tPrev);
      const xCurr = xAnchor + tCurr * attemptOpts.h;
      const diffX = Math.abs(attemptOpts.h) * diffT;
      const residual = Math.abs(evalPoly(tCurr) - attemptOpts.ybar);
      const inBracket = xCurr >= bracket.x0 && xCurr <= bracket.x1;

      let conclusion = 'continue';
      if (diffT < attemptOpts.tolT && residual < attemptOpts.tolRes && inBracket) {
        if (attemptOpts.qCoeff != null) {
          boundX = (attemptOpts.qCoeff / (1 - attemptOpts.qCoeff)) * diffX;
          const lo = xCurr - boundX;
          const hi = xCurr + boundX;
          if (boundX < attemptOpts.tolX && lo >= bracket.x0 && hi <= bracket.x1) {
            status = 'CONVERGED_CERTIFIED';
            conclusion = 'certified';
            history.push({ step, t: tCurr, x: xCurr, diffT, diffX, residual, inBracket, conclusion });
            break;
          }
        }
        status = 'CANDIDATE_NUMERIC';
        conclusion = 'candidate';
        history.push({ step, t: tCurr, x: xCurr, diffT, diffX, residual, inBracket, conclusion });
        break;
      }

      history.push({ step, t: tCurr, x: xCurr, diffT, diffX, residual, inBracket, conclusion });
    }

    const finiteFinal = isFinite(tCurr);
    const finalX = finiteFinal ? xAnchor + tCurr * attemptOpts.h : null;
    const finalInBracket = finalX != null && finalX >= bracket.x0 && finalX <= bracket.x1;
    if (status === 'MAX_ITER' && !finalInBracket) status = 'OUT_OF_BRACKET';

    const finalResidual = finiteFinal ? Math.abs(evalPoly(tCurr) - attemptOpts.ybar) : null;
    const checkY = finiteFinal ? evalPoly(tCurr) : null;
    const converged = status === 'CONVERGED_CERTIFIED' || status === 'CANDIDATE_NUMERIC';
    const residualTol = Math.max(10 * attemptOpts.tolRes, 1e-10) * yScale;
    const valid = converged && finalX != null && finalInBracket &&
      finalResidual != null && finalResidual <= residualTol;

    return {
      method: method,
      valid: valid,
      reason: valid ? null : status,
      status: status,
      converged: converged,
      xbar: finalX,
      t: finiteFinal ? tCurr : null,
      t0: t0,
      h: attemptOpts.h,
      history: history,
      bracket: bracket,
      anchorIndex: anchorIndex,
      anchorX: xAnchor,
      residual: finalResidual,
      checkY: checkY,
      iterations: iterations,
      stencil: { start: stencil.start, end: stencil.end, xs: subX, ys: subY },
      diffTable: D,
      c1: c1,
      baseY: baseY,
      m: m,
      illConditioned: illConditioned,
      cDiff: cDiff,
      cNormalized: cNormalized,
      dividedNodes: dividedNodes,
      cDivided: cDivided,
      productTable: productTable,
      polyT: polyT
    };
  }

  function chooseAttempt(attempts) {
    const finite = attempts.filter(function (a) { return a.xbar != null && isFinite(a.xbar); });
    const valid = finite.filter(function (a) { return a.valid; });
    const converged = finite.filter(function (a) { return a.converged; });
    const pool = valid.length ? valid : converged;
    if (pool.length) {
      return pool.slice().sort(function (a, b) {
        const ra = a.residual == null ? Infinity : a.residual;
        const rb = b.residual == null ? Infinity : b.residual;
        if (ra !== rb) return ra - rb;
        return a.iterations - b.iterations;
      })[0];
    }
    if (finite.length) {
      return finite.slice().sort(function (a, b) { return a.iterations - b.iterations; })[0];
    }
    return attempts[0];
  }

  function runDirections(xs, ys, bracket, stencil, directions, attemptOpts) {
    const attempts = [];
    for (let i = 0; i < directions.length; i++) {
      const isForward = directions[i];
      try {
        attempts.push(runIterationAttempt(xs, ys, bracket, stencil, isForward, attemptOpts));
      } catch (error) {
        attempts.push(failedAttempt(
          isForward,
          error && error.message ? error.message : 'numeric',
          attemptOpts,
          bracket
        ));
      }
    }
    return attempts;
  }

  function yScaleOf(ys, ybar) {
    let maxAbsY = 0;
    for (let i = 0; i < ys.length; i++) {
      if (Math.abs(ys[i]) > maxAbsY) maxAbsY = Math.abs(ys[i]);
    }
    return Math.max(1, maxAbsY, Math.abs(ybar));
  }

  function directionsFor(method) {
    if (method === 'forward') return [true];
    if (method === 'backward') return [false];
    return [true, false];
  }

  /**
   * CÁCH 2: Phương pháp Lặp Nội suy Ngược (Điểm bất động t = phi(t)) (Mục 7.3 - Algorithm.md)
   * Bảng Trace chi tiết từng bước, kiểm tra điều kiện dừng và trạng thái hội tụ.
   */
  function solveIteration(xs, ys, ybar, opts) {
    opts = opts || {};
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 2) throw new Error('tooFew');
    if (!isFinite(ybar)) throw new Error('numeric');

    const spacing = PPS.newton.checkEqualSpacing(xs);
    const h = spacing.h;
    const yScale = yScaleOf(ys, ybar);

    const brackets = findBrackets(xs, ys, ybar);
    if (!brackets.length) throw new Error('noBracket');
    const bracket = brackets[0];

    if (bracket.exact) {
      return {
        mode: 'iteration',
        exact: true,
        xbar: bracket.x,
        nodeIndex: bracket.index,
        residual: 0,
        checkY: ybar,
        history: [],
        status: 'CONVERGED_CERTIFIED',
        converged: true,
        iterations: 0,
        h: h,
        method: null,
        bracket: bracket,
        stencil: null,
        attempts: []
      };
    }

    const nPoints = normalizePoints(opts.nPoints);
    const run = monotoneRun(xs, ys, bracket.k);
    const stencil = selectStencil(run.start, run.end, bracket.k, nPoints);
    const attemptOpts = buildAttemptOptions(opts, ybar, h, yScale);
    const attempts = runDirections(
      xs, ys, bracket, stencil,
      directionsFor(opts.method || 'auto'),
      attemptOpts
    );
    const best = chooseAttempt(attempts);

    return Object.assign({
      mode: 'iteration',
      exact: false,
      nPoints: nPoints,
      run: run,
      stencil: stencil,
      attempts: attempts
    }, best);
  }

  // Wrapper hàm solve tương thích ngược
  function solve(xs, ys, ybar, opts) {
    opts = opts || {};
    if (opts.mode === 'function') {
      return solveInverseFunction(xs, ys, ybar, opts);
    }
    return solveIteration(xs, ys, ybar, opts);
  }

  // ============================================================
  // TÍNH TẤT CẢ CÁC KHOẢNG CÁCH LY NGHIỆM (MỖI KHOẢNG MỘT NGHIỆM)
  // ============================================================

  function solveAllInverseFunction(xs, ys, ybar, opts) {
    opts = opts || {};
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 2) throw new Error('tooFew');
    if (!isFinite(ybar)) throw new Error('numeric');

    const intervals = scanMonotonicIntervals(xs, ys);
    const brackets = findBrackets(xs, ys, ybar);
    if (!brackets.length) throw new Error('noBracket');

    const nPoints = normalizePoints(opts.nPoints);
    const results = [];

    for (const bracket of brackets) {
      if (bracket.exact) {
        results.push({
          mode: 'function',
          exact: true,
          bracket: bracket,
          xbar: bracket.x,
          ybar: ybar,
          residual: 0,
          derivative: null,
          inBracket: true,
          intervals: intervals,
          nPoints: nPoints,
          stencil: null,
          fallback: false
        });
        continue;
      }

      const window = inverseWindow(xs, ys, bracket, nPoints);
      const evaluated = evaluateInverseWindow(window.rawX, window.rawY, bracket, ybar);
      const residual = forwardResidual(xs, ys, evaluated.xbar, ybar);

      const bracketStr = '[' + PPS.format.fmt(bracket.x0, 3) + ', ' + PPS.format.fmt(bracket.x1, 3) + ']';
      const stencilStr = '[' + PPS.format.fmt(window.rawX[0], 3) + ', ' + PPS.format.fmt(window.rawX[window.rawX.length - 1], 3) + ']';
      const comparisons = [{
        rootName: 'x*',
        bracketStr: bracketStr,
        stencilStr: stencilStr,
        stencil: window.stencil,
        methodName: 'Nội suy hàm ngược x = P(y)',
        xbar: evaluated.xbar,
        result: evaluated
      }];

      results.push(Object.assign({
        mode: 'function',
        exact: false,
        bracket: bracket,
        ybar: ybar,
        nPoints: nPoints,
        run: window.run,
        stencil: {
          start: window.stencil.start,
          end: window.stencil.end,
          count: window.stencil.count,
          xs: window.rawX,
          ys: window.rawY
        },
        intervals: intervals,
        comparisons: comparisons
      }, evaluated, { residual: residual }));
    }

    return results;
  }

  function buildStencilComparison(xs, ys, bracket, run, nPoints, attemptOpts) {
    if (!bracket || bracket.exact) return [];
    const k = bracket.k;
    const m = (nPoints === 'all' || typeof nPoints !== 'number')
      ? (run.end - run.start + 1)
      : Math.min(nPoints, run.end - run.start + 1);
    if (m < 2) return [];

    const comparisons = [];
    const formatRange = function (a, b) {
      return '[' + Number(a).toFixed(3) + ', ' + Number(b).toFixed(3) + ']';
    };
    const bracketStr = formatRange(bracket.x0, bracket.x1);

    // 1. Central (quanh tâm)
    const half = Math.floor((m - 2) / 2);
    let s1_start = Math.max(run.start, k - half);
    let s1_end = s1_start + m - 1;
    if (s1_end > run.end) {
      s1_end = run.end;
      s1_start = Math.max(run.start, s1_end - (m - 1));
    }
    const s1 = { start: s1_start, end: s1_end };
    try {
      const r1 = runIterationAttempt(xs, ys, bracket, s1, true, attemptOpts);
      if (r1 && r1.xbar != null && isFinite(r1.xbar)) {
        comparisons.push({
          rootName: 'x*',
          bracketStr: bracketStr,
          stencilStr: formatRange(xs[s1.start], xs[s1.end]),
          stencil: s1,
          methodName: 'Newton tiến / lùi (quanh tâm)',
          xbar: r1.xbar,
          result: r1
        });
      }
    } catch (e) {}

    // 2. Backward from x_{k+1}
    const s2_end = k + 1;
    const s2_start = s2_end - (m - 1);
    if (s2_start >= run.start && s2_end <= run.end) {
      const s2 = { start: s2_start, end: s2_end };
      try {
        const r2 = runIterationAttempt(xs, ys, bracket, s2, false, attemptOpts);
        if (r2 && r2.xbar != null && isFinite(r2.xbar)) {
          const subK1 = PPS.central ? PPS.central.sub(k + 1) : String(k + 1);
          comparisons.push({
            rootName: 'x*',
            bracketStr: bracketStr,
            stencilStr: formatRange(xs[s2.start], xs[s2.end]),
            stencil: s2,
            methodName: 'Newton lùi (từ x' + subK1 + ' = ' + Number(xs[k + 1]).toFixed(3) + ')',
            xbar: r2.xbar,
            result: r2
          });
        }
      } catch (e) {}
    }

    // 3. Forward from x_k
    const s3_start = k;
    const s3_end = s3_start + m - 1;
    if (s3_start >= run.start && s3_end <= run.end) {
      const s3 = { start: s3_start, end: s3_end };
      try {
        const r3 = runIterationAttempt(xs, ys, bracket, s3, true, attemptOpts);
        if (r3 && r3.xbar != null && isFinite(r3.xbar)) {
          const subK = PPS.central ? PPS.central.sub(k) : String(k);
          comparisons.push({
            rootName: 'x*',
            bracketStr: bracketStr,
            stencilStr: formatRange(xs[s3.start], xs[s3.end]),
            stencil: s3,
            methodName: 'Newton tiến (từ x' + subK + ' = ' + Number(xs[k]).toFixed(3) + ')',
            xbar: r3.xbar,
            result: r3
          });
        }
      } catch (e) {}
    }

    return comparisons;
  }

  function solveAllIteration(xs, ys, ybar, opts) {
    opts = opts || {};
    if (!xs.length || xs.length !== ys.length) throw new Error('needPoints');
    if (xs.length < 2) throw new Error('tooFew');
    if (!isFinite(ybar)) throw new Error('numeric');

    const spacing = PPS.newton.checkEqualSpacing(xs);
    const h = spacing.h;
    const yScale = yScaleOf(ys, ybar);

    const brackets = findBrackets(xs, ys, ybar);
    if (!brackets.length) throw new Error('noBracket');

    const nPoints = normalizePoints(opts.nPoints);
    const attemptOpts = buildAttemptOptions(opts, ybar, h, yScale);
    const directions = directionsFor(opts.method || 'auto');
    const results = [];

    for (let i = 0; i < brackets.length; i++) {
      const bracket = brackets[i];

      if (bracket.exact) {
        results.push({
          mode: 'iteration',
          exact: true,
          bracket: bracket,
          xbar: bracket.x,
          nodeIndex: bracket.index,
          residual: 0,
          checkY: ybar,
          history: [],
          status: 'CONVERGED_CERTIFIED',
          converged: true,
          iterations: 0,
          h: h,
          method: null,
          stencil: null,
          attempts: []
        });
        continue;
      }

      const run = monotoneRun(xs, ys, bracket.k);
      const stencil = selectStencil(run.start, run.end, bracket.k, nPoints);
      const attempts = runDirections(xs, ys, bracket, stencil, directions, attemptOpts);
      const best = chooseAttempt(attempts);
      const comparisons = buildStencilComparison(xs, ys, bracket, run, nPoints, attemptOpts);

      results.push(Object.assign({
        mode: 'iteration',
        exact: false,
        nPoints: nPoints,
        run: run,
        stencil: stencil,
        attempts: attempts,
        comparisons: comparisons
      }, best));
    }

    return results;
  }

  function solveAll(xs, ys, ybar, opts) {
    opts = opts || {};
    if (opts.mode === 'function') {
      return solveAllInverseFunction(xs, ys, ybar, opts);
    }
    return solveAllIteration(xs, ys, ybar, opts);
  }

  PPS.inverse = {
    solve: solve,
    solveAll: solveAll,
    solveInverseFunction: solveInverseFunction,
    solveIteration: solveIteration,
    solveAllInverseFunction: solveAllInverseFunction,
    solveAllIteration: solveAllIteration,
    findBracket: findBracket,
    findBrackets: findBrackets,
    scanMonotonicIntervals: scanMonotonicIntervals,
    runIterationAttempt: runIterationAttempt,
    buildStencilComparison: buildStencilComparison
  };
})();
