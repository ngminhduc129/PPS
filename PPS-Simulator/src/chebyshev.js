/* ============================================================
   chebyshev.js — Chebyshev optimal interpolation nodes,
                  theoretical error bounds, and Runge phenomenon demo
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const { factorial } = PPS.format;

  /**
   * Sinh mốc nội suy Chebyshev trên đoạn [a, b], bậc n (n + 1 mốc).
   * k = 0..n: theta_k = (2k + 1)*PI / (2*(n + 1)), t_k = cos(theta_k)
   * x_k = (a + b)/2 + ((b - a)/2) * t_k
   */
  function generateNodes(a, b, n) {
    a = Number(a);
    b = Number(b);
    n = Math.floor(Number(n));
    if (!isFinite(a) || !isFinite(b)) throw new Error('numeric');
    if (!(b > a)) throw new Error('invalidInterval');
    if (n < 1 || !isFinite(n)) throw new Error('invalidDegree');

    const xMid = (a + b) / 2;
    const halfLen = (b - a) / 2;
    const totalNodes = n + 1;
    const nodes = [];

    for (let k = 0; k < totalNodes; k++) {
      const theta = ((2 * k + 1) * Math.PI) / (2 * totalNodes);
      const t = Math.cos(theta);
      const x = xMid + halfLen * t;
      nodes.push({
        k: k,
        theta: theta,
        t: t,
        x: x
      });
    }

    // Sắp xếp tăng dần theo hoành độ x để thuận tiện biểu diễn và nội suy
    const sorted = nodes.slice().sort(function (p, q) { return p.x - q.x; });
    return {
      a: a,
      b: b,
      n: n,
      totalNodes: totalNodes,
      rawNodes: nodes,
      sortedNodes: sorted,
      xs: sorted.map(function (p) { return p.x; })
    };
  }

  /**
   * Sinh mốc cách đều trên đoạn [a, b] với n + 1 mốc để so sánh.
   */
  function equispacedNodes(a, b, n) {
    a = Number(a);
    b = Number(b);
    n = Math.floor(Number(n));
    const h = (b - a) / n;
    const xs = [];
    for (let i = 0; i <= n; i++) {
      xs.push(a + i * h);
    }
    return xs;
  }

  /**
   * Đánh giá sai số chặn trên lý thuyết Chebyshev:
   * max |f(x) - P_n(x)| <= [M_{n+1} / ((n + 1)! * 2^n)] * ((b - a) / 2)^(n + 1)
   */
  function errorBound(a, b, n, M) {
    M = Number(M);
    if (!isFinite(M) || M <= 0) return null;
    const halfLen = (b - a) / 2;
    const denom = factorial(n + 1) * Math.pow(2, n);
    const num = M * Math.pow(halfLen, n + 1);
    return num / denom;
  }

  /**
   * Trình thông dịch hàm f(x) an toàn cho các hàm mẫu và biểu thức toán học.
   */
  function parseFunction(expr) {
    expr = (expr || '').trim();
    if (!expr) expr = '1 / (1 + 25 * x * x)'; // Hàm Runge kinh điển

    // Các hàm toán học thông dụng
    return function (x) {
      return (function (x, sin, cos, tan, exp, log, abs, sqrt, pow, PI, E) {
        try {
          return Number(eval(expr));
        } catch (e) {
          return NaN;
        }
      })(
        x,
        Math.sin,
        Math.cos,
        Math.tan,
        Math.exp,
        Math.log,
        Math.abs,
        Math.sqrt,
        Math.pow,
        Math.PI,
        Math.E
      );
    };
  }

  /**
   * So sánh nội suy Lagrange mốc cách đều vs mốc Chebyshev trên [a, b].
   */
  function compareRunge(fnExpr, a, b, n, sampleCount) {
    sampleCount = sampleCount || 200;
    const fn = parseFunction(fnExpr);
    const cheb = generateNodes(a, b, n);
    const equiX = equispacedNodes(a, b, n);

    const chebY = cheb.xs.map(fn);
    const equiY = equiX.map(fn);

    const chebInterp = PPS.lagrange.interpolate(cheb.xs, chebY);
    const equiInterp = PPS.lagrange.interpolate(equiX, equiY);

    // Lấy mẫu mịn trên [a, b]
    const step = (b - a) / (sampleCount - 1);
    const grid = [];
    let maxErrCheb = 0;
    let maxErrEqui = 0;

    for (let i = 0; i < sampleCount; i++) {
      const x = a + i * step;
      const yTrue = fn(x);
      const yCheb = chebInterp.evaluate(x);
      const yEqui = equiInterp.evaluate(x);

      const errCheb = Math.abs(yTrue - yCheb);
      const errEqui = Math.abs(yTrue - yEqui);

      if (isFinite(errCheb)) maxErrCheb = Math.max(maxErrCheb, errCheb);
      if (isFinite(errEqui)) maxErrEqui = Math.max(maxErrEqui, errEqui);

      grid.push({
        x: x,
        yTrue: yTrue,
        yCheb: yCheb,
        yEqui: yEqui,
        errCheb: errCheb,
        errEqui: errEqui
      });
    }

    return {
      a: a,
      b: b,
      n: n,
      fnExpr: fnExpr,
      chebNodes: cheb.sortedNodes,
      chebX: cheb.xs,
      chebY: chebY,
      equiX: equiX,
      equiY: equiY,
      grid: grid,
      maxErrCheb: maxErrCheb,
      maxErrEqui: maxErrEqui,
      chebInterp: chebInterp,
      equiInterp: equiInterp
    };
  }

  PPS.chebyshev = {
    generateNodes: generateNodes,
    equispacedNodes: equispacedNodes,
    errorBound: errorBound,
    parseFunction: parseFunction,
    compareRunge: compareRunge
  };
})();
