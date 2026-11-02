/* ============================================================
   run-tests.mjs — engine tests (Node, no framework)
   Usage: node tests/run-tests.mjs
   ============================================================ */

import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { fileURLToPath } from 'node:url';
import assert from 'node:assert/strict';

const here = path.dirname(fileURLToPath(import.meta.url));
const srcDir = path.resolve(here, '..', 'src');

globalThis.window = globalThis;

const files = [
  'format.js', 'i18n.js', 'linalg.js', 'polynomial.js', 'horner.js', 'chebyshev.js', 'lagrange.js',
  'newton.js', 'central.js', 'inverse.js', 'spline.js', 'leastSquares.js'
];

for (const file of files) {
  const code = fs.readFileSync(path.join(srcDir, file), 'utf8');
  vm.runInThisContext(code, { filename: file });
}

const PPS = globalThis.PPS;
let passed = 0;
let failed = 0;

function test(name, fn) {
  try {
    fn();
    passed++;
    console.log('  ✓ ' + name);
  } catch (error) {
    failed++;
    console.error('  ✗ ' + name);
    console.error('    ' + (error && error.stack ? error.stack.split('\n').slice(0, 3).join('\n    ') : error));
  }
}

function close(actual, expected, eps = 1e-8) {
  assert.ok(Math.abs(actual - expected) < eps, `expected ${expected}, got ${actual}`);
}

console.log('\nHorner');
test('P(2) and derivatives of 2x^3 - 3x^2 + 4x - 5', () => {
  const r = PPS.horner.evaluate([2, -3, 4, -5], 2, 3);
  close(r.value, 7);
  close(r.derivatives[1], 16);
  close(r.derivatives[2], 18);
  close(r.derivatives[3], 12);
  close(r.derivative, 12);
});

test('k beyond degree gives 0', () => {
  const r = PPS.horner.evaluate([4, -3, 2, -5], 2, 7);
  close(r.derivative, 0);
  assert.equal(r.kTooBig, true);
});

console.log('\nPolynomial');
test('multiply and evaluate', () => {
  const p = PPS.poly.mul([1, 1], [-1, 1]); // x^2 - 1
  close(PPS.poly.evaluate(p, 3), 8);
  close(PPS.poly.evaluate(p, -2), 3);
});

console.log('\nLagrange');
test('interpolates all nodes', () => {
  const xs = [-1, 0, 1, 2];
  const ys = [2, -1, 3, 1];
  const r = PPS.lagrange.interpolate(xs, ys);
  xs.forEach((x, i) => close(r.evaluate(x), ys[i], 1e-9));
  assert.equal(r.degree, 3);
});

test('computes product table omega(x) and synthetic division table', () => {
  const xs = [-1, 0, 1, 2];
  const ys = [2, -1, 3, 1];
  const r = PPS.lagrange.interpolate(xs, ys);
  assert.ok(Array.isArray(r.productTable));
  assert.equal(r.productTable.length, 5); // steps 0..4
  assert.ok(Array.isArray(r.divisionTable));
  assert.equal(r.divisionTable.length, 4); // nodes 0..3

  // omega(x) = (x+1)x(x-1)(x-2) = x^4 - 2x^3 - x^2 + 2x
  close(r.omega[0], 0, 1e-9);
  close(r.omega[1], 2, 1e-9);
  close(r.omega[2], -1, 1e-9);
  close(r.omega[3], -2, 1e-9);
  close(r.omega[4], 1, 1e-9);

  // divisionTable for x0 = -1: remainder is 0, D0 is -6
  close(r.divisionTable[0].remainder, 0, 1e-9);
  close(r.divisionTable[0].Di, -6, 1e-9);
  close(r.divisionTable[0].ci, 2 / -6, 1e-9);

  r.divisionTable.forEach(row => close(row.remainder, 0, 1e-9));
});

test('rejects duplicate nodes', () => {
  assert.throws(() => PPS.lagrange.interpolate([0, 0, 1], [1, 2, 3]), /dupX/);
});

console.log('\nNewton divided differences');
test('matches Lagrange at an interior point', () => {
  const xs = [-1, 0, 1, 2];
  const ys = [2, -1, 3, 1];
  const lag = PPS.lagrange.interpolate(xs, ys);
  const net = PPS.newton.interpolate(xs, ys);
  close(net.evaluate(0.37), lag.evaluate(0.37), 1e-9);
  xs.forEach((x, i) => close(net.evaluate(x), ys[i], 1e-9));
});

console.log('\nNewton equal spacing');
test('forward reproduces a cubic', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = xs.map(x => x * x * x);
  const r = PPS.newton.equalSpacing(xs, ys, 1.5, 'forward');
  close(r.value, 3.375, 1e-10);
  close(r.derivative, 6.75, 1e-10);
  assert.equal(r.mode, 'forward');
  assert.ok(Array.isArray(r.cForward) && r.cForward.length === 5);
  assert.ok(Array.isArray(r.cBackward) && r.cBackward.length === 5);
  assert.ok(Array.isArray(r.productTable) && r.productTable.length === 5);
  assert.ok(Array.isArray(r.polyX) && r.polyX.length > 0);
});

test('backward reproduces a cubic', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = xs.map(x => x * x * x);
  const r = PPS.newton.equalSpacing(xs, ys, 3.4, 'backward');
  close(r.value, Math.pow(3.4, 3), 1e-10);
  close(r.derivative, 3 * 3.4 * 3.4, 1e-9);
});

test('rejects non-uniform spacing', () => {
  assert.throws(() => PPS.newton.equalSpacing([0, 1, 3], [0, 1, 9], 1.5, 'auto'), /spacing/);
});

console.log('\nCentral formulas');
const cxs = [-2, -1, 0, 1, 2];
const cys = cxs.map(x => x * x * x);
for (const formula of ['gauss1', 'gauss2', 'stirling']) {
  test(formula + ' reproduces a cubic at x = 0.3', () => {
    const r = PPS.central.interpolate(cxs, cys, 0.3, { formula, N: 2 });
    close(r.value, 0.027, 1e-9);
    close(r.derivative, 0.27, 1e-9);
    assert.ok(Array.isArray(r.cVector));
    assert.ok(Array.isArray(r.dVector));
    assert.ok(Array.isArray(r.productTable));
    assert.ok(Array.isArray(r.polyX));
  });
}

test('Bessel reproduces a cubic', () => {
  const xs = [-1, 0, 1, 2];
  const ys = xs.map(x => x * x * x);
  const r = PPS.central.interpolate(xs, ys, 0.4, { formula: 'bessel', N: 2 });
  close(r.value, 0.064, 1e-9);
  close(r.derivative, 0.48, 1e-9);
});

test('auto formula chooses Stirling on a central point', () => {
  const r = PPS.central.interpolate(cxs, cys, 0.1, { formula: 'auto', N: 2 });
  assert.equal(r.formula, 'stirling');
});

test('central auto selects Bessel off-centre', () => {
  const r = PPS.central.interpolate(cxs, cys, 0.6, { formula: 'auto', N: 'auto' });
  assert.equal(r.formula, 'bessel');
});

console.log('\nInverse interpolation');
test('solves x^2 = 6.25', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = xs.map(x => x * x);
  const r = PPS.inverse.solve(xs, ys, 6.25, { method: 'auto', eps: 1e-10, maxIter: 60 });
  close(r.xbar, 2.5, 1e-7);
  assert.equal(r.converged, true);
});

test('returns exact node when ybar matches', () => {
  const r = PPS.inverse.solve([0, 1, 2], [1, 2, 5], 2, {});
  assert.equal(r.exact, true);
  close(r.xbar, 1);
});

test('throws when ybar is outside the data range', () => {
  assert.throws(() => PPS.inverse.solve([0, 1, 2], [0, 1, 2], 5, {}), /noBracket/);
});

test('findBrackets returns all bracketing intervals', () => {
  const brackets = PPS.inverse.findBrackets([0, 1, 2, 3], [0, 2, 0, 2], 1);
  assert.equal(brackets.length, 3);
  assert.equal(brackets[0].k, 0);
  assert.equal(brackets[1].k, 1);
  assert.equal(brackets[2].k, 2);
});

test('solveAll iteration returns all roots for oscillating data', () => {
  const res = PPS.inverse.solveAll([0, 1, 2, 3], [0, 2, 0, 2], 1, { mode: 'iteration', maxIter: 30 });
  assert.equal(res.length, 3);
  res.forEach(function (r) {
    assert.equal(r.mode, 'iteration');
    assert.equal(r.converged, true);
    assert.equal(typeof r.xbar, 'number');
    close(r.xbar, r.bracket.k + 0.5, 1e-8);
  });
});

test('solveAll function returns all roots for oscillating data', () => {
  const res = PPS.inverse.solveAll([0, 1, 2, 3], [0, 2, 1, 3], 1.5, { mode: 'function' });
  assert.equal(res.length, 3);
  res.forEach(function (r) {
    assert.equal(r.mode, 'function');
  });
});

test('iteration matches Python 11_Loop_Inverse sample with all points', () => {
  const xs = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1];
  const ys = [-1.121982, -1.214072, -1.256192, -1.254842, -1.214912, -1.139932, -1.032282, -0.523962, -0.292592];
  const res = PPS.inverse.solveAll(xs, ys, -1.25, {
    mode: 'iteration',
    nPoints: 'all',
    tolT: 1e-8,
    tolRes: 1e-8,
    tolX: 1e-8,
    maxIter: 100
  });
  assert.equal(res.length, 2);
  close(res[0].xbar, -0.5653971, 1e-5);
  close(res[1].xbar, -0.0995931, 1e-5);
  assert.equal(res[0].converged, true);
  assert.equal(res[1].converged, true);
  assert.equal(res[1].method, 'backward');
  const forward = res[1].attempts.filter(function (a) { return a.method === 'forward'; })[0];
  assert.equal(forward.converged, false);
  assert.equal(forward.status, 'DIVERGED');
});

test('iteration clips the stencil to the monotone branch', () => {
  const xs = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1];
  const ys = [-1.121982, -1.214072, -1.256192, -1.254842, -1.214912, -1.139932, -1.032282, -0.523962, -0.292592];
  const res = PPS.inverse.solveAll(xs, ys, -1.25, { mode: 'iteration', nPoints: 'all' });
  assert.equal(res[0].stencil.start, 0);
  assert.equal(res[0].stencil.end, 2);
  assert.equal(res[1].stencil.start, 2);
  assert.equal(res[1].stencil.end, 8);
});

test('iteration centers a 4-point stencil around the bracket', () => {
  const xs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
  const ys = xs.map(function (x) { return x * x; });
  const res = PPS.inverse.solveIteration(xs, ys, 20.25, { nPoints: 4, method: 'forward' });
  assert.equal(res.stencil.start, 3);
  assert.equal(res.stencil.end, 6);
  close(res.xbar, 4.5, 1e-6);
});

test('iteration reproduces Algorithm.md 7.3.2 trace', () => {
  const res = PPS.inverse.solveIteration(
    [2.380, 2.495, 2.610, 2.725],
    [-1.5819, -1.5127, -1.3977, -1.3814],
    -1.43,
    { nPoints: 4, method: 'forward', eps: 1e-6, tolRes: 1e-6, tolX: 1e-6, maxIter: 50, qCoeff: 0.75 }
  );
  close(res.xbar, 2.5738863, 1e-6);
  assert.equal(res.iterations, 31);
  assert.equal(res.status, 'CONVERGED_CERTIFIED');
});

const inverseSampleX = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1];
const inverseSampleY = [-1.121982, -1.214072, -1.256192, -1.254842, -1.214912, -1.139932, -1.032282, -0.523962, -0.292592];

test('function method: n_points changes the stencil and the result', () => {
  const n2 = PPS.inverse.solveInverseFunction(inverseSampleX, inverseSampleY, -1.25, { nPoints: 2 });
  const n3 = PPS.inverse.solveInverseFunction(inverseSampleX, inverseSampleY, -1.25, { nPoints: 3 });

  assert.equal(n2.stencil.start, 1);
  assert.equal(n2.stencil.end, 2);
  assert.equal(n3.stencil.start, 0);
  assert.equal(n3.stencil.end, 2);

  close(n2.xbar, -0.5367521, 1e-6);
  close(n3.xbar, -0.5420907, 1e-6);
  assert.ok(Math.abs(n2.xbar - n3.xbar) > 1e-4, 'n_points must change the root');
});

test('function method clips n_points to the monotone branch', () => {
  const res = PPS.inverse.solveInverseFunction(inverseSampleX, inverseSampleY, -1.25, { nPoints: 20 });
  assert.equal(res.stencil.start, 0);
  assert.equal(res.stencil.end, 2);
  assert.equal(res.stencil.count, 3);
  assert.equal(res.nPoints, 20);
});

test('function method falls back to linear inverse interpolation outside the bracket', () => {
  const res = PPS.inverse.solveAllInverseFunction(inverseSampleX, inverseSampleY, -1.25, { nPoints: 4 });
  assert.equal(res.length, 2);

  const second = res[1];
  assert.equal(second.fallback, true);
  assert.equal(second.stencil.start, 2);
  assert.equal(second.stencil.end, 5);
  close(second.xbar, -0.2196844, 1e-6);
  assert.ok(second.polyX < second.bracket.x0 || second.polyX > second.bracket.x1);
  assert.equal(second.inBracket, true);
});

test('function method keeps every root inside its bracket', () => {
  const res = PPS.inverse.solveAllInverseFunction(inverseSampleX, inverseSampleY, -1.25, { nPoints: 'all' });
  assert.equal(res.length, 2);
  res.forEach(function (r) {
    assert.ok(r.xbar >= r.bracket.x0 && r.xbar <= r.bracket.x1);
    assert.equal(r.inBracket, true);
  });
  assert.equal(res[0].stencil.start, 0);
  assert.equal(res[0].stencil.end, 2);
  assert.equal(res[1].stencil.start, 2);
  assert.equal(res[1].stencil.end, 8);
});

console.log('\nSplines');
const sxs = [0, 1, 2, 3, 4];
const sys = [1, 2, 0, 2, 1];

test('linear spline passes through nodes and is piecewise linear', () => {
  const r = PPS.spline.linear(sxs, sys);
  PPS.spline.verify(sxs, sys, r).forEach(v => close(v.error, 0, 1e-12));
  close(r.segments[0].a, 1);
  close(r.segments[0].b, 1);
});

test('quadratic spline passes through nodes and has continuous slope', () => {
  const r = PPS.spline.quadratic(sxs, sys, 0);
  PPS.spline.verify(sxs, sys, r).forEach(v => close(v.error, 0, 1e-9));
  for (let k = 0; k < r.segments.length - 1; k++) {
    const left = 2 * r.segments[k].a * r.segments[k].x1 + r.segments[k].b;
    const right = 2 * r.segments[k + 1].a * r.segments[k + 1].x0 + r.segments[k + 1].b;
    close(left, right, 1e-9);
  }
});

test('natural cubic spline has zero second derivative at ends', () => {
  const r = PPS.spline.cubic(sxs, sys);
  close(r.alpha[0], 0, 1e-12);
  close(r.alpha[r.alpha.length - 1], 0, 1e-12);
  PPS.spline.verify(sxs, sys, r).forEach(v => close(v.error, 0, 1e-8));
});

test('cubic spline handles two nodes', () => {
  const r = PPS.spline.cubic([0, 1], [1, 3]);
  close(r.segments[0].evaluate(0.5), 2, 1e-12);
});

console.log('\nLeast squares');
test('recovers exact linear data', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = xs.map(x => 2 + 3 * x);
  const r = PPS.leastSquares.linear(xs, ys, ['one', 'x']);
  close(r.coef[0], 2, 1e-9);
  close(r.coef[1], 3, 1e-9);
  close(r.sse, 0, 1e-15);
  assert.ok(r.r2 > 0.999999);
});

test('recovers exact quadratic data', () => {
  const xs = [-2, -1, 0, 1, 2, 3];
  const ys = xs.map(x => 1 - 2 * x + 0.5 * x * x);
  const r = PPS.leastSquares.linear(xs, ys, ['one', 'x', 'x2']);
  close(r.coef[0], 1, 1e-9);
  close(r.coef[1], -2, 1e-9);
  close(r.coef[2], 0.5, 1e-9);
  close(r.sse, 0, 1e-12);
});

test('log-linear recovers y = 2 e^(0.5 x) + 1', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = xs.map(x => 2 * Math.exp(0.5 * x) + 1);
  const r = PPS.leastSquares.logLinear(xs, ys, ['x'], 1);
  close(r.a, 2, 1e-9);
  close(r.b[0], 0.5, 1e-9);
  close(r.sse, 0, 1e-10);
});

test('log-linear rejects mixed signs', () => {
  assert.throws(() => PPS.leastSquares.logLinear([0, 1, 2], [1, -1, 2], ['x'], 0), /sign|sameSign/);
});

console.log('\nLinear algebra');
test('Gaussian elimination with pivoting', () => {
  const A = [[0, 2, 1], [1, -2, -3], [-1, 1, 2]];
  const b = [5, -3, 1];
  const x = PPS.linalg.solve(A, b).solution;
  close(x[0], 0, 1e-12);
  close(x[1], 3, 1e-12);
  close(x[2], -1, 1e-12);
});

test('Thomas algorithm', () => {
  const x = PPS.linalg.thomas([1, 1], [4, 4, 4], [1, 1], [5, 5, 5]);
  close(x[0], 15 / 14, 1e-12);
  close(x[1], 5 / 7, 1e-12);
});

console.log('\nChebyshev (Section 2)');
test('generates optimal Chebyshev nodes on [-1, 1]', () => {
  const res = PPS.chebyshev.generateNodes(-1, 1, 3); // 4 nodes
  assert.equal(res.totalNodes, 4);
  assert.equal(res.sortedNodes.length, 4);
  // cos(7pi/8), cos(5pi/8), cos(3pi/8), cos(pi/8)
  close(res.xs[0], Math.cos(7 * Math.PI / 8), 1e-12);
  close(res.xs[3], Math.cos(Math.PI / 8), 1e-12);
});

test('theoretical error bound on Chebyshev nodes', () => {
  const bound = PPS.chebyshev.errorBound(-1, 1, 3, 24);
  // 24 * (1)^4 / (4! * 2^3) = 24 / (24 * 8) = 1/8 = 0.125
  close(bound, 0.125, 1e-12);
});

test('Runge comparison runs without errors', () => {
  const comp = PPS.chebyshev.compareRunge('1 / (1 + 25 * x * x)', -1, 1, 6, 50);
  assert.ok(comp.maxErrCheb >= 0);
  assert.ok(comp.maxErrEqui >= 0);
});

console.log('\nHorner enhancements (Section 1)');
test('reverse Horner multiplication', () => {
  // P(x) = (x - 2)*(3x + 4) + 5 = 3x^2 - 2x - 3
  // qAsc = [4, 3], c = 2, r = 5 => aAsc = [-3, -2, 3]
  const aAsc = PPS.horner.reverseMultiplyAsc([4, 3], 2, 5);
  close(aAsc[0], -3);
  close(aAsc[1], -2);
  close(aAsc[2], 3);
});

console.log('\nNewton enhancements (Section 4)');
test('divided differences returns forward & backward vectors and derivatives', () => {
  const xs = [1.57, 1.59, 1.62, 1.64];
  const ys = [2.4625, 2.5857, 2.772, 2.897];
  const r = PPS.newton.interpolate(xs, ys);
  assert.equal(r.cForward.length, 4);
  assert.equal(r.cBackward.length, 4);
  assert.equal(r.productTable.length, 4);
  const derivs = r.evalDerivatives(1.64, 2);
  close(derivs[0], 2.897, 1e-8);
});

console.log('\nInverse function interpolation (Section 7)');
test('solves inverse function method x = P(y)', () => {
  const xs = [1.40, 1.50, 1.60, 1.70];
  const ys = [2.44018, 2.55290, 2.65511, 2.74677];
  const res = PPS.inverse.solveInverseFunction(xs, ys, 2.60);
  assert.equal(res.exact, false);
  close(res.xbar, 1.54477, 1e-4); // Matches 1.5447682 in Algorithm.md
  assert.ok(Array.isArray(res.cForward) && res.cForward.length === 4);
  assert.ok(Array.isArray(res.cBackward) && res.cBackward.length === 4);
  assert.ok(Array.isArray(res.productTable) && res.productTable.length === 4);
  assert.ok(Array.isArray(res.polyCanonical) && res.polyCanonical.length === 4);
  close(res.cForward[0], 1.40, 1e-8);
  close(res.cForward[1], 0.887154, 1e-4);
});

test('iteration method computes product table and difference coefficient vectors', () => {
  const xs = [2.380, 2.495, 2.610, 2.725];
  const ys = [-1.5819, -1.5127, -1.3977, -1.3814];
  const res = PPS.inverse.solveIteration(xs, ys, -1.43, { nPoints: 4, method: 'forward' });
  assert.ok(Array.isArray(res.cDiff) && res.cDiff.length === 4);
  assert.ok(Array.isArray(res.cNormalized) && res.cNormalized.length === 4);
  assert.ok(Array.isArray(res.cDivided) && res.cDivided.length === 4);
  assert.ok(Array.isArray(res.productTable) && res.productTable.length === 4);
  assert.ok(Array.isArray(res.polyT) && res.polyT.length > 0);
  close(res.cDiff[0], -1.5819, 1e-8);
  close(res.cDiff[1], 0.0692, 1e-8);
  close(res.cDiff[2], 0.0458, 1e-8);
  close(res.cDiff[3], -0.1445, 1e-8);
  close(res.cNormalized[1], 0.0692, 1e-8);
  close(res.cNormalized[2], 0.0458 / 2, 1e-8);
  close(res.cNormalized[3], -0.1445 / 6, 1e-8);
});

console.log('\nSpline boundary conditions (Section 8)');
test('clamped cubic spline satisfies derivative boundary conditions', () => {
  const sxs = [0, 1, 2];
  const sys = [0, 1, 8]; // x^3
  const r = PPS.spline.cubic(sxs, sys, { bcType: 'clamped', s0: 0, sn: 12 });
  close(r.segments[0].evalDerivative(0), 0, 1e-9);
  close(r.segments[1].evalDerivative(2), 12, 1e-9);
});

console.log('\nLeast squares power law (Section 9)');
test('recovers power law y = 2 * x^2', () => {
  const pxs = [1, 2, 3, 4, 5];
  const pys = pxs.map(x => 2 * x * x);
  const r = PPS.leastSquares.fitPower(pxs, pys);
  close(r.a, 2, 1e-8);
  close(r.b, 2, 1e-8);
  close(r.sse, 0, 1e-8);
});

console.log('\nAudit verification & UI/UX enhancements');
test('Stirling factorString uses integer roots (t² - 1), not fraction (t² - 0.25)', () => {
  const xs = [-2, -1, 0, 1, 2];
  const ys = [16, 1, 0, 1, 16]; // x^4
  const r = PPS.central.interpolate(xs, ys, 0.1, { formula: 'stirling', N: 2 });
  const term3 = r.terms.find(t => t.j === 3);
  const term4 = r.terms.find(t => t.j === 4);
  assert.ok(term3, 'term3 should exist');
  assert.ok(term3.factorString.includes('t² − 1'), 'term3 factorString should have (t² − 1)');
  assert.ok(!term3.factorString.includes('0.25'), 'term3 factorString should not have 0.25');
  assert.ok(term4, 'term4 should exist');
  assert.ok(term4.factorString.includes('t² − 1'), 'term4 factorString should have (t² − 1)');
  assert.ok(!term4.factorString.includes('0.25'), 'term4 factorString should not have 0.25');
});

test('all algorithm error codes are defined in both vi and en dictionaries', () => {
  const errCodes = [
    'basis', 'central', 'domain', 'dupX', 'dupY', 'empty', 'invalidDegree',
    'invalidInterval', 'leastSquares', 'needPoints', 'noBracket', 'numeric',
    'order', 'parse', 'positiveOnly', 'range', 'sameSign', 'sign', 'singular',
    'spacing', 'tooFew', 'xbar'
  ];
  for (const code of errCodes) {
    assert.ok(PPS.i18n.dict.vi['errors.' + code], 'Missing vi translation for errors.' + code);
    assert.ok(PPS.i18n.dict.en['errors.' + code], 'Missing en translation for errors.' + code);
  }
});

test('least squares fitPower throws positiveOnly on negative or zero data', () => {
  assert.throws(() => PPS.leastSquares.fitPower([0, 1, 2], [1, 2, 4]), /positiveOnly/);
  assert.throws(() => PPS.leastSquares.fitPower([1, 2, 3], [-1, 2, 4]), /positiveOnly/);
});

test('Chebyshev generateNodes validates interval and degree', () => {
  assert.throws(() => PPS.chebyshev.generateNodes(2, 1, 4), /invalidInterval/);
  assert.throws(() => PPS.chebyshev.generateNodes(0, 1, 0), /invalidDegree/);
});

test('Spline validates strictly increasing nodes', () => {
  assert.throws(() => PPS.spline.linear([1, 0.5, 2], [1, 2, 3]), /order/);
});

test('default rounding is 7 decimal places and can be updated dynamically', () => {
  PPS.format.setDecimals(7);
  assert.equal(PPS.format.getDecimals(), 7);
  assert.equal(PPS.format.fmt(1 / 3), '0.3333333');
  assert.equal(PPS.format.fmt(1.23456789), '1.2345679');

  PPS.format.setDecimals(4);
  assert.equal(PPS.format.getDecimals(), 4);
  assert.equal(PPS.format.fmt(1 / 3), '0.3333');
  assert.equal(PPS.format.fmt(1.23456789), '1.2346');

  // restore default 7
  PPS.format.setDecimals(7);
  assert.equal(PPS.format.getDecimals(), 7);
});

test('PPS.format.subscript and superscript exist and format properly', () => {
  assert.equal(typeof PPS.format.subscript, 'function');
  assert.equal(typeof PPS.format.superscript, 'function');
  assert.equal(PPS.format.subscript(0), '₀');
  assert.equal(PPS.format.subscript(12), '₁₂');
  assert.equal(PPS.format.subscript(-3), '₋₃');
  assert.equal(PPS.format.superscript(3), '³');
});

console.log('\nCross-validation: Mathematical equivalence across interpolation methods');
test('Lagrange, Newton Divided, Newton Forward, and Newton Backward produce identical results', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = [2, 4.5, 11, 23.5, 44]; // Cubic polynomial 0.5 x^3 + x^2 + x + 2
  const xStar = 2.4;

  const lagVal = PPS.lagrange.interpolate(xs, ys).evaluate(xStar);
  const netDivVal = PPS.newton.interpolate(xs, ys).evaluate(xStar);
  const netFwdVal = PPS.newton.equalSpacing(xs, ys, xStar, 'forward').value;
  const netBwdVal = PPS.newton.equalSpacing(xs, ys, xStar, 'backward').value;

  close(lagVal, netDivVal, 1e-9);
  close(lagVal, netFwdVal, 1e-9);
  close(lagVal, netBwdVal, 1e-9);
  close(lagVal, 15.168, 1e-3);
});

test('Central formulas (Gauss I, Gauss II, Stirling) match Lagrange on the same central window', () => {
  const cxs = [-2, -1, 0, 1, 2];
  const cys = [-7, 0, 1, 2, 9]; // x^3 + 1
  const xStar = 0.35;

  const lagVal = PPS.lagrange.interpolate(cxs, cys).evaluate(xStar);
  const g1 = PPS.central.interpolate(cxs, cys, xStar, { formula: 'gauss1', N: 2 });
  const g2 = PPS.central.interpolate(cxs, cys, xStar, { formula: 'gauss2', N: 2 });
  const stirling = PPS.central.interpolate(cxs, cys, xStar, { formula: 'stirling', N: 2 });

  close(g1.value, lagVal, 1e-9);
  close(g2.value, lagVal, 1e-9);
  close(stirling.value, lagVal, 1e-9);
});

console.log('\nDerivative validation: Analytical vs Finite Difference');
test('analytical derivatives match central finite difference approximation', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = [1, 3, 9, 27, 81];
  const xStar = 2.3;
  const delta = 1e-6;

  // Newton equal spacing forward
  const fwd = PPS.newton.equalSpacing(xs, ys, xStar, 'forward');
  const fwdPlus = PPS.newton.equalSpacing(xs, ys, xStar + delta, 'forward').value;
  const fwdMinus = PPS.newton.equalSpacing(xs, ys, xStar - delta, 'forward').value;
  const numDeriv = (fwdPlus - fwdMinus) / (2 * delta);
  close(fwd.derivative, numDeriv, 1e-5);

  // Central Stirling derivative
  const cxs = [0, 1, 2, 3, 4];
  const cys = [1, 4, 11, 26, 53];
  const stir = PPS.central.interpolate(cxs, cys, 2.2, { formula: 'stirling', N: 2 });
  const stirPlus = PPS.central.interpolate(cxs, cys, 2.2 + delta, { formula: 'stirling', N: 2 }).value;
  const stirMinus = PPS.central.interpolate(cxs, cys, 2.2 - delta, { formula: 'stirling', N: 2 }).value;
  const numStirDeriv = (stirPlus - stirMinus) / (2 * delta);
  close(stir.derivative, numStirDeriv, 1e-5);
});

console.log('\nSpline continuity: C0, C1, and C2 at interior nodes');
test('cubic spline satisfies C0, C1, C2 continuity across segments', () => {
  const xs = [0, 1.5, 3.2, 5.0, 7.5];
  const ys = [1.2, 3.8, 2.1, 6.4, 4.0];
  const cs = PPS.spline.cubic(xs, ys, { bcType: 'natural' });

  for (let i = 0; i < cs.segments.length - 1; i++) {
    const leftSeg = cs.segments[i];
    const rightSeg = cs.segments[i + 1];
    const xNode = xs[i + 1];

    // C0 continuity: S_i(x_{i+1}) == S_{i+1}(x_{i+1}) == y_{i+1}
    const yLeft = leftSeg.evaluate(xNode);
    const yRight = rightSeg.evaluate(xNode);
    close(yLeft, ys[i + 1], 1e-9);
    close(yRight, ys[i + 1], 1e-9);

    // C1 continuity: S'_i(x_{i+1}) == S'_{i+1}(x_{i+1})
    const dyLeft = leftSeg.evalDerivative(xNode);
    const dyRight = rightSeg.evalDerivative(xNode);
    close(dyLeft, dyRight, 1e-8);

    // C2 continuity: S''_i(x_{i+1}) == S''_{i+1}(x_{i+1}) == M_{i+1}
    const d2yLeft = 2 * leftSeg.local.C + 6 * leftSeg.local.D * (xNode - leftSeg.local.xi);
    const d2yRight = 2 * rightSeg.local.C;
    close(d2yLeft, d2yRight, 1e-8);
    close(d2yLeft, cs.M[i + 1], 1e-8);
  }

  // Natural boundary conditions: M_0 = 0 and M_n = 0
  close(cs.M[0], 0, 1e-9);
  close(cs.M[xs.length - 1], 0, 1e-9);
});

console.log('\nLeast squares: Multi-basis & Normal Equations');
test('linear least squares satisfies normal equations and recovers trigonometric data', () => {
  // y = 2 + 3 sin(x) - 1.5 cos(x)
  const xs = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
  const ys = xs.map(x => 2 + 3 * Math.sin(x) - 1.5 * Math.cos(x));
  const res = PPS.leastSquares.linear(xs, ys, ['one', 'sin', 'cos']);

  close(res.coef[0], 2, 1e-8);
  close(res.coef[1], 3, 1e-8);
  close(res.coef[2], -1.5, 1e-8);
  close(res.sse, 0, 1e-8);
  close(res.r2, 1, 1e-8);

  // Normal equations residual ||M * c - rhs|| ~ 0
  const Mc = PPS.linalg.matvec(res.M, res.coef);
  for (let i = 0; i < Mc.length; i++) {
    close(Mc[i], res.rhs[i], 1e-8);
  }
});

test('least squares computes design matrix Phi, normal matrix M, inverse invM, and pseudoinverse pinv', () => {
  const xs = [1, 2, 3, 4, 5];
  const ys = [2, 5, 10, 17, 26]; // 1 + x^2
  const res = PPS.leastSquares.linear(xs, ys, ['one', 'x', 'x2']);

  // Dimensions of Phi (5 x 3)
  assert.equal(res.phi.length, 5);
  assert.equal(res.phi[0].length, 3);

  // M = Phi^T * Phi (3 x 3)
  assert.equal(res.M.length, 3);
  assert.equal(res.M[0].length, 3);

  // invM = (Phi^T * Phi)^-1
  assert.ok(res.invM, 'invM must be computed');
  assert.equal(res.invM.length, 3);
  assert.equal(res.invM[0].length, 3);

  // M * invM should be Identity
  const I = PPS.linalg.matmul(res.M, res.invM);
  for (let r = 0; r < 3; r++) {
    for (let c = 0; c < 3; c++) {
      close(I[r][c], r === c ? 1 : 0, 1e-10);
    }
  }

  // Solution via invM * rhs matches coef
  const cFromInv = PPS.linalg.matvec(res.invM, res.rhs);
  for (let i = 0; i < 3; i++) {
    close(cFromInv[i], res.coef[i], 1e-10);
  }

  // Pseudoinverse pinv * ys matches coef
  assert.ok(res.pinv, 'pinv must be computed');
  const cFromPinv = PPS.linalg.matvec(res.pinv, ys);
  for (let i = 0; i < 3; i++) {
    close(cFromPinv[i], res.coef[i], 1e-10);
  }

  // Square system (N = m = 2) computes invPhi directly
  const squareRes = PPS.leastSquares.linear([1, 3], [3, 7], ['one', 'x']);
  assert.ok(squareRes.invPhi, 'invPhi must be computed for square system');
  const I_phi = PPS.linalg.matmul(squareRes.phi, squareRes.invPhi);
  for (let r = 0; r < 2; r++) {
    for (let c = 0; c < 2; c++) {
      close(I_phi[r][c], r === c ? 1 : 0, 1e-10);
    }
  }

  // Summation table matches normal equation sums
  const basis = PPS.leastSquares.getBasis(['one', 'x', 'x2']);
  const st = PPS.leastSquares.buildSumTable(xs, ys, basis);
  close(st.sums.x, 15, 1e-10);
  close(st.sums.y, 60, 1e-10);
  close(st.sums.prod_one, res.rhs[0], 1e-10);
  close(st.sums.prod_x, res.rhs[1], 1e-10);
  close(st.sums.prod_x2, res.rhs[2], 1e-10);
  assert.ok(res.r2adj <= res.r2, 'Adjusted R^2 must be <= R^2');
});

test('least squares logLinear supports negative y with C0=0 and basis [x, 1/x]', () => {
  // Synthetic data y = -2 * exp(0.5*x - 1/x)
  const testXs = [1, 2, 3, 4, 5];
  const testYs = testXs.map(x => -2 * Math.exp(0.5 * x - 1 / x));
  const res = PPS.leastSquares.logLinear(testXs, testYs, ['x', 'inv'], 0);
  close(res.a, -2, 1e-7);
  close(res.b[0], 0.5, 1e-7);
  close(res.b[1], -1.0, 1e-7);
  close(res.sse, 0, 1e-10);
  assert.equal(res.sign, -1);
  assert.equal(res.rhs.length, 3);
});

console.log('\nLinear algebra: Thomas TDMA matches Gaussian elimination');
test('Thomas algorithm matches Gaussian elimination on tridiagonal systems', () => {
  const lower = [1, 2, 1];
  const diag = [4, 5, 6, 7];
  const upper = [2, 1, 3];
  const rhs = [5, 12, 18, 25];

  const thomasSol = PPS.linalg.thomas(lower, diag, upper, rhs);

  // Full matrix form for Gaussian elimination
  const A = [
    [4, 2, 0, 0],
    [1, 5, 1, 0],
    [0, 2, 6, 3],
    [0, 0, 1, 7]
  ];
  const gaussSol = PPS.linalg.solve(A, rhs).solution;

  for (let i = 0; i < 4; i++) {
    close(thomasSol[i], gaussSol[i], 1e-9);
  }
});

console.log('\n' + '─'.repeat(52));
console.log(`  ${passed} passed, ${failed} failed`);
if (failed) process.exitCode = 1;
