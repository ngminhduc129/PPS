/* ============================================================
   linalg.js — Gaussian elimination with partial pivoting,
              Thomas algorithm for tridiagonal systems
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});

  function cloneMatrix(A) {
    return A.map(function (row) { return row.slice(); });
  }

  /**
   * Solve A x = b using Gaussian elimination with partial pivoting.
   * A: n x n array of arrays, b: length n array.
   * Returns { solution, U, pivots }.
   */
  function solve(A, b) {
    const n = A.length;
    if (!n || b.length !== n) throw new Error('range');
    const M = A.map(function (row, i) { return row.slice().concat([b[i]]); });

    for (let col = 0; col < n; col++) {
      let pivotRow = col;
      let pivotValue = Math.abs(M[col][col]);
      for (let r = col + 1; r < n; r++) {
        const value = Math.abs(M[r][col]);
        if (value > pivotValue) { pivotValue = value; pivotRow = r; }
      }
      if (pivotValue < 1e-13) throw new Error('singular');
      if (pivotRow !== col) {
        const tmp = M[col]; M[col] = M[pivotRow]; M[pivotRow] = tmp;
      }
      const pivot = M[col][col];
      for (let r = col + 1; r < n; r++) {
        const factor = M[r][col] / pivot;
        if (Math.abs(factor) < 1e-18) continue;
        M[r][col] = 0;
        for (let c = col + 1; c <= n; c++) {
          M[r][c] -= factor * M[col][c];
        }
      }
    }

    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      let sum = M[i][n];
      for (let j = i + 1; j < n; j++) sum -= M[i][j] * x[j];
      if (Math.abs(M[i][i]) < 1e-13) throw new Error('singular');
      x[i] = sum / M[i][i];
    }

    return { solution: x, augmented: M };
  }

  /**
   * Thomas algorithm for tridiagonal systems.
   * lower[i] = A[i+1][i] (length n-1), diag[i] = A[i][i],
   * upper[i] = A[i][i+1] (length n-1), d = right-hand side.
   */
  function thomas(lower, diag, upper, d) {
    const n = diag.length;
    if (n === 1) return [d[0] / diag[0]];
    const c = new Array(n - 1).fill(0);
    const dd = new Array(n).fill(0);
    const x = new Array(n).fill(0);

    dd[0] = diag[0];
    if (Math.abs(dd[0]) < 1e-13) throw new Error('singular');
    c[0] = upper[0] / dd[0];
    dd[0] = d[0] / dd[0];

    for (let i = 1; i < n; i++) {
      const denom = diag[i] - lower[i - 1] * c[i - 1];
      if (Math.abs(denom) < 1e-13) throw new Error('singular');
      if (i < n - 1) c[i] = upper[i] / denom;
      dd[i] = (d[i] - lower[i - 1] * dd[i - 1]) / denom;
    }

    x[n - 1] = dd[n - 1];
    for (let i = n - 2; i >= 0; i--) {
      x[i] = dd[i] - c[i] * x[i + 1];
    }
    return x;
  }

  function transpose(A) {
    if (!A.length) return [];
    return A[0].map(function (_, j) { return A.map(function (row) { return row[j]; }); });
  }

  function matmul(A, B) {
    if (!A.length || !B.length) return [];
    const n = A.length, m = B[0].length, k = B.length;
    const out = Array.from({ length: n }, function () { return new Array(m).fill(0); });
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < m; j++) {
        let sum = 0;
        for (let t = 0; t < k; t++) sum += A[i][t] * B[t][j];
        out[i][j] = sum;
      }
    }
    return out;
  }

  function matvec(A, v) {
    return A.map(function (row) {
      let sum = 0;
      for (let j = 0; j < row.length; j++) sum += row[j] * v[j];
      return sum;
    });
  }

  /**
   * Compute matrix inverse A^-1 using Gauss-Jordan elimination with partial pivoting.
   * A: n x n array of arrays.
   * Returns n x n inverse array of arrays.
   */
  function invert(A) {
    const n = A.length;
    if (!n || !A[0] || A[0].length !== n) throw new Error('range');

    const aug = new Array(n);
    for (let i = 0; i < n; i++) {
      aug[i] = new Array(2 * n).fill(0);
      for (let j = 0; j < n; j++) aug[i][j] = A[i][j];
      aug[i][n + i] = 1;
    }

    for (let col = 0; col < n; col++) {
      let pivotRow = col;
      let pivotVal = Math.abs(aug[col][col]);
      for (let r = col + 1; r < n; r++) {
        const val = Math.abs(aug[r][col]);
        if (val > pivotVal) {
          pivotVal = val;
          pivotRow = r;
        }
      }
      if (pivotVal < 1e-13) throw new Error('singular');

      if (pivotRow !== col) {
        const tmp = aug[col];
        aug[col] = aug[pivotRow];
        aug[pivotRow] = tmp;
      }

      const diag = aug[col][col];
      for (let c = 0; c < 2 * n; c++) {
        aug[col][c] /= diag;
      }

      for (let r = 0; r < n; r++) {
        if (r === col) continue;
        const factor = aug[r][col];
        if (Math.abs(factor) < 1e-18) continue;
        for (let c = 0; c < 2 * n; c++) {
          aug[r][c] -= factor * aug[col][c];
        }
      }
    }

    const inv = new Array(n);
    for (let i = 0; i < n; i++) {
      inv[i] = aug[i].slice(n);
    }
    return inv;
  }

  PPS.linalg = {
    solve: solve,
    thomas: thomas,
    transpose: transpose,
    matmul: matmul,
    matvec: matvec,
    invert: invert,
    cloneMatrix: cloneMatrix
  };
})();
