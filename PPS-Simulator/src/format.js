/* ============================================================
   format.js — number formatting, factorials, polynomial strings
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});

  const SUPER = {
    '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
    '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'
  };

  const SUB = {
    '-': '₋', '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
    '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'
  };

  function superscript(n) {
    return String(n).split('').map(function (ch) { return SUPER[ch] || ch; }).join('');
  }

  function subscript(n) {
    return String(n).split('').map(function (ch) { return SUB[ch] || ch; }).join('');
  }

  function factorial(n) {
    n = Math.round(n);
    if (n < 0) return NaN;
    let f = 1;
    for (let i = 2; i <= n; i++) f *= i;
    return f;
  }

  function isNum(x) {
    return typeof x === 'number' && isFinite(x);
  }

  let currentDecimals = 7;
  let hasUserCustomDecimals = false;
  try {
    if (typeof localStorage !== 'undefined') {
      const saved = localStorage.getItem('pps_decimals');
      if (saved != null && saved !== '') {
        const parsed = parseInt(saved, 10);
        if (!isNaN(parsed) && parsed >= 0 && parsed <= 15) {
          currentDecimals = parsed;
          hasUserCustomDecimals = true;
        }
      }
    }
  } catch (e) {}

  function getDecimals() {
    return currentDecimals;
  }

  function setDecimals(n) {
    const d = Math.max(0, Math.min(15, Math.round(Number(n) || 0)));
    currentDecimals = d;
    hasUserCustomDecimals = true;
    try {
      if (typeof localStorage !== 'undefined') {
        localStorage.setItem('pps_decimals', String(d));
      }
    } catch (e) {}
    return currentDecimals;
  }

  function fmt(x, digits, opts) {
    opts = opts || {};
    // Mặc định làm tròn 7 chữ số sau dấu phẩy (hoặc theo currentDecimals)
    const d = (digits == null || (hasUserCustomDecimals && !opts.exactPrecision)) ? currentDecimals : digits;
    if (typeof x !== 'number') return String(x);
    if (!isFinite(x)) return x > 0 ? '+∞' : (x < 0 ? '−∞' : 'NaN');
    if (Math.abs(x) < 1e-14) return '0';

    const a = Math.abs(x);
    if (a >= 1e12 || (a < Math.pow(10, -(d + 1)) && a >= 1e-14)) {
      return x.toExponential(4).replace('e+', 'e').replace(/\.?0+e/, 'e');
    }

    let s = x.toFixed(d);
    if (!opts.fixed && s.indexOf('.') !== -1) {
      s = s.replace(/0+$/, '').replace(/\.$/, '');
    }
    if (s === '-0' || s === '−0') return '0';
    return s;
  }

  function coefFmt(c, isConstant) {
    if (Math.abs(c - 1) < 1e-12 && !isConstant) return '';
    if (Math.abs(c + 1) < 1e-12 && !isConstant) return '−';
    return fmt(c, currentDecimals);
  }

  function polyToString(coefs, varName, opts) {
    varName = varName || 'x';
    opts = opts || {};
    const digits = opts.digits == null ? currentDecimals : opts.digits;
    const terms = [];
    for (let i = coefs.length - 1; i >= 0; i--) {
      let c = coefs[i];
      if (!isNum(c)) continue;
      if (Math.abs(c) < 1e-11) continue;
      const isConst = i === 0;
      const sign = c < 0 ? '−' : '+';
      const abs = Math.abs(c);
      let body;
      if (isConst) {
        body = fmt(abs, digits);
      } else {
        const coefStr = Math.abs(abs - 1) < 1e-12 ? '' : fmt(abs, digits);
        const power = i === 1 ? varName : varName + superscript(i);
        body = coefStr + power;
      }
      terms.push({ sign: sign, body: body });
    }
    if (!terms.length) return '0';
    let out = '';
    terms.forEach(function (term, index) {
      if (index === 0) {
        out += (term.sign === '−' ? '−' : '') + term.body;
      } else {
        out += ' ' + (term.sign === '−' ? '−' : '+') + ' ' + term.body;
      }
    });
    return out;
  }

  function trimCoefs(coefs, eps) {
    eps = eps == null ? 1e-12 : eps;
    let n = coefs.length;
    while (n > 1 && Math.abs(coefs[n - 1]) < eps) n--;
    return coefs.slice(0, n);
  }

  function parseList(str) {
    if (str == null) return [];
    const cleaned = String(str).trim();
    if (!cleaned) return [];
    const parts = cleaned.split(/[\s,;]+/).filter(function (p) { return p !== ''; });
    const out = [];
    for (let i = 0; i < parts.length; i++) {
      const part = parts[i];
      const normalized = part.replace(/−/g, '-');
      const value = Number(normalized);
      if (!isFinite(value)) return { error: 'parse' };
      out.push(value);
    }
    return { values: out };
  }

  function parseNumber(str) {
    if (str == null) return NaN;
    const normalized = String(str).trim().replace(/−/g, '-');
    if (!normalized) return NaN;
    return Number(normalized);
  }

  function detectDelimiter(text) {
    if (!text) return ',';
    const sample = text.slice(0, 4096);
    const comma = (sample.match(/,/g) || []).length;
    const semicolon = (sample.match(/;/g) || []).length;
    const tab = (sample.match(/\t/g) || []).length;
    if (tab > comma && tab > semicolon) return '\t';
    if (semicolon > comma && semicolon > tab) return ';';
    return ',';
  }

  function parseCsvText(text, opts) {
    opts = opts || {};
    if (text == null) return { rows: [] };
    const raw = String(text).replace(/\r\n/g, '\n').replace(/\r/g, '\n');
    const delimiter = opts.delimiter || detectDelimiter(raw);
    const tempRows = [];
    let lineCells = [];
    let cur = '';
    let inQuotes = false;
    for (let i = 0; i < raw.length; i++) {
      const ch = raw[i];
      if (ch === '"') {
        if (inQuotes && raw[i + 1] === '"') {
          cur += '"';
          i++;
          continue;
        }
        inQuotes = !inQuotes;
        continue;
      }
      if (ch === '\n' && !inQuotes) {
        lineCells.push(cur);
        tempRows.push(lineCells);
        lineCells = [];
        cur = '';
        continue;
      }
      if (ch === delimiter && !inQuotes) {
        lineCells.push(cur);
        cur = '';
        continue;
      }
      cur += ch;
    }
    if (cur.length || raw[raw.length - 1] === '\n' || raw[raw.length - 1] === delimiter) {
      lineCells.push(cur);
      tempRows.push(lineCells);
    }
    const outRows = tempRows.filter(function (r) {
      if (!r.length) return false;
      const joined = r.join('').trim();
      if (!joined) return false;
      if (joined.indexOf('#') === 0) return false;
      return true;
    }).map(function (r) {
      return r.map(function (c) {
        return String(c).trim();
      });
    });
    return { rows: outRows, delimiter: delimiter };
  }

  function rowsToXY(rows, opts) {
    opts = opts || {};
    if (!rows || !rows.length) return { error: 'importFormat' };
    let dataRows = rows.slice();
    const first = dataRows[0];
    if (first && first.length >= 2) {
      const h0 = String(first[0]).toLowerCase().replace(/\s/g, '');
      const h1 = String(first[1]).toLowerCase().replace(/\s/g, '');
      if ((h0 === 'x' || h0.indexOf('x') === 0) && (h1 === 'y' || h1.indexOf('y') === 0)) {
        dataRows = dataRows.slice(1);
      } else if ((h0 === 'y' || h0.indexOf('y') === 0) && (h1 === 'x' || h1.indexOf('x') === 0)) {
        dataRows = dataRows.slice(1);
        const swapped = [];
        for (let i = 0; i < dataRows.length; i++) {
          if (dataRows[i] && dataRows[i].length >= 2) {
            swapped.push([dataRows[i][1], dataRows[i][0]]);
          }
        }
        dataRows = swapped;
      }
    }
    const xs = [];
    const ys = [];
    for (let i = 0; i < dataRows.length; i++) {
      const r = dataRows[i];
      if (!r || r.length < 2) continue;
      const sx = String(r[0]).trim().replace(/,/g, '.');
      const sy = String(r[1]).trim().replace(/,/g, '.');
      const nx = Number(sx.replace(/\s/g, ''));
      const ny = Number(sy.replace(/\s/g, ''));
      if (!isFinite(nx) || !isFinite(ny)) continue;
      xs.push(nx);
      ys.push(ny);
    }
    if (!xs.length) return { error: 'importFormat' };
    return { xs: xs, ys: ys };
  }

  function formatVector(arr, digits) {
    return arr.map(function (x) { return fmt(x, digits); }).join(', ');
  }

  function signed(x, digits) {
    if (!isNum(x)) return '—';
    return (x >= 0 ? '+' : '−') + fmt(Math.abs(x), digits == null ? 4 : digits);
  }

  function ordinal(n) {
    return String(n);
  }

  PPS.format = {
    fmt: fmt,
    getDecimals: getDecimals,
    setDecimals: setDecimals,
    superscript: superscript,
    subscript: subscript,
    factorial: factorial,
    polyToString: polyToString,
    trimCoefs: trimCoefs,
    parseList: parseList,
    parseNumber: parseNumber,
    formatVector: formatVector,
    signed: signed,
    isNum: isNum,
    parseCsvText: parseCsvText,
    detectDelimiter: detectDelimiter,
    rowsToXY: rowsToXY
  };
})();
