/* ============================================================
   View: Inverse Interpolation (Section 7 - Algorithm.md)
   Method 1: Inverse Function Interpolation x = P(y)
   Method 2: Fixed-point iteration loop t = phi(t)
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;

  const STATUS_CLASS = {
    CONVERGED_CERTIFIED: 'green',
    CANDIDATE_NUMERIC: 'cyan',
    OUT_OF_BRACKET: 'orange',
    MAX_ITER: 'orange',
    DIVERGED: 'red',
    ZERO_DIFF: 'red'
  };

  const STATUS_KEY = {
    CONVERGED_CERTIFIED: 'inverse.statusCertified',
    CANDIDATE_NUMERIC: 'inverse.statusCandidate',
    OUT_OF_BRACKET: 'inverse.statusOutOfBracket',
    MAX_ITER: 'inverse.statusMaxIter',
    DIVERGED: 'inverse.statusDiverged',
    ZERO_DIFF: 'inverse.statusZeroDiff'
  };

  const CONCLUSION_KEY = {
    init: 'inverse.stepInit',
    continue: 'inverse.stepContinue',
    certified: 'inverse.stepCertified',
    candidate: 'inverse.stepCandidate',
    diverged: 'inverse.stepDiverged'
  };

  function statusClass(status) {
    return STATUS_CLASS[status] || 'dim';
  }

  function statusLabel(status) {
    return STATUS_KEY[status] ? t(STATUS_KEY[status]) : String(status || '—');
  }

  function methodLabel(method) {
    if (method === 'forward') return t('inverse.methodForward');
    if (method === 'backward') return t('inverse.methodBackward');
    return '—';
  }

  function conclusionLabel(key) {
    return CONCLUSION_KEY[key] ? t(CONCLUSION_KEY[key]) : String(key || '');
  }

  function valueOrDash(value, digits) {
    return value == null || !isFinite(value) ? '—' : fmt(value, digits);
  }

  function renderIterationDetail(res) {
    const nodes = [];

    if (res.exact) {
      nodes.push(ui.note(t('inverse.exactNode')));
      return nodes;
    }

    if (res.bracket && res.bracket.x0 != null) {
      nodes.push(ui.note(t('inverse.bracketFound', {
        x0: fmt(res.bracket.x0, 4),
        x1: fmt(res.bracket.x1, 4),
        y0: fmt(res.bracket.y0, 4),
        y1: fmt(res.bracket.y1, 4)
      })));
    }

    const kpis = [
      { key: 'x*', value: valueOrDash(res.xbar, 8), className: 'accent' },
      { key: 't*', value: valueOrDash(res.t, 8) },
      { key: 'h', value: valueOrDash(res.h, 6) },
      { key: t('inverse.direction'), value: methodLabel(res.method) },
      { key: t('inverse.iterations'), value: String(res.iterations) },
      { key: t('inverse.checkY'), value: valueOrDash(res.checkY, 8) },
      { key: t('inverse.residual'), value: valueOrDash(res.residual, 8), className: 'green' },
      { key: t('inverse.status'), value: statusLabel(res.status), className: statusClass(res.status) }
    ];
    if (res.boundX != null) {
      kpis.push({ key: t('inverse.boundX'), value: '≤ ' + valueOrDash(res.boundX, 8), className: 'green' });
    }
    nodes.push(ui.kv(kpis));

    if (res.attempts && res.attempts.length > 1) {
      const attemptRows = res.attempts.map(function (a) {
        return [
          { text: methodLabel(a.method), class: 'rowhead' },
          { text: valueOrDash(a.xbar, 8) },
          { text: String(a.iterations) },
          { text: valueOrDash(a.residual, 8) },
          { text: statusLabel(a.status), class: statusClass(a.status) },
          { text: a.valid ? t('inverse.validYes') : t('inverse.validNo'), class: a.valid ? 'green' : 'orange' }
        ];
      });
      nodes.push(ui.section(t('inverse.attemptsTitle'), ui.table([
        t('inverse.direction'),
        t('inverse.colRoot'),
        t('inverse.iterations'),
        t('inverse.colResidual'),
        t('inverse.colStatus'),
        t('inverse.valid')
      ], attemptRows)));
    }

    if (res.stencil) {
      const st = res.stencil;
      const pointRows = st.xs.map(function (x, i) {
        const isAnchor = (res.method === 'forward' && i === 0) ||
          (res.method === 'backward' && i === st.xs.length - 1);
        return [
          { text: String(i), class: 'rowhead' },
          { text: fmt(x, 8), class: isAnchor ? 'hl' : '' },
          { text: fmt(st.ys[i], 8) }
        ];
      });
      nodes.push(ui.section(
        t('inverse.stencilTitle'),
        ui.table([t('common.node'), 'x_i', 'y_i'], pointRows),
        t('inverse.stencilRange', { start: st.start, end: st.end })
      ));
    }

    if (res.diffTable) {
      const D = res.diffTable;
      const m = D[0].length;
      const headers = ['i', 'y_i'];
      for (let j = 1; j < m; j++) headers.push('Δ^' + j + ' y');
      const rows = [];
      for (let i = 0; i < m; i++) {
        const cells = [
          { text: String(i), class: 'rowhead' },
          { text: fmt(D[0][i], 8) }
        ];
        for (let j = 1; j < m; j++) {
          const v = D[j][i];
          const isC1 = j === 1 && (res.method === 'forward' ? i === 0 : i === m - 2);
          cells.push({
            text: v == null ? '—' : fmt(v, 8),
            class: isC1 ? 'hl' : ''
          });
        }
        rows.push(cells);
      }
      nodes.push(ui.section(t('inverse.diffTableTitle'), ui.table(headers, rows)));
    }

    if (res.cDiff && res.cDivided) {
      const isFwd = res.method === 'forward';
      const diffLabel = isFwd ? t('inverse.vectorFiniteForward') : t('inverse.vectorFiniteBackward');
      const diffName = isFwd ? 'C_tiến' : 'C_lùi';
      const diffStr = diffName + ' = [' + res.cDiff.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
      const normStr = 'D = [' + res.cNormalized.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
      const divNodesStr = '[' + res.dividedNodes.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
      const divCoeffsStr = '[' + res.cDivided.map(function (v) { return fmt(v, 6); }).join(', ') + ']';

      nodes.push(ui.section(t('inverse.diffVectorsTitle'), [
        ui.el('p', { class: 'math-label', text: diffLabel }),
        ui.formula(diffStr, true),
        ui.el('p', { class: 'math-label', text: t('inverse.vectorNormalized'), style: 'margin-top: 8px;' }),
        ui.formula(normStr, true),
        ui.el('p', { class: 'math-label', text: t('inverse.vectorDivided'), style: 'margin-top: 8px;' }),
        ui.el('p', { class: 'math-label dim', text: t('inverse.dividedNodesOrder') + ' ' + divNodesStr }),
        ui.formula(divCoeffsStr, true)
      ]));
    }

    if (res.productTable && res.productTable.length) {
      const m = res.productTable.length;
      const prodHeaders = [t('inverse.factorCol')];
      for (let d = m - 1; d >= 0; d--) {
        prodHeaders.push(d === 0 ? 't⁰' : (d === 1 ? 't' : 't' + PPS.format.superscript(d)));
      }
      const prodRows = res.productTable.map(function (row) {
        const cells = [
          { text: row.factor || '1', class: 'rowhead left' }
        ];
        row.polyDesc.forEach(function (coefVal) {
          cells.push({ text: fmt(coefVal, 6) });
        });
        return cells;
      });
      const prodTableEl = ui.table(prodHeaders, prodRows);

      const stepHeaders = [
        t('inverse.termCol'),
        t('inverse.diffCol'),
        t('inverse.dCol'),
        t('inverse.biCol'),
        t('inverse.niCol')
      ];
      const stepRows = res.productTable.map(function (row, idx) {
        return [
          { text: String(idx), class: 'rowhead' },
          { text: fmt(row.diff, 6) },
          { text: fmt(row.coef, 6), class: 'hl' },
          { text: row.factorStr || '1', class: 'left' },
          { text: PPS.format.polyToString(row.termPoly, 't'), class: 'left' }
        ];
      });
      const stepTableEl = ui.table(stepHeaders, stepRows);

      const prodSectionChildren = [
        prodTableEl,
        ui.el('p', { class: 'math-label', text: t('inverse.productTableStepsTitle'), style: 'margin-top: 14px;' }),
        stepTableEl
      ];
      if (res.polyT) {
        prodSectionChildren.push(
          ui.el('p', { class: 'math-label', text: t('inverse.canonicalPolyT'), style: 'margin-top: 14px;' }),
          ui.formula('P(t) = ' + PPS.format.polyToString(res.polyT, 't'), true)
        );
      }
      nodes.push(ui.section(t('inverse.productTableTitleT'), prodSectionChildren));
    }

    if (res.stencil && res.diffTable && res.c1 != null) {
      const D = res.diffTable;
      const m = res.stencil.xs.length;
      const parts = [];
      for (let r = 2; r < m; r++) {
        const delta = res.method === 'forward' ? D[r][0] : D[r][m - 1 - r];
        const coef = '(' + fmt(delta, 8) + ' / (' + r + '! · ' + fmt(res.c1, 8) + '))';
        const factors = [];
        for (let j = 0; j < r; j++) {
          if (j === 0) factors.push('t');
          else if (res.method === 'forward') factors.push('(t − ' + j + ')');
          else factors.push('(t + ' + j + ')');
        }
        parts.push(coef + ' · ' + factors.join(''));
      }
      const anchorName = res.method === 'forward' ? 'y_0' : 'y_' + (m - 1);
      const t0Text = 't⁰ = (y* − ' + anchorName + ') / c₁ = ' + fmt(res.t0, 10);
      const phiText = parts.length
        ? 'φ(t) = ' + fmt(res.t0, 10) + ' − [ ' + parts.join(' + ') + ' ]'
        : 'φ(t) = ' + fmt(res.t0, 10);
      const wrap = ui.el('div');
      wrap.appendChild(ui.formula(t0Text, true));
      wrap.appendChild(ui.formula(phiText, true));
      nodes.push(ui.section(t('inverse.phiTitle'), wrap));
    }

    if (res.history && res.history.length) {
      const headers = [
        t('inverse.colStep'),
        t('inverse.colT'),
        t('inverse.colX'),
        '|Δt|',
        '|Δx|',
        t('inverse.colResidual'),
        t('inverse.colConclusion')
      ];
      const traceRows = res.history.map(function (row) {
        return [
          { text: String(row.step), class: 'rowhead' },
          { text: valueOrDash(row.t, 8), class: 'hl' },
          { text: valueOrDash(row.x, 8) },
          { text: valueOrDash(row.diffT, 8) },
          { text: valueOrDash(row.diffX, 8) },
          { text: valueOrDash(row.residual, 8) },
          { text: conclusionLabel(row.conclusion) }
        ];
      });
      nodes.push(ui.section(
        t('inverse.traceTitle'),
        ui.table(headers, traceRows),
        res.history.length + ' ' + t('common.steps').toLowerCase()
      ));
    }

    return nodes;
  }

  function renderFunctionDetail(res) {
    const nodes = [];

    if (res.exact) {
      nodes.push(ui.note(t('inverse.exactNode')));
      return nodes;
    }

    nodes.push(ui.kv([
      { key: 'x* = P(ȳ)', value: valueOrDash(res.xbar, 8), className: 'accent' },
      { key: 'Residual |P(x*) - ȳ|', value: valueOrDash(res.residual, 8), className: 'green' },
      { key: "(f⁻¹)'(ȳ) xấp xỉ", value: res.derivative != null ? fmt(res.derivative, 6) : '—' },
      { key: t('inverse.status'), value: res.inBracket ? t('inverse.inBracketYes') : t('inverse.inBracketNo'), className: res.inBracket ? 'green' : 'orange' }
    ]));

    if (res.bracket) {
      nodes.push(ui.note(
        t('inverse.bracketFound', {
          x0: fmt(res.bracket.x0, 4),
          x1: fmt(res.bracket.x1, 4),
          y0: fmt(res.bracket.y0, 4),
          y1: fmt(res.bracket.y1, 4)
        })
      ));
    }

    if (!res.exact) {
      if (res.stencil) {
        const st = res.stencil;
        if (res.nPoints !== 'all' && st.count < res.nPoints) {
          nodes.push(ui.note(
            t('inverse.stencilClipped', { asked: res.nPoints, used: st.count })
          ));
        }
        const pointRows = st.xs.map(function (x, i) {
          return [
            { text: String(st.start + i), class: 'rowhead' },
            { text: fmt(x, 8) },
            { text: fmt(st.ys[i], 8) }
          ];
        });
        nodes.push(ui.section(
          t('inverse.stencilTitle'),
          ui.table([t('common.node'), 'x_i', 'y_i'], pointRows),
          t('inverse.stencilRange', { start: st.start, end: st.end }) +
            ' · ' + t('inverse.stencilUsed') + ': ' + st.count
        ));
      }

      if (res.fallback) {
        nodes.push(ui.note(t('inverse.fallbackLinear', {
          poly: valueOrDash(res.polyX, 8),
          x0: fmt(res.bracket.x0, 6),
          x1: fmt(res.bracket.x1, 6)
        })));
      }

      // Bảng tỷ sai phân theo biến Y
      if (res.tableY && res.subX) {
        const m = res.subX.length;
        const headers = [t('common.node'), 'Y_i (Biến)', 'X_i = F[i, 0] (Hàm)'];
        for (let j = 1; j < m; j++) headers.push('F[i, ' + j + ']');
        const rows = [];
        for (let i = 0; i < m; i++) {
          const cells = [
            { text: String(i), class: 'rowhead' },
            { text: fmt(res.subY[i], 6), class: 'hl' },
            { text: fmt(res.subX[i], 6) }
          ];
          for (let j = 1; j < m; j++) {
            const val = res.tableY[i][j];
            cells.push({
              text: val != null ? fmt(val, 6) : '—',
              class: i === 0 && val != null ? 'hl' : ''
            });
          }
          rows.push(cells);
        }
        nodes.push(ui.section(t('inverse.tableYTitle'), ui.table(headers, rows)));
      }

      // Vector hệ số tỷ sai phân theo biến Y (Mục 7.2 Algorithm.md)
      if (res.cForward && res.cForward.length) {
        const fwdStr = 'd_tiến = [' + res.cForward.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const vectorChildren = [
          ui.el('p', { class: 'math-label', text: t('inverse.vectorDividedY') }),
          ui.formula(fwdStr, true)
        ];
        if (res.cBackward && res.cBackward.length) {
          const bwdStr = 'd_lùi = [' + res.cBackward.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
          vectorChildren.push(
            ui.el('p', { class: 'math-label', text: t('inverse.vectorDividedYBwd'), style: 'margin-top: 8px;' }),
            ui.formula(bwdStr, true)
          );
        }
        nodes.push(ui.section(t('inverse.diffVectorsTitle'), vectorChildren));
      }

      // Bảng tích lũy thừa B_i(y) (Product Table) theo biến y
      if (res.productTable && res.productTable.length) {
        const m = res.productTable.length;
        const prodHeaders = [t('common.node') + ' Y_k'];
        for (let d = m - 1; d >= 0; d--) {
          prodHeaders.push(d === 0 ? 'y⁰' : (d === 1 ? 'y' : 'y' + PPS.format.superscript(d)));
        }
        const prodRows = res.productTable.map(function (row) {
          const cells = [
            { text: row.node === null ? '—' : fmt(row.node, 4), class: 'rowhead' }
          ];
          row.polyDesc.forEach(function (coefVal) {
            cells.push({ text: fmt(coefVal, 6) });
          });
          return cells;
        });
        const prodTableEl = ui.table(prodHeaders, prodRows);
        const prodSectionChildren = [prodTableEl];
        if (res.polyCanonical) {
          prodSectionChildren.push(
            ui.el('p', { class: 'math-label', text: t('inverse.canonicalPolyY'), style: 'margin-top: 14px;' }),
            ui.formula('x = P(y) = ' + PPS.format.polyToString(res.polyCanonical, 'y'), true)
          );
        }
        nodes.push(ui.section(t('inverse.productTableTitleY'), prodSectionChildren));
      }
    }

    return nodes;
  }

  const PRESETS = [
    // Ví dụ mẫu Mục 7.3.2 (Lặp điểm bất động)
    {
      mode: 'iteration',
      xs: [2.380, 2.495, 2.610, 2.725],
      ys: [-1.5819, -1.5127, -1.3977, -1.3814], // Δ³y₀ = -0.1445 như Algorithm.md
      ybar: -1.43,
      eps: 1e-6,
      maxIter: 40,
      qCoeff: 0.75,
      nPoints: 4
    },
    // Ví dụ mẫu 11_Loop_Inverse/input.txt (nhiều khoảng cách ly)
    {
      mode: 'iteration',
      xs: [-1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00],
      ys: [-1.121982, -1.214072, -1.256192, -1.254842, -1.214912, -1.139932, -1.032282, -0.523962, -0.292592],
      ybar: -1.25,
      eps: 1e-8,
      maxIter: 100,
      qCoeff: '',
      nPoints: 'all'
    },
    // Ví dụ mẫu Mục 7.2.3 (Hàm ngược)
    {
      mode: 'function',
      xs: [1.40, 1.50, 1.60, 1.70],
      ys: [2.44018, 2.55290, 2.65511, 2.74677],
      ybar: 2.60,
      eps: 1e-6,
      maxIter: 50,
      qCoeff: '',
      nPoints: 4
    },
    {
      mode: 'iteration',
      xs: [0, 1, 2, 3, 4],
      ys: [0, 1, 4, 9, 16],
      ybar: 6.25,
      eps: 1e-7,
      maxIter: 50,
      qCoeff: '',
      nPoints: 4
    }
  ];

  function render(root) {
    let presetIndex = 0;

    const modeSel = ui.select([
      { value: 'iteration', label: t('inverse.modeIter') },
      { value: 'function', label: t('inverse.modeFunc') }
    ], PRESETS[0].mode);

    const xy = ui.xyControl({
      xs: PRESETS[0].xs,
      ys: PRESETS[0].ys,
      xLabel: t('inverse.xs'),
      yLabel: t('inverse.ys'),
      xHint: t('lagrange.xsHint'),
      yHint: t('lagrange.xsHint')
    });
    const ybarIn = ui.input({ value: String(PRESETS[0].ybar) });
    const epsIn = ui.input({ value: String(PRESETS[0].eps) });
    const iterIn = ui.input({ value: String(PRESETS[0].maxIter) });
    const qIn = ui.input({ value: String(PRESETS[0].qCoeff || '') });
    const nPointsIn = ui.input({ value: String(PRESETS[0].nPoints || 4) });

    const methodSel = ui.select([
      { value: 'auto', label: t('common.auto') },
      { value: 'forward', label: t('common.forward') },
      { value: 'backward', label: t('common.backward') }
    ], 'auto');

    const out = ui.el('div');

    const form = ui.card(t('common.input'),
      xy.node,
      ui.grid(
        ui.field(t('inverse.approach'), modeSel),
        ui.field(t('inverse.ybar'), ybarIn),
        ui.field(t('inverse.eps'), epsIn),
        ui.field(t('inverse.maxIter'), iterIn),
        ui.field(t('inverse.nPoints'), nPointsIn, t('inverse.nPointsHint')),
        ui.field(t('inverse.qCoeff'), qIn, t('inverse.qCoeffHint')),
        ui.field(t('inverse.method'), methodSel)
      ),
      ui.actions(
        ui.button('▶ ' + t('common.run'), run, { variant: 'primary' }),
        ui.button('✨ ' + t('common.example'), function () {
          presetIndex = (presetIndex + 1) % PRESETS.length;
          const p = PRESETS[presetIndex];
          modeSel.value = p.mode;
          xy.setValues(p.xs, p.ys);
          ybarIn.value = String(p.ybar);
          epsIn.value = String(p.eps);
          iterIn.value = String(p.maxIter);
          qIn.value = String(p.qCoeff || '');
          nPointsIn.value = String(p.nPoints || 4);
          run();
        }, { variant: 'ghost' }),
        (function(){
          const imp = ui.importDataButton({ onParsed: function(res){ xy.setValues(res.xs||[], res.ys||[]); run(); } });
          imp.node.style.marginLeft = 'auto';
          return imp.node;
        })()
      )
    );

    root.appendChild(form);
    root.appendChild(out);

    function run() {
      out.innerHTML = '';
      try {
        const { xs, ys } = xy.getValues();
        const ybar = ui.parseNumber(ybarIn.value);
        const mode = modeSel.value;
        const eps = epsIn.value.trim() ? ui.parseNumber(epsIn.value) : 1e-8;
        const maxIter = iterIn.value.trim() ? Math.floor(ui.parseNumber(iterIn.value)) : 100;
        const qCoeff = qIn.value.trim() ? ui.parseNumber(qIn.value) : null;
        const method = methodSel.value;
        const nPointsRaw = nPointsIn.value.trim();
        let nPoints = 4;
        if (nPointsRaw.toLowerCase() === 'all') nPoints = 'all';
        else if (nPointsRaw !== '') nPoints = Math.max(2, Math.floor(ui.parseNumber(nPointsRaw)));

        // Tính TẤT CẢ các khoảng cách ly nghiệm
        const all = PPS.inverse.solveAll(xs, ys, ybar, {
          mode: mode,
          method: method,
          nPoints: nPoints,
          tolT: eps,
          tolRes: eps,
          tolX: eps,
          maxIter: maxIter,
          qCoeff: qCoeff
        });

        const totalRoots = all.length;

        // BẢNG TỔNG HỢP MỐC SỬ DỤNG VÀ KHOẢNG CÁCH LY (CHUẨN ĐỀ THI)
        const allExamRows = [];
        all.forEach(function (res, resIdx) {
          const rootLabel = all.length > 1 ? ('x*' + (PPS.central ? PPS.central.sub(resIdx + 1) : ('_' + (resIdx + 1)))) : 'x*';
          if (res.exact) {
            allExamRows.push([
              { text: rootLabel, class: 'rowhead hl' },
              { text: '[' + fmt(res.bracket.x, 3) + ', ' + fmt(res.bracket.x, 3) + ']' },
              { text: '[' + fmt(res.bracket.x, 3) + ']' },
              { text: 'Trùng mốc x' + (res.nodeIndex != null && PPS.central ? PPS.central.sub(res.nodeIndex) : '') },
              { text: fmt(res.xbar, 6), class: 'accent hl', style: 'font-family: monospace; font-weight: bold;' }
            ]);
          } else if (res.comparisons && res.comparisons.length) {
            res.comparisons.forEach(function (cmp) {
              allExamRows.push([
                { text: rootLabel, class: 'rowhead hl' },
                { text: cmp.bracketStr },
                { text: cmp.stencilStr },
                { text: cmp.methodName },
                { text: fmt(cmp.xbar, 6), class: 'accent hl', style: 'font-family: monospace; font-weight: bold;' }
              ]);
            });
          } else {
            const bStr = res.bracket ? ('[' + fmt(res.bracket.x0, 3) + ', ' + fmt(res.bracket.x1, 3) + ']') : '—';
            const sStr = res.stencil && res.stencil.xs ? ('[' + fmt(res.stencil.xs[0], 3) + ', ' + fmt(res.stencil.xs[res.stencil.xs.length - 1], 3) + ']') : '—';
            allExamRows.push([
              { text: rootLabel, class: 'rowhead hl' },
              { text: bStr },
              { text: sStr },
              { text: res.mode === 'function' ? 'Nội suy hàm ngược x = P(y)' : methodLabel(res.method) },
              { text: fmt(res.xbar, 6), class: 'accent hl', style: 'font-family: monospace; font-weight: bold;' }
            ]);
          }
        });

        if (allExamRows.length) {
          const examTable = ui.table([
            t('inverse.colRootSymbol'),
            t('inverse.colBracket'),
            t('inverse.colStencil'),
            t('inverse.colMethodApplied'),
            t('inverse.colResult')
          ], allExamRows);
          out.appendChild(ui.section(t('inverse.examTableTitle'), examTable));
        }

        if (totalRoots === 1 && mode === 'function') {
          // CÁCH 1: Nội suy hàm ngược x = P(y)
          renderFunctionDetail(all[0]).forEach(function (node) {
            out.appendChild(node);
          });
        } else if (totalRoots === 1) {
          // CÁCH 2: đúng một nghiệm — hiển thị đầy đủ mốc, sai phân, φ(t), trace
          renderIterationDetail(all[0]).forEach(function (node) {
            out.appendChild(node);
          });
        } else {
          // Nhiều nghiệm: bảng tổng hợp + chi tiết từng khoảng
          out.appendChild(ui.note(t('inverse.multipleRoots', { n: totalRoots })));

          const tableRows = all.map(function (r, i) {
            let statusText;
            let cls;
            if (r.exact) {
              statusText = t('inverse.statusExact');
              cls = 'green';
            } else if (r.mode === 'function') {
              statusText = r.inBracket ? t('inverse.inBracketYes') : t('inverse.inBracketNo');
              if (r.fallback) statusText += ' · ' + t('inverse.fallbackShort');
              cls = r.inBracket ? 'green' : 'orange';
            } else {
              statusText = statusLabel(r.status);
              cls = statusClass(r.status);
            }
            const bracketX = r.exact
              ? fmt(r.bracket.x, 4)
              : fmt(r.bracket.x0, 4) + ' → ' + fmt(r.bracket.x1, 4);
            const bracketY = r.exact
              ? fmt(r.bracket.y, 4)
              : fmt(r.bracket.y0, 4) + ' → ' + fmt(r.bracket.y1, 4);
            return [
              { text: String(i + 1), class: 'rowhead' },
              { text: bracketX, class: 'left' },
              { text: bracketY, class: 'left' },
              { text: r.mode === 'function' ? 'P(y)' : methodLabel(r.method) },
              { text: r.stencil ? (r.stencil.start + '→' + r.stencil.end) : '—' },
              { text: valueOrDash(r.xbar, 8), class: 'hl accent' },
              { text: valueOrDash(r.checkY, 8) },
              { text: valueOrDash(r.residual, 8) },
              { text: r.iterations == null ? '—' : String(r.iterations) },
              { text: statusText, class: cls }
            ];
          });
          out.appendChild(ui.section(t('inverse.rootTable'),
            ui.table([
              t('inverse.colNo'),
              t('inverse.colBracketX'),
              t('inverse.colBracketY'),
              t('inverse.direction'),
              t('inverse.stencilShort'),
              t('inverse.colRoot'),
              t('inverse.checkY'),
              t('inverse.colResidual'),
              t('inverse.iterations'),
              t('inverse.colStatus')
            ], tableRows)));

          if (mode === 'iteration') {
            all.forEach(function (r, i) {
              if (r.exact) return;
              out.appendChild(ui.section(
                t('inverse.bracketDetail', { n: i + 1, method: methodLabel(r.method) }),
                renderIterationDetail(r)
              ));
            });
          } else {
            all.forEach(function (r, i) {
              if (r.exact) return;
              out.appendChild(ui.section(
                t('inverse.bracketDetail', { n: i + 1, method: 'P(y)' }),
                renderFunctionDetail(r)
              ));
            });
          }
        }

        // --- Đồ thị trực quan: đánh dấu TẤT CẢ các nghiệm ---
        const chart = ui.chart(t('common.chart'), 380);
        const minX = Math.min.apply(null, xs);
        const maxX = Math.max.apply(null, xs);
        const margin = Math.max(0.5, (maxX - minX) * 0.15);

        let fwdInterp = null;
        try { fwdInterp = PPS.newton.interpolate(xs, ys); } catch (e) { /* ignore */ }

        const curve = fwdInterp ? PPS.poly.sample(fwdInterp.coef, minX - margin, maxX + margin, 250) : [];
        const nodeMarkers = xs.map(function (x, i) {
          return { x: x, y: ys[i], color: 'var(--accent)', r: 5, label: '(' + fmt(x, 2) + ', ' + fmt(ys[i], 2) + ')' };
        });

        // Các điểm nghiệm x* nằm trên đường y = ybar (bỏ qua ứng viên không hữu hạn)
        const rootMarkers = all.filter(function (r) {
          return r.xbar != null && isFinite(r.xbar);
        }).map(function (r, i) {
          return {
            x: r.xbar,
            y: ybar,
            color: 'var(--green)',
            r: 6,
            label: 'x*' + (i + 1) + ' = ' + fmt(r.xbar, 5)
          };
        });

        chart.draw({
          series: [
            { points: curve, color: 'var(--accent)', width: 2.2, label: 'P(x) nội suy' },
            {
              points: [{ x: minX - margin, y: ybar }, { x: maxX + margin, y: ybar }],
              color: 'var(--orange)',
              width: 1.8,
              dash: '6 4',
              label: 'ȳ = ' + fmt(ybar, 4)
            }
          ],
          points: nodeMarkers.concat(rootMarkers),
          highlightPoints: rootMarkers,
          yLabel: 'y',
          xLabel: 'x'
        });

        out.appendChild(chart.node);
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
  }

  PPS.views.inverse = { render: render };
})();
