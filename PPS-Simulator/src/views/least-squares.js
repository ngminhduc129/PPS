/* ============================================================
   View: least squares — linear-in-parameters, log-linear, and power law
   Updated based on Section 9 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;

  const LINEAR_PRESET = {
    xs: [0, 1, 2, 3, 4, 5],
    ys: [1.1, 2.9, 4.8, 7.2, 8.9, 11.2],
    basis: ['one', 'x']
  };
  const LOG_PRESET = {
    xs: [0, 1, 2, 3, 4],
    ys: [3, 4.297, 6.436, 10.482, 18.182],
    basis: ['x'],
    c0: 1
  };
  const POWER_PRESET = {
    xs: [1, 2, 3, 4, 5],
    ys: [2.1, 7.9, 18.2, 31.8, 50.1] // y ≈ 2 * x^2
  };

  function basisOptions() {
    return PPS.leastSquares.BASIS.map(function (b) {
      return { value: b.key, label: t(b.labelKey) };
    });
  }

  function joinTerms(terms) {
    if (!terms.length) return '0';
    let out = '';
    terms.forEach(function (term, index) {
      const negative = term.charAt(0) === '−';
      const body = negative ? term.slice(1) : term;
      if (index === 0) out += (negative ? '−' : '') + body;
      else out += negative ? ' − ' + body : ' + ' + body;
    });
    return out;
  }

  function linearModel(result) {
    const terms = [];
    result.basis.forEach(function (b, j) {
      const c = result.coef[j];
      if (Math.abs(c) < 1e-12) return;
      const label = t(b.labelKey);
      const abs = Math.abs(c);
      const sign = c < 0 ? '−' : '';
      if (label === '1') terms.push(sign + fmt(abs, 6));
      else if (Math.abs(abs - 1) < 1e-12) terms.push(sign + label);
      else terms.push(sign + fmt(abs, 6) + '·' + label);
    });
    return 'g(x) = ' + joinTerms(terms);
  }

  function logModel(result) {
    const terms = [];
    result.basis.forEach(function (b, j) {
      const c = result.b[j];
      if (Math.abs(c) < 1e-12) return;
      const label = t(b.labelKey);
      const abs = Math.abs(c);
      const sign = c < 0 ? '−' : '';
      if (Math.abs(abs - 1) < 1e-12) terms.push(sign + label);
      else terms.push(sign + fmt(abs, 6) + '·' + label);
    });
    const expPart = terms.length ? 'exp(' + joinTerms(terms) + ')' : '1';
    let out = fmt(result.a, 6) + '·' + expPart;
    if (Math.abs(result.c0) > 1e-12) {
      out += (result.c0 > 0 ? ' + ' : ' − ') + fmt(Math.abs(result.c0), 6);
    }
    return 'g(x) = ' + out;
  }

  function matrixTable(mat, rowLabels, colLabels) {
    const headers = [''].concat(colLabels);
    const rows = mat.map(function (row, i) {
      return [{ text: rowLabels[i], class: 'rowhead' }].concat(row.map(function (v) {
        return { text: fmt(v, 6) };
      }));
    });
    return ui.table(headers, rows);
  }

  function buildMatricesSection(result, xs, labels, yLabel) {
    yLabel = yLabel || 'y';
    const children = [];

    // 1. Table Phi (n x m)
    const phiHeaders = ['i', 'xᵢ'].concat(labels);
    const phiRows = result.phi.map(function (row, i) {
      return [{ text: String(i), class: 'rowhead' }, { text: fmt(xs[i], 5) }].concat(row.map(function (v) {
        return { text: fmt(v, 6) };
      }));
    });
    const phiTable = ui.table(phiHeaders, phiRows);

    // 2. Table M = Phi^T Phi (m x m)
    const mTable = matrixTable(result.M, labels, labels);

    // Grid for Phi and M
    const grid1 = ui.el('div', { class: 'two-col' },
      ui.el('div', {}, ui.el('div', { class: 'section-title', text: t('lsq.phi') }), phiTable),
      ui.el('div', {}, ui.el('div', { class: 'section-title', text: t('lsq.M') }), mTable)
    );
    children.push(grid1);

    // 3. RHS vector B = Phi^T y and detailed components
    const rhsHeaders = ['j', t('lsq.basisCol'), t('lsq.sumFormulaCol'), t('lsq.rhsValCol'), t('lsq.paramSolCol')];
    const rhsRows = labels.map(function (label, j) {
      const formulaStr = label === '1' ? 'Σ ' + yLabel + 'ᵢ' : 'Σ ' + label + ' · ' + yLabel + 'ᵢ';
      const cName = 'c' + PPS.central.sub(j + 1);
      return [
        { text: String(j + 1), class: 'rowhead' },
        { text: label, class: 'hl' },
        { text: formulaStr },
        { text: fmt(result.rhs[j], 8), class: 'accent', style: 'font-family: monospace; font-weight: bold;' },
        { text: cName + ' = ' + fmt(result.coef[j], 8), style: 'font-family: monospace;' }
      ];
    });
    const rhsTable = ui.table(rhsHeaders, rhsRows);

    // 4. Explicit Algebraic Normal Equations System
    const eqLines = labels.map(function (_, j) {
      const terms = [];
      labels.forEach(function (lbl, k) {
        const val = result.M[j][k];
        const cName = 'c' + PPS.central.sub(k + 1);
        terms.push(fmt(val, 6) + ' · ' + cName);
      });
      return terms.join(' + ') + ' = ' + fmt(result.rhs[j], 6);
    });
    const systemWrap = ui.el('div', { class: 'card', style: 'margin-top: 12px; background: var(--panel-2);' },
      ui.el('p', { class: 'math-label', text: t('lsq.normSystem'), style: 'font-weight: 700; margin-bottom: 6px;' }),
      ui.formula(eqLines.join('\n'), true)
    );

    children.push(ui.el('div', { style: 'margin-top: 14px;' },
      ui.el('div', { class: 'section-title', text: t('lsq.rhs') }),
      rhsTable,
      systemWrap
    ));

    children.push(ui.note(t('lsq.linearNote')));

    return ui.section(t('lsq.matrices'), children);
  }

  function makePresetBar(items) {
    const bar = ui.el('div', { class: 'lsq-preset-bar' });
    bar.appendChild(ui.el('span', { class: 'lsq-preset-label', text: t('lsq.quickPresets') }));
    items.forEach(function (it) {
      const btn = ui.el('button', { class: 'lsq-preset-btn', type: 'button', text: it.label });
      btn.onclick = function (e) {
        e.preventDefault();
        it.onClick();
      };
      bar.appendChild(btn);
    });
    return bar;
  }

  function buildSumTableSection(result, xs, ys, labels, yLabel) {
    yLabel = yLabel || 'y';
    const basisList = result.extended || result.basis;
    const st = PPS.leastSquares.buildSumTable(xs, ys, basisList);
    const headers = st.cols.map(function (c) {
      return c.label.replace('y', yLabel);
    });
    const rows = st.rows.map(function (r) {
      return st.cols.map(function (c, idx) {
        if (idx === 0) return { text: String(r[c.key]), class: 'rowhead' };
        return { text: fmt(r[c.key], 5) };
      });
    });
    const sumRow = st.cols.map(function (c, idx) {
      if (idx === 0) return { text: 'Σ', class: 'rowhead hl', style: 'font-weight: 800;' };
      return { text: fmt(st.sums[c.key], 6), class: 'accent', style: 'font-weight: 700;' };
    });
    rows.push(sumRow);
    return ui.section(t('lsq.sumTableTitle'), ui.table(headers, rows));
  }

  function buildPredictionCard(result, xStar) {
    if (xStar == null || !isFinite(xStar)) return null;
    let yStar;
    try {
      yStar = result.evaluate(xStar);
    } catch (e) {
      return null;
    }
    if (!isFinite(yStar)) return null;
    return {
      node: ui.el('div', { class: 'lsq-pred-card' },
        ui.el('div', { style: 'font-weight: 700; font-size: 13.5px; color: var(--accent); margin-bottom: 4px;' },
          '★ ' + t('lsq.predictionTitle') + ': x* = ' + fmt(xStar, 6)
        ),
        ui.el('div', { style: 'font-size: 16px; font-family: var(--font-mono); font-weight: 800;' },
          'g(' + fmt(xStar, 6) + ') = ' + fmt(yStar, 8)
        )
      ),
      yStar: yStar
    };
  }

  function chartFit(result, nodes, title, xStar, yStar) {
    const xs = nodes.map(function (p) { return p.x; });
    if (xStar != null && isFinite(xStar)) xs.push(xStar);
    const minX = Math.min.apply(null, xs);
    const maxX = Math.max.apply(null, xs);
    const span = Math.max(1, maxX - minX);
    const low = minX - span * 0.12;
    const high = maxX + span * 0.12;

    const curve = [];
    const step = (high - low) / 200;
    for (let i = 0; i <= 200; i++) {
      const x = low + i * step;
      try {
        const y = result.evaluate(x);
        if (isFinite(y)) curve.push({ x: x, y: y });
      } catch (e) { /* domain */ }
    }

    const nodeMarkers = nodes.map(function (p) {
      return { x: p.x, y: p.y, color: 'var(--accent)', r: 5, label: '(' + fmt(p.x, 2) + ', ' + fmt(p.y, 2) + ')' };
    });

    const highlightPoints = [];
    if (xStar != null && isFinite(yStar)) {
      highlightPoints.push({
        x: xStar,
        y: yStar,
        color: 'var(--yellow)',
        r: 7,
        label: 'x* = ' + fmt(xStar, 4) + ' → g(x*) = ' + fmt(yStar, 4)
      });
    }

    const residuals = (result.residuals || []).map(function (r) {
      return { x: r.x, y0: r.y, y1: r.fitted };
    });

    const chart = ui.chart(title, 360);
    const handle = chart.draw({
      series: [{
        points: curve,
        color: 'var(--accent)',
        width: 2.6,
        label: t('lsq.fit')
      }],
      points: nodeMarkers,
      highlightPoints: highlightPoints,
      residuals: residuals,
      xLabel: 'x',
      yLabel: 'y'
    });
    return { node: chart.node, handle: handle, curve: curve };
  }

  function renderLinear(tabRoot) {
    const xy = ui.xyControl({
      xs: LINEAR_PRESET.xs,
      ys: LINEAR_PRESET.ys,
      xLabel: t('lsq.xs'),
      yLabel: t('lsq.ys'),
      xHint: t('lagrange.xsHint'),
      yHint: t('lagrange.xsHint')
    });
    const checks = ui.checkGroup(basisOptions(), LINEAR_PRESET.basis);
    const xStarIn = ui.input({ placeholder: 'Ví dụ: 2.5' });
    const out = ui.el('div');

    const presets = makePresetBar([
      { label: t('lsq.presetLinear'), onClick: function () { checks.setValues(['one', 'x']); run(); } },
      { label: t('lsq.presetParabol'), onClick: function () { checks.setValues(['one', 'x', 'x2']); run(); } },
      { label: t('lsq.presetCubic'), onClick: function () { checks.setValues(['one', 'x', 'x2', 'x3']); run(); } },
      { label: t('lsq.presetHyperbolic'), onClick: function () { checks.setValues(['one', 'inv']); run(); } },
      { label: t('lsq.presetLog'), onClick: function () { checks.setValues(['one', 'ln']); run(); } }
    ]);

    const form = ui.card(t('common.input'),
      presets,
      xy.node,
      ui.grid(
        ui.field(t('lsq.chooseBasis'), checks.node, t('lsq.chooseBasisHint'), { span2: true }),
        ui.field(t('lsq.xStar'), xStarIn, t('lsq.xStarHint'))
      ),
      ui.actions(
        ui.button('▶ ' + t('common.run'), run, { variant: 'primary' }),
        (function(){
          const imp = ui.importDataButton({ onParsed: function(res){ xy.setValues(res.xs||[], res.ys||[]); run(); } });
          imp.node.style.marginLeft = 'auto';
          return imp.node;
        })()
      )
    );

    tabRoot.appendChild(form);
    tabRoot.appendChild(out);

    function run() {
      out.innerHTML = '';
      try {
        const { xs, ys } = xy.getValues();
        const basis = checks.values();
        const result = PPS.leastSquares.linear(xs, ys, basis);
        const labels = result.basis.map(function (b) { return t(b.labelKey); });

        let xStar = null;
        if (xStarIn.value.trim() !== '') {
          xStar = Number(xStarIn.value.trim().replace(',', '.'));
          if (!isFinite(xStar)) xStar = null;
        }

        // Bộ chỉ số sai số đầy đủ (Mục 9.2 - Algorithm.md)
        out.appendChild(ui.kv([
          { key: t('common.sse'), value: fmt(result.sse, 6), className: result.sse < 1e-9 ? 'green' : 'accent' },
          { key: 'MAE', value: fmt(result.mae, 6) },
          { key: 'RMSE', value: fmt(result.rmse, 6) },
          { key: 'Max Error', value: fmt(result.maxError, 6) },
          { key: 'R²', value: fmt(result.r2, 6), className: result.r2 > 0.99 ? 'green' : '' },
          { key: t('lsq.r2adj'), value: fmt(result.r2adj, 6) },
          { key: 'SST', value: fmt(result.sst, 6) }
        ]));

        out.appendChild(ui.section(t('lsq.model'), ui.formula(linearModel(result), true)));

        // Đánh giá dự báo tại điểm x* (nếu có)
        const pred = buildPredictionCard(result, xStar);
        if (pred) {
          out.appendChild(pred.node);
        }

        out.appendChild(ui.note(t('lsq.linearNote')));

        const paramRows = result.coef.map(function (c, j) {
          return [
            { text: 'c' + PPS.central.sub(j + 1), class: 'rowhead' },
            { text: labels[j], class: 'left' },
            { text: fmt(c, 8), class: 'hl' }
          ];
        });
        out.appendChild(ui.section(t('lsq.params'),
          ui.table([t('lsq.paramName'), t('common.basis'), t('common.coefficient')], paramRows)));

        // Bảng tính tổng tích lũy (Summation table phục vụ bài thi giấy)
        out.appendChild(buildSumTableSection(result, xs, ys, labels, 'y'));

        // Trích xuất các Vector kết quả Bình phương tối thiểu
        const cVectorStr = 'c = [' + result.coef.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const bVectorStr = 'B = Φᵀy = [' + result.rhs.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const rVectorStr = 'r = [' + result.residuals.map(function (r) { return fmt(r.residual, 6); }).join(', ') + ']';

        out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), [
          ui.el('p', { class: 'math-label', text: 'Vector nghiệm hệ phương trình chuẩn c (các tham số mô hình):' }),
          ui.formula(cVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector vế phải B = Φᵀy (với Bⱼ = Σ φⱼ(xᵢ)yᵢ):', style: 'margin-top: 8px;' }),
          ui.formula(bVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector phần dư r_i = y_i − g(x_i):', style: 'margin-top: 8px;' }),
          ui.formula(rVectorStr, true)
        ]));

        // Ma trận thiết kế Phi, ma trận hệ chuẩn M, ma trận nghịch đảo và vế phải
        out.appendChild(buildMatricesSection(result, xs, labels, 'y'));

        const residualRows = result.residuals.map(function (r, i) {
          return [
            { text: String(i), class: 'rowhead' },
            { text: fmt(r.x, 6) },
            { text: fmt(r.y, 6) },
            { text: fmt(r.fitted, 6) },
            { text: fmt(r.residual, 4), class: Math.abs(r.residual) < 1e-8 ? 'hl-green' : 'hl-orange' }
          ];
        });
        out.appendChild(ui.section(t('lsq.residuals'),
          ui.table(['i', 'x', 'y', 'g(x)', t('common.residual')], residualRows)));

        const fitInfo = chartFit(result, result.nodes, t('common.chart'), xStar, pred ? pred.yStar : null);
        out.appendChild(ui.section(t('common.chart'), fitInfo.node));
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
    return tabRoot;
  }

  function renderLog(tabRoot) {
    const xy = ui.xyControl({
      xs: LOG_PRESET.xs,
      ys: LOG_PRESET.ys,
      xLabel: t('lsq.xs'),
      yLabel: t('lsq.ys'),
      xHint: t('lagrange.xsHint'),
      yHint: t('lagrange.xsHint')
    });
    const c0In = ui.input({ value: String(LOG_PRESET.c0) });
    const checks = ui.checkGroup(basisOptions().filter(function (b) {
      return b.value !== 'one';
    }), LOG_PRESET.basis);
    const xStarIn = ui.input({ placeholder: 'Ví dụ: 10.0' });
    const out = ui.el('div');

    const presets = makePresetBar([
      { label: t('lsq.presetExpStd'), onClick: function () { c0In.value = '0'; checks.setValues(['x']); run(); } },
      { label: t('lsq.presetExpShift'), onClick: function () { c0In.value = '1'; checks.setValues(['x']); run(); } },
      { label: t('lsq.presetMidterm20251'), onClick: function () { c0In.value = '0'; checks.setValues(['x', 'inv']); run(); } },
      { label: t('lsq.presetMidterm20251_x2'), onClick: function () { c0In.value = '0'; checks.setValues(['x2', 'inv']); run(); } }
    ]);

    const form = ui.card(t('common.input'),
      presets,
      xy.node,
      ui.grid(
        ui.field(t('lsq.c0'), c0In, t('lsq.c0Hint')),
        ui.field(t('lsq.chooseBasis'), checks.node, t('lsq.chooseBasisHint'), { span2: true }),
        ui.field(t('lsq.xStar'), xStarIn, t('lsq.xStarHint'))
      ),
      ui.actions(
        ui.button('▶ ' + t('common.run'), run, { variant: 'primary' }),
        (function(){
          const imp = ui.importDataButton({ onParsed: function(res){ xy.setValues(res.xs||[], res.ys||[]); run(); } });
          imp.node.style.marginLeft = 'auto';
          return imp.node;
        })()
      )
    );

    tabRoot.appendChild(form);
    tabRoot.appendChild(out);

    function run() {
      out.innerHTML = '';
      try {
        const { xs, ys } = xy.getValues();
        const c0 = ui.parseNumber(c0In.value, 'c0');
        const basis = checks.values();
        const result = PPS.leastSquares.logLinear(xs, ys, basis, c0);
        const labels = result.basis.map(function (b) { return t(b.labelKey); });

        let xStar = null;
        if (xStarIn.value.trim() !== '') {
          xStar = Number(xStarIn.value.trim().replace(',', '.'));
          if (!isFinite(xStar)) xStar = null;
        }

        out.appendChild(ui.kv([
          { key: t('lsq.aValue'), value: fmt(result.a, 8), className: 'accent' },
          { key: t('common.sign'), value: result.sign > 0 ? t('common.positive') : t('common.negative'), small: true },
          { key: t('common.sse'), value: fmt(result.sse, 6), className: result.sse < 1e-9 ? 'green' : 'orange' },
          { key: 'RMSE', value: fmt(result.rmse, 6) },
          { key: 'MAE', value: fmt(result.mae, 6) },
          { key: 'R²', value: fmt(result.r2, 6), className: result.r2 > 0.99 ? 'green' : '' },
          { key: t('lsq.r2adj'), value: fmt(result.r2adj, 6) }
        ]));

        out.appendChild(ui.section(t('lsq.originalModel'),
          ui.formula(logModel(result), true)));
        out.appendChild(ui.section(t('lsq.logModel'),
          ui.formula('Y = ' + fmt(result.A, 6) + result.b.map(function (b, j) {
            if (Math.abs(b) < 1e-12) return '';
            return (b >= 0 ? ' + ' : ' − ') + fmt(Math.abs(b), 6) + '·' + labels[j];
          }).join('') + '   (Y = ln|y − C₀|)', true)));

        // Đánh giá dự báo tại điểm x* (nếu có)
        const pred = buildPredictionCard(result, xStar);
        if (pred) {
          out.appendChild(pred.node);
        }

        const extLabels = result.extended.map(function (b) { return t(b.labelKey); });

        // Bảng tính tổng tích lũy của mô hình tuyến tính hóa Y = ln|y - C0|
        const Y = ys.map(function (y) { return Math.log(Math.abs(y - c0)); });
        out.appendChild(buildSumTableSection(result, xs, Y, extLabels, 'Y'));

        out.appendChild(buildMatricesSection(result, xs, extLabels, 'Y'));

        const fitInfo = chartFit(result, result.nodes, t('common.chart'), xStar, pred ? pred.yStar : null);
        out.appendChild(ui.section(t('common.chart'), fitInfo.node));
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
    return tabRoot;
  }

  function renderPower(tabRoot) {
    const xy = ui.xyControl({
      xs: POWER_PRESET.xs,
      ys: POWER_PRESET.ys,
      xLabel: t('lsq.xs'),
      yLabel: t('lsq.ys'),
      xHint: 'Các giá trị x > 0',
      yHint: 'Các giá trị y > 0'
    });
    const xStarIn = ui.input({ placeholder: 'Ví dụ: 3.5' });
    const out = ui.el('div');

    const form = ui.card(t('common.input'),
      xy.node,
      ui.grid(
        ui.field(t('lsq.xStar'), xStarIn, t('lsq.xStarHint'))
      ),
      ui.actions(
        ui.button('▶ ' + t('common.run'), run, { variant: 'primary' }),
        (function(){
          const imp = ui.importDataButton({ onParsed: function(res){ xy.setValues(res.xs||[], res.ys||[]); run(); } });
          imp.node.style.marginLeft = 'auto';
          return imp.node;
        })()
      )
    );

    tabRoot.appendChild(form);
    tabRoot.appendChild(out);

    function run() {
      out.innerHTML = '';
      try {
        const { xs, ys } = xy.getValues();
        const result = PPS.leastSquares.fitPower(xs, ys);

        let xStar = null;
        if (xStarIn.value.trim() !== '') {
          xStar = Number(xStarIn.value.trim().replace(',', '.'));
          if (!isFinite(xStar)) xStar = null;
        }

        out.appendChild(ui.kv([
          { key: 'Hệ số a', value: fmt(result.a, 8), className: 'accent' },
          { key: 'Số mũ b', value: fmt(result.b, 8), className: 'hl' },
          { key: t('common.sse'), value: fmt(result.sse, 6), className: 'green' },
          { key: 'RMSE', value: fmt(result.rmse, 6) },
          { key: 'R²', value: fmt(result.r2, 6), className: result.r2 > 0.99 ? 'green' : '' },
          { key: t('lsq.r2adj'), value: fmt(result.r2adj, 6) }
        ]));

        const modelStr = 'y = ' + fmt(result.a, 6) + ' · x' + PPS.format.superscript(Number(result.b.toFixed(2)));
        out.appendChild(ui.section(t('lsq.model'), ui.formula(modelStr, true)));

        // Đánh giá dự báo tại điểm x* (nếu có)
        const pred = buildPredictionCard(result, xStar);
        if (pred) {
          out.appendChild(pred.node);
        }

        out.appendChild(ui.note('Tuyến tính hóa qua log-log: ln(y) = ln(a) + b·ln(x). Yêu cầu dữ liệu x > 0 và y > 0.'));

        out.appendChild(buildMatricesSection(result, xs, ['1', 'ln(x)'], 'ln(y)'));

        const fitInfo = chartFit(result, result.nodes, t('common.chart'), xStar, pred ? pred.yStar : null);
        out.appendChild(ui.section(t('common.chart'), fitInfo.node));
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
    return tabRoot;
  }

  function render(root) {
    const tabControl = ui.tabs([
      { key: 'linear', label: t('lsq.tabLinear'), render: function () { const node = ui.el('div'); return renderLinear(node); } },
      { key: 'log', label: t('lsq.tabLog'), render: function () { const node = ui.el('div'); return renderLog(node); } },
      { key: 'power', label: t('lsq.tabPower'), render: function () { const node = ui.el('div'); return renderPower(node); } }
    ], 'linear');
    root.appendChild(tabControl.node);
  }

  PPS.views['least-squares'] = { render: render };
})();
