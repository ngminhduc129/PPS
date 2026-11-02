/* ============================================================
   View: Chebyshev optimal nodes & Runge phenomenon demo
   Based on Section 2 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;

  const PRESETS = [
    { a: -1, b: 1, n: 8, fn: '1 / (1 + 25 * x * x)', M: '' },
    { a: -1, b: 1, n: 10, fn: '1 / (1 + 25 * x * x)', M: '' },
    { a: 0, b: 3, n: 5, fn: 'Math.exp(x)', M: '20.1' },
    { a: -2, b: 2, n: 6, fn: 'Math.cos(x)', M: '1' }
  ];

  function render(root) {
    let presetIdx = 0;

    const aIn = ui.input({ value: String(PRESETS[0].a) });
    const bIn = ui.input({ value: String(PRESETS[0].b) });
    const nIn = ui.input({ value: String(PRESETS[0].n) });
    const mIn = ui.input({ value: String(PRESETS[0].M) });
    const fnIn = ui.input({ value: PRESETS[0].fn });
    const out = ui.el('div');

    const runBtn = ui.button('▶ ' + t('common.run'), function () { run(); }, { variant: 'primary' });
    const exampleBtn = ui.button('✨ ' + t('common.example'), function () {
      presetIdx = (presetIdx + 1) % PRESETS.length;
      const p = PRESETS[presetIdx];
      aIn.value = String(p.a);
      bIn.value = String(p.b);
      nIn.value = String(p.n);
      mIn.value = String(p.M);
      fnIn.value = p.fn;
      run();
    }, { variant: 'ghost' });

    const form = ui.card(t('common.input'),
      ui.grid(
        ui.field(t('cheb.a'), aIn, t('cheb.aHint')),
        ui.field(t('cheb.b'), bIn, t('cheb.bHint')),
        ui.field(t('cheb.n'), nIn, t('cheb.nHint')),
        ui.field(t('cheb.mBound'), mIn, t('cheb.mBoundHint')),
        ui.field(t('cheb.fn'), fnIn, t('cheb.fnHint'), { span2: true })
      ),
      ui.actions(runBtn, exampleBtn)
    );

    root.appendChild(form);
    root.appendChild(out);

    function run() {
      out.innerHTML = '';
      try {
        const a = ui.parseNumber(aIn.value);
        const b = ui.parseNumber(bIn.value);
        const n = Math.floor(ui.parseNumber(nIn.value));
        const M = mIn.value.trim() ? ui.parseNumber(mIn.value) : null;
        const fnExpr = fnIn.value.trim();

        if (!(b > a)) throw new Error('invalidInterval');
        if (n < 1 || n > 25) throw new Error('invalidDegree');

        const sim = PPS.chebyshev.compareRunge(fnExpr, a, b, n, 250);
        const boundVal = M ? PPS.chebyshev.errorBound(a, b, n, M) : null;

        // KPI summary
        const kpis = [
          { key: t('cheb.totalNodes'), value: String(n + 1), className: 'accent' },
          { key: t('cheb.maxErrCheb'), value: fmt(sim.maxErrCheb, 6), className: 'green' },
          { key: t('cheb.maxErrEqui'), value: fmt(sim.maxErrEqui, 6), className: sim.maxErrEqui > sim.maxErrCheb * 2 ? 'red' : 'orange' }
        ];
        if (boundVal != null) {
          kpis.push({ key: t('cheb.boundTheory'), value: fmt(boundVal, 6), className: 'cyan' });
        }
        out.appendChild(ui.kv(kpis));

        // Note about Runge phenomenon
        if (sim.maxErrEqui > sim.maxErrCheb * 1.5) {
          out.appendChild(ui.note(t('cheb.rungeObserved', {
            ratio: (sim.maxErrEqui / Math.max(1e-12, sim.maxErrCheb)).toFixed(1)
          }), 'warn'));
        }

        // Table of Chebyshev nodes
        const tableRows = sim.chebNodes.map(function (node, idx) {
          return [
            { text: 'x' + PPS.format.subscript(idx), class: 'rowhead' },
            { text: String(node.k) },
            { text: fmt(node.theta, 6) },
            { text: fmt(node.t, 8) },
            { text: fmt(node.x, 8), class: 'hl' },
            { text: fmt(PPS.chebyshev.equispacedNodes(a, b, n)[idx], 8) }
          ];
        });

        const table = ui.table([
          t('common.node'),
          'k',
          'θ_k (rad)',
          't_k = cos(θ_k)',
          t('cheb.colChebX'),
          t('cheb.colEquiX')
        ], tableRows);

        out.appendChild(ui.section(t('cheb.nodesTableTitle'), table, (n + 1) + ' ' + t('common.node').toLowerCase()));

        // Trích xuất các Vector kết quả Chebyshev
        const equiNodes = PPS.chebyshev.equispacedNodes(a, b, n);
        const chebXStr = 'X_cheb = [' + sim.chebX.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const equiXStr = 'X_equi = [' + equiNodes.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const thetaStr = 'θ = [' + sim.chebNodes.map(function (nd) { return fmt(nd.theta, 6); }).join(', ') + ']';
        const tStr = 't_k = cos(θ_k) = [' + sim.chebNodes.map(function (nd) { return fmt(nd.t, 6); }).join(', ') + ']';

        out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), [
          ui.el('p', { class: 'math-label', text: 'Vector mốc tối ưu Chebyshev X_cheb trên đoạn [' + a + ', ' + b + ']:' }),
          ui.formula(chebXStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector mốc cách đều đối chứng X_equi:', style: 'margin-top: 8px;' }),
          ui.formula(equiXStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector góc lượng giác θ_k = (2k + 1)π / (2n + 2):', style: 'margin-top: 8px;' }),
          ui.formula(thetaStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector mốc chuẩn hóa trên [-1, 1]: t_k = cos(θ_k):', style: 'margin-top: 8px;' }),
          ui.formula(tStr, true)
        ]));

        // Chart comparison
        const chart = ui.chart(t('common.chart') + ' — ' + t('cheb.chartTitle'), 400);

        const truePoints = sim.grid.map(function (p) { return { x: p.x, y: p.yTrue }; });
        const chebPoints = sim.grid.map(function (p) { return { x: p.x, y: p.yCheb }; });
        const equiPoints = sim.grid.map(function (p) { return { x: p.x, y: p.yEqui }; });

        const chebNodeMarkers = sim.chebX.map(function (x, i) {
          return { x: x, y: sim.chebY[i], color: 'var(--green)', r: 4.5, label: 'Cheb ' + fmt(x, 3) };
        });

        chart.draw({
          series: [
            { points: truePoints, color: 'var(--accent)', width: 2.5, label: 'f(x) chuẩn' },
            { points: chebPoints, color: 'var(--green)', width: 2.2, label: 'Nội suy Chebyshev (P_Cheb)' },
            { points: equiPoints, color: 'var(--red)', width: 1.8, dash: '6 4', label: 'Nội suy Mốc cách đều (P_Equi)' }
          ],
          points: chebNodeMarkers,
          xLabel: 'x',
          yLabel: 'y'
        });

        out.appendChild(chart.node);
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
  }

  PPS.views.chebyshev = { render: render };
})();
