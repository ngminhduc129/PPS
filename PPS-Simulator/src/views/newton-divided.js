/* ============================================================
   View: Newton interpolation with divided differences
   Updated based on Section 4 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;
  const polyToString = PPS.format.polyToString;

  const PRESETS = [
    // Ví dụ chuẩn 7 mốc từ Mục 4.3 - Algorithm.md
    {
      xs: [1.57, 1.59, 1.62, 1.64, 1.67, 1.71, 1.74],
      ys: [2.4625, 2.5857, 2.772, 2.897, 3.086, 3.3404, 3.5328],
      xbar: 1.64
    },
    { xs: [-1, 0, 1, 2], ys: [2, -1, 3, 1], xbar: 0.5 },
    { xs: [0, 2, 3, 5], ys: [1, 3, 2, 6], xbar: 1 }
  ];

  function render(root) {
    let presetIndex = 0;
    const xy = ui.xyControl({
      xs: PRESETS[0].xs,
      ys: PRESETS[0].ys,
      xLabel: t('newtonDiv.xs'),
      yLabel: t('newtonDiv.ys'),
      xHint: t('newtonDiv.xsHint') || t('lagrange.xsHint'),
      yHint: t('newtonDiv.ysHint') || t('lagrange.ysHint')
    });
    const xbarIn = ui.input({ value: String(PRESETS[0].xbar) });
    const out = ui.el('div');

    const form = ui.card(t('common.input'),
      xy.node,
      ui.grid(
        ui.field(t('newtonDiv.xEval'), xbarIn, t('newtonDiv.xEvalHint') || t('lagrange.xEvalHint'))
      ),
      ui.actions(
        ui.button('▶ ' + t('common.run'), run, { variant: 'primary' }),
        ui.button('✨ ' + t('common.example'), function () {
          presetIndex = (presetIndex + 1) % PRESETS.length;
          const p = PRESETS[presetIndex];
          xy.setValues(p.xs, p.ys);
          xbarIn.value = String(p.xbar);
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
        const result = PPS.newton.interpolate(xs, ys);
        const xbarText = xbarIn.value.trim();
        const xbar = xbarText ? ui.parseNumber(xbarText, 'xbar') : null;

        const summary = [
          { key: t('common.degreeN'), value: String(result.degree) },
          { key: t('common.nodes'), value: String(xs.length) }
        ];

        let derivVals = null;
        if (xbar != null) {
          const pVal = result.evaluate(xbar);
          derivVals = result.evalDerivatives(xbar, 2);
          summary.push({ key: 'P(x*)', value: fmt(pVal, 8), className: 'accent' });
          if (derivVals.length > 1) {
            summary.push({ key: "P'(x*)", value: fmt(derivVals[1], 8), className: 'green' });
          }
          if (derivVals.length > 2) {
            summary.push({ key: "P''(x*)", value: fmt(derivVals[2], 8), className: 'orange' });
          }
        }
        out.appendChild(ui.kv(summary));

        // Trích xuất Vector hệ số tiến & lùi (Mục 4.3 A - Algorithm.md)
        const forwardVectorStr = 'C_tiến = [' + result.cForward.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const backwardVectorStr = 'C_lùi = [' + result.cBackward.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), [
          ui.el('p', { class: 'math-label', text: 'Vector hệ số nội suy tiến (đường chéo trên):' }),
          ui.formula(forwardVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector hệ số nội suy lùi (đường chéo dưới / hàng cuối):', style: 'margin-top: 8px;' }),
          ui.formula(backwardVectorStr, true)
        ]));

        // Đa thức chính tắc
        const polyFormula = polyToString(result.coef, 'x');
        out.appendChild(ui.section(t('common.canonicalPoly'), ui.formula('P(x) = ' + polyFormula, true)));

        // Bảng tỷ sai phân 2D
        const n = xs.length;
        const headers = [t('common.node'), 'x_i', 'y_i = f[x_i]'];
        for (let j = 1; j < n; j++) headers.push('f[' + j + ']');
        const rows = [];
        for (let i = 0; i < n; i++) {
          const cells = [
            { text: String(i), class: 'rowhead' },
            { text: fmt(xs[i], 6), class: 'hl' },
            { text: fmt(ys[i], 6) }
          ];
          for (let j = 1; j < n; j++) {
            const val = result.table[i][j];
            if (val === null || val === undefined) {
              cells.push({ text: '—', class: 'muted' });
            } else {
              const isFwd = i === j;
              const isBwd = i === n - 1;
              cells.push({
                text: fmt(val, 6),
                class: isFwd ? 'hl' : (isBwd ? 'accent' : '')
              });
            }
          }
          rows.push(cells);
        }
        const table = ui.table(headers, rows);

        // Bảng tích lũy thừa (Product Table - Mục 4.3 B)
        const prodHeaders = [t('common.node') + ' x_k'];
        for (let d = n - 1; d >= 0; d--) {
          prodHeaders.push(d === 0 ? 'x⁰' : (d === 1 ? 'x' : 'x' + PPS.format.superscript(d)));
        }
        const prodRows = result.productTable.map(function (row) {
          const cells = [
            { text: row.node === null ? '—' : fmt(row.node, 4), class: 'rowhead' }
          ];
          row.polyDesc.forEach(function (coefVal) {
            cells.push({ text: fmt(coefVal, 6) });
          });
          return cells;
        });
        const prodTableEl = ui.table(prodHeaders, prodRows);

        // Animation & Highlights
        const anim = ui.animBar({
          count: result.terms.length,
          delay: 720,
          onStep: function (index) {
            ui.highlightRow(table, index);
            const subTerms = result.terms.slice(0, index + 1);
            let partial = [subTerms[0].coef];
            for (let k = 1; k < subTerms.length; k++) {
              partial = PPS.poly.add(partial, PPS.poly.scale(subTerms[k].omega, subTerms[k].coef));
            }
            const curve = PPS.poly.sample(partial, minX - margin, maxX + margin, 250);
            handle.update({
              series: [{
                points: curve,
                color: 'var(--accent)',
                width: 2.2,
                label: 'P_' + index + '(x)'
              }],
              marker: xbar != null ? { x: xbar, y: PPS.poly.evaluate(partial, xbar), color: 'var(--green)' } : null
            });
          },
          onReset: function () {
            ui.clearAnim(table);
            handle.update({ series: baseSeries, marker: null });
          }
        });

        out.appendChild(ui.section(t('newtonDiv.dividedTable'), [
          table,
          anim.node
        ], result.terms.length + ' ' + t('common.terms').toLowerCase()));

        out.appendChild(ui.section(t('newtonDiv.productTableTitle'), prodTableEl));

        // Chart
        const minX = Math.min.apply(null, xs);
        const maxX = Math.max.apply(null, xs);
        const margin = Math.max(0.8, (maxX - minX) * 0.18);
        const fullCurve = PPS.poly.sample(result.coef, minX - margin, maxX + margin, 320);

        const nodeMarkers = xs.map(function (x, i) {
          return { x: x, y: ys[i], color: 'var(--accent)', r: 5, label: '(' + fmt(x, 2) + ', ' + fmt(ys[i], 2) + ')' };
        });

        const baseSeries = [{
          points: fullCurve,
          color: 'var(--accent)',
          width: 2.4,
          label: 'P(x)'
        }];

        const chart = ui.chart(t('common.chart'), 380);
        const handle = chart.draw({
          series: baseSeries,
          points: nodeMarkers,
          vline: xbar != null ? { x: xbar, label: 'x* = ' + fmt(xbar, 4), color: 'var(--green)' } : null,
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

  PPS.views['newton-div'] = { render: render };
})();
