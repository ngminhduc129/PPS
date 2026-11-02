/* ============================================================
   View: Newton forward/backward interpolation (equal spacing)
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;
  const polyToString = PPS.format.polyToString;
  const poly = PPS.poly;

  const PRESETS = [
    { xs: [0, 1, 2, 3, 4], ys: [1, 2, 5, 10, 17], xbar: 1.5 },
    { xs: [0, 0.5, 1, 1.5, 2, 2.5, 3], ys: [0, 0.4794, 0.8415, 0.9975, 0.9093, 0.5985, 0.1411], xbar: 1.2 },
    { xs: [1, 2, 3, 4, 5, 6], ys: [2, 3, 10, 29, 66, 127], xbar: 4.7 }
  ];

  function render(root) {
    let presetIndex = 0;
    const xy = ui.xyControl({
      xs: PRESETS[0].xs,
      ys: PRESETS[0].ys,
      xLabel: t('newtonEq.xs'),
      yLabel: t('newtonEq.ys'),
      xHint: t('lagrange.xsHint'),
      yHint: t('lagrange.ysHint')
    });
    const xbarIn = ui.input({ value: String(PRESETS[0].xbar) });
    const modeSel = ui.select([
      { value: 'auto', label: t('common.auto') },
      { value: 'forward', label: t('common.forward') },
      { value: 'backward', label: t('common.backward') }
    ], 'auto');
    const out = ui.el('div');

    const form = ui.card(t('common.input'),
      xy.node,
      ui.grid(
        ui.field(t('newtonEq.xbar'), xbarIn),
        ui.field(t('newtonEq.mode'), modeSel, t('newtonEq.modeHint'))
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
        const xbar = ui.parseNumber(xbarIn.value, 'xbar');
        const result = PPS.newton.equalSpacing(xs, ys, xbar, modeSel.value);
        const n = xs.length;

        const modeLabel = result.mode === 'forward' ? t('common.forward') : t('common.backward');
        out.appendChild(ui.kv([
          { key: t('newtonEq.h'), value: fmt(result.h, 8) },
          { key: 't', value: fmt(result.t, 8), className: 'accent' },
          { key: t('common.method'), value: modeLabel, small: true },
          { key: 'P(x̄)', value: fmt(result.value, 8), className: 'green' },
          { key: "P′(x̄)", value: fmt(result.derivative, 8), className: 'orange' }
        ]));

        out.appendChild(ui.note(
          result.mode === 'forward' ? t('newtonEq.forwardNote') : t('newtonEq.backwardNote'),
          'good'
        ));

        out.appendChild(ui.section(t('common.formula'),
          ui.formula(result.mode === 'forward' ? t('newtonEq.formulaForward') : t('newtonEq.formulaBackward'), true)));

        // Difference table
        const headers = [t('newtonEq.diff')];
        for (let i = 0; i < n; i++) headers.push('i = ' + i);
        const rows = [];
        for (let j = 0; j < n; j++) {
          const row = [{ text: 'Δ' + PPS.format.superscript(j), class: 'rowhead' }];
          for (let i = 0; i < n; i++) {
            const value = result.table[j][i];
            if (value == null) {
              row.push('');
            } else {
              const used = result.mode === 'forward'
                ? (i === 0)
                : (i === n - 1 - j);
              row.push({ text: fmt(value, 7), class: used ? 'hl' : '' });
            }
          }
          rows.push(row);
        }
        const diffTable = ui.table(headers, rows);
        const diffSection = ui.section(t('newtonEq.table'), diffTable);
        out.appendChild(diffSection);

        // Trích xuất Vector hệ số sai phân & tỷ sai phân
        const fwdVectorStr = 'C_tiến = [' + result.cForward.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const bwdVectorStr = 'C_lùi = [' + result.cBackward.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const normVectorStr = 'D = [' + result.cNormalized.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const divNodes = result.mode === 'forward' ? xs : xs.slice().reverse();
        const divNodesStr = '[' + divNodes.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const divVectorStr = '[' + result.cDivided.map(function (v) { return fmt(v, 6); }).join(', ') + ']';

        out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), [
          ui.el('p', { class: 'math-label', text: 'Vector hệ số sai phân tiến (đường chéo trên):' }),
          ui.formula(fwdVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector hệ số sai phân lùi (đường chéo dưới / hàng cuối):', style: 'margin-top: 8px;' }),
          ui.formula(bwdVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector hệ số chuẩn hóa D_i = C_i / i!:', style: 'margin-top: 8px;' }),
          ui.formula(normVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector hệ số tỷ sai phân Newton [f[x₀], f[x₀,x₁], …]:', style: 'margin-top: 8px;' }),
          ui.el('p', { class: 'math-label dim', text: 'Thứ tự mốc: ' + divNodesStr }),
          ui.formula(divVectorStr, true)
        ]));

        // Bảng tích lũy thừa B_i(t) (Product Table)
        if (result.productTable && result.productTable.length) {
          const prodHeaders = ['Nhân tử'];
          for (let d = n - 1; d >= 0; d--) {
            prodHeaders.push(d === 0 ? 't⁰' : (d === 1 ? 't' : 't' + PPS.format.superscript(d)));
          }
          const prodRows = result.productTable.map(function (row) {
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
            'Bậc i',
            'Sai phân C_i',
            'Hệ số D_i = C_i/i!',
            'Biểu thức tích B_i(t)',
            'Đa thức số hạng N_i(t)'
          ];
          const stepRows = result.productTable.map(function (row, idx) {
            return [
              { text: String(idx), class: 'rowhead' },
              { text: fmt(row.diff, 6) },
              { text: fmt(row.coef, 6), class: 'hl' },
              { text: row.factorStr || '1', class: 'left' },
              { text: polyToString(row.termPoly, 't'), class: 'left' }
            ];
          });
          const stepTableEl = ui.table(stepHeaders, stepRows);

          out.appendChild(ui.section(t('newtonDiv.productTableTitle'), [
            prodTableEl,
            ui.el('p', { class: 'math-label', text: 'Bảng các bước dựng đa thức số hạng N_i(t) = D_i · B_i(t):', style: 'margin-top: 14px;' }),
            stepTableEl
          ]));
        }

        // Đa thức chính tắc theo biến t và biến thực x
        const polyNodes = [
          ui.el('p', { class: 'math-label', text: 'Theo biến chuẩn hóa t:' }),
          ui.formula('P(t) = ' + polyToString(result.polyT, 't'), true)
        ];
        if (result.polyX) {
          polyNodes.push(
            ui.el('p', { class: 'math-label', text: 'Quy đổi theo biến thực x:', style: 'margin-top: 8px;' }),
            ui.formula('P(x) = ' + polyToString(result.polyX, 'x'), true)
          );
        }
        out.appendChild(ui.section(t('common.canonicalPoly'), polyNodes));

        // Coefficients in t
        const termRows = result.terms.map(function (term) {
          return [
            { text: String(term.j), class: 'rowhead' },
            { text: term.sourceLabel, class: 'left' },
            { text: fmt(term.diff, 7) },
            { text: fmt(term.coef, 7), class: 'hl' },
            { text: term.factor, class: 'left' },
            { text: polyToString(term.poly, 't'), class: 'left dim' }
          ];
        });
        out.appendChild(ui.section(t('newtonEq.coefficients'),
          ui.table([
            t('newtonEq.term'),
            t('newtonEq.diff'),
            t('common.value'),
            t('common.coefficient'),
            t('central.factor'),
            'ωⱼ(t)'
          ], termRows)));

        out.appendChild(ui.section(t('newtonEq.valueDeriv'),
          ui.formula('f(x̄) ≈ ' + fmt(result.value, 8) + '   ·   f′(x̄) ≈ ' + fmt(result.derivative, 8), true)));
        const chart = ui.chart(t('common.chart'), 370);
        out.appendChild(chart.node);
        const xmin = xs[0];
        const xmax = xs[n - 1];
        const span = Math.max(0.4, (xmax - xmin) * 0.08);
        const evalAt = function (x) {
          return poly.evaluate(result.polyT, (x - result.anchorX) / result.h);
        };
        const curve = PPS.chart.sample(evalAt, xmin - span, xmax + span, 340);
        const handle = chart.draw({
          series: [{ points: curve, color: 'var(--accent)', width: 2.4, label: 'Pₙ(x)' }],
          points: xs.map(function (x, i) {
            return { x: x, y: ys[i], color: 'var(--text-dim)', r: 4, label: 'x = ' + fmt(x, 5) };
          }),
          highlightPoints: [
            { x: result.anchorX, y: ys[result.anchorIndex], color: 'var(--green)', label: result.mode === 'forward' ? 'x₀' : 'xₙ' },
            { x: xbar, y: result.value, color: 'var(--orange)', label: 'P(x̄) = ' + fmt(result.value, 5) }
          ],
          vline: { x: xbar, label: 'x̄ = ' + fmt(xbar, 6), color: 'var(--orange)' },
          xLabel: 'x',
          yLabel: 'y'
        });

        // Partial sums of the Newton terms
        const partial = [];
        let acc = [0];
        result.terms.forEach(function (term) {
          acc = PPS.poly.add(acc, PPS.poly.scale(term.poly, term.coef));
          partial.push(acc);
        });

        function sourceCell(j) {
          const i = result.mode === 'forward' ? 0 : n - 1 - j;
          const row = diffTable.querySelectorAll('tbody tr')[j];
          return row ? row.children[i + 1] : null;
        }

        const anim = ui.animBar({
          count: result.terms.length,
          delay: 680,
          onStep: function (index) {
            ui.clearAnim(diffTable);
            ui.highlightRow(diffTable, index);
            for (let j = 0; j <= index; j++) {
              const cell = sourceCell(j);
              if (cell) cell.classList.add('anim-cell');
            }
            const points = PPS.chart.sample(function (x) {
              return PPS.poly.evaluate(partial[index], (x - result.anchorX) / result.h);
            }, xmin - span, xmax + span, 340);
            handle.update({
              series: [{ points: points, color: 'var(--accent)', width: 2.4, label: 'Pₙ(x)' }],
              marker: {
                x: xbar,
                y: PPS.poly.evaluate(partial[index], result.t),
                color: 'var(--orange)',
                label: 'P' + PPS.central.sub(index) + '(x̄)'
              }
            });
          },
          onReset: function () {
            ui.clearAnim(diffTable);
            handle.update({
              series: [{ points: curve, color: 'var(--accent)', width: 2.4, label: 'Pₙ(x)' }],
              marker: null
            });
          }
        });
        diffSection.appendChild(anim.node);

        // small verification table of used points
        const checkRows = result.usedIndices.slice(0, Math.min(6, n)).map(function (i) {
          return [
            { text: String(i), class: 'rowhead' },
            { text: fmt(xs[i], 6) },
            { text: fmt(ys[i], 6) },
            { text: fmt(PPS.poly.evaluate(result.polyT, (xs[i] - result.anchorX) / result.h), 7), class: 'hl-green' }
          ];
        });
        out.appendChild(ui.section(t('common.verification'),
          ui.table([t('common.node'), 'xᵢ', 'yᵢ', 'P(xᵢ)'], checkRows)));
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
  }

  PPS.views['newton-eq'] = { render: render };
})();
