/* ============================================================
   View: central interpolation — Gauss I/II, Stirling, Bessel
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;
  const polyToString = PPS.format.polyToString;

  const PRESETS = [
    { xs: [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4], ys: [1, 1.3125, 1.5, 1.5625, 1.5, 1.3125, 1, 0.5625, 0], xbar: 1.9 },
    { xs: [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4], ys: [1, 0.875, 1, 1.375, 2, 2.875, 4, 5.375, 7], xbar: 2.2 },
    { xs: [1, 1.2, 1.4, 1.6, 1.8, 2, 2.2], ys: [2, 3.78, 6.96, 12.1, 20.4, 33, 51.5], xbar: 1.55 }
  ];

  function formulaLabel(key) {
    const map = {
      gauss1: t('central.gauss1'),
      gauss2: t('central.gauss2'),
      stirling: t('central.stirling'),
      bessel: t('central.bessel')
    };
    return map[key] || key;
  }

  function render(root) {
    let presetIndex = 0;
    const xy = ui.xyControl({
      xs: PRESETS[0].xs,
      ys: PRESETS[0].ys,
      xLabel: t('central.xs'),
      yLabel: t('central.ys'),
      xHint: t('lagrange.xsHint'),
      yHint: t('lagrange.xsHint')
    });
    const xbarIn = ui.input({ value: String(PRESETS[0].xbar) });
    const formulaSel = ui.select([
      { value: 'auto', label: t('central.auto') },
      { value: 'gauss1', label: t('central.gauss1') },
      { value: 'gauss2', label: t('central.gauss2') },
      { value: 'stirling', label: t('central.stirling') },
      { value: 'bessel', label: t('central.bessel') }
    ], 'auto');
    const nIn = ui.input({ value: 'auto' });
    const kIn = ui.input({ value: '1', placeholder: '1, 2, ...' });
    const out = ui.el('div');

    const form = ui.card(t('common.input'),
      xy.node,
      ui.grid(
        ui.field(t('central.xbar'), xbarIn),
        ui.field(t('central.formula'), formulaSel),
        ui.field('N', nIn, t('central.NHint')),
        ui.field(t('central.orderK'), kIn, t('central.orderKHint'))
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
        let N = 'auto';
        const nText = nIn.value.trim();
        if (nText && nText.toLowerCase() !== 'auto') {
          const parsed = ui.parseNumber(nText);
          N = Math.max(0, Math.round(parsed));
        }
        const result = PPS.central.interpolate(xs, ys, xbar, {
          formula: formulaSel.value,
          N: N
        });

        const varName = result.variable;

        let kOrder = 1;
        const kText = kIn.value.trim();
        if (kText) {
          const parsedK = ui.parseNumber(kText);
          kOrder = Math.max(1, Math.round(parsedK));
        }

        let selectedDeriv = null;
        if (result.derivatives && result.derivatives.length) {
          if (kOrder <= result.derivatives.length) {
            selectedDeriv = result.derivatives[kOrder - 1];
          } else {
            selectedDeriv = {
              k: kOrder,
              symbol: "f^(" + kOrder + ")(x̄)",
              symbolT: "P^(" + kOrder + ")(" + varName + "*)",
              valT: 0,
              valX: 0,
              polyT: [0],
              polyX: [0],
              formulaStr: "P^(" + kOrder + ")(" + varName + "*) / h" + (kOrder > 1 ? PPS.format.superscript(kOrder) : "")
            };
          }
        }

        const summary = [
          { key: t('central.formulaLabel'), value: formulaLabel(result.formula), small: true, className: 'accent' },
          { key: 'N', value: String(result.N) },
          { key: t('central.center'), value: fmt(result.centerX, 6) },
          { key: 'f(x̄)', value: fmt(result.value, 8), className: 'green' }
        ];
        if (selectedDeriv) {
          summary.push({
            key: selectedDeriv.symbol,
            value: fmt(selectedDeriv.valX, 8),
            className: 'orange'
          });
        } else {
          summary.push({ key: 'f′(x̄)', value: fmt(result.derivative, 8), className: 'orange' });
        }
        out.appendChild(ui.kv(summary));

        let intervalInfo = null;
        for (let i = 0; i < xs.length - 1; i++) {
          if (xbar >= xs[i] - 1e-10 && xbar <= xs[i + 1] + 1e-10) {
            const k1 = i;
            const k2 = i + 1;
            const sub1 = PPS.central ? PPS.central.sub(k1) : String(k1);
            const sub2 = PPS.central ? PPS.central.sub(k2) : String(k2);
            const startX = xs[result.blockStart];
            const endX = xs[result.blockStart + result.blockSize - 1];
            intervalInfo = {
              sub1: sub1,
              sub2: sub2,
              x1: xs[k1],
              x2: xs[k2],
              startX: startX,
              endX: endX
            };
            break;
          }
        }
        if (intervalInfo) {
          out.appendChild(ui.note(
            'Điểm x̄ = ' + fmt(xbar, 4) + ' nằm trong khoảng trung tâm [x' + intervalInfo.sub1 + ', x' + intervalInfo.sub2 + '] = [' + fmt(intervalInfo.x1, 4) + ', ' + fmt(intervalInfo.x2, 4) + '] (với bước h = ' + fmt(result.h, 4) + '). Mốc trung tâm x₀ = ' + fmt(result.centerX, 4) + ' · Mốc sử dụng chi tiết: [' + fmt(intervalInfo.startX, 4) + ', ' + fmt(intervalInfo.endX, 4) + '] (gồm ' + result.blockSize + ' mốc).',
            'good'
          ));
        }

        if (result.autoReason) {
          const reasonKey = result.autoReason === 'stirling' ? 'central.reasonStirling'
            : result.autoReason === 'bessel' ? 'central.reasonBessel' : 'central.reasonGauss';
          out.appendChild(ui.note(t('central.formulaPickAuto', {
            formula: formulaLabel(result.formula),
            reason: t(reasonKey)
          }), 'good'));
        }
        if (result.clamped) {
          out.appendChild(ui.note(t('central.clamped', { max: result.maxN, n: result.N }), 'warn'));
        }

        out.appendChild(ui.section(t('common.formula'), [
          ui.formula(polyToString(result.polyT, varName), true),
          ui.formula(
            (varName === 't'
              ? 't = (x̄ − x₀)/h = ' + fmt(result.t, 8)
              : 't = (x̄ − x₀)/h = ' + fmt(result.t, 8) + '   ·   u = t − 1/2 = ' + fmt(result.u, 8))
            + '   ·   f⁽ᵏ⁾(x̄) ≈ P⁽ᵏ⁾(' + varName + ')/hᵏ',
            true
          )
        ]));

        // terms
        const termRows = result.terms.map(function (term) {
          const sourceLabels = term.sources.map(function (s) { return s.label; }).join(' + ');
          return [
            { text: String(term.j), class: 'rowhead' },
            { text: sourceLabels, class: 'left' },
            { text: fmt(term.coef, 8) },
            { text: term.factorString, class: 'left' },
            { text: polyToString(term.poly, varName), class: 'left dim' }
          ];
        });
        out.appendChild(ui.section(t('central.termTable'),
          ui.table([
            t('central.order'),
            t('central.source'),
            t('common.coefficient'),
            t('central.factor'),
            'Pⱼ(' + varName + ')'
          ], termRows),
          result.terms.length + ' ' + t('common.steps').toLowerCase()));

        // central difference table
        const highlight = {};
        result.highlight.forEach(function (cell) {
          highlight[cell.j + ':' + cell.local] = true;
        });
        const headers = ['Δ' + PPS.format.superscript(0) + ' / ∇'];
        result.columns.forEach(function (col) {
          headers.push('y' + PPS.central.sub(col.rel));
        });
        const rows = [];
        const size = result.columns.length;
        for (let j = 0; j < size; j++) {
          const row = [{ text: 'Δ' + PPS.format.superscript(j), class: 'rowhead' }];
          for (let local = 0; local < size; local++) {
            const value = result.table[j][local];
            if (value == null) {
              row.push('');
            } else {
              row.push({
                text: fmt(value, 7),
                class: highlight[j + ':' + local] ? 'hl' : '',
                title: 'Δ' + PPS.format.superscript(j) + 'y' + PPS.central.sub(result.columns[local].rel)
              });
            }
          }
          rows.push(row);
        }
        const diffTable = ui.table(headers, rows);
        const diffSection = ui.section(t('central.table'), [
          diffTable,
          ui.note(t('central.derivativeNote'))
        ], result.blockSize + ' ' + t('common.nodes'));
        out.appendChild(diffSection);

        // Trích xuất Vector hệ số sai phân & hệ số nội suy
        if (result.cVector && result.dVector) {
          const zStr = 'Z = [' + result.nodesZ.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
          const cStr = 'C = [' + result.cVector.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
          const dStr = 'D = [' + result.dVector.map(function (v) { return fmt(v, 6); }).join(', ') + ']';

          out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), [
            ui.el('p', { class: 'math-label', text: 'Vector mốc dữ liệu trích xuất Z:' }),
            ui.formula(zStr, true),
            ui.el('p', { class: 'math-label', text: 'Vector hệ số sai phân trích xuất C (theo công thức ' + formulaLabel(result.formula) + '):', style: 'margin-top: 8px;' }),
            ui.formula(cStr, true),
            ui.el('p', { class: 'math-label', text: 'Vector hệ số chuẩn hóa đa thức D = [D₀, D₁, …]:', style: 'margin-top: 8px;' }),
            ui.formula(dStr, true)
          ]));
        }

        // Bảng tích lũy thừa (Product Table)
        if (result.productTable && result.productTable.length) {
          const m = result.productTable.length;
          const prodHeaders = ['Nhân tử'];
          for (let d = m - 1; d >= 0; d--) {
            prodHeaders.push(d === 0 ? varName + '⁰' : (d === 1 ? varName : varName + PPS.format.superscript(d)));
          }
          const prodRows = result.productTable.map(function (row) {
            const cells = [
              { text: row.factorStr || '1', class: 'rowhead left' }
            ];
            row.polyDesc.forEach(function (coefVal) {
              cells.push({ text: fmt(coefVal, 6) });
            });
            return cells;
          });
          const prodTableEl = ui.table(prodHeaders, prodRows);

          const stepHeaders = [
            'Bậc j',
            'Sai phân C_j',
            'Hệ số D_j',
            'Biểu thức tích B_j(' + varName + ')',
            'Đa thức số hạng N_j(' + varName + ')'
          ];
          const stepRows = result.productTable.map(function (row) {
            return [
              { text: String(row.j), class: 'rowhead' },
              { text: fmt(row.diff, 6) },
              { text: fmt(row.coef, 6), class: 'hl' },
              { text: row.factorStr || '1', class: 'left' },
              { text: polyToString(row.termPoly, varName), class: 'left' }
            ];
          });
          const stepTableEl = ui.table(stepHeaders, stepRows);

          out.appendChild(ui.section(t('newtonDiv.productTableTitle'), [
            prodTableEl,
            ui.el('p', { class: 'math-label', text: 'Bảng các bước dựng đa thức số hạng N_j(' + varName + ') = D_j · B_j(' + varName + '):', style: 'margin-top: 14px;' }),
            stepTableEl
          ]));
        }

        // Đa thức chính tắc theo biến chuẩn hóa và biến thực x
        const polyNodes = [
          ui.el('p', { class: 'math-label', text: 'Theo biến chuẩn hóa ' + varName + ':' }),
          ui.formula('P(' + varName + ') = ' + polyToString(result.polyT, varName), true)
        ];
        if (result.polyX) {
          polyNodes.push(
            ui.el('p', { class: 'math-label', text: 'Quy đổi theo biến thực x:', style: 'margin-top: 8px;' }),
            ui.formula('P(x) = ' + polyToString(result.polyX, 'x'), true)
          );
        }
        out.appendChild(ui.section(t('common.canonicalPoly'), polyNodes));

        // Bảng đạo hàm các cấp tại x̄
        if (result.derivatives && result.derivatives.length) {
          const derivRows = result.derivatives.map(function (d) {
            const isSelected = d.k === kOrder;
            return [
              { text: String(d.k), class: 'rowhead' + (isSelected ? ' hl' : '') },
              { text: d.symbol, class: 'left' + (isSelected ? ' hl' : '') },
              { text: fmt(d.valX, 8), class: (isSelected ? 'hl bold' : '') },
              { text: fmt(d.valT, 8) },
              { text: d.formulaStr, class: 'left dim' },
              { text: polyToString(d.polyT, varName), class: 'left dim' }
            ];
          });
          const derivTableEl = ui.table([
            t('central.orderK') || 'Cấp k',
            t('central.derivSymbol') || 'Ký hiệu',
            t('central.derivValX') || 'Giá trị f^(k)(x̄)',
            (t('central.derivValT') || 'Giá trị') + ' (P⁽ᵏ⁾(' + varName + '*))',
            t('central.derivFormula') || 'Công thức quy đổi',
            (t('central.derivPoly') || 'Đa thức đạo hàm') + ' P⁽ᵏ⁾(' + varName + ')'
          ], derivRows);

          const derivNotes = [
            derivTableEl,
            ui.note(
              'Quy đổi đạo hàm: Với biến chuẩn hóa ' + varName + ' (bước h = ' + fmt(result.h, 6) + '), ta có d' + varName + '/dx = 1/h. ' +
              'Do đó: dᵏf/dxᵏ(x̄) ≈ P⁽ᵏ⁾(' + varName + '*) / hᵏ' +
              (selectedDeriv ? ' ⇒ Đạo hàm cấp k = ' + kOrder + ': ' + selectedDeriv.symbol + ' ≈ ' + fmt(selectedDeriv.valX, 8) : '') + '.',
              'good'
            )
          ];
          out.appendChild(ui.section(t('central.derivTableTitle') || 'Bảng đạo hàm các cấp tại x̄ (Stirling & Bessel)', derivNotes));
        }

        // chart
        const xmin = xs[0];
        const xmax = xs[xs.length - 1];
        const span = Math.max(0.3, (xmax - xmin) * 0.04);
        const lo = xmin - span;
        const hi = xmax + span;
        const curve = PPS.chart.sample(result.evalAt, lo, hi, 360);
        const chart = ui.chart(t('common.chart'), 380);
        out.appendChild(chart.node);

        const baseSeries = [{ points: curve, color: 'var(--accent)', width: 2.4, label: 'P(' + varName + ')' }];
        const cfg = {
          series: baseSeries,
          points: xs.map(function (x, i) {
            return { x: x, y: ys[i], color: 'var(--text-dim)', r: 4, label: 'x = ' + fmt(x, 5) };
          }),
          highlightPoints: result.usedIndices.map(function (i) {
            return { x: xs[i], y: ys[i], color: 'var(--green)', r: 4 };
          }).concat([{ x: xbar, y: result.value, color: 'var(--orange)', label: 'f(x̄) = ' + fmt(result.value, 5) }]),
          vline: { x: xbar, label: 'x̄ = ' + fmt(xbar, 6), color: 'var(--orange)' },
          xLabel: 'x',
          yLabel: 'y'
        };

        // Partial sums t ↦ c₀ + c₁ω₁(t) + … + cⱼωⱼ(t)
        function varOf(x) {
          return varName === 'u' ? (x - result.centerX) / result.h - 0.5 : (x - result.centerX) / result.h;
        }
        function partialSampler(coefs) {
          return function (x) { return PPS.poly.evaluate(coefs, varOf(x)); };
        }
        const partials = [];
        let accPoly = [0];
        result.terms.forEach(function (term) {
          accPoly = PPS.poly.add(accPoly, term.poly);
          partials.push(accPoly);
        });
        const domain = [];
        partials.forEach(function (coefs, i) {
          if (i === partials.length - 1) return; // full curve already included
          domain.push.apply(domain, PPS.chart.sample(partialSampler(coefs), lo, hi, 90));
        });
        cfg.domain = domain;

        const handle = chart.draw(cfg);

        function highlightTermSources(step) {
          ui.clearAnim(diffTable);
          const trs = diffTable.querySelectorAll('tbody tr');
          for (let k = 0; k <= step; k++) {
            const term = result.terms[k];
            if (!term) continue;
            term.sources.forEach(function (s) {
              const row = trs[s.j];
              if (!row) return;
              const cell = row.children[s.local + 1];
              if (cell) cell.classList.add('anim-cell');
            });
            if (trs[term.j]) trs[term.j].classList.add('anim-flash');
          }
        }

        const steps = result.terms.length;
        const anim = ui.animBar({
          count: steps,
          delay: 700,
          onStep: function (index) {
            highlightTermSources(index);
            const points = PPS.chart.sample(partialSampler(partials[index]), lo, hi, 360);
            handle.update({
              series: [{ points: points, color: 'var(--accent)', width: 2.4, label: 'P(' + varName + ')' }],
              marker: {
                x: xbar,
                y: PPS.poly.evaluate(partials[index], varOf(xbar)),
                color: 'var(--orange)',
                label: 'P' + PPS.central.sub(index) + '(x̄)'
              }
            });
          },
          onReset: function () {
            ui.clearAnim(diffTable);
            handle.update({ series: baseSeries, marker: null });
          }
        });
        diffSection.appendChild(anim.node);

        const usedTable = result.usedIndices.map(function (i) {
          return [
            { text: 'x' + PPS.central.sub(result.columns[i - result.blockStart].rel), class: 'rowhead' },
            { text: fmt(xs[i], 6) },
            { text: fmt(ys[i], 6) }
          ];
        });
        out.appendChild(ui.section(t('central.block'),
          ui.table([t('common.node'), 'x', 'y'], usedTable),
          result.blockSize + ' ' + t('common.nodes')));

        if (result.formula === 'stirling') out.appendChild(ui.note(t('central.stirlingNote')));
        else if (result.formula === 'bessel') out.appendChild(ui.note(t('central.besselNote')));
        else out.appendChild(ui.note(t('central.gaussNote')));
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
  }

  PPS.views.central = { render: render };
})();
