/* ============================================================
   View: Lagrange interpolation
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;
  const { polyToString } = PPS.format;

  const PRESETS = [
    { xs: [-1, 0, 1, 2], ys: [2, -1, 3, 1], xbar: 0.5 },
    { xs: [0, 1, 2, 3], ys: [1, 2, 0, 5], xbar: 1.5 },
    { xs: [-2, -1, 1, 3], ys: [4, 1, 1, 9], xbar: 0 }
  ];

  function padPolyDesc(desc, totalDegree) {
    const res = new Array(totalDegree).fill(0);
    const start = totalDegree - desc.length;
    for (let i = 0; i < desc.length; i++) {
      res[start + i] = desc[i];
    }
    return res;
  }

  function render(root) {
    let presetIndex = 0;
    const xy = ui.xyControl({
      xs: PRESETS[0].xs,
      ys: PRESETS[0].ys,
      xLabel: t('lagrange.xs'),
      yLabel: t('lagrange.ys'),
      xHint: t('lagrange.xsHint'),
      yHint: t('lagrange.ysHint')
    });
    const xbarIn = ui.input({ value: String(PRESETS[0].xbar) });
    const out = ui.el('div');

    const form = ui.card(t('common.input'),
      xy.node,
      ui.grid(
        ui.field(t('lagrange.xEval'), xbarIn, t('lagrange.xEvalHint'))
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
          const imp = ui.importDataButton({
            onParsed: function(res){
              xy.setValues(res.xs || [], res.ys || []);
              run();
            }
          });
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
        const result = PPS.lagrange.interpolate(xs, ys);
        const xbarText = xbarIn.value.trim();
        let xbar = null;
        if (xbarText) {
          xbar = ui.parseNumber(xbarText, 'xbar');
        }

        const summary = [
          { key: t('common.degreeN'), value: String(result.degree) },
          { key: t('common.nodes'), value: String(xs.length) }
        ];
        if (xbar != null) {
          summary.push({ key: 'Lₙ(x̄)', value: fmt(result.evaluate(xbar), 8), className: 'accent' });
        }
        out.appendChild(ui.kv(summary));

        out.appendChild(ui.section(t('common.polynomial'),
          ui.formula(polyToString(result.coef, 'x'), true)));

        // Trích xuất các Vector kết quả Lagrange
        const n = xs.length;
        const xVectorStr = 'X = [' + xs.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const yVectorStr = 'Y = [' + ys.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const dVectorStr = 'D = [' + result.basis.map(function (b) { return fmt(b.Di, 6); }).join(', ') + ']';
        const wVectorStr = 'w = [' + result.basis.map(function (b) { return fmt(1 / b.Di, 6); }).join(', ') + ']';
        const coefDesc = result.coef.slice().reverse();
        const pVectorStr = 'P_coef (xⁿ … x⁰) = [' + coefDesc.map(function (c) { return fmt(c, 6); }).join(', ') + ']';

        out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), [
          ui.el('p', { class: 'math-label', text: 'Vector mốc dữ liệu X:' }),
          ui.formula(xVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector giá trị hàm Y:', style: 'margin-top: 8px;' }),
          ui.formula(yVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector mẫu số Lagrange D_i = ∏_{j≠i} (x_i − x_j):', style: 'margin-top: 8px;' }),
          ui.formula(dVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector trọng số Barycentric w_i = 1 / D_i:', style: 'margin-top: 8px;' }),
          ui.formula(wVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector hệ số đa thức chính tắc (lũy thừa giảm dần):', style: 'margin-top: 8px;' }),
          ui.formula(pVectorStr, true)
        ]));

        const basisRows = result.basis.map(function (b) {
          return [
            { text: String(b.i), class: 'rowhead' },
            { text: fmt(b.xi, 6) },
            { text: fmt(b.yi, 6) },
            { text: fmt(b.Di, 6) },
            { text: polyToString(b.wi, 'x'), class: 'left' },
            { text: polyToString(b.li, 'x'), class: 'left' }
          ];
        });
        const basisTable = ui.table([
          t('lagrange.node'),
          'xᵢ',
          'yᵢ',
          t('lagrange.di'),
          t('lagrange.wi'),
          t('lagrange.li')
        ], basisRows);
        const basisSection = ui.section(t('lagrange.basisTable'), basisTable);
        out.appendChild(basisSection);

        // 1. BẢNG TÍCH (Product Table): Tính omega(x) = prod(x - x_j)
        const omegaFormulaStr = 'ω(x) = ' + polyToString(result.omega, 'x');
        const prodHeaders = ['Bước', t('lagrange.factor')];
        for (let d = n; d >= 0; d--) {
          prodHeaders.push(d === 0 ? 'x⁰' : (d === 1 ? 'x' : 'x' + PPS.format.superscript(d)));
        }
        const prodRows = result.productTable.map(function (stepRow) {
          const cells = [
            { text: String(stepRow.step), class: 'rowhead' },
            { text: stepRow.factor, class: 'left' }
          ];
          stepRow.polyDesc.forEach(function (coefVal) {
            cells.push({ text: fmt(coefVal) });
          });
          return cells;
        });
        out.appendChild(ui.section(t('lagrange.productTable'), [
          ui.el('p', { class: 'math-label', text: t('lagrange.omegaFormula') + ':' }),
          ui.formula(omegaFormulaStr, true),
          ui.table(prodHeaders, prodRows)
        ]));

        // 2. BẢNG CHIA (Division Table): Chia Horner omega(x) cho (x - x_i) để tìm w_i(x)
        const divHeaders = [
          t('lagrange.node'),
          'xᵢ',
          t('lagrange.quotient'),
          t('lagrange.remainder'),
          t('lagrange.di'),
          'yᵢ',
          'cᵢ = yᵢ/Dᵢ',
          t('lagrange.contribPoly')
        ];
        const divRows = result.divisionTable.map(function (divRow) {
          return [
            { text: String(divRow.i), class: 'rowhead' },
            { text: fmt(divRow.xi) },
            { text: polyToString(divRow.wiAsc, 'x'), class: 'left hl' },
            { text: fmt(divRow.remainder), class: 'center' },
            { text: fmt(divRow.Di) },
            { text: fmt(divRow.yi) },
            { text: fmt(divRow.ci) },
            { text: polyToString(divRow.termPoly, 'x'), class: 'left' }
          ];
        });
        out.appendChild(ui.section(t('lagrange.divisionTable'), ui.table(divHeaders, divRows)));

        // 3. MA TRẬN HỆ SỐ TỬ SỐ w_i(x) (Hệ số theo lũy thừa giảm dần)
        const wiHeaders = ['Mốc xᵢ'];
        for (let d = n - 1; d >= 0; d--) {
          wiHeaders.push(d === 0 ? 'x⁰' : (d === 1 ? 'x' : 'x' + PPS.format.superscript(d)));
        }
        const wiRows = result.divisionTable.map(function (b) {
          const cells = [
            { text: fmt(b.xi), class: 'rowhead' }
          ];
          b.wiDesc.forEach(function (coefVal) {
            cells.push({ text: fmt(coefVal) });
          });
          return cells;
        });
        out.appendChild(ui.section(t('lagrange.wi') + ' — wᵢ(x)', ui.table(wiHeaders, wiRows)));

        const coefRows = result.coef.map(function (c, power) {
          return [
            { text: String(power), class: 'rowhead' },
            { text: power === 0 ? '1' : (power === 1 ? 'x' : 'x' + PPS.format.superscript(power)), class: 'left' },
            { text: fmt(c, 8) }
          ];
        });
        out.appendChild(ui.section(t('lagrange.finalCoefs'),
          ui.table([t('common.degree'), t('lagrange.coefDegree'), t('common.coefficient')], coefRows)));

        const checks = result.basis.map(function (b) {
          const value = result.evaluate(b.xi);
          return [
            { text: String(b.i), class: 'rowhead' },
            { text: fmt(b.xi, 6) },
            { text: fmt(b.yi, 6) },
            { text: fmt(value, 8), class: Math.abs(value - b.yi) < 1e-9 ? 'hl-green' : 'hl-orange' },
            { text: fmt(value - b.yi, 3), class: 'dim' }
          ];
        });
        out.appendChild(ui.section(t('lagrange.nodeCheck'),
          ui.table([t('lagrange.node'), 'xᵢ', t('lagrange.trueY'), t('lagrange.computedY'), t('lagrange.error')], checks)));

        const chart = ui.chart(t('common.chart'), 370);
        out.appendChild(chart.node);

        const xmin = Math.min.apply(null, xs);
        const xmax = Math.max.apply(null, xs);
        const span = Math.max(0.5, (xmax - xmin) * 0.12);
        const lo = xmin - span;
        const hi = xmax + span;
        const curve = PPS.chart.sample(result.evaluate, lo, hi, 320);

        // Partial sums of yᵢ·lᵢ(x) for the step animation
        const partialCoefs = [];
        let running = [0];
        result.basis.forEach(function (b) {
          running = PPS.poly.add(running, PPS.poly.scale(b.li, b.yi));
          partialCoefs.push(running);
        });

        const cfg = {
          series: [{ points: curve, color: 'var(--accent)', width: 2.4, label: 'Lₙ(x)' }],
          points: xs.map(function (x, i) {
            return { x: x, y: ys[i], color: 'var(--green)', label: 'x = ' + fmt(x, 5) + ', y = ' + fmt(ys[i], 5) };
          }),
          xLabel: 'x',
          yLabel: 'y'
        };
        if (xbar != null) {
          const ybar = result.evaluate(xbar);
          cfg.vline = { x: xbar, label: 'x̄ = ' + fmt(xbar, 6), color: 'var(--orange)' };
          cfg.highlightPoints = [{ x: xbar, y: ybar, color: 'var(--orange)', label: 'Lₙ(x̄) = ' + fmt(ybar, 6) }];
        }
        const handle = chart.draw(cfg);

        const steps = result.basis.length;
        const anim = ui.animBar({
          count: steps,
          delay: 700,
          onStep: function (index) {
            ui.highlightRow(basisTable, index);
            const revealed = PPS.chart.sample(function (x) {
              return PPS.poly.evaluate(partialCoefs[index], x);
            }, lo, hi, 320);
            const update = {
              series: [{ points: revealed, color: 'var(--accent)', width: 2.4, label: 'Lₙ(x)' }]
            };
            if (xbar != null) {
              update.marker = {
                x: xbar,
                y: PPS.poly.evaluate(partialCoefs[index], xbar),
                color: 'var(--orange)',
                label: 'Σ l₀…l' + PPS.central.sub(index)
              };
            }
            handle.update(update);
          },
          onReset: function () {
            ui.clearAnim(basisTable);
            handle.update({
              series: [{ points: curve, color: 'var(--accent)', width: 2.4, label: 'Lₙ(x)' }],
              marker: null
            });
          }
        });
        basisSection.appendChild(anim.node);
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
  }

  PPS.views.lagrange = { render: render };
})();
