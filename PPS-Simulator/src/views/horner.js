/* ============================================================
   View: Horner scheme
   Updated based on Section 1 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;
  const poly = PPS.poly;

  const PRESETS = [
    { coef: [2, -4, 5, -3], c: 2, k: 3 }, // Ví dụ chuẩn Section 1.3
    { coef: [2, -3, 0, 4, -5], c: 1.5, k: 2 },
    { coef: [1, 0, -2, 0, 3, -1], c: -1, k: 1 }
  ];

  function arrayNode(values, highlightIndex) {
    const span = ui.el('span');
    values.forEach(function (value, index) {
      if (index > 0) span.appendChild(document.createTextNode(', '));
      const part = ui.el('span', { text: fmt(value, 6) });
      if (index === highlightIndex) {
        part.className = 'accent';
        part.style.fontWeight = '800';
      }
      span.appendChild(part);
    });
    return span;
  }

  function render(root) {
    let presetIndex = 0;

    const coefIn = ui.input({ value: PRESETS[0].coef.join(', ') });
    const cIn = ui.input({ value: String(PRESETS[0].c) });
    const kIn = ui.input({ value: String(PRESETS[0].k) });
    const out = ui.el('div');

    const runBtn = ui.button('▶ ' + t('common.run'), function () { run(); }, { variant: 'primary' });
    const exampleBtn = ui.button('✨ ' + t('common.example'), function () {
      presetIndex = (presetIndex + 1) % PRESETS.length;
      const preset = PRESETS[presetIndex];
      coefIn.value = preset.coef.join(', ');
      cIn.value = String(preset.c);
      kIn.value = String(preset.k);
      run();
    }, { variant: 'ghost' });

    const form = ui.card(t('common.input'),
      ui.grid(
        ui.field(t('horner.coeffs'), coefIn, t('horner.coeffsHint'), { span2: true }),
        ui.field(t('horner.c'), cIn, t('horner.cHint')),
        ui.field(t('horner.k'), kIn, t('horner.kHint'))
      ),
      ui.actions(runBtn, exampleBtn)
    );

    root.appendChild(form);
    root.appendChild(out);

    function run() {
      out.innerHTML = '';
      try {
        const coef = ui.parseList(coefIn.value);
        const c = ui.parseNumber(cIn.value);
        const kRaw = ui.parseNumber(kIn.value, 'invalidK');
        if (kRaw < 0 || Math.abs(kRaw - Math.round(kRaw)) > 1e-9) throw new Error('invalidK');
        const k = Math.round(kRaw);
        const result = PPS.horner.evaluate(coef, c, k);

        const kLabel = 'P' + (k === 0 ? '' : '(' + PPS.format.superscript(k) + ')') + '(c)';
        out.appendChild(ui.kv([
          { key: 'P(c)', value: fmt(result.value, 8), className: 'accent' },
          { key: kLabel, value: fmt(result.derivative, 8), className: result.kTooBig ? 'orange' : 'green' },
          { key: t('horner.degree'), value: String(result.N) },
          { key: t('horner.passes'), value: String(result.passes.length) }
        ]));

        if (result.kTooBig) {
          out.appendChild(ui.note(t('horner.kTooBig'), 'warn'));
        }

        // Đa thức ban đầu
        const polyOriginal = PPS.format.polyToString(coef.slice().reverse(), 'x');
        const polyCard = ui.section(t('common.polynomial'),
          ui.formula(polyOriginal, true));
        out.appendChild(polyCard);

        // Đa thức thương Q(x) khi chia cho (x - c)
        if (result.quotientDesc && result.quotientDesc.length) {
          const qStr = PPS.format.polyToString(result.quotientDesc.slice().reverse(), 'x');
          const remStr = fmt(result.value, 6);
          const cSign = c >= 0 ? ' - ' + fmt(c, 4) : ' + ' + fmt(-c, 4);
          const divStr = 'P(x) = (x' + cSign + ') · [' + qStr + '] + ' + remStr;
          out.appendChild(ui.section(t('horner.divisionTitle'), ui.formula(divStr, true)));
        }

        // Khai triển Taylor tại lân cận điểm c
        const taylorStrParts = [];
        result.taylorTerms.forEach(function (term) {
          const m = term.m;
          const rem = term.rem;
          if (Math.abs(rem) < 1e-12) return;
          const remStr = fmt(rem, 6);
          if (m === 0) {
            taylorStrParts.push(remStr);
          } else {
            const cStr = c >= 0 ? '(x - ' + fmt(c, 4) + ')' : '(x + ' + fmt(-c, 4) + ')';
            const powStr = m === 1 ? '' : PPS.format.superscript(m);
            taylorStrParts.push(remStr + ' ' + cStr + powStr);
          }
        });
        const taylorFull = 'P(x) = ' + (taylorStrParts.join(' + ').replace(/\+\s*-/g, '- ') || '0');
        out.appendChild(ui.section(t('horner.taylorTitle'), ui.formula(taylorFull, true)));

        // Trích xuất các Vector kết quả Horner
        const aVectorStr = 'a = [' + coef.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const qVectorStr = result.quotientDesc && result.quotientDesc.length
          ? 'Q = [' + result.quotientDesc.map(function (v) { return fmt(v, 6); }).join(', ') + ']'
          : null;
        const bVectorStr = 'b_Taylor = [' + result.passes.map(function (p) { return fmt(p.remainder, 6); }).join(', ') + ']';
        const dVectorStr = 'P⁽ᵐ⁾(c) = [' + result.passes.map(function (p) { return fmt(p.derivativeValue, 6); }).join(', ') + ']';

        const vectorChildren = [
          ui.el('p', { class: 'math-label', text: 'Vector hệ số đa thức ban đầu (lũy thừa giảm dần a_n … a_0):' }),
          ui.formula(aVectorStr, true)
        ];
        if (qVectorStr) {
          vectorChildren.push(
            ui.el('p', { class: 'math-label', text: 'Vector hệ số đa thức thương Q(x) khi chia cho (x − c):', style: 'margin-top: 8px;' }),
            ui.formula(qVectorStr, true)
          );
        }
        vectorChildren.push(
          ui.el('p', { class: 'math-label', text: 'Vector hệ số Taylor b_m = P⁽ᵐ⁾(c) / m! (các phần dư):', style: 'margin-top: 8px;' }),
          ui.formula(bVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector giá trị đạo hàm các cấp tại c: P⁽ᵐ⁾(c):', style: 'margin-top: 8px;' }),
          ui.formula(dVectorStr, true)
        );
        out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), vectorChildren));

        // Passes table
        const rows = result.passes.map(function (pass) {
          return [
            { text: String(pass.m), class: 'rowhead' },
            { node: arrayNode(pass.before.slice(0, result.N - pass.m + 1)), class: 'left' },
            { node: arrayNode(pass.after.slice(0, result.N - pass.m + 1), pass.remainderIndex), class: 'left' },
            { text: fmt(pass.remainder, 8), class: 'hl' },
            { text: fmt(pass.derivativeValue, 8) }
          ];
        });
        const table = ui.table([
          t('horner.pass'),
          t('horner.before'),
          t('horner.after'),
          t('horner.remainder'),
          t('horner.derivativeValue')
        ], rows);

        // Chart + animation
        const span = Math.max(3, Math.abs(c) * 0.6);
        const curve = poly.sample(result.passes[0].before.slice().reverse(), c - span, c + span, 300);
        const tangentSeries = result.derivatives.length > 1 ? {
          points: [
            { x: c - span, y: result.value - result.derivatives[1] * span },
            { x: c + span, y: result.value + result.derivatives[1] * span }
          ],
          color: 'var(--orange)',
          width: 1.7,
          dash: '7 5',
          label: t('horner.tangent')
        } : null;

        const chart = ui.chart(t('common.chart'), 360);
        const baseSeries = [{
          points: curve,
          color: 'var(--accent)',
          width: 2.4,
          label: 'P(x)'
        }].concat(tangentSeries ? [tangentSeries] : []);

        const handle = chart.draw({
          series: baseSeries,
          points: [{ x: c, y: result.value, color: 'var(--green)', r: 5.5, label: 'c = ' + fmt(c, 6) }],
          highlightPoints: [{ x: c, y: result.value, color: 'var(--green)' }],
          vline: { x: c, label: 'c = ' + fmt(c, 6), color: 'var(--green)' },
          xLabel: 'x',
          yLabel: 'P(x)'
        });

        const steps = result.passes.length;
        const anim = ui.animBar({
          count: steps,
          delay: 720,
          onStep: function (index) {
            ui.highlightRow(table, index);
            const reveal = (index + 1) / steps;
            const series = [{
              points: curve,
              color: 'var(--accent)',
              width: 2.4,
              label: 'P(x)',
              reveal: reveal
            }];
            if (tangentSeries && index === steps - 1) series.push(tangentSeries);
            handle.update({
              series: series,
              marker: { x: c, y: result.value, color: 'var(--green)', label: 'P(c) = ' + fmt(result.value, 6) }
            });
          },
          onReset: function () {
            ui.clearAnim(table);
            handle.update({ series: baseSeries, marker: null });
          }
        });

        out.appendChild(ui.section(t('horner.passTable'), [
          table,
          anim.node
        ], steps + ' ' + t('common.steps').toLowerCase()));

        out.appendChild(ui.note(t('horner.noteFormula')));
        out.appendChild(chart.node);
      } catch (error) {
        out.appendChild(ui.errorBox(ui.errText(error)));
      }
    }

    run();
  }

  PPS.views.horner = { render: render };
})();
