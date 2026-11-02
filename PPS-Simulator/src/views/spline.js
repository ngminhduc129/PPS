/* ============================================================
   View: Linear, Quadratic and Cubic Splines
   Supports Natural, Given Moments, and Clamped Boundary Conditions
   Updated based on Section 8 of Algorithm.md
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };
  const fmt = PPS.format.fmt;
  const polyToString = PPS.format.polyToString;

  const PRESETS = [
    { xs: [0, 1, 2, 3, 4], ys: [1, 2, 0, 2, 1], z0: 0, bc: 'natural', bVal0: '0', bVal1: '0' },
    { xs: [1, 2, 3, 4, 5], ys: [2, 3, 5, 4, 6], z0: 1, bc: 'clamped', bVal0: '1', bVal1: '0.5' },
    { xs: [0, 0.5, 1, 1.5, 2], ys: [0, 0.4794, 0.8415, 0.9975, 0.9093], z0: 1, bc: 'moments', bVal0: '0', bVal1: '-0.8' }
  ];

  const COLORS = ['#5b9dff', '#3ddc97', '#ffb454', '#c08cff', '#4fd6e3', '#ff6b81', '#7c6cff'];

  function render(root) {
    let presetIndex = 0;
    const xy = ui.xyControl({
      xs: PRESETS[0].xs,
      ys: PRESETS[0].ys,
      xLabel: t('spline.xs'),
      yLabel: t('spline.ys'),
      xHint: t('lagrange.xsHint'),
      yHint: t('lagrange.xsHint')
    });

    const degSel = ui.select([
      { value: '3', label: t('spline.deg3') },
      { value: '2', label: t('spline.deg2') },
      { value: '1', label: t('spline.deg1') }
    ], '3');

    const bcSel = ui.select([
      { value: 'natural', label: t('spline.bcNatural') },
      { value: 'moments', label: t('spline.bcMoments') },
      { value: 'clamped', label: t('spline.bcClamped') }
    ], 'natural');

    const bVal0In = ui.input({ value: '0' });
    const bVal1In = ui.input({ value: '0' });
    const z0In = ui.input({ value: '0' });

    const out = ui.el('div');

    const form = ui.card(t('common.input'),
      xy.node,
      ui.grid(
        ui.field(t('spline.degree'), degSel),
        ui.field(t('spline.boundaryType'), bcSel),
        ui.field(t('spline.bcVal0'), bVal0In, t('spline.bcVal0Hint')),
        ui.field(t('spline.bcVal1'), bVal1In, t('spline.bcVal1Hint')),
        ui.field(t('spline.boundaryZ0'), z0In, t('spline.boundaryZ0Hint'))
      ),
      ui.actions(
        ui.button('▶ ' + t('common.run'), run, { variant: 'primary' }),
        ui.button('✨ ' + t('common.example'), function () {
          presetIndex = (presetIndex + 1) % PRESETS.length;
          const p = PRESETS[presetIndex];
          xy.setValues(p.xs, p.ys);
          degSel.value = '3';
          bcSel.value = p.bc;
          bVal0In.value = p.bVal0;
          bVal1In.value = p.bVal1;
          z0In.value = String(p.z0);
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
        if (xs.length !== ys.length) throw new Error('needPoints');

        const deg = degSel.value;
        const bc = bcSel.value;
        const b0 = bVal0In.value.trim() ? ui.parseNumber(bVal0In.value) : 0;
        const b1 = bVal1In.value.trim() ? ui.parseNumber(bVal1In.value) : 0;
        const z0 = z0In.value.trim() ? ui.parseNumber(z0In.value) : 0;

        let spline;
        if (deg === '1') {
          spline = PPS.spline.linear(xs, ys);
        } else if (deg === '2') {
          spline = PPS.spline.quadratic(xs, ys, z0);
        } else {
          spline = PPS.spline.cubic(xs, ys, {
            bcType: bc,
            M0: b0,
            Mn: b1,
            s0: b0,
            sn: b1
          });
        }

        const summary = [
          { key: t('common.degreeN'), value: String(spline.degree) },
          { key: t('common.segments'), value: String(spline.segments.length) },
          { key: t('common.nodes'), value: String(xs.length) }
        ];
        if (deg === '3') {
          summary.push({ key: t('spline.bcApplied'), value: bc.toUpperCase(), className: 'accent' });
        }
        out.appendChild(ui.kv(summary));

        // Bảng đa thức từng đoạn S_i(x)
        const segHeaders = [t('common.segment'), t('common.interval'), t('spline.localFormula'), t('common.canonicalPoly')];
        const segRows = spline.segments.map(function (seg, idx) {
          const rangeText = '[' + fmt(seg.x0, 4) + ', ' + fmt(seg.x1, 4) + ']';
          let localStr = '—';
          if (seg.local) {
            const loc = seg.local;
            const xi = fmt(loc.xi, 4);
            const dx = '(x - ' + xi + ')';
            const parts = [fmt(loc.A, 6)];
            if (Math.abs(loc.B) > 1e-12) parts.push(fmt(loc.B, 6) + ' ' + dx);
            if (Math.abs(loc.C) > 1e-12) parts.push(fmt(loc.C, 6) + ' ' + dx + '²');
            if (loc.D && Math.abs(loc.D) > 1e-12) parts.push(fmt(loc.D, 6) + ' ' + dx + '³');
            localStr = parts.join(' + ').replace(/\+\s*-/g, '- ');
          }
          return [
            { text: 'S' + PPS.format.subscript(idx), class: 'rowhead' },
            { text: rangeText },
            { text: localStr, class: 'hl' },
            { text: polyToString(seg.coef, 'x') }
          ];
        });
        out.appendChild(ui.section(t('spline.segTableTitle'), ui.table(segHeaders, segRows)));

        // Trích xuất các Vector kết quả Spline
        const hArr = [];
        for (let i = 0; i < xs.length - 1; i++) hArr.push(xs[i + 1] - xs[i]);
        const xVectorStr = 'X = [' + xs.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const yVectorStr = 'Y = [' + ys.map(function (v) { return fmt(v, 6); }).join(', ') + ']';
        const hVectorStr = 'h = [' + hArr.map(function (v) { return fmt(v, 6); }).join(', ') + ']';

        const splineVectorNodes = [
          ui.el('p', { class: 'math-label', text: 'Vector mốc ghép X:' }),
          ui.formula(xVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector giá trị hàm Y:', style: 'margin-top: 8px;' }),
          ui.formula(yVectorStr, true),
          ui.el('p', { class: 'math-label', text: 'Vector độ dài các khoảng h_k = x_k − x_{k−1}:', style: 'margin-top: 8px;' }),
          ui.formula(hVectorStr, true)
        ];
        if (spline.M) {
          const mVectorStr = 'M = [' + spline.M.map(function (mVal) { return fmt(mVal, 6); }).join(', ') + ']';
          splineVectorNodes.push(
            ui.el('p', { class: 'math-label', text: 'Vector moment uốn M_i (đạo hàm cấp 2 tại các nút):', style: 'margin-top: 8px;' }),
            ui.formula(mVectorStr, true)
          );
        }
        if (spline.z) {
          const zVectorStr = 'z = [' + spline.z.map(function (zVal) { return fmt(zVal, 6); }).join(', ') + ']';
          splineVectorNodes.push(
            ui.el('p', { class: 'math-label', text: 'Vector hệ số góc z_i = S′(x_i):', style: 'margin-top: 8px;' }),
            ui.formula(zVectorStr, true)
          );
        }
        out.appendChild(ui.section(t('newtonDiv.vectorsTitle'), splineVectorNodes));

        // Hệ phương trình 3 đường chéo và nghiệm TDMA nếu là bậc 3
        if (deg === '3' && spline.system && spline.system.length) {
          const sysHeaders = ['Phương trình', 'a_i (đường chéo dưới)', 'b_i (đường chéo chính)', 'c_i (đường chéo trên)', 'd_i (vế phải)'];
          const sysRows = spline.system.map(function (eq) {
            return [
              { text: 'Eq ' + eq.i, class: 'rowhead' },
              { text: fmt(eq.left, 6) },
              { text: fmt(eq.mid, 6), class: 'hl' },
              { text: fmt(eq.right, 6) },
              { text: fmt(eq.rhs, 6) }
            ];
          });
          const mVectorStr = 'M = [' + spline.M.map(function (mVal) { return fmt(mVal, 6); }).join(', ') + ']';
          out.appendChild(ui.section(t('spline.tridiagonalTitle'), [
            ui.table(sysHeaders, sysRows),
            ui.el('p', { class: 'math-label', text: 'Vector moment uốn M_i giải được bằng Thomas TDMA:' }),
            ui.formula(mVectorStr, true)
          ]));
        }

        // Đồ thị ghép trơn
        const chart = ui.chart(t('common.chart'), 380);
        const series = [];
        spline.segments.forEach(function (seg, idx) {
          const pts = [];
          const step = (seg.x1 - seg.x0) / 40;
          for (let s = 0; s <= 40; s++) {
            const x = seg.x0 + s * step;
            pts.push({ x: x, y: seg.evaluate(x) });
          }
          series.push({
            points: pts,
            color: COLORS[idx % COLORS.length],
            width: 2.5,
            label: 'S' + PPS.format.subscript(idx)
          });
        });

        const nodeMarkers = xs.map(function (x, i) {
          return { x: x, y: ys[i], color: 'var(--accent)', r: 5, label: '(' + fmt(x, 2) + ', ' + fmt(ys[i], 2) + ')' };
        });

        chart.draw({
          series: series,
          points: nodeMarkers,
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

  PPS.views.spline = { render: render };
})();
