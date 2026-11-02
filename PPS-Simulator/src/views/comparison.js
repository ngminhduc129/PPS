/* ============================================================
   View: Algorithm comparison table (Section 10 - Algorithm.md)
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  PPS.views = PPS.views || {};
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };

  const ROWS_DATA = [
    {
      nameKey: 'comp.horner',
      nodesKey: 'comp.hornerNodes',
      purposeKey: 'comp.hornerPurpose',
      complexity: 'O((K+1)N)',
      loopKey: 'comp.hornerLoop'
    },
    {
      nameKey: 'comp.cheb',
      nodesKey: 'comp.chebNodes',
      purposeKey: 'comp.chebPurpose',
      complexity: 'O(N)',
      loopKey: 'comp.chebLoop'
    },
    {
      nameKey: 'comp.lagrange',
      nodesKey: 'comp.lagrangeNodes',
      purposeKey: 'comp.lagrangePurpose',
      complexity: 'O(N²)',
      loopKey: 'comp.lagrangeLoop'
    },
    {
      nameKey: 'comp.newtonDiv',
      nodesKey: 'comp.newtonDivNodes',
      purposeKey: 'comp.newtonDivPurpose',
      complexity: 'O(N²)',
      loopKey: 'comp.newtonDivLoop'
    },
    {
      nameKey: 'comp.newtonEq',
      nodesKey: 'comp.newtonEqNodes',
      purposeKey: 'comp.newtonEqPurpose',
      complexity: 'O(N²)',
      loopKey: 'comp.newtonEqLoop'
    },
    {
      nameKey: 'comp.central',
      nodesKey: 'comp.centralNodes',
      purposeKey: 'comp.centralPurpose',
      complexity: 'O(N²)',
      loopKey: 'comp.centralLoop'
    },
    {
      nameKey: 'comp.invFunc',
      nodesKey: 'comp.invFuncNodes',
      purposeKey: 'comp.invFuncPurpose',
      complexity: 'O(N²)',
      loopKey: 'comp.invFuncLoop'
    },
    {
      nameKey: 'comp.invIter',
      nodesKey: 'comp.invIterNodes',
      purposeKey: 'comp.invIterPurpose',
      complexity: 'O(N² + KN)',
      loopKey: 'comp.invIterLoop'
    },
    {
      nameKey: 'comp.spline',
      nodesKey: 'comp.splineNodes',
      purposeKey: 'comp.splinePurpose',
      complexity: 'O(N)',
      loopKey: 'comp.splineLoop'
    },
    {
      nameKey: 'comp.leastSquares',
      nodesKey: 'comp.leastSquaresNodes',
      purposeKey: 'comp.leastSquaresPurpose',
      complexity: 'O(Nm² + m³)',
      loopKey: 'comp.leastSquaresLoop'
    }
  ];

  function render(root) {
    const card = ui.card(t('comp.title'),
      ui.el('p', { class: 'page-desc', text: t('comp.subtitle') })
    );

    const headers = [
      t('comp.colAlg'),
      t('comp.colNodes'),
      t('comp.colPurpose'),
      t('comp.colComplexity'),
      t('comp.colLoop')
    ];

    const rows = ROWS_DATA.map(function (item) {
      return [
        { text: t(item.nameKey), class: 'hl' },
        { text: t(item.nodesKey) },
        { text: t(item.purposeKey) },
        { text: item.complexity, class: 'accent', style: 'font-family: monospace; font-weight: bold;' },
        { text: t(item.loopKey) }
      ];
    });

    const table = ui.table(headers, rows);
    card.appendChild(table);
    root.appendChild(card);
  }

  PPS.views.comparison = { render: render };
})();
