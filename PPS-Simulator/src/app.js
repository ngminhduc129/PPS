/* ============================================================
   app.js — navigation, theme, language switching
   ============================================================ */

(function () {
  const PPS = window.PPS;
  const ui = PPS.ui;
  const t = function (k, v) { return PPS.i18n.t(k, v); };

  const ALGORITHMS = [
    {
      id: 'horner', num: 1,
      titleKey: 'alg.horner.title', tagKey: 'alg.horner.tag', descKey: 'alg.horner.desc'
    },
    {
      id: 'chebyshev', num: 2,
      titleKey: 'alg.chebyshev.title', tagKey: 'alg.chebyshev.tag', descKey: 'alg.chebyshev.desc'
    },
    {
      id: 'lagrange', num: 3,
      titleKey: 'alg.lagrange.title', tagKey: 'alg.lagrange.tag', descKey: 'alg.lagrange.desc'
    },
    {
      id: 'newton-div', num: 4,
      titleKey: 'alg.newtonDiv.title', tagKey: 'alg.newtonDiv.tag', descKey: 'alg.newtonDiv.desc'
    },
    {
      id: 'newton-eq', num: 5,
      titleKey: 'alg.newtonEq.title', tagKey: 'alg.newtonEq.tag', descKey: 'alg.newtonEq.desc'
    },
    {
      id: 'central', num: 6,
      titleKey: 'alg.central.title', tagKey: 'alg.central.tag', descKey: 'alg.central.desc'
    },
    {
      id: 'inverse', num: 7,
      titleKey: 'alg.inverse.title', tagKey: 'alg.inverse.tag', descKey: 'alg.inverse.desc'
    },
    {
      id: 'spline', num: 8,
      titleKey: 'alg.spline.title', tagKey: 'alg.spline.tag', descKey: 'alg.spline.desc'
    },
    {
      id: 'least-squares', num: 9,
      titleKey: 'alg.leastSquares.title', tagKey: 'alg.leastSquares.tag', descKey: 'alg.leastSquares.desc'
    },
    {
      id: 'comparison', num: 10,
      titleKey: 'alg.comparison.title', tagKey: 'alg.comparison.tag', descKey: 'alg.comparison.desc'
    }
  ];

  const GROUPS = [
    { key: 'nav.interp', ids: ['horner', 'chebyshev', 'lagrange', 'newton-div', 'newton-eq'] },
    { key: 'nav.center', ids: ['central', 'inverse'] },
    { key: 'nav.fit', ids: ['spline', 'least-squares'] },
    { key: 'nav.comp', ids: ['comparison'] }
  ];

  const sidebar = document.getElementById('sidebar');
  const main = document.getElementById('main');
  const themeBtn = document.getElementById('theme-toggle');
  const precisionInput = document.getElementById('precision-input');
  const copyReportBtn = document.getElementById('copy-report-btn');
  const exportBtn = document.getElementById('export-btn');
  const langButtons = document.querySelectorAll('.lang-btn');

  function currentId() {
    const hash = (location.hash || '').replace('#', '');
    if (ALGORITHMS.some(function (a) { return a.id === hash; })) return hash;
    return 'horner';
  }

  function algorithm(id) {
    return ALGORITHMS.find(function (a) { return a.id === id; }) || ALGORITHMS[0];
  }

  function renderSidebar() {
    sidebar.innerHTML = '';
    const active = currentId();
    GROUPS.forEach(function (group) {
      sidebar.appendChild(ui.el('div', { class: 'nav-group-title', text: t(group.key) }));
      group.ids.forEach(function (id) {
        const meta = algorithm(id);
        const btn = ui.el('button', {
          type: 'button',
          class: 'nav-item' + (id === active ? ' active' : ''),
          'aria-current': id === active ? 'page' : 'false'
        });
        btn.appendChild(ui.el('span', { class: 'nav-num', text: String(meta.num) }));
        const textWrap = ui.el('span', { class: 'nav-item-text' });
        textWrap.appendChild(ui.el('span', { class: 'nav-item-title', text: t(meta.titleKey) }));
        textWrap.appendChild(ui.el('span', { class: 'nav-item-tag', text: t(meta.tagKey) }));
        btn.appendChild(textWrap);
        btn.addEventListener('click', function () {
          if (currentId() === id) {
            renderCurrent();
          } else {
            location.hash = id;
          }
        });
        sidebar.appendChild(btn);
      });
    });
  }

  function renderCurrent() {
    const id = currentId();
    const meta = algorithm(id);
    if (PPS.anim && PPS.anim.stopAll) PPS.anim.stopAll();
    main.innerHTML = '';
    main.appendChild(ui.pageHeader(
      t(meta.titleKey),
      t(meta.descKey),
      '#' + meta.num + ' · ' + t(meta.tagKey)
    ));

    const root = ui.el('div', { id: 'view-root' });
    main.appendChild(root);

    const view = PPS.views[id];
    if (view && view.render) {
      try {
        view.render(root);
      } catch (error) {
        console.error(error);
        root.appendChild(ui.errorBox(error && error.message ? error.message : String(error)));
      }
    } else {
      root.appendChild(ui.errorBox('View not found: ' + id));
    }

    document.title = t(meta.titleKey) + ' — PPS Simulator';
    window.scrollTo({ top: 0, behavior: 'smooth' });
  }

  function applyLanguage(lang) {
    PPS.i18n.setLang(lang);
    PPS.i18n.apply();
    langButtons.forEach(function (btn) {
      btn.classList.toggle('active', btn.getAttribute('data-lang') === lang);
    });
    renderSidebar();
    renderCurrent();
    applyTheme(document.documentElement.getAttribute('data-theme') || 'dark');
  }

  function applyTheme(theme) {
    document.documentElement.setAttribute('data-theme', theme);
    themeBtn.textContent = theme === 'dark' ? '☀️' : '🌙';
    themeBtn.setAttribute('title', theme === 'dark' ? t('app.themeToLight') : t('app.themeToDark'));
    try { localStorage.setItem('pps.theme', theme); } catch (e) { /* ignore */ }
  }

  function initTheme() {
    let theme = 'dark';
    try {
      const stored = localStorage.getItem('pps.theme');
      if (stored === 'light' || stored === 'dark') theme = stored;
      else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) theme = 'light';
    } catch (e) { /* ignore */ }
    applyTheme(theme);
  }

  function init() {
    initTheme();
    PPS.i18n.apply();
    langButtons.forEach(function (btn) {
      btn.classList.toggle('active', btn.getAttribute('data-lang') === PPS.i18n.lang);
      btn.addEventListener('click', function () {
        const lang = btn.getAttribute('data-lang');
        if (lang !== PPS.i18n.lang) applyLanguage(lang);
      });
    });
    themeBtn.addEventListener('click', function () {
      const current = document.documentElement.getAttribute('data-theme');
      applyTheme(current === 'dark' ? 'light' : 'dark');
    });
    if (precisionInput) {
      if (PPS.format && PPS.format.getDecimals) {
        precisionInput.value = PPS.format.getDecimals();
      }
      precisionInput.addEventListener('input', function () {
        const val = parseInt(this.value, 10);
        if (!isNaN(val) && val >= 0 && val <= 15) {
          PPS.format.setDecimals(val);
          const runBtn = document.querySelector('#view-root button.btn-primary');
          if (runBtn) {
            runBtn.click();
          }
        }
      });
      precisionInput.addEventListener('change', function () {
        let val = parseInt(this.value, 10);
        if (isNaN(val) || val < 0) val = 0;
        else if (val > 15) val = 15;
        this.value = val;
        PPS.format.setDecimals(val);
        const runBtn = document.querySelector('#view-root button.btn-primary');
        if (runBtn) {
          runBtn.click();
        }
      });
    }
    if (copyReportBtn) {
      copyReportBtn.addEventListener('click', function () {
        const root = document.getElementById('view-root');
        if (!PPS.export || !PPS.export.copyCurrent(root, copyReportBtn)) {
          alert(t('common.noResult'));
        }
      });
    }
    exportBtn.addEventListener('click', function () {
      const root = document.getElementById('view-root');
      const id = currentId();
      if (!PPS.export || !PPS.export.exportCurrent(root, id)) {
        alert(t('common.noResult'));
      }
    });
    window.addEventListener('hashchange', function () {
      renderSidebar();
      renderCurrent();
    });
    renderSidebar();
    renderCurrent();
  }

  PPS.app = {
    algorithms: ALGORITHMS,
    refresh: function () {
      renderSidebar();
      renderCurrent();
    }
  };

  // UI is loaded when app.js runs; wait for DOM just in case.
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
