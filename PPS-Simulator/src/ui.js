/* ============================================================
   ui.js — small DOM helpers for the simulator views
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const t = function (key, vars) { return PPS.i18n.t(key, vars); };

  function el(tag, attrs) {
    const node = document.createElement(tag);
    attrs = attrs || {};
    Object.keys(attrs).forEach(function (key) {
      const value = attrs[key];
      if (value == null) return;
      if (key === 'class') node.className = value;
      else if (key === 'text') node.textContent = value;
      else if (key === 'style' && typeof value === 'object') {
        Object.keys(value).forEach(function (prop) { node.style[prop] = value[prop]; });
      } else if (key === 'dataset' && typeof value === 'object') {
        Object.keys(value).forEach(function (prop) { node.dataset[prop] = value[prop]; });
      } else if (key === 'onclick' && typeof value === 'function') {
        node.addEventListener('click', value);
      } else if (key === 'html') {
        node.innerHTML = value;
      } else {
        node.setAttribute(key, value);
      }
    });
    const children = Array.prototype.slice.call(arguments, 2);
    appendChildren(node, children);
    return node;
  }

  function appendChildren(node, children) {
    children.forEach(function (child) {
      if (child == null || child === false) return;
      if (Array.isArray(child)) appendChildren(node, child);
      else if (child instanceof Node) node.appendChild(child);
      else node.appendChild(document.createTextNode(String(child)));
    });
  }

  function field(label, control, hint, options) {
    options = options || {};
    const wrap = el('div', { class: 'field' + (options.span2 ? ' span-2' : '') });
    if (label) wrap.appendChild(el('label', { text: label }));
    wrap.appendChild(control);
    if (hint) wrap.appendChild(el('span', { class: 'hint', text: hint }));
    return wrap;
  }

  function input(opts) {
    opts = opts || {};
    const node = el('input', {
      class: 'input' + (opts.mono === false ? '' : ' mono'),
      type: 'text',
      value: opts.value == null ? '' : String(opts.value),
      placeholder: opts.placeholder || ''
    });
    if (opts.autocomplete) node.setAttribute('autocomplete', opts.autocomplete);
    return node;
  }

  function textarea(opts) {
    opts = opts || {};
    const node = el('textarea', {
      class: 'textarea' + (opts.mono === false ? '' : ' mono'),
      rows: opts.rows || 3,
      placeholder: opts.placeholder || ''
    });
    node.value = opts.value == null ? '' : String(opts.value);
    return node;
  }

  function select(options, value) {
    const node = el('select', { class: 'select' });
    options.forEach(function (opt) {
      const option = el('option', { value: opt.value, text: opt.label });
      if (String(opt.value) === String(value)) option.selected = true;
      node.appendChild(option);
    });
    return node;
  }

  function button(label, onClick, opts) {
    opts = opts || {};
    const node = el('button', {
      type: 'button',
      class: 'btn ' + (opts.variant === 'ghost' ? 'btn-ghost' : 'btn-primary') + (opts.small ? ' btn-sm' : ''),
      text: label
    });
    if (onClick) node.addEventListener('click', onClick);
    return node;
  }

  function card(title, children) {
    const wrap = el('section', { class: 'card' });
    if (title) wrap.appendChild(el('div', { class: 'card-title', text: title }));
    const rest = Array.prototype.slice.call(arguments, 2);
    appendChildren(wrap, [children].concat(rest));
    return wrap;
  }

  function grid(children) {
    const wrap = el('div', { class: 'form-grid' });
    appendChildren(wrap, Array.prototype.slice.call(arguments));
    return wrap;
  }

  function actions(children) {
    const wrap = el('div', { class: 'actions' });
    appendChildren(wrap, Array.prototype.slice.call(arguments));
    return wrap;
  }

  function note(text, kind, icon) {
    const wrap = el('div', { class: 'note' + (kind ? ' ' + kind : '') });
    wrap.appendChild(el('span', { class: 'icon', text: icon || (kind === 'warn' ? '⚠' : (kind === 'good' ? '✓' : 'ℹ')) }));
    wrap.appendChild(el('span', { text: text }));
    return wrap;
  }

  function errorBox(message) {
    return el('div', { class: 'error-box', text: message });
  }

  function errText(error) {
    const code = error && error.message ? error.message : 'range';
    const dict = PPS.i18n.dict[PPS.i18n.lang] || {};
    if (dict['errors.' + code]) return dict['errors.' + code];
    if (error && error.message && !/^[a-zA-Z]+$/.test(error.message)) return error.message;
    return t('errors.' + code);
  }

  function kv(items) {
    const wrap = el('div', { class: 'kv-grid' });
    items.forEach(function (item) {
      const cardEl = el('div', { class: 'kv' + (item.className ? ' ' + item.className : '') });
      cardEl.appendChild(el('span', { class: 'k', text: item.key }));
      cardEl.appendChild(el('span', { class: 'v' + (item.small ? ' small' : ''), text: item.value }));
      if (item.sub) cardEl.appendChild(el('span', { class: 'sub', text: item.sub }));
      wrap.appendChild(cardEl);
    });
    return wrap;
  }

  function table(headers, rows, opts) {
    opts = opts || {};
    const wrap = el('div', { class: 'table-wrap' });
    const tbl = el('table', { class: 'data-table' });
    const thead = el('thead');
    const headRow = el('tr');
    headers.forEach(function (header) {
      const th = el('th');
      if (header != null && typeof header === 'object') {
        if (header.text != null) th.textContent = header.text;
        if (header.colspan) th.setAttribute('colspan', header.colspan);
        if (header.rowspan) th.setAttribute('rowspan', header.rowspan);
        if (header.class) th.className = header.class;
      } else {
        th.textContent = header == null ? '' : String(header);
      }
      headRow.appendChild(th);
    });
    thead.appendChild(headRow);
    tbl.appendChild(thead);
    const tbody = el('tbody');
    rows.forEach(function (row, rowIndex) {
      const tr = el('tr');
      if (opts.rowClass) {
        const cls = typeof opts.rowClass === 'function' ? opts.rowClass(row, rowIndex) : opts.rowClass;
        if (cls) tr.className = cls;
      }
      row.forEach(function (cell) {
        const td = el('td');
        if (cell instanceof Node) {
          td.appendChild(cell);
        } else if (cell != null && typeof cell === 'object' && !Array.isArray(cell)) {
          if (cell.node instanceof Node) td.appendChild(cell.node);
          else if (cell.text != null) td.textContent = String(cell.text);
          if (cell.class) td.className = cell.class;
          if (cell.title) td.setAttribute('title', cell.title);
          if (cell.colspan) td.setAttribute('colspan', cell.colspan);
          if (cell.rowspan) td.setAttribute('rowspan', cell.rowspan);
        } else {
          td.textContent = cell == null ? '' : String(cell);
        }
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
    tbl.appendChild(tbody);
    wrap.appendChild(tbl);
    return wrap;
  }

  function tabs(items, initialKey) {
    const wrap = el('div');
    const bar = el('div', { class: 'tabs', role: 'tablist' });
    const panel = el('div', { class: 'tab-panel' });
    const buttons = {};
    let active = initialKey || (items[0] && items[0].key);

    function renderPanel() {
      const item = items.find(function (it) { return it.key === active; });
      panel.innerHTML = '';
      if (item && item.render) {
        const content = item.render();
        if (Array.isArray(content)) appendChildren(panel, content);
        else if (content instanceof Node) panel.appendChild(content);
      }
    }

    function select(key) {
      active = key;
      Object.keys(buttons).forEach(function (k) {
        buttons[k].classList.toggle('active', k === key);
      });
      renderPanel();
    }

    items.forEach(function (item) {
      const btn = el('button', {
        type: 'button',
        class: 'tab-btn' + (item.key === active ? ' active' : ''),
        text: item.label,
        role: 'tab'
      });
      btn.addEventListener('click', function () { select(item.key); });
      buttons[item.key] = btn;
      bar.appendChild(btn);
    });

    wrap.appendChild(bar);
    wrap.appendChild(panel);
    renderPanel();
    return { node: wrap, select: select, panel: panel, get active() { return active; } };
  }

  function chart(title, height) {
    const holder = el('div', { class: 'chart-holder' });
    const cardEl = el('section', { class: 'card chart-card' });
    if (title) cardEl.appendChild(el('div', { class: 'card-title', text: title }));
    cardEl.appendChild(holder);
    let handle = null;
    return {
      node: cardEl,
      holder: holder,
      draw: function (cfg) {
        cfg = cfg || {};
        cfg.height = cfg.height || height || 360;
        handle = PPS.chart.create(holder, cfg);
        return handle;
      },
      update: function (cfg) {
        if (!handle) return null;
        cfg = cfg || {};
        return handle.update(cfg);
      },
      get handle() { return handle; }
    };
  }

  const SPEEDS = [
    { value: '0.5', label: '0.5×' },
    { value: '1', label: '1×' },
    { value: '2', label: '2×' }
  ];

  /**
   * Animation control bar.
   * opts: { count, delay, label, onStep(i), onDone(), onReset() }
   */
  function animBar(opts) {
    opts = opts || {};
    const count = Math.max(1, Math.floor(opts.count || 1));
    const wrap = el('div', { class: 'anim-bar' });

    const playBtn = button('▶ ' + t('common.animation'), null, { variant: 'ghost', small: true });
    const nextBtn = button('⏭ ' + t('common.step'), null, { variant: 'ghost', small: true });
    const resetBtn = button('↺ ' + t('common.reset'), null, { variant: 'ghost', small: true });
    const speedSel = select(SPEEDS, '1');
    speedSel.classList.add('select-sm');
    speedSel.setAttribute('aria-label', t('common.speed'));

    const status = el('span', { class: 'anim-status' });
    const track = el('span', { class: 'anim-track' });
    const fill = el('span', { class: 'anim-fill' });
    track.appendChild(fill);

    const controller = PPS.anim.create({
      count: count,
      delay: opts.delay || 620,
      onStep: function (index, state) {
        if (opts.onStep) opts.onStep(index, state);
      },
      onDone: function () {
        if (opts.onDone) opts.onDone();
      },
      onReset: function () {
        if (opts.onReset) opts.onReset();
      },
      onState: function (state) {
        playBtn.textContent = state.playing
          ? '⏸ ' + t('common.pause')
          : (state.done ? '↺ ' + t('common.replay') : '▶ ' + t('common.animation'));
        playBtn.classList.toggle('active', state.playing);
        nextBtn.disabled = state.playing;
        const progress = state.index < 0 ? 0 : (state.index + 1) / state.count;
        fill.style.width = Math.round(progress * 100) + '%';
        status.textContent = state.index < 0
          ? t('common.stepsCount', { n: state.count })
          : t('common.stepOf', { i: state.index + 1, n: state.count });
      }
    });

    playBtn.addEventListener('click', function () { controller.toggle(); });
    nextBtn.addEventListener('click', function () { controller.next(); });
    resetBtn.addEventListener('click', function () { controller.reset(); });
    speedSel.addEventListener('change', function () { controller.setSpeed(Number(speedSel.value)); });

    wrap.appendChild(playBtn);
    wrap.appendChild(nextBtn);
    wrap.appendChild(resetBtn);
    if (count > 2) wrap.appendChild(speedSel);
    wrap.appendChild(track);
    wrap.appendChild(status);

    return { node: wrap, controller: controller };
  }

  function clearAnim(tableEl) {
    if (!tableEl) return;
    tableEl.querySelectorAll('.anim-target, .anim-cell, .anim-flash').forEach(function (node) {
      node.classList.remove('anim-target', 'anim-cell', 'anim-flash');
    });
  }

  function highlightRow(tableEl, index) {
    clearAnim(tableEl);
    if (!tableEl || index == null || index < 0) return;
    const row = tableEl.querySelectorAll('tbody tr')[index];
    if (row) row.classList.add('anim-target');
  }

  function flashCells(cells) {
    (cells || []).forEach(function (cell) {
      if (cell) cell.classList.add('anim-cell');
    });
  }

  function readFileAsText(file) {
    return new Promise(function (resolve, reject) {
      const reader = new FileReader();
      reader.onload = function (e) {
        resolve(String(e.target.result || ''));
      };
      reader.onerror = function () { reject(reader.error); };
      reader.readAsText(file, 'utf-8');
    });
  }

  function readFileAsArrayBuffer(file) {
    return new Promise(function (resolve, reject) {
      const reader = new FileReader();
      reader.onload = function (e) {
        resolve(e.target.result);
      };
      reader.onerror = function () { reject(reader.error); };
      reader.readAsArrayBuffer(file);
    });
  }

  function importDataButton(opts) {
    opts = opts || {};
    const input = el('input', { type: 'file', accept: opts.accept || '.csv,.txt,.tsv,.xlsx,.xls' });
    input.style.display = 'none';
    const btn = button('📁 ' + t('common.importData'), function () {
      input.click();
    }, { variant: 'ghost' });
    btn.style.order = '99';
    input.addEventListener('change', function () {
      const file = input.files && input.files[0];
      if (!file) return;
      const name = (file.name || '').toLowerCase();
      const isXlsx = name.indexOf('.xlsx') !== -1 || name.indexOf('.xls') !== -1;
      if (isXlsx && (typeof window.XLSX !== 'undefined' || typeof XLSX !== 'undefined')) {
        readFileAsArrayBuffer(file).then(function (buf) {
          try {
            const X = window.XLSX || XLSX;
            const wb = X.read(buf, { type: 'array' });
            const sheetName = wb.SheetNames && wb.SheetNames[0];
            const sheet = wb.Sheets[sheetName];
            const csv = X.utils.sheet_to_csv(sheet);
            const parsed = PPS.format.parseCsvText(csv);
            const xy = PPS.format.rowsToXY(parsed.rows);
            if (xy.error) throw new Error(xy.error);
            if (opts.onParsed) opts.onParsed({ xs: xy.xs, ys: xy.ys, meta: { file: file.name, sheet: sheetName } });
          } catch (e) {
            if (opts.onError) opts.onError(e); else alert(t('common.importError'));
          }
        }).catch(function (e) {
          if (opts.onError) opts.onError(e); else alert(t('common.importError'));
        });
      } else {
        readFileAsText(file).then(function (text) {
          try {
            const parsed = PPS.format.parseCsvText(text);
            const xy = PPS.format.rowsToXY(parsed.rows);
            if (xy.error) throw new Error(xy.error);
            if (opts.onParsed) opts.onParsed({ xs: xy.xs, ys: xy.ys, meta: { file: file.name } });
          } catch (e) {
            if (opts.onError) opts.onError(e); else alert(t('common.importError'));
          }
        }).catch(function (e) {
          if (opts.onError) opts.onError(e); else alert(t('common.importError'));
        });
      }
      input.value = '';
    });
    return { node: btn, input: input };
  }

  function parseList(str) {
    const result = PPS.format.parseList(str);
    if (result.error) throw new Error('parse');
    if (!result.values.length) throw new Error('empty');
    return result.values;
  }

  function parseNumber(str, code) {
    const value = PPS.format.parseNumber(str);
    if (!isFinite(value)) throw new Error(code || 'numeric');
    return value;
  }

  function xyControl(opts) {
    opts = opts || {};
    let currentMode = opts.initialMode || 'horizontal';
    const xs = opts.xs || [];
    const ys = opts.ys || [];
    const xLabel = opts.xLabel || t('lagrange.xs') || 'Mốc x';
    const yLabel = opts.yLabel || t('lagrange.ys') || 'Giá trị y';
    const xHint = opts.xHint || t('lagrange.xsHint') || '';
    const yHint = opts.yHint || t('lagrange.ysHint') || '';
    const onChange = opts.onChange;

    const xsIn = input({ value: Array.isArray(xs) ? xs.join(', ') : String(xs) });
    const ysIn = input({ value: Array.isArray(ys) ? ys.join(', ') : String(ys) });

    function formatVerticalText(arrX, arrY) {
      const lines = [];
      const n = Math.max(arrX.length, arrY.length);
      for (let i = 0; i < n; i++) {
        const xStr = i < arrX.length ? String(arrX[i]) : '';
        const yStr = i < arrY.length ? String(arrY[i]) : '';
        if (xStr !== '' || yStr !== '') {
          lines.push(xStr.padEnd(8, ' ') + '\t' + yStr);
        }
      }
      return lines.join('\n');
    }

    function parseCellNumber(str) {
      if (str == null) return NaN;
      const s = String(str).trim().replace(/−/g, '-').replace(/,/g, '.');
      if (!s) return NaN;
      return Number(s);
    }

    function splitLineXY(line) {
      line = line.trim();
      if (!line || line.indexOf('#') === 0) return null;

      // 1. If tab exists, split by tab
      if (line.indexOf('\t') !== -1) {
        const parts = line.split('\t').map(function (p) { return p.trim(); }).filter(Boolean);
        if (parts.length >= 2) return [parts[0], parts[1]];
      }

      // 2. If semicolon exists, split by semicolon
      if (line.indexOf(';') !== -1) {
        const parts = line.split(';').map(function (p) { return p.trim(); }).filter(Boolean);
        if (parts.length >= 2) return [parts[0], parts[1]];
      }

      // 3. If pipe '|' exists (markdown table row)
      if (line.indexOf('|') !== -1) {
        const parts = line.split('|').map(function (p) { return p.trim(); }).filter(Boolean);
        if (parts.length >= 2) return [parts[0], parts[1]];
      }

      // 4. Try whitespace split (handles "1,5 2,7", "1.5 2.7", "1,5, 2,7", "1, 2")
      const whiteParts = line.split(/\s+/).map(function (p) { return p.trim(); }).filter(Boolean);
      if (whiteParts.length >= 2) {
        let p0 = whiteParts[0];
        let p1 = whiteParts[1];
        if (p0.charAt(p0.length - 1) === ',') p0 = p0.slice(0, -1);
        const x = parseCellNumber(p0);
        const y = parseCellNumber(p1);
        if (isFinite(x) && isFinite(y)) {
          return [p0, p1];
        }
      }

      // 5. If comma exists without whitespace (e.g. "1.5,2.7" or "1,2")
      const commaCount = (line.match(/,/g) || []).length;
      if (commaCount === 1) {
        const parts = line.split(',').map(function (p) { return p.trim(); }).filter(Boolean);
        if (parts.length >= 2) {
          const x = parseCellNumber(parts[0]);
          const y = parseCellNumber(parts[1]);
          if (isFinite(x) && isFinite(y)) {
            return [parts[0], parts[1]];
          }
        }
      }

      return null;
    }

    function parseVerticalText(text) {
      if (!text) return { xs: [], ys: [] };
      const lines = text.split(/\r?\n/).map(function (l) { return l.trim(); }).filter(function (l) {
        return l && l.indexOf('#') !== 0;
      });
      const parsedX = [];
      const parsedY = [];
      for (let i = 0; i < lines.length; i++) {
        const pair = splitLineXY(lines[i]);
        if (pair) {
          const vx = parseCellNumber(pair[0]);
          const vy = parseCellNumber(pair[1]);
          if (isFinite(vx) && isFinite(vy)) {
            parsedX.push(vx);
            parsedY.push(vy);
          }
        }
      }
      return { xs: parsedX, ys: parsedY };
    }

    const initXArr = Array.isArray(xs) ? xs : (function () { try { return parseList(String(xs)); } catch (e) { return []; } })();
    const initYArr = Array.isArray(ys) ? ys : (function () { try { return parseList(String(ys)); } catch (e) { return []; } })();

    const verticalIn = textarea({
      rows: opts.rows || 6,
      mono: true,
      placeholder: '# x\ty\n1.0\t2.5\n2.0\t5.0\n3.0\t7.5',
      value: formatVerticalText(initXArr, initYArr)
    });

    const countBadge = el('span', { class: 'xy-badge-count', text: '' });
    function updateBadge(cnt) {
      countBadge.textContent = t('common.pointsCount', { count: cnt });
    }
    updateBadge(Math.min(initXArr.length, initYArr.length));

    const horizBtn = el('button', {
      class: 'xy-mode-btn' + (currentMode === 'horizontal' ? ' active' : ''),
      type: 'button',
      text: t('common.horizMode')
    });
    const vertBtn = el('button', {
      class: 'xy-mode-btn' + (currentMode === 'vertical' ? ' active' : ''),
      type: 'button',
      text: t('common.vertMode')
    });

    const horizContainer = el('div', { class: 'form-grid', style: currentMode === 'horizontal' ? '' : 'display:none;' },
      field(xLabel, xsIn, xHint),
      field(yLabel, ysIn, yHint)
    );

    const vertContainer = el('div', { class: 'xy-vertical-wrap', style: currentMode === 'vertical' ? '' : 'display:none;' },
      el('div', { class: 'xy-vertical-header' },
        el('span', { text: 'Cột 1: ' + xLabel }),
        el('span', { text: 'Cột 2: ' + yLabel })
      ),
      verticalIn,
      el('div', { class: 'xy-vertical-meta' },
        el('span', { text: t('common.vertHint') }),
        countBadge
      )
    );

    function syncFromHoriz() {
      try {
        const curX = parseList(xsIn.value);
        const curY = parseList(ysIn.value);
        verticalIn.value = formatVerticalText(curX, curY);
        updateBadge(Math.min(curX.length, curY.length));
      } catch (e) {}
    }

    function syncFromVert() {
      const res = parseVerticalText(verticalIn.value);
      xsIn.value = res.xs.join(', ');
      ysIn.value = res.ys.join(', ');
      updateBadge(res.xs.length);
    }

    horizBtn.onclick = function () {
      if (currentMode === 'horizontal') return;
      currentMode = 'horizontal';
      horizBtn.classList.add('active');
      vertBtn.classList.remove('active');
      syncFromVert();
      horizContainer.style.display = '';
      vertContainer.style.display = 'none';
      if (onChange) onChange();
    };

    vertBtn.onclick = function () {
      if (currentMode === 'vertical') return;
      currentMode = 'vertical';
      vertBtn.classList.add('active');
      horizBtn.classList.remove('active');
      syncFromHoriz();
      horizContainer.style.display = 'none';
      vertContainer.style.display = '';
      if (onChange) onChange();
    };

    xsIn.addEventListener('input', function () {
      syncFromHoriz();
      if (onChange) onChange();
    });
    ysIn.addEventListener('input', function () {
      syncFromHoriz();
      if (onChange) onChange();
    });
    verticalIn.addEventListener('input', function () {
      syncFromVert();
      if (onChange) onChange();
    });

    const modeToolbar = el('div', { class: 'xy-mode-bar' },
      el('div', { class: 'xy-mode-title', text: t('common.xyInputMode') }),
      el('div', { class: 'xy-mode-switch' }, horizBtn, vertBtn)
    );

    const wrap = el('div', { class: 'xy-control-group' + (opts.span2 ? ' span-2' : '') },
      modeToolbar,
      horizContainer,
      vertContainer
    );

    function getValues() {
      if (currentMode === 'vertical') {
        const res = parseVerticalText(verticalIn.value);
        if (!res.xs.length) throw new Error('empty');
        if (res.xs.length !== res.ys.length) throw new Error('needPoints');
        return res;
      }
      const vx = parseList(xsIn.value);
      const vy = parseList(ysIn.value);
      if (vx.length !== vy.length) throw new Error('needPoints');
      return { xs: vx, ys: vy };
    }

    function setValues(newXs, newYs) {
      const arrX = newXs || [];
      const arrY = newYs || [];
      xsIn.value = arrX.join(', ');
      ysIn.value = arrY.join(', ');
      verticalIn.value = formatVerticalText(arrX, arrY);
      updateBadge(Math.min(arrX.length, arrY.length));
    }

    return {
      node: wrap,
      xsIn: xsIn,
      ysIn: ysIn,
      verticalIn: verticalIn,
      getValues: getValues,
      setValues: setValues,
      getMode: function () { return currentMode; },
      setMode: function (m) {
        if (m === 'vertical') vertBtn.onclick();
        else horizBtn.onclick();
      }
    };
  }

  function checkGroup(options, checked) {
    const wrap = el('div', { class: 'check-group' });
    const inputs = [];
    options.forEach(function (opt) {
      const label = el('label', { class: 'check-chip' });
      const box = el('input', { type: 'checkbox', value: opt.value });
      box.checked = (checked || []).indexOf(opt.value) !== -1;
      label.appendChild(box);
      label.appendChild(el('span', { text: opt.label }));
      wrap.appendChild(label);
      inputs.push(box);
    });
    return {
      node: wrap,
      values: function () {
        return inputs.filter(function (b) { return b.checked; }).map(function (b) { return b.value; });
      },
      setValues: function (vals) {
        vals = vals || [];
        inputs.forEach(function (b) {
          b.checked = vals.indexOf(b.value) !== -1;
        });
      }
    };
  }

  function section(title, content, extra) {
    const wrap = el('div', { class: 'result-section' });
    const head = el('div', { class: 'section-title' });
    head.appendChild(el('span', { text: title }));
    if (extra) head.appendChild(el('span', { class: 'count', text: extra }));
    wrap.appendChild(head);
    appendChildren(wrap, [content]);
    return wrap;
  }

  function badge(text, kind) {
    return el('span', { class: 'badge' + (kind ? ' ' + kind : ''), text: text });
  }

  function fallbackCopy(text, cb) {
    if (typeof document === 'undefined') return;
    const ta = document.createElement('textarea');
    ta.value = text;
    ta.style.position = 'fixed';
    ta.style.top = '-9999px';
    ta.style.opacity = '0';
    document.body.appendChild(ta);
    ta.focus();
    ta.select();
    try {
      document.execCommand('copy');
      if (cb) cb();
    } catch (e) {
      console.warn('Copy failed', e);
    }
    document.body.removeChild(ta);
  }

  function copyToClipboard(text, btn) {
    if (!text) return;
    function showSuccess() {
      if (!btn) return;
      const origText = btn.getAttribute('data-orig-title') || btn.innerHTML;
      btn.setAttribute('data-orig-title', origText);
      const copiedText = (PPS.i18n && typeof PPS.i18n.t === 'function') ? PPS.i18n.t('common.copied') : 'Đã chép!';
      btn.textContent = '✓ ' + copiedText;
      btn.classList.add('copied');
      setTimeout(function () {
        btn.innerHTML = origText;
        btn.classList.remove('copied');
      }, 1600);
    }
    if (typeof navigator !== 'undefined' && navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(showSuccess).catch(function () {
        fallbackCopy(text, showSuccess);
      });
    } else {
      fallbackCopy(text, showSuccess);
    }
  }

  function copyButton(getTextFn, label) {
    const defaultLabel = (PPS.i18n && typeof PPS.i18n.t === 'function') ? PPS.i18n.t('common.copy') : 'Sao chép';
    const btnText = label || defaultLabel;
    const btn = el('button', {
      type: 'button',
      class: 'copy-btn',
      title: btnText
    });
    btn.innerHTML = '📋 ' + btnText;
    btn.addEventListener('click', function (e) {
      e.stopPropagation();
      const val = typeof getTextFn === 'function' ? getTextFn() : getTextFn;
      copyToClipboard(val, btn);
    });
    return btn;
  }

  function formula(text, wrapLines, opts) {
    opts = opts || {};
    const div = el('div', { class: 'formula' + (wrapLines ? ' wrap' : '') });
    const content = el('span', { class: 'formula-text', text: text });
    div.appendChild(content);

    if (!opts.noCopy) {
      const copyBtn = el('button', {
        type: 'button',
        class: 'copy-icon-btn',
        title: (PPS.i18n && typeof PPS.i18n.t === 'function') ? PPS.i18n.t('common.copy') : 'Sao chép'
      });
      copyBtn.setAttribute('aria-label', 'Copy formula');
      copyBtn.innerHTML = '📋';
      copyBtn.addEventListener('click', function (e) {
        e.stopPropagation();
        copyToClipboard(text, copyBtn);
      });
      div.appendChild(copyBtn);
    }
    return div;
  }

  function pageHeader(title, desc, eyebrow) {
    const head = el('div', { class: 'page-head' });
    if (eyebrow) head.appendChild(el('div', { class: 'page-eyebrow', text: eyebrow }));
    head.appendChild(el('h2', { text: title }));
    if (desc) head.appendChild(el('p', { class: 'page-desc', text: desc }));
    return head;
  }

  function resultHead(title, badgeEl) {
    const head = el('div', { class: 'result-head' });
    head.appendChild(el('h3', { text: title }));
    if (badgeEl) head.appendChild(badgeEl);
    return head;
  }

  PPS.ui = {
    el: el,
    field: field,
    input: input,
    textarea: textarea,
    select: select,
    button: button,
    card: card,
    grid: grid,
    actions: actions,
    note: note,
    errorBox: errorBox,
    errText: errText,
    kv: kv,
    table: table,
    tabs: tabs,
    chart: chart,
    animBar: animBar,
    clearAnim: clearAnim,
    highlightRow: highlightRow,
    flashCells: flashCells,
    parseList: parseList,
    parseNumber: parseNumber,
    xyControl: xyControl,
    checkGroup: checkGroup,
    section: section,
    badge: badge,
    formula: formula,
    copyToClipboard: copyToClipboard,
    copyButton: copyButton,
    pageHeader: pageHeader,
    resultHead: resultHead,
    importDataButton: importDataButton,
    readFileAsText: readFileAsText,
    readFileAsArrayBuffer: readFileAsArrayBuffer
  };
})();
