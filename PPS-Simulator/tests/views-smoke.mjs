/* ============================================================
   views-smoke.mjs — headless DOM smoke tests for all 10 views
   Usage: node tests/views-smoke.mjs
   ============================================================ */

import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { fileURLToPath } from 'node:url';
import assert from 'node:assert/strict';

const here = path.dirname(fileURLToPath(import.meta.url));
const srcDir = path.resolve(here, '..', 'src');
const viewsDir = path.join(srcDir, 'views');

// ------------------------------------------------------------
// Minimal DOM Mock
// ------------------------------------------------------------
class Node {
  constructor() {
    this.childNodes = [];
    this.parentNode = null;
  }
  appendChild(child) {
    if (child.parentNode) {
      child.parentNode.removeChild(child);
    }
    child.parentNode = this;
    this.childNodes.push(child);
    return child;
  }
  removeChild(child) {
    const idx = this.childNodes.indexOf(child);
    if (idx !== -1) {
      this.childNodes.splice(idx, 1);
      child.parentNode = null;
    }
    return child;
  }
  get children() {
    return this.childNodes.filter(n => n instanceof Element);
  }
}

class TextNode extends Node {
  constructor(text) {
    super();
    this.textContent = String(text);
  }
}

class Element extends Node {
  constructor(tagName, namespaceURI = null) {
    super();
    this.tagName = tagName.toUpperCase();
    this.namespaceURI = namespaceURI;
    this.attributes = {};
    this.style = {};
    this.dataset = {};
    this.classList = {
      _classes: new Set(),
      add: (...cls) => cls.forEach(c => this.classList._classes.add(c)),
      remove: (...cls) => cls.forEach(c => this.classList._classes.delete(c)),
      contains: (c) => this.classList._classes.has(c),
      toggle: (c) => this.classList._classes.has(c) ? this.classList._classes.delete(c) : this.classList._classes.add(c)
    };
    this._listeners = {};
    this._value = '';
    this.id = '';
  }

  get className() {
    return Array.from(this.classList._classes).join(' ');
  }
  set className(val) {
    this.classList._classes.clear();
    if (val) val.split(/\s+/).filter(Boolean).forEach(c => this.classList._classes.add(c));
  }

  get textContent() {
    return this.childNodes.map(n => n.textContent || '').join('');
  }
  set textContent(val) {
    this.childNodes = [];
    if (val) this.appendChild(new TextNode(val));
  }

  get innerHTML() {
    return this.textContent;
  }
  set innerHTML(val) {
    this.textContent = val.replace(/<[^>]*>/g, '');
  }

  get value() {
    return this._value;
  }
  set value(val) {
    this._value = String(val);
  }

  setAttribute(k, v) {
    this.attributes[k] = String(v);
    if (k === 'class') this.className = String(v);
    if (k === 'id') this.id = String(v);
    if (k === 'value') this._value = String(v);
  }
  getAttribute(k) {
    return this.attributes[k] != null ? this.attributes[k] : null;
  }
  removeAttribute(k) {
    delete this.attributes[k];
  }

  addEventListener(type, fn) {
    this._listeners[type] = this._listeners[type] || [];
    this._listeners[type].push(fn);
  }
  removeEventListener(type, fn) {
    if (!this._listeners[type]) return;
    this._listeners[type] = this._listeners[type].filter(l => l !== fn);
  }
  dispatchEvent(evt) {
    const list = this._listeners[evt.type] || [];
    list.forEach(fn => fn.call(this, evt));
  }

  click() {
    this.dispatchEvent({ type: 'click', target: this, preventDefault: () => {}, stopPropagation: () => {} });
  }

  getBoundingClientRect() {
    return { left: 0, top: 0, right: 600, bottom: 400, width: 600, height: 400 };
  }

  querySelector(sel) {
    return this.querySelectorAll(sel)[0] || null;
  }
  querySelectorAll(sel) {
    const results = [];
    function walk(node) {
      if (node instanceof Element) {
        if (sel.startsWith('.') && node.classList.contains(sel.slice(1))) results.push(node);
        else if (sel.startsWith('#') && node.id === sel.slice(1)) results.push(node);
        else if (node.tagName.toLowerCase() === sel.toLowerCase()) results.push(node);
      }
      for (const child of node.childNodes) {
        walk(child);
      }
    }
    for (const child of this.childNodes) {
      walk(child);
    }
    return results;
  }
}

const mockDoc = {
  createElement: (tag) => new Element(tag),
  createElementNS: (ns, tag) => new Element(tag, ns),
  createTextNode: (text) => new TextNode(text),
  getElementById: (id) => null,
  body: new Element('body'),
  addEventListener: () => {},
  removeEventListener: () => {}
};

// Setup global context
globalThis.window = globalThis;
globalThis.document = mockDoc;
globalThis.Node = Node;
globalThis.Element = Element;
globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 16);
globalThis.cancelAnimationFrame = (id) => clearTimeout(id);
if (typeof globalThis.navigator === 'undefined') {
  globalThis.navigator = { clipboard: { writeText: async () => {} } };
} else {
  try {
    Object.defineProperty(globalThis.navigator, 'clipboard', {
      value: { writeText: async () => {} },
      configurable: true,
      writable: true
    });
  } catch (e) {}
}

// Load engine scripts
const coreFiles = [
  'format.js', 'i18n.js', 'linalg.js', 'polynomial.js', 'horner.js',
  'chebyshev.js', 'lagrange.js', 'newton.js', 'central.js',
  'inverse.js', 'spline.js', 'leastSquares.js',
  'ui.js', 'chart.js', 'export.js'
];

for (const file of coreFiles) {
  const code = fs.readFileSync(path.join(srcDir, file), 'utf8');
  vm.runInThisContext(code, { filename: file });
}

// Load view scripts
const viewFiles = [
  { id: 'horner', file: 'horner.js' },
  { id: 'chebyshev', file: 'chebyshev.js' },
  { id: 'lagrange', file: 'lagrange.js' },
  { id: 'newton-div', file: 'newton-divided.js' },
  { id: 'newton-eq', file: 'newton-equal.js' },
  { id: 'central', file: 'central.js' },
  { id: 'inverse', file: 'inverse.js' },
  { id: 'spline', file: 'spline.js' },
  { id: 'least-squares', file: 'least-squares.js' },
  { id: 'comparison', file: 'comparison.js' }
];

for (const v of viewFiles) {
  const code = fs.readFileSync(path.join(viewsDir, v.file), 'utf8');
  vm.runInThisContext(code, { filename: v.file });
}

const PPS = globalThis.PPS;
let passed = 0;
let failed = 0;

function test(name, fn) {
  try {
    fn();
    passed++;
    console.log('  ✓ ' + name);
  } catch (error) {
    failed++;
    console.error('  ✗ ' + name);
    console.error('    ' + (error && error.stack ? error.stack.split('\n').slice(0, 4).join('\n    ') : error));
  }
}

console.log('\n====================================================');
console.log('Running UI View Smoke Tests (Headless DOM)');
console.log('====================================================\n');

for (const v of viewFiles) {
  test(`View: ${v.id} (${v.file}) renders and calculates correctly`, () => {
    const view = PPS.views[v.id];
    assert.ok(view, `PPS.views["${v.id}"] must be defined`);
    assert.strictEqual(typeof view.render, 'function', `PPS.views["${v.id}"].render must be a function`);

    const root = mockDoc.createElement('div');
    view.render(root);

    // Root must have content
    assert.ok(root.childNodes.length > 0, `Root element must contain rendered children`);

    // Check if any error boxes were generated during default preset run
    const errorBoxes = root.querySelectorAll('.card-error');
    assert.strictEqual(errorBoxes.length, 0, `View rendered with unexpected error boxes: ${errorBoxes.map(b => b.textContent).join('; ')}`);

    // Find buttons and click Example / Run if available
    const buttons = root.querySelectorAll('button');
    if (buttons.length > 0) {
      // Click each button to ensure event handlers don't crash
      buttons.forEach(btn => {
        btn.click();
      });

      // Re-verify no error boxes appeared after button triggers
      const postClickErrors = root.querySelectorAll('.card-error');
      assert.strictEqual(postClickErrors.length, 0, `Button click produced errors: ${postClickErrors.map(b => b.textContent).join('; ')}`);
    } else {
      // Static reference views (e.g. comparison table) must render tables or cards
      const tables = root.querySelectorAll('table');
      assert.ok(tables.length >= 1, `Static view must render at least one table`);
    }
  });
}

test('UI component: xyControl syncs between horizontal and vertical modes', () => {
  const ctrl = PPS.ui.xyControl({
    xs: [1, 2, 3],
    ys: [10, 20, 30]
  });

  // Default is horizontal
  assert.equal(ctrl.getMode(), 'horizontal');
  let vals = ctrl.getValues();
  assert.deepEqual(vals.xs, [1, 2, 3]);
  assert.deepEqual(vals.ys, [10, 20, 30]);

  // Switch to vertical mode
  ctrl.setMode('vertical');
  assert.equal(ctrl.getMode(), 'vertical');
  vals = ctrl.getValues();
  assert.deepEqual(vals.xs, [1, 2, 3]);
  assert.deepEqual(vals.ys, [10, 20, 30]);

  // Edit in vertical mode (paste vertical data)
  ctrl.verticalIn.value = '4.0\t40.0\n5.0\t50.0\n6.0\t60.0';
  ctrl.verticalIn.dispatchEvent({ type: 'input' });
  vals = ctrl.getValues();
  assert.deepEqual(vals.xs, [4, 5, 6]);
  assert.deepEqual(vals.ys, [40, 50, 60]);

  // Switch back to horizontal mode - inputs should be synchronized
  ctrl.setMode('horizontal');
  assert.equal(ctrl.xsIn.value, '4, 5, 6');
  assert.equal(ctrl.ysIn.value, '40, 50, 60');

  // Programmatic setValues
  ctrl.setValues([7, 8], [70, 80]);
  vals = ctrl.getValues();
  assert.deepEqual(vals.xs, [7, 8]);
  assert.deepEqual(vals.ys, [70, 80]);

  // Test equivalence of comma ',' and dot '.' decimals in vertical mode
  ctrl.setMode('vertical');
  // 1. Whitespace separated with comma decimals
  ctrl.verticalIn.value = '1,5   2,7\n3,0   4,25';
  ctrl.verticalIn.dispatchEvent({ type: 'input' });
  vals = ctrl.getValues();
  assert.deepEqual(vals.xs, [1.5, 3.0]);
  assert.deepEqual(vals.ys, [2.7, 4.25]);

  // 2. Tab separated with comma decimals (Excel style)
  ctrl.verticalIn.value = '1,5\t2,7\n3,0\t4,25';
  ctrl.verticalIn.dispatchEvent({ type: 'input' });
  vals = ctrl.getValues();
  assert.deepEqual(vals.xs, [1.5, 3.0]);
  assert.deepEqual(vals.ys, [2.7, 4.25]);

  // 3. Dot decimals
  ctrl.verticalIn.value = '1.5   2.7\n3.0   4.25';
  ctrl.verticalIn.dispatchEvent({ type: 'input' });
  const valsDot = ctrl.getValues();
  assert.deepEqual(valsDot.xs, vals.xs);
  assert.deepEqual(valsDot.ys, vals.ys);

  // 4. Negative numbers with comma decimals
  ctrl.verticalIn.value = '-1,25   0,5\n−2,75   -3,125';
  ctrl.verticalIn.dispatchEvent({ type: 'input' });
  vals = ctrl.getValues();
  assert.deepEqual(vals.xs, [-1.25, -2.75]);
  assert.deepEqual(vals.ys, [0.5, -3.125]);
});

console.log('\n────────────────────────────────────────────────────');
console.log(`  ${passed} passed, ${failed} failed`);
console.log('────────────────────────────────────────────────────\n');

if (failed > 0) {
  process.exit(1);
}
