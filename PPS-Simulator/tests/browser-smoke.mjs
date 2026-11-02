/* ============================================================
   browser-smoke.mjs — headless Edge smoke test over CDP.
   Usage: node tests/browser-smoke.mjs
   ============================================================ */

import { spawn } from 'node:child_process';
import { pathToFileURL } from 'node:url';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const indexUrl = pathToFileURL(path.resolve(here, '..', 'index.html')).href;
const EDGE = 'C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe';
const PORT = 9333;
const PROFILE = 'C:\\Users\\minhd\\AppData\\Local\\Temp\\opencode\\pps-edge-smoke';

const edge = spawn(EDGE, [
  '--headless',
  '--disable-gpu',
  '--no-first-run',
  '--disable-extensions',
  '--allow-file-access-from-files',
  '--remote-debugging-port=' + PORT,
  '--user-data-dir=' + PROFILE,
  'about:blank'
], { stdio: 'ignore' });

function sleep(ms) { return new Promise(resolve => setTimeout(resolve, ms)); }

async function waitForDevtools() {
  for (let i = 0; i < 50; i++) {
    try {
      const res = await fetch('http://127.0.0.1:' + PORT + '/json/version');
      if (res.ok) return;
    } catch (e) { /* retry */ }
    await sleep(200);
  }
  throw new Error('DevTools endpoint did not start');
}

let ws;
let seq = 0;
const pending = new Map();
const pageErrors = [];
let loadResolver = null;

function send(method, params) {
  return new Promise((resolve, reject) => {
    const id = ++seq;
    pending.set(id, { resolve, reject });
    ws.send(JSON.stringify({ id, method, params: params || {} }));
  });
}

async function evalJS(expression) {
  const result = await send('Runtime.evaluate', {
    expression,
    returnByValue: true,
    awaitPromise: true
  });
  if (result.exceptionDetails) {
    throw new Error('Evaluation failed: ' + (result.exceptionDetails.text || '') + ' ' +
      (result.exceptionDetails.exception && result.exceptionDetails.exception.description || ''));
  }
  return result.result.value;
}

async function waitLoad() {
  await new Promise((resolve, reject) => {
    const timer = setTimeout(() => reject(new Error('Page load timeout')), 8000);
    loadResolver = () => { clearTimeout(timer); resolve(); };
  });
}

async function main() {
  await waitForDevtools();
  const list = await (await fetch('http://127.0.0.1:' + PORT + '/json/list')).json();
  const page = list.find(t => t.type === 'page');
  if (!page) throw new Error('No page target');

  ws = new WebSocket(page.webSocketDebuggerUrl);
  await new Promise((resolve, reject) => {
    ws.onopen = resolve;
    ws.onerror = reject;
  });

  ws.onmessage = (event) => {
    const msg = JSON.parse(event.data);
    if (msg.id && pending.has(msg.id)) {
      const p = pending.get(msg.id);
      pending.delete(msg.id);
      if (msg.error) p.reject(new Error(msg.error.message));
      else p.resolve(msg.result);
      return;
    }
    if (msg.method === 'Runtime.exceptionThrown') {
      const details = msg.params && msg.params.exceptionDetails;
      pageErrors.push('EXCEPTION: ' + (details && details.exception && details.exception.description || details.text));
    }
    if (msg.method === 'Runtime.consoleAPICalled' && msg.params.type === 'error') {
      pageErrors.push('CONSOLE: ' + msg.params.args.map(a => a.value || a.description).join(' '));
    }
    if (msg.method === 'Page.loadEventFired' && loadResolver) {
      const resolve = loadResolver;
      loadResolver = null;
      resolve();
    }
  };

  await send('Page.enable');
  await send('Runtime.enable');

  const loadPromise = waitLoad();
  await send('Page.navigate', { url: indexUrl });
  await loadPromise;
  await sleep(400);

  await evalJS(`
    window.__errs = [];
    window.addEventListener('error', e => window.__errs.push(String(e.message)));
    window.addEventListener('unhandledrejection', e => window.__errs.push(String(e.reason)));
    true
  `);

  const ids = ['horner', 'lagrange', 'newton-div', 'newton-eq', 'central', 'inverse', 'spline', 'least-squares'];
  const results = [];
  let failures = 0;

  for (const id of ids) {
    await evalJS(`location.hash = '${id}'`);
    await sleep(500);
    let info = await evalJS(`(() => {
      const root = document.querySelector('#view-root');
      const errs = document.querySelectorAll('.error-box').length;
      const svgs = document.querySelectorAll('.chart-holder svg').length;
      const tabs = document.querySelectorAll('.tab-btn').length;
      return { children: root ? root.children.length : 0, errorBoxes: errs, svgs, tabs,
               text: root ? root.innerText.slice(0, 90) : '' };
    })()`);

    // click the sample-data button
    await evalJS(`(() => {
      const btns = Array.from(document.querySelectorAll('.actions button'));
      const sample = btns.find(b => b.textContent.includes('Ví dụ') || b.textContent.includes('example'));
      if (sample) sample.click();
      return true;
    })()`);
    await sleep(650);

    // exercise tabs
    const tabCount = await evalJS(`document.querySelectorAll('.tab-btn').length`);
    for (let i = 0; i < tabCount; i++) {
      await evalJS(`(() => { const t = document.querySelectorAll('.tab-btn')[${i}]; if (t) t.click(); return true; })()`);
      await sleep(260);
    }
    await sleep(400);

    // test animation: click play button, wait for completion
    await evalJS(`(() => {
      const btns = Array.from(document.querySelectorAll('.anim-bar button'));
      const play = btns.find(b => b.textContent.includes('Chạy') || b.textContent.includes('Play') || b.textContent.includes('animation') || b.textContent.includes('▶'));
      if (play) { play.click(); return 'clicked'; }
      return 'none';
    })()`);
    await sleep(2500);
    const animAfter = await evalJS(`(() => ({
      animTargets: document.querySelectorAll('.anim-target').length,
      animCells: document.querySelectorAll('.anim-cell').length,
      errorBoxes: document.querySelectorAll('.error-box').length,
      windowErrors: window.__errs.length
    }))()`);

    // test pause mid-run
    await evalJS(`(() => {
      const btns = Array.from(document.querySelectorAll('.anim-bar button'));
      const reset = btns.find(b => b.textContent.includes('Đặt') || b.textContent.includes('Reset'));
      if (reset) { reset.click(); return true; }
      return false;
    })()`);
    await sleep(200);
    await evalJS(`(() => {
      const btns = Array.from(document.querySelectorAll('.anim-bar button'));
      const play = btns.find(b => b.textContent.includes('Chạy') || b.textContent.includes('Play') || b.textContent.includes('animation') || b.textContent.includes('▶'));
      if (play) { play.click(); return true; }
      return false;
    })()`);
    await sleep(400);
    await evalJS(`(() => {
      const btns = Array.from(document.querySelectorAll('.anim-bar button'));
      const pause = btns.find(b => b.textContent.includes('Tạm') || b.textContent.includes('Pause') || b.textContent.includes('⏸'));
      if (pause) { pause.click(); return true; }
      return false;
    })()`);
    await sleep(300);

    const after = await evalJS(`(() => ({
      errorBoxes: document.querySelectorAll('.error-box').length,
      svgs: document.querySelectorAll('.chart-holder svg').length,
      windowErrors: window.__errs.length,
      text: (document.querySelector('#view-root') || {}).innerText ? document.querySelector('#view-root').innerText.slice(0, 140) : ''
    }))()`);

    const bad = info.children === 0 || after.errorBoxes > 0 || after.windowErrors > 0 || info.svgs === 0 || animAfter.errorBoxes > 0 || animAfter.animTargets > 0;
    if (bad) failures++;
    results.push({ id, ...info, after, bad, animAfter });
  }

  // language switch
  await evalJS(`location.hash = 'horner'`);
  await sleep(350);
  await evalJS(`document.getElementById('lang-en').click()`);
  await sleep(450);
  const langCheck = await evalJS(`(() => ({
    heading: (document.querySelector('.page-head h2') || {}).textContent,
    nav: (document.querySelector('.nav-item-title') || {}).textContent,
    card: (document.querySelector('#view-root .card-title') || {}).textContent,
    errors: document.querySelectorAll('.error-box').length,
    windowErrors: window.__errs.length
  }))()`);
  if (langCheck.errors || langCheck.windowErrors || langCheck.heading !== 'Horner scheme' || langCheck.card !== 'Input data') failures++;

  // English rendering of every view
  const enResults = [];
  for (const id of ids) {
    await evalJS(`location.hash = '${id}'`);
    await sleep(320);
    const info = await evalJS(`(() => {
      const text = (document.querySelector('#view-root') || {}).innerText || '';
      return {
        errors: document.querySelectorAll('.error-box').length,
        text: text.slice(0, 260),
        bad: /undefined|NaN|common\\.|alg\\.|errors\\./.test(text) ||
             /(Dữ liệu|Chạy|Đồ thị|Nội suy|Giá trị|Hệ số|Bậc|Mốc|Sai số|Công thức|Điểm|Nghiệm|Phương trình)/.test(text)
      };
    })()`);
    if (info.bad) failures++;
    enResults.push({ id, ...info });
  }

  await evalJS(`document.getElementById('lang-vi').click()`);
  await sleep(350);

  // theme switch
  await evalJS(`location.hash = 'horner'`);
  await sleep(250);
  const themeBefore = await evalJS(`document.documentElement.getAttribute('data-theme')`);
  await evalJS(`document.getElementById('theme-toggle').click()`);
  const themeAfter = await evalJS(`document.documentElement.getAttribute('data-theme')`);
  if (themeBefore === themeAfter) failures++;

  console.log('\nBrowser smoke results');
  for (const r of results) {
    const status = r.bad ? '✗' : '✓';
    console.log(`  ${status} ${r.id.padEnd(14)} cards=${r.children} svgs=${r.after.svgs} tabs=${r.after.svgs ? r.tabs : r.tabs} errors=${r.after.errorBoxes + r.after.windowErrors}`);
    if (r.bad) {
      console.log('      before: ' + JSON.stringify(r).slice(0, 300));
      console.log('      text: ' + r.after.text.replace(/\n/g, ' | ').slice(0, 240));
    }
  }
  console.log('  language:', JSON.stringify(langCheck));
  for (const r of enResults) {
    console.log(`  ${r.bad ? '✗' : '✓'} EN ${r.id.padEnd(14)} errors=${r.errors}`);
    if (r.bad) console.log('      ' + r.text.replace(/\n/g, ' | ').slice(0, 260));
  }
  console.log('  theme:', themeBefore, '->', themeAfter);
  if (pageErrors.length) {
    console.log('\n  Page exceptions:');
    pageErrors.slice(0, 10).forEach(e => console.log('   ' + e));
  }
  console.log('\n  ' + (failures ? failures + ' failure(s)' : 'all UI checks passed'));
  if (failures) process.exitCode = 1;

  ws.close();
  edge.kill();
  setTimeout(() => process.exit(process.exitCode || 0), 300);
}

main().catch(error => {
  console.error(error);
  edge.kill();
  process.exit(1);
});
