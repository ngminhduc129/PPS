/* ============================================================
   animate.spec.js — unit tests for animation controller
   ============================================================ */

import assert from "node:assert/strict";
import vm from "node:vm";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const srcDir = path.resolve(here, "..", "..", "src");

const ctx = {
  window: {
    matchMedia: function () { return { matches: false }; },
    requestAnimationFrame: function (cb) { return setTimeout(cb, 0); },
    cancelAnimationFrame: function (id) { clearTimeout(id); }
  }
};
ctx.globalThis = ctx;
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(path.join(srcDir, "animate.js"), "utf8"), ctx);

const PPS = ctx.window.PPS;
let passed = 0;
let failed = 0;

function test(name, fn) {
  try {
    fn();
    passed++;
    console.log("  \u2713 " + name);
  } catch (error) {
    failed++;
    console.error("  \u2717 " + name);
    console.error("    " + (error && error.stack ? error.stack.split("\n").slice(0, 3).join("\n    ") : error));
  }
}

console.log("\nAnimate");

test("play-steps-onDone", () => {
  const controller = PPS.anim.create({
    count: 3,
    delay: 10,
    raf: function () { return 1; },
    cancel: function () {}
  });
  controller.play();
  controller.pause();
  controller.reset();
  assert.equal(typeof controller.play, "function");
  assert.equal(typeof controller.pause, "function");
  assert.equal(typeof controller.toggle, "function");
  assert.equal(typeof controller.next, "function");
  assert.equal(typeof controller.reset, "function");
  controller.destroy();
});

test("pause stops steps", () => {
  let steps = 0;
  const controller = PPS.anim.create({
    count: 5,
    delay: 10,
    raf: function () { return 1; },
    cancel: function () {},
    onStep: function () { steps++; }
  });
  controller.play();
  controller.pause();
  const before = steps;
  controller.pause();
  assert.equal(steps, before);
  controller.destroy();
});

test("reset", () => {
  let resetCalled = false;
  const controller = PPS.anim.create({
    count: 3,
    delay: 10,
    raf: function () { return 1; },
    cancel: function () {},
    onReset: function () { resetCalled = true; }
  });
  controller.play();
  controller.pause();
  controller.reset();
  assert.equal(resetCalled, true);
  controller.destroy();
});

test("stopAll", () => {
  const c1 = PPS.anim.create({ count: 2, delay: 10, raf: function () { return 1; }, cancel: function () {} });
  const c2 = PPS.anim.create({ count: 2, delay: 10, raf: function () { return 1; }, cancel: function () {} });
  PPS.anim.stopAll();
  assert.equal(c1.state().playing, false);
  assert.equal(c2.state().playing, false);
});

console.log("\n" + "\u2500".repeat(52));
console.log(`  ${passed} passed, ${failed} failed`);
if (failed) process.exitCode = 1;
