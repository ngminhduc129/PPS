/* ============================================================
   animate.js — tiny step animation controller (rAF based)
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});

  function reduceMotion() {
    return !!(window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches);
  }

  const registry = [];

  /**
   * opts:
   *   count      number of steps
   *   delay      milliseconds between steps (default 700)
   *   onStep(i)  invoked for every step (0-based)
   *   onDone()   invoked when the last step has been shown
   *   onReset()  invoked when the sequence is rewound
   *   onState(s) invoked whenever index / playing / speed changes
   *   raf/cancel optional scheduler injection (for tests)
   */
  function create(opts) {
    opts = opts || {};
    const count = Math.max(0, Math.floor(opts.count || 0));
    const baseDelay = opts.delay || 700;
    const rafFn = opts.raf || function (cb) { return window.requestAnimationFrame(cb); };
    const cancelFn = opts.cancel || function (id) { window.cancelAnimationFrame(id); };

    let index = -1;
    let playing = false;
    let speed = 1;
    let acc = 0;
    let last = null;
    let raf = null;

    function state() {
      return {
        index: index,
        count: count,
        playing: playing,
        speed: speed,
        done: count > 0 && index >= count - 1
      };
    }

    function notify() {
      if (opts.onState) opts.onState(state());
    }

    function emit() {
      if (opts.onStep) opts.onStep(index, state());
    }

    function stopFrame() {
      if (raf != null) { cancelFn(raf); raf = null; }
      last = null;
      acc = 0;
    }

    function finish() {
      playing = false;
      stopFrame();
      notify();
      if (opts.onDone) opts.onDone(state());
    }

    function tick(ts) {
      if (!playing) return;
      if (last == null) last = ts;
      acc += Math.max(0, ts - last);
      last = ts;
      const interval = baseDelay / speed;
      while (acc >= interval && index < count - 1) {
        acc -= interval;
        index++;
        emit();
      }
      if (index >= count - 1) {
        finish();
        return;
      }
      raf = rafFn(tick);
    }

    function play() {
      if (!count) return;
      if (index >= count - 1) doReset();
      if (reduceMotion()) {
        while (index < count - 1) { index++; emit(); }
        finish();
        return;
      }
      playing = true;
      stopFrame();
      index++;
      emit();
      notify();
      if (index >= count - 1) { finish(); return; }
      raf = rafFn(tick);
    }

    function pause() {
      if (!playing) return;
      playing = false;
      stopFrame();
      notify();
    }

    function toggle() {
      if (playing) pause();
      else play();
    }

    function next() {
      if (!count) return;
      if (index >= count - 1) doReset();
      pause();
      index++;
      emit();
      if (index >= count - 1) finish();
      else notify();
    }

    function doReset() {
      index = -1;
      acc = 0;
      if (opts.onReset) opts.onReset();
    }

    function reset() {
      pause();
      doReset();
      notify();
    }

    function setSpeed(value) {
      speed = value === 0.5 || value === 2 ? value : 1;
      notify();
    }

    function destroy() {
      pause();
      const position = registry.indexOf(controller);
      if (position !== -1) registry.splice(position, 1);
    }

    const controller = {
      play: play,
      pause: pause,
      toggle: toggle,
      next: next,
      reset: reset,
      setSpeed: setSpeed,
      destroy: destroy,
      state: state,
      get running() { return playing; }
    };

    registry.push(controller);
    return controller;
  }

  function stopAll() {
    registry.slice().forEach(function (controller) { controller.destroy(); });
  }

  PPS.anim = {
    create: create,
    stopAll: stopAll,
    reduceMotion: reduceMotion
  };
})();
