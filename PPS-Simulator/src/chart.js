/* ============================================================
   chart.js — dependency-free SVG plotting with animatable layers
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});
  const { fmt } = PPS.format;
  const NS = 'http://www.w3.org/2000/svg';

  function svgEl(tag, attrs) {
    const node = document.createElementNS(NS, tag);
    if (attrs) {
      Object.keys(attrs).forEach(function (key) {
        node.setAttribute(key, attrs[key]);
      });
    }
    return node;
  }

  function niceNum(range, round) {
    const exp = Math.floor(Math.log10(range));
    const frac = range / Math.pow(10, exp);
    let nice;
    if (round) {
      if (frac < 1.5) nice = 1;
      else if (frac < 3) nice = 2;
      else if (frac < 7) nice = 5;
      else nice = 10;
    } else {
      if (frac <= 1) nice = 1;
      else if (frac <= 2) nice = 2;
      else if (frac <= 5) nice = 5;
      else nice = 10;
    }
    return nice * Math.pow(10, exp);
  }

  function ticks(min, max, count) {
    if (!isFinite(min) || !isFinite(max)) return [0, 1];
    if (min === max) {
      const delta = Math.abs(min) > 1 ? Math.abs(min) * 0.2 : 1;
      min -= delta; max += delta;
    }
    const range = niceNum(max - min, false);
    const step = niceNum(range / Math.max(1, count - 1), true);
    const start = Math.floor(min / step) * step;
    const end = Math.ceil(max / step) * step;
    const out = [];
    for (let v = start; v <= end + step * 1e-9; v += step) {
      out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
      if (out.length > 200) break;
    }
    return out;
  }

  function pathData(points, X, Y, reveal) {
    const pts = (points || []).filter(function (p) {
      return isFinite(p.x) && isFinite(p.y);
    }).sort(function (a, b) { return a.x - b.x; });
    if (!pts.length) return '';

    let limit = pts.length;
    if (reveal != null && reveal < 1) {
      limit = Math.max(1, Math.ceil(pts.length * Math.max(0, reveal)));
    }

    let d = '';
    let penDown = false;
    for (let i = 0; i < limit; i++) {
      const p = pts[i];
      if (!isFinite(p.y)) { penDown = false; continue; }
      d += (penDown ? 'L' : 'M') + X(p.x).toFixed(2) + ' ' + Y(p.y).toFixed(2) + ' ';
      penDown = true;
    }
    return d.trim();
  }

  function create(holder, cfg) {
    cfg = cfg || {};
    holder.innerHTML = '';

    const width = Math.max(280, Math.round(holder.clientWidth || cfg.width || 760));
    const height = cfg.height || 360;
    const compact = width < 470;

    // ---- collect ranges (domain hints never rendered) ----
    const xsAll = [];
    const ysAll = [];
    function addPoint(x, y) {
      if (isFinite(x) && isFinite(y)) { xsAll.push(x); ysAll.push(y); }
    }
    function addList(list) {
      (list || []).forEach(function (p) { addPoint(p.x, p.y); });
    }
    addList(cfg.domain);
    addList(cfg.points);
    addList(cfg.highlightPoints);
    (cfg.series || []).forEach(function (s) { addList(s.points); });
    (cfg.residuals || []).forEach(function (r) {
      addPoint(r.x, r.y0);
      addPoint(r.x, r.y1);
    });
    if (cfg.hline) addPoint(cfg.hline.x != null ? cfg.hline.x : 0, cfg.hline.y);
    if (cfg.vline) addPoint(cfg.vline.x, cfg.vline.y != null ? cfg.vline.y : 0);
    if (cfg.marker) addPoint(cfg.marker.x, cfg.marker.y);

    const svg = svgEl('svg', {
      viewBox: '0 0 ' + width + ' ' + height,
      preserveAspectRatio: 'xMidYMid meet',
      role: 'img'
    });
    svg.style.width = '100%';
    svg.style.height = 'auto';

    if (!xsAll.length) {
      const empty = svgEl('text', {
        x: width / 2, y: height / 2, 'text-anchor': 'middle',
        fill: 'var(--text-dim)', 'font-size': 13
      });
      empty.textContent = cfg.emptyText || '—';
      svg.appendChild(empty);
      holder.appendChild(svg);
      return {
        svg: svg,
        update: function () { return this; },
        destroy: function () {}
      };
    }

    let xmin = cfg.xRange ? cfg.xRange[0] : Math.min.apply(null, xsAll);
    let xmax = cfg.xRange ? cfg.xRange[1] : Math.max.apply(null, xsAll);
    let ymin = cfg.yRange ? cfg.yRange[0] : Math.min.apply(null, ysAll);
    let ymax = cfg.yRange ? cfg.yRange[1] : Math.max.apply(null, ysAll);

    if (!isFinite(xmin) || !isFinite(xmax)) { xmin = 0; xmax = 1; }
    if (!isFinite(ymin) || !isFinite(ymax)) { ymin = 0; ymax = 1; }
    if (xmin === xmax) { xmin -= 1; xmax += 1; }
    if (ymin === ymax) {
      const pad = Math.abs(ymin) > 1 ? Math.abs(ymin) * 0.25 : 1;
      ymin -= pad; ymax += pad;
    }
    const xPad = (xmax - xmin) * 0.06;
    const yPad = (ymax - ymin) * 0.1;
    xmin -= xPad; xmax += xPad;
    ymin -= yPad; ymax += yPad;

    const m = compact
      ? { l: 44, r: 12, t: 18, b: 38 }
      : { l: 64, r: 18, t: 26, b: 44 };
    const plotW = width - m.l - m.r;
    const plotH = height - m.t - m.b;
    const X = function (x) { return m.l + ((x - xmin) / (xmax - xmin)) * plotW; };
    const Y = function (y) { return m.t + plotH - ((y - ymin) / (ymax - ymin)) * plotH; };

    const xTicks = ticks(xmin, xmax, Math.max(3, Math.floor(plotW / (compact ? 62 : 90))));
    const yTicks = ticks(ymin, ymax, Math.max(3, Math.floor(plotH / 55)));

    // ---- static frame ----
    svg.appendChild(svgEl('rect', {
      x: m.l, y: m.t, width: plotW, height: plotH,
      fill: 'var(--chart-bg)', stroke: 'var(--border-soft)', rx: 6
    }));

    xTicks.forEach(function (value) {
      const px = X(value);
      svg.appendChild(svgEl('line', {
        x1: px, y1: m.t, x2: px, y2: m.t + plotH,
        stroke: 'var(--grid)', 'stroke-width': 1
      }));
      const label = svgEl('text', {
        x: px, y: height - m.b + 18, 'text-anchor': 'middle',
        fill: 'var(--text-dim)', 'font-size': 11, 'font-family': 'var(--font-mono)'
      });
      label.textContent = fmt(value, 4);
      svg.appendChild(label);
    });
    yTicks.forEach(function (value) {
      const py = Y(value);
      svg.appendChild(svgEl('line', {
        x1: m.l, y1: py, x2: m.l + plotW, y2: py,
        stroke: 'var(--grid)', 'stroke-width': 1
      }));
      const label = svgEl('text', {
        x: m.l - 10, y: py + 4, 'text-anchor': 'end',
        fill: 'var(--text-dim)', 'font-size': 11, 'font-family': 'var(--font-mono)'
      });
      label.textContent = fmt(value, 4);
      svg.appendChild(label);
    });

    if (ymin < 0 && ymax > 0) {
      svg.appendChild(svgEl('line', {
        x1: m.l, y1: Y(0), x2: m.l + plotW, y2: Y(0),
        stroke: 'var(--border)', 'stroke-width': 1.4
      }));
    }
    if (xmin < 0 && xmax > 0) {
      svg.appendChild(svgEl('line', {
        x1: X(0), y1: m.t, x2: X(0), y2: m.t + plotH,
        stroke: 'var(--border)', 'stroke-width': 1.4
      }));
    }

    if (cfg.xLabel) {
      const lx = svgEl('text', {
        x: m.l + plotW / 2, y: height - 8, 'text-anchor': 'middle',
        fill: 'var(--text-soft)', 'font-size': 12, 'font-weight': 600
      });
      lx.textContent = cfg.xLabel;
      svg.appendChild(lx);
    }
    if (cfg.yLabel) {
      if (compact) {
        const ly = svgEl('text', {
          x: m.l, y: m.t - 6, 'text-anchor': 'start',
          fill: 'var(--text-soft)', 'font-size': 11, 'font-weight': 600
        });
        ly.textContent = cfg.yLabel;
        svg.appendChild(ly);
      } else {
        const ly = svgEl('text', {
          x: 14, y: m.t + plotH / 2, 'text-anchor': 'middle',
          fill: 'var(--text-soft)', 'font-size': 12, 'font-weight': 600,
          transform: 'rotate(-90 14 ' + (m.t + plotH / 2) + ')'
        });
        ly.textContent = cfg.yLabel;
        svg.appendChild(ly);
      }
    }

    // ---- dynamic layers (order matters) ----
    const layers = {
      residuals: svgEl('g', { class: 'layer-residuals' }),
      guides: svgEl('g', { class: 'layer-guides' }),
      series: svgEl('g', { class: 'layer-series' }),
      points: svgEl('g', { class: 'layer-points' }),
      highlights: svgEl('g', { class: 'layer-highlights' }),
      marker: svgEl('g', { class: 'layer-marker' }),
      tooltip: svgEl('g', { class: 'layer-tooltip', style: 'pointer-events: none; opacity: 0; transition: opacity 100ms ease;' })
    };
    ['residuals', 'guides', 'series', 'points', 'highlights', 'marker', 'tooltip'].forEach(function (key) {
      svg.appendChild(layers[key]);
    });

    const legend = document.createElement('div');
    legend.className = 'chart-legend';
    holder.appendChild(svg);
    holder.appendChild(legend);

    let current = cfg;

    function clearLayer(layer) {
      while (layer.firstChild) layer.removeChild(layer.firstChild);
    }

    function renderDynamic() {
      const state = current;
      Object.keys(layers).forEach(function (key) { clearLayer(layers[key]); });

      // residuals
      (state.residuals || []).forEach(function (r) {
        if (!isFinite(r.x) || !isFinite(r.y0) || !isFinite(r.y1)) return;
        layers.residuals.appendChild(svgEl('line', {
          x1: X(r.x), y1: Y(r.y0), x2: X(r.x), y2: Y(r.y1),
          stroke: r.color || 'var(--red)', 'stroke-width': 1.4,
          'stroke-dasharray': '3 3', opacity: 0.75
        }));
      });

      // guides
      if (state.hline && isFinite(state.hline.y)) {
        const hy = Y(state.hline.y);
        layers.guides.appendChild(svgEl('line', {
          x1: m.l, y1: hy, x2: m.l + plotW, y2: hy,
          stroke: state.hline.color || 'var(--orange)', 'stroke-width': 1.6,
          'stroke-dasharray': '7 5'
        }));
        if (state.hline.label) {
          const hl = svgEl('text', {
            x: m.l + plotW - 4, y: hy - 6, 'text-anchor': 'end',
            fill: state.hline.color || 'var(--orange)', 'font-size': 11, 'font-weight': 700
          });
          hl.textContent = state.hline.label;
          layers.guides.appendChild(hl);
        }
      }
      if (state.vline && isFinite(state.vline.x)) {
        const vx = X(state.vline.x);
        layers.guides.appendChild(svgEl('line', {
          x1: vx, y1: m.t, x2: vx, y2: m.t + plotH,
          stroke: state.vline.color || 'var(--green)', 'stroke-width': 1.6,
          'stroke-dasharray': '7 5'
        }));
        if (state.vline.label) {
          const vl = svgEl('text', {
            x: vx + 5, y: m.t + 13, 'text-anchor': 'start',
            fill: state.vline.color || 'var(--green)', 'font-size': 11, 'font-weight': 700
          });
          vl.textContent = state.vline.label;
          layers.guides.appendChild(vl);
        }
      }

      // series
      const legendItems = [];
      (state.series || []).forEach(function (series) {
        if (series.visible === false) return;
        const d = pathData(series.points, X, Y, series.reveal);
        if (d) {
          layers.series.appendChild(svgEl('path', {
            d: d,
            fill: 'none',
            stroke: series.color || 'var(--accent)',
            'stroke-width': series.width || 2.2,
            'stroke-linecap': 'round',
            'stroke-linejoin': 'round',
            'stroke-dasharray': series.dash || 'none',
            opacity: series.opacity == null ? 1 : series.opacity
          }));
        }
        if (series.markers) {
          (series.points || []).forEach(function (p) {
            if (!isFinite(p.x) || !isFinite(p.y)) return;
            layers.series.appendChild(svgEl('circle', {
              cx: X(p.x), cy: Y(p.y), r: 2.4,
              fill: series.color || 'var(--accent)',
              stroke: 'var(--chart-bg)', 'stroke-width': 1
            }));
          });
        }
        if (series.label && series.color) {
          legendItems.push({ label: series.label, color: series.color });
        }
      });

      // data points
      (state.points || []).forEach(function (p) {
        if (!isFinite(p.x) || !isFinite(p.y)) return;
        const circle = svgEl('circle', {
          cx: X(p.x), cy: Y(p.y), r: p.r || 4.6,
          fill: p.color || 'var(--accent)',
          stroke: 'var(--chart-bg)', 'stroke-width': 2
        });
        const title = svgEl('title');
        title.textContent = (p.label ? p.label + ': ' : '') + '(' + fmt(p.x, 4) + ', ' + fmt(p.y, 4) + ')';
        circle.appendChild(title);
        layers.points.appendChild(circle);
      });

      // highlighted points
      (state.highlightPoints || []).forEach(function (p) {
        if (!isFinite(p.x) || !isFinite(p.y)) return;
        const color = p.color || 'var(--green)';
        const outer = svgEl('circle', {
          cx: X(p.x), cy: Y(p.y), r: 8.5,
          fill: 'none', stroke: color, 'stroke-width': 2, opacity: 0.9
        });
        const inner = svgEl('circle', {
          cx: X(p.x), cy: Y(p.y), r: 4,
          fill: color, stroke: 'var(--chart-bg)', 'stroke-width': 1.5
        });
        const title = svgEl('title');
        title.textContent = (p.label ? p.label + ': ' : '') + '(' + fmt(p.x, 4) + ', ' + fmt(p.y, 4) + ')';
        inner.appendChild(title);
        layers.highlights.appendChild(outer);
        layers.highlights.appendChild(inner);
        if (p.label) {
          const t = svgEl('text', {
            x: X(p.x) + 11, y: Y(p.y) - 9,
            fill: color, 'font-size': 10.5, 'font-weight': 700
          });
          t.textContent = p.label;
          layers.highlights.appendChild(t);
        }
      });

      // moving marker
      if (state.marker && isFinite(state.marker.x) && isFinite(state.marker.y)) {
        const color = state.marker.color || 'var(--accent)';
        layers.marker.appendChild(svgEl('circle', {
          cx: X(state.marker.x), cy: Y(state.marker.y), r: 9,
          fill: 'none', stroke: color, 'stroke-width': 1.6, opacity: 0.45
        }));
        layers.marker.appendChild(svgEl('circle', {
          cx: X(state.marker.x), cy: Y(state.marker.y), r: 5,
          fill: color, stroke: 'var(--chart-bg)', 'stroke-width': 1.5
        }));
        if (state.marker.label) {
          const label = svgEl('text', {
            x: X(state.marker.x) + 12, y: Y(state.marker.y) - 10,
            fill: color, 'font-size': 10.5, 'font-weight': 700
          });
          label.textContent = state.marker.label;
          layers.marker.appendChild(label);
        }
      }

      // legend (HTML, rebuilt to follow visible series)
      legend.innerHTML = '';
      legendItems.forEach(function (item) {
        const span = document.createElement('span');
        span.className = 'legend-item';
        const dot = document.createElement('span');
        dot.className = 'legend-dot';
        dot.style.background = item.color;
        span.appendChild(dot);
        span.appendChild(document.createTextNode(item.label));
        legend.appendChild(span);
      });
      legend.style.display = legendItems.length ? 'flex' : 'none';
    }

    // Interactive mouse hover crosshair & tooltip
    const crosshairV = svgEl('line', {
      y1: m.t, y2: m.t + plotH,
      stroke: 'var(--accent)', 'stroke-width': 1, 'stroke-dasharray': '3 3', opacity: 0.6
    });
    const crosshairH = svgEl('line', {
      x1: m.l, x2: m.l + plotW,
      stroke: 'var(--accent)', 'stroke-width': 1, 'stroke-dasharray': '3 3', opacity: 0.6
    });
    const tipBg = svgEl('rect', {
      rx: 5, ry: 5, fill: 'var(--panel-2)', stroke: 'var(--border)', 'stroke-width': 1,
      filter: 'drop-shadow(0 2px 6px rgba(0,0,0,0.28))'
    });
    const tipText = svgEl('text', {
      fill: 'var(--text)', 'font-size': 11, 'font-family': 'var(--font-mono)', 'font-weight': 600
    });
    layers.tooltip.appendChild(crosshairV);
    layers.tooltip.appendChild(crosshairH);
    layers.tooltip.appendChild(tipBg);
    layers.tooltip.appendChild(tipText);

    svg.addEventListener('mousemove', function (e) {
      const rect = svg.getBoundingClientRect();
      if (!rect.width || !rect.height) return;
      const scaleX = width / rect.width;
      const scaleY = height / rect.height;
      const px = (e.clientX - rect.left) * scaleX;
      const py = (e.clientY - rect.top) * scaleY;

      if (px < m.l || px > m.l + plotW || py < m.t || py > m.t + plotH) {
        layers.tooltip.style.opacity = '0';
        return;
      }

      // Check for snapping to nearest data point within 16px
      let snapped = null;
      let minD = 16;
      const allPts = (current.points || []).concat(current.highlightPoints || []);
      for (let i = 0; i < allPts.length; i++) {
        const pt = allPts[i];
        if (!isFinite(pt.x) || !isFinite(pt.y)) continue;
        const dx = X(pt.x) - px;
        const dy = Y(pt.y) - py;
        const d = Math.sqrt(dx * dx + dy * dy);
        if (d < minD) {
          minD = d;
          snapped = pt;
        }
      }

      const curX = snapped ? X(snapped.x) : px;
      const curY = snapped ? Y(snapped.y) : py;
      const valX = snapped ? snapped.x : (xmin + ((px - m.l) / plotW) * (xmax - xmin));
      const valY = snapped ? snapped.y : (ymin + ((m.t + plotH - py) / plotH) * (ymax - ymin));

      crosshairV.setAttribute('x1', curX);
      crosshairV.setAttribute('x2', curX);
      crosshairH.setAttribute('y1', curY);
      crosshairH.setAttribute('y2', curY);

      const labelStr = (snapped && snapped.label ? snapped.label + ' ' : '') + '(' + fmt(valX, 4) + ', ' + fmt(valY, 4) + ')';
      tipText.textContent = labelStr;
      const textLen = labelStr.length * 7 + 16;
      const bgW = Math.max(80, textLen);
      const bgH = 22;
      let bgX = curX + 10;
      let bgY = curY - 26;
      if (bgX + bgW > m.l + plotW) bgX = curX - bgW - 10;
      if (bgY < m.t) bgY = curY + 10;

      tipBg.setAttribute('x', bgX);
      tipBg.setAttribute('y', bgY);
      tipBg.setAttribute('width', bgW);
      tipBg.setAttribute('height', bgH);
      tipText.setAttribute('x', bgX + 8);
      tipText.setAttribute('y', bgY + 15);

      layers.tooltip.style.opacity = '1';
    });

    svg.addEventListener('mouseleave', function () {
      layers.tooltip.style.opacity = '0';
    });

    const handle = {
      svg: svg,
      scales: { X: X, Y: Y, xmin: xmin, xmax: xmax, ymin: ymin, ymax: ymax },
      update: function (patch) {
        if (patch) current = Object.assign({}, current, patch);
        renderDynamic();
        return handle;
      },
      destroy: function () {
        holder.innerHTML = '';
      }
    };

    renderDynamic();
    return handle;
  }

  function draw(holder, cfg) {
    return create(holder, cfg).svg;
  }

  function sampleFunction(fn, xmin, xmax, count) {
    const points = [];
    const n = count || 260;
    for (let i = 0; i <= n; i++) {
      const x = xmin + (xmax - xmin) * i / n;
      let y;
      try { y = fn(x); } catch (e) { y = NaN; }
      points.push({ x: x, y: y });
    }
    return points;
  }

  PPS.chart = {
    create: create,
    draw: draw,
    sample: sampleFunction,
    ticks: ticks
  };
})();
