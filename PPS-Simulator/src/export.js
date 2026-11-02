/* ============================================================
   export.js — xuất báo cáo kết quả hiện tại ra file .txt
   ============================================================ */

(function () {
  const PPS = (window.PPS = window.PPS || {});

  // Trích xuất text dễ đọc từ một phần tử DOM (dùng cho báo cáo text).
  function elementText(el) {
    if (!el || !el.tagName) return [];
    const tag = el.tagName.toLowerCase();
    const cls = el.className || '';

    // Bỏ qua các điều khiển / animation / biểu đồ (không mang nghĩa trong text)
    if (cls.indexOf('anim-bar') !== -1 ||
        cls.indexOf('btn') !== -1 ||
        cls.indexOf('select') !== -1 ||
        cls.indexOf('input') !== -1 ||
        cls.indexOf('textarea') !== -1 ||
        cls.indexOf('chart-holder') !== -1 ||
        cls.indexOf('chart-card') !== -1 ||
        cls.indexOf('chart-legend') !== -1 ||
        cls.indexOf('legend-item') !== -1) {
      return [];
    }
    if (tag === 'svg' || tag === 'canvas') return [];

    const lines = [];

    if (tag === 'h1' || tag === 'h2' || tag === 'h3') {
      const text = el.textContent.trim();
      if (text) lines.push(text);
    } else if (tag === 'table') {
      // Lặp từng hàng (tr), bao gồm cả thead và tbody
      const rows = el.querySelectorAll('tr');
      rows.forEach(function (row) {
        const cells = Array.from(row.querySelectorAll('th, td')).map(function (c) {
          return c.textContent.trim();
        });
        if (cells.some(function (c) { return c.length > 0; })) {
          lines.push(cells.join(' | '));
        }
      });
    } else if (cls.indexOf('kv-grid') !== -1) {
      Array.from(el.children).forEach(function (kv) {
        const k = kv.querySelector('.k');
        const v = kv.querySelector('.v');
        const sub = kv.querySelector('.sub');
        let line = '';
        if (k) line += k.textContent.trim() + ': ';
        if (v) line += v.textContent.trim();
        if (sub) line += ' (' + sub.textContent.trim() + ')';
        if (line.trim()) lines.push(line);
      });
    } else if (cls.indexOf('formula') !== -1) {
      const text = el.textContent.trim();
      if (text) lines.push('Công thức: ' + text);
    } else if (cls.indexOf('section-title') !== -1) {
      const text = el.textContent.trim();
      if (text) lines.push('--- ' + text + ' ---');
    } else if (cls.indexOf('note') !== -1 || cls.indexOf('error-box') !== -1 || cls.indexOf('result-head') !== -1) {
      const text = el.textContent.trim();
      if (text) lines.push(text);
    } else if (cls.indexOf('page-head') !== -1) {
      // Tiêu đề trang được xử lý riêng ở captureText
    } else if (cls.indexOf('card') !== -1 ||
               cls.indexOf('field') !== -1 ||
               cls.indexOf('form-grid') !== -1 ||
               cls.indexOf('actions') !== -1 ||
               cls.indexOf('result-section') !== -1) {
      // Bỏ qua khung form, chỉ quét sâu vào các phần nội dung đã được xử lý ở trên
    } else {
      // Phần tử không xác định: quét con để lấy text
      const text = el.textContent.trim();
      if (text && el.children.length === 0) {
        lines.push(text);
      } else {
        Array.from(el.children).forEach(function (child) {
          Array.prototype.push.apply(lines, elementText(child));
        });
      }
    }

    return lines.filter(function (l) { return l && l.trim().length > 0; });
  }

  // Gom báo cáo từ container (thường là #view-root).
  function captureText(container) {
    if (!container) return '';
    const lines = [];

    // Tiêu đề trang
    const pageHead = container.querySelector('.page-head');
    if (pageHead) {
      const eyebrow = pageHead.querySelector('.page-eyebrow');
      const h2 = pageHead.querySelector('h2');
      const p = pageHead.querySelector('.page-desc');
      if (eyebrow) lines.push(eyebrow.textContent.trim());
      if (h2) lines.push(h2.textContent.trim());
      if (p) lines.push(p.textContent.trim());
      lines.push('');
    }

    // Nội dung các section
    Array.from(container.children).forEach(function (child) {
      const part = elementText(child);
      if (part.length) {
        Array.prototype.push.apply(lines, part);
        lines.push('');
      }
    });

    return lines.join('\n').trim();
  }

  // Tạo file .txt và triggers download.
  function downloadText(filename, text) {
    const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  PPS.export = {
    captureText: captureText,
    downloadText: downloadText,
    // Xuất báo cáo của thuật toán hiện tại.
    exportCurrent: function (container, algoId) {
      const text = captureText(container);
      if (!text || !text.trim()) return false;
      const now = new Date();
      const ts = now.getFullYear() +
        String(now.getMonth() + 1).padStart(2, '0') +
        String(now.getDate()).padStart(2, '0') + '_' +
        String(now.getHours()).padStart(2, '0') +
        String(now.getMinutes()).padStart(2, '0');
      const name = (algoId || 'pps') + '_' + ts + '.txt';
      downloadText(name, text);
      return true;
    },
    // Sao chép toàn bộ báo cáo vào clipboard
    copyCurrent: function (container, btn) {
      const text = captureText(container);
      if (!text || !text.trim()) return false;
      if (PPS.ui && PPS.ui.copyToClipboard) {
        PPS.ui.copyToClipboard(text, btn);
      } else if (navigator.clipboard) {
        navigator.clipboard.writeText(text);
      }
      return true;
    }
  };
})();
