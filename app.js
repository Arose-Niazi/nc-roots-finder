/* ============================================================
   NC Root Finder — Root Finding Calculator
   Complete rewrite: modern architecture, step-by-step solutions
   ============================================================ */

(function () {
  'use strict';

  // ===== MATH EVALUATION =====
  var compiledExpr = null;
  var compiledDeriv = null;

  function compileExpression(exprStr) {
    try {
      var node = math.parse(exprStr);
      compiledExpr = node.compile();
      // Compute symbolic derivative for Newton-Raphson
      try {
        var derivNode = math.derivative(node, 'x');
        compiledDeriv = derivNode.compile();
      } catch (e) {
        compiledDeriv = null;
      }
      // Test evaluation
      compiledExpr.evaluate({ x: 1 });
      return true;
    } catch (e) {
      compiledExpr = null;
      compiledDeriv = null;
      return false;
    }
  }

  function f(x) {
    return compiledExpr.evaluate({ x: x });
  }

  function fPrime(x) {
    if (compiledDeriv) {
      return compiledDeriv.evaluate({ x: x });
    }
    // Numerical derivative fallback
    var h = 1e-8;
    return (f(x + h) - f(x - h)) / (2 * h);
  }

  function fmt(val, precision) {
    if (!isFinite(val)) return String(val);
    return Number(val.toFixed(precision));
  }

  // ===== ROOT FINDING METHODS =====

  function bisection(a, b, tol, maxIter, precision) {
    var fa = f(a), fb = f(b);
    if (fa * fb > 0) {
      throw new Error('f(a) and f(b) must have opposite signs. f(' + a + ') = ' + fmt(fa, precision) + ', f(' + b + ') = ' + fmt(fb, precision));
    }

    var iterations = [];
    var converged = false;
    var root = null;

    for (var i = 1; i <= maxIter; i++) {
      var c = (a + b) / 2;
      var fc = f(c);
      var err = Math.abs(b - a) / 2;

      iterations.push({
        n: i, a: a, b: b, c: c,
        fa: fa, fb: fb, fc: fc,
        error: err
      });

      if (Math.abs(fc) < tol || err < tol) {
        converged = true;
        root = c;
        break;
      }

      if (fa * fc < 0) {
        b = c;
        fb = fc;
      } else {
        a = c;
        fa = fc;
      }

      if (i === maxIter) {
        root = c;
      }
    }

    return { method: 'Bisection', root: root, iterations: iterations, converged: converged };
  }

  function regulaFalsi(a, b, tol, maxIter, precision) {
    var fa = f(a), fb = f(b);
    if (fa * fb > 0) {
      throw new Error('f(a) and f(b) must have opposite signs. f(' + a + ') = ' + fmt(fa, precision) + ', f(' + b + ') = ' + fmt(fb, precision));
    }

    var iterations = [];
    var converged = false;
    var root = null;
    var prevC = null;

    for (var i = 1; i <= maxIter; i++) {
      var c = (a * fb - b * fa) / (fb - fa);
      var fc = f(c);
      var err = prevC !== null ? Math.abs(c - prevC) : Math.abs(b - a);

      iterations.push({
        n: i, a: a, b: b, c: c,
        fa: fa, fb: fb, fc: fc,
        error: err
      });

      if (Math.abs(fc) < tol || (prevC !== null && Math.abs(c - prevC) < tol)) {
        converged = true;
        root = c;
        break;
      }

      prevC = c;

      if (fa * fc < 0) {
        b = c;
        fb = fc;
      } else {
        a = c;
        fa = fc;
      }

      if (i === maxIter) {
        root = c;
      }
    }

    return { method: 'Regula Falsi', root: root, iterations: iterations, converged: converged };
  }

  function newtonRaphson(x0, tol, maxIter, precision) {
    var iterations = [];
    var converged = false;
    var root = null;
    var x = x0;

    for (var i = 1; i <= maxIter; i++) {
      var fx = f(x);
      var fpx = fPrime(x);

      if (Math.abs(fpx) < 1e-14) {
        throw new Error('Derivative is zero at x = ' + fmt(x, precision) + '. Newton-Raphson method cannot continue.');
      }

      var xNew = x - fx / fpx;
      var err = Math.abs(xNew - x);

      iterations.push({
        n: i, x: x, fx: fx, fpx: fpx,
        xNew: xNew, error: err
      });

      if (Math.abs(f(xNew)) < tol || err < tol) {
        converged = true;
        root = xNew;
        break;
      }

      if (!isFinite(xNew)) {
        throw new Error('Method diverged at iteration ' + i + '. Try a different initial guess.');
      }

      x = xNew;

      if (i === maxIter) {
        root = xNew;
      }
    }

    return { method: 'Newton-Raphson', root: root, iterations: iterations, converged: converged };
  }

  function secant(x0, x1, tol, maxIter, precision) {
    var iterations = [];
    var converged = false;
    var root = null;

    var xPrev = x0, xCurr = x1;
    var fPrev = f(xPrev), fCurr = f(xCurr);

    for (var i = 1; i <= maxIter; i++) {
      if (Math.abs(fCurr - fPrev) < 1e-14) {
        throw new Error('Division by zero: f(x' + (i) + ') ≈ f(x' + (i - 1) + '). Try different initial values.');
      }

      var xNew = xCurr - fCurr * (xCurr - xPrev) / (fCurr - fPrev);
      var fNew = f(xNew);
      var err = Math.abs(xNew - xCurr);

      iterations.push({
        n: i, xPrev: xPrev, xCurr: xCurr,
        fPrev: fPrev, fCurr: fCurr,
        xNew: xNew, fNew: fNew, error: err
      });

      if (Math.abs(fNew) < tol || err < tol) {
        converged = true;
        root = xNew;
        break;
      }

      if (!isFinite(xNew)) {
        throw new Error('Method diverged at iteration ' + i + '. Try different initial values.');
      }

      xPrev = xCurr;
      fPrev = fCurr;
      xCurr = xNew;
      fCurr = fNew;

      if (i === maxIter) {
        root = xNew;
      }
    }

    return { method: 'Secant', root: root, iterations: iterations, converged: converged };
  }

  // ===== RENDER ITERATION TABLE =====
  function renderIterationTable(result, precision) {
    var method = result.method;
    var iters = result.iterations;
    var html = '<div class="iter-table-wrapper"><table class="iter-table"><thead><tr>';

    if (method === 'Bisection' || method === 'Regula Falsi') {
      html += '<th>n</th><th>a</th><th>b</th><th>c</th><th>f(a)</th><th>f(b)</th><th>f(c)</th><th>|Error|</th>';
    } else if (method === 'Newton-Raphson') {
      html += '<th>n</th><th>x<sub>n</sub></th><th>f(x<sub>n</sub>)</th><th>f\'(x<sub>n</sub>)</th><th>x<sub>n+1</sub></th><th>|Error|</th>';
    } else if (method === 'Secant') {
      html += '<th>n</th><th>x<sub>n-1</sub></th><th>x<sub>n</sub></th><th>f(x<sub>n-1</sub>)</th><th>f(x<sub>n</sub>)</th><th>x<sub>n+1</sub></th><th>|Error|</th>';
    }

    html += '</tr></thead><tbody>';

    for (var i = 0; i < iters.length; i++) {
      var it = iters[i];
      var isLast = i === iters.length - 1 && result.converged;
      html += '<tr' + (isLast ? ' class="converged"' : '') + '>';

      if (method === 'Bisection' || method === 'Regula Falsi') {
        html += '<td>' + it.n + '</td>';
        html += '<td>' + fmt(it.a, precision) + '</td>';
        html += '<td>' + fmt(it.b, precision) + '</td>';
        html += '<td>' + fmt(it.c, precision) + '</td>';
        html += '<td>' + fmt(it.fa, precision) + '</td>';
        html += '<td>' + fmt(it.fb, precision) + '</td>';
        html += '<td>' + fmt(it.fc, precision) + '</td>';
        html += '<td>' + fmt(it.error, precision) + '</td>';
      } else if (method === 'Newton-Raphson') {
        html += '<td>' + it.n + '</td>';
        html += '<td>' + fmt(it.x, precision) + '</td>';
        html += '<td>' + fmt(it.fx, precision) + '</td>';
        html += '<td>' + fmt(it.fpx, precision) + '</td>';
        html += '<td>' + fmt(it.xNew, precision) + '</td>';
        html += '<td>' + fmt(it.error, precision) + '</td>';
      } else if (method === 'Secant') {
        html += '<td>' + it.n + '</td>';
        html += '<td>' + fmt(it.xPrev, precision) + '</td>';
        html += '<td>' + fmt(it.xCurr, precision) + '</td>';
        html += '<td>' + fmt(it.fPrev, precision) + '</td>';
        html += '<td>' + fmt(it.fCurr, precision) + '</td>';
        html += '<td>' + fmt(it.xNew, precision) + '</td>';
        html += '<td>' + fmt(it.error, precision) + '</td>';
      }

      html += '</tr>';
    }

    html += '</tbody></table></div>';
    return html;
  }

  // ===== RENDER STEP-BY-STEP =====
  function renderSteps(result, precision, exprStr) {
    var iters = result.iterations;
    var method = result.method;
    var html = '';

    html += '<p>Finding roots of <strong>f(x) = ' + escapeHTML(exprStr) + '</strong> using the <strong>' + method + '</strong> method.</p>';

    if (method === 'Bisection') {
      html += '<p>The Bisection method repeatedly halves the interval [a, b] where f(a) and f(b) have opposite signs, choosing the subinterval that brackets the root.</p>';
    } else if (method === 'Regula Falsi') {
      html += '<p>The Regula Falsi (False Position) method uses linear interpolation between (a, f(a)) and (b, f(b)) to find a better approximation of the root.</p>';
    } else if (method === 'Newton-Raphson') {
      html += '<p>Newton-Raphson uses the tangent line at x<sub>n</sub> to approximate the root: x<sub>n+1</sub> = x<sub>n</sub> − f(x<sub>n</sub>) / f\'(x<sub>n</sub>).</p>';
    } else if (method === 'Secant') {
      html += '<p>The Secant method approximates the derivative using two previous points: x<sub>n+1</sub> = x<sub>n</sub> − f(x<sub>n</sub>)(x<sub>n</sub> − x<sub>n-1</sub>) / (f(x<sub>n</sub>) − f(x<sub>n-1</sub>)).</p>';
    }

    for (var i = 0; i < iters.length; i++) {
      var it = iters[i];
      html += '<h3 class="step-heading">Iteration ' + it.n + '</h3>';
      html += '<div class="step-block">';

      if (method === 'Bisection') {
        html += '<div class="step-line">Interval: [' + fmt(it.a, precision) + ', ' + fmt(it.b, precision) + ']</div>';
        html += '<div class="step-line">c = (a + b) / 2 = (' + fmt(it.a, precision) + ' + ' + fmt(it.b, precision) + ') / 2 = <strong>' + fmt(it.c, precision) + '</strong></div>';
        html += '<div class="step-line">f(a) = ' + fmt(it.fa, precision) + ', f(c) = ' + fmt(it.fc, precision) + '</div>';
        if (i < iters.length - 1 || !result.converged) {
          if (it.fa * it.fc < 0) {
            html += '<div class="step-line">f(a)·f(c) < 0 → root is in [' + fmt(it.a, precision) + ', ' + fmt(it.c, precision) + ']</div>';
          } else {
            html += '<div class="step-line">f(a)·f(c) > 0 → root is in [' + fmt(it.c, precision) + ', ' + fmt(it.b, precision) + ']</div>';
          }
        }
      } else if (method === 'Regula Falsi') {
        html += '<div class="step-line">Interval: [' + fmt(it.a, precision) + ', ' + fmt(it.b, precision) + ']</div>';
        html += '<div class="step-line">f(a) = ' + fmt(it.fa, precision) + ', f(b) = ' + fmt(it.fb, precision) + '</div>';
        html += '<div class="step-line">c = (a·f(b) − b·f(a)) / (f(b) − f(a))</div>';
        html += '<div class="step-line">c = (' + fmt(it.a, precision) + '·' + fmt(it.fb, precision) + ' − ' + fmt(it.b, precision) + '·' + fmt(it.fa, precision) + ') / (' + fmt(it.fb, precision) + ' − ' + fmt(it.fa, precision) + ')</div>';
        html += '<div class="step-line">c = <strong>' + fmt(it.c, precision) + '</strong></div>';
        html += '<div class="step-line">f(c) = ' + fmt(it.fc, precision) + '</div>';
      } else if (method === 'Newton-Raphson') {
        html += '<div class="step-line">x<sub>' + (it.n - 1) + '</sub> = ' + fmt(it.x, precision) + '</div>';
        html += '<div class="step-line">f(x<sub>' + (it.n - 1) + '</sub>) = ' + fmt(it.fx, precision) + '</div>';
        html += '<div class="step-line">f\'(x<sub>' + (it.n - 1) + '</sub>) = ' + fmt(it.fpx, precision) + '</div>';
        html += '<div class="step-line">x<sub>' + it.n + '</sub> = ' + fmt(it.x, precision) + ' − (' + fmt(it.fx, precision) + ') / (' + fmt(it.fpx, precision) + ') = <strong>' + fmt(it.xNew, precision) + '</strong></div>';
      } else if (method === 'Secant') {
        html += '<div class="step-line">x<sub>' + (it.n - 1) + '</sub> = ' + fmt(it.xPrev, precision) + ', x<sub>' + it.n + '</sub> = ' + fmt(it.xCurr, precision) + '</div>';
        html += '<div class="step-line">f(x<sub>' + (it.n - 1) + '</sub>) = ' + fmt(it.fPrev, precision) + ', f(x<sub>' + it.n + '</sub>) = ' + fmt(it.fCurr, precision) + '</div>';
        html += '<div class="step-line">x<sub>' + (it.n + 1) + '</sub> = ' + fmt(it.xCurr, precision) + ' − ' + fmt(it.fCurr, precision) + ' × (' + fmt(it.xCurr, precision) + ' − ' + fmt(it.xPrev, precision) + ') / (' + fmt(it.fCurr, precision) + ' − ' + fmt(it.fPrev, precision) + ')</div>';
        html += '<div class="step-line">x<sub>' + (it.n + 1) + '</sub> = <strong>' + fmt(it.xNew, precision) + '</strong></div>';
      }

      html += '<div class="step-line">|Error| = ' + fmt(it.error, precision) + '</div>';
      html += '</div>';
    }

    return html;
  }

  // ===== CONVERGENCE GRAPH =====
  function drawConvergenceGraph(canvas, result, precision) {
    var ctx = canvas.getContext('2d');
    var dpr = window.devicePixelRatio || 1;
    var rect = canvas.getBoundingClientRect();
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
    ctx.scale(dpr, dpr);
    var W = rect.width, H = rect.height;

    // Get theme colors
    var style = getComputedStyle(document.documentElement);
    var bgColor = style.getPropertyValue('--bg-secondary').trim();
    var textColor = style.getPropertyValue('--text-muted').trim();
    var borderColor = style.getPropertyValue('--border-color').trim();
    var accentColor = style.getPropertyValue('--accent').trim();
    var dangerColor = style.getPropertyValue('--danger').trim();

    // Clear
    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, W, H);

    var iters = result.iterations;
    if (iters.length < 2) return;

    // Extract approximation values and errors
    var approxValues = [];
    var errorValues = [];
    for (var i = 0; i < iters.length; i++) {
      var it = iters[i];
      if (result.method === 'Newton-Raphson') {
        approxValues.push(it.xNew);
      } else if (result.method === 'Secant') {
        approxValues.push(it.xNew);
      } else {
        approxValues.push(it.c);
      }
      errorValues.push(Math.max(it.error, 1e-16));
    }

    // Padding
    var pad = { top: 30, right: 30, bottom: 40, left: 60 };
    var plotW = W - pad.left - pad.right;
    var plotH = H - pad.top - pad.bottom;

    // Use log scale for errors
    var logErrors = errorValues.map(function (e) { return Math.log10(e); });
    var minLog = Math.min.apply(null, logErrors);
    var maxLog = Math.max.apply(null, logErrors);
    if (maxLog - minLog < 1) { minLog -= 0.5; maxLog += 0.5; }
    var logRange = maxLog - minLog;

    // Axes
    ctx.strokeStyle = borderColor;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(pad.left, pad.top);
    ctx.lineTo(pad.left, pad.top + plotH);
    ctx.lineTo(pad.left + plotW, pad.top + plotH);
    ctx.stroke();

    // Grid lines & labels
    ctx.fillStyle = textColor;
    ctx.font = '11px -apple-system, sans-serif';
    ctx.textAlign = 'center';

    // X axis labels
    var xStep = Math.max(1, Math.floor(iters.length / 8));
    for (var i2 = 0; i2 < iters.length; i2 += xStep) {
      var x = pad.left + (i2 / (iters.length - 1)) * plotW;
      ctx.fillText(i2 + 1, x, pad.top + plotH + 18);

      ctx.strokeStyle = borderColor;
      ctx.globalAlpha = 0.3;
      ctx.beginPath();
      ctx.moveTo(x, pad.top);
      ctx.lineTo(x, pad.top + plotH);
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    // Y axis labels (log scale)
    ctx.textAlign = 'right';
    var yTicks = 5;
    for (var t = 0; t <= yTicks; t++) {
      var logVal = maxLog - (t / yTicks) * logRange;
      var y = pad.top + (t / yTicks) * plotH;
      ctx.fillText('10^' + logVal.toFixed(1), pad.left - 8, y + 4);

      ctx.strokeStyle = borderColor;
      ctx.globalAlpha = 0.3;
      ctx.beginPath();
      ctx.moveTo(pad.left, y);
      ctx.lineTo(pad.left + plotW, y);
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    // Axis titles
    ctx.fillStyle = textColor;
    ctx.textAlign = 'center';
    ctx.font = '12px -apple-system, sans-serif';
    ctx.fillText('Iteration', pad.left + plotW / 2, H - 5);

    ctx.save();
    ctx.translate(14, pad.top + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Error (log₁₀)', 0, 0);
    ctx.restore();

    // Plot error line
    ctx.strokeStyle = dangerColor;
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (var i3 = 0; i3 < logErrors.length; i3++) {
      var px = pad.left + (i3 / (iters.length - 1)) * plotW;
      var py = pad.top + ((maxLog - logErrors[i3]) / logRange) * plotH;
      if (i3 === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Plot dots
    for (var i4 = 0; i4 < logErrors.length; i4++) {
      var dx = pad.left + (i4 / (iters.length - 1)) * plotW;
      var dy = pad.top + ((maxLog - logErrors[i4]) / logRange) * plotH;
      ctx.fillStyle = i4 === logErrors.length - 1 && result.converged ? accentColor : dangerColor;
      ctx.beginPath();
      ctx.arc(dx, dy, 4, 0, 2 * Math.PI);
      ctx.fill();
    }

    // Plot root approximation line
    var minApprox = Math.min.apply(null, approxValues);
    var maxApprox = Math.max.apply(null, approxValues);
    var approxRange = maxApprox - minApprox;
    if (approxRange < 1e-10) approxRange = 1;

    ctx.strokeStyle = accentColor;
    ctx.lineWidth = 2;
    ctx.setLineDash([5, 3]);
    ctx.beginPath();
    for (var i5 = 0; i5 < approxValues.length; i5++) {
      var ax = pad.left + (i5 / (iters.length - 1)) * plotW;
      var ay = pad.top + ((maxApprox - approxValues[i5]) / approxRange) * plotH;
      if (i5 === 0) ctx.moveTo(ax, ay);
      else ctx.lineTo(ax, ay);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // ===== LATEX EXPORT =====
  function generateLaTeX(result, precision, exprStr) {
    var tex = '% Root Finding: ' + result.method + ' Method\n';
    tex += '% f(x) = ' + exprStr + '\n\n';

    if (result.converged) {
      tex += '\\text{Root found: } x \\approx ' + fmt(result.root, precision) + ' \\\\\n';
      tex += '\\text{After } ' + result.iterations.length + ' \\text{ iterations} \\\\\n\n';
    }

    // Iteration table
    tex += '\\begin{array}{';
    var method = result.method;
    if (method === 'Bisection' || method === 'Regula Falsi') {
      tex += '|c|c|c|c|c|c|c|c|}\n\\hline\n';
      tex += 'n & a & b & c & f(a) & f(b) & f(c) & |\\text{Error}| \\\\\n\\hline\n';
      for (var i = 0; i < result.iterations.length; i++) {
        var it = result.iterations[i];
        tex += it.n + ' & ' + fmt(it.a, precision) + ' & ' + fmt(it.b, precision) + ' & ' + fmt(it.c, precision) + ' & ' + fmt(it.fa, precision) + ' & ' + fmt(it.fb, precision) + ' & ' + fmt(it.fc, precision) + ' & ' + fmt(it.error, precision) + ' \\\\\n';
      }
    } else if (method === 'Newton-Raphson') {
      tex += '|c|c|c|c|c|c|}\n\\hline\n';
      tex += 'n & x_n & f(x_n) & f\'(x_n) & x_{n+1} & |\\text{Error}| \\\\\n\\hline\n';
      for (var j = 0; j < result.iterations.length; j++) {
        var jt = result.iterations[j];
        tex += jt.n + ' & ' + fmt(jt.x, precision) + ' & ' + fmt(jt.fx, precision) + ' & ' + fmt(jt.fpx, precision) + ' & ' + fmt(jt.xNew, precision) + ' & ' + fmt(jt.error, precision) + ' \\\\\n';
      }
    } else if (method === 'Secant') {
      tex += '|c|c|c|c|c|c|c|}\n\\hline\n';
      tex += 'n & x_{n-1} & x_n & f(x_{n-1}) & f(x_n) & x_{n+1} & |\\text{Error}| \\\\\n\\hline\n';
      for (var k = 0; k < result.iterations.length; k++) {
        var kt = result.iterations[k];
        tex += kt.n + ' & ' + fmt(kt.xPrev, precision) + ' & ' + fmt(kt.xCurr, precision) + ' & ' + fmt(kt.fPrev, precision) + ' & ' + fmt(kt.fCurr, precision) + ' & ' + fmt(kt.xNew, precision) + ' & ' + fmt(kt.error, precision) + ' \\\\\n';
      }
    }

    tex += '\\hline\n\\end{array}\n';
    return tex;
  }

  // ===== HISTORY =====
  var HISTORY_KEY = 'nc_roots_history';

  function loadHistory() {
    try { return JSON.parse(localStorage.getItem(HISTORY_KEY)) || []; }
    catch (e) { return []; }
  }

  function saveHistory(entry) {
    var history = loadHistory();
    history.unshift(entry);
    if (history.length > 50) history.length = 50;
    localStorage.setItem(HISTORY_KEY, JSON.stringify(history));
  }

  function clearHistory() {
    localStorage.removeItem(HISTORY_KEY);
    renderHistoryPanel();
  }

  function renderHistoryPanel() {
    var list = document.getElementById('history-list');
    var history = loadHistory();
    if (history.length === 0) {
      list.innerHTML = '<p class="history-empty">No calculations yet</p>';
      return;
    }
    var html = '';
    for (var i = 0; i < history.length; i++) {
      var h = history[i];
      html += '<div class="history-item" data-index="' + i + '">';
      html += '<div class="history-meta">';
      html += '<span class="history-method">' + h.method + '</span>';
      html += '<span class="history-badge">' + h.iterations + ' iters</span>';
      html += '<span class="history-date">' + new Date(h.date).toLocaleDateString() + '</span>';
      html += '</div>';
      html += '<div class="history-preview">f(x) = ' + escapeHTML(h.equation) + '</div>';
      html += '</div>';
    }
    list.innerHTML = html;

    list.addEventListener('click', function (e) {
      var item = e.target.closest('.history-item');
      if (!item) return;
      var idx = parseInt(item.getAttribute('data-index'));
      if (!isNaN(idx) && history[idx]) {
        loadFromHistory(history[idx]);
        document.getElementById('history-panel').classList.remove('open');
      }
    });
  }

  function loadFromHistory(h) {
    document.getElementById('equation-input').value = h.equation;
    document.getElementById('method-select').value = h.methodKey;
    updateMethodUI();
    if (h.a !== undefined) document.getElementById('interval-a').value = h.a;
    if (h.b !== undefined) document.getElementById('interval-b').value = h.b;
    if (h.x0 !== undefined) document.getElementById('x0-input').value = h.x0;
    if (h.x1 !== undefined) document.getElementById('x1-input').value = h.x1;
    if (h.precision) document.getElementById('precision-input').value = h.precision;
    if (h.maxIter) document.getElementById('max-iter-input').value = h.maxIter;
    if (h.tolerance) document.getElementById('tolerance-input').value = h.tolerance;
  }

  // ===== URL SHARING =====
  function encodeToURL(params) {
    var urlParams = new URLSearchParams();
    urlParams.set('eq', params.equation);
    urlParams.set('method', params.methodKey);
    if (params.a !== undefined) urlParams.set('a', params.a);
    if (params.b !== undefined) urlParams.set('b', params.b);
    if (params.x0 !== undefined) urlParams.set('x0', params.x0);
    if (params.x1 !== undefined) urlParams.set('x1', params.x1);
    urlParams.set('p', params.precision);
    urlParams.set('max', params.maxIter);
    urlParams.set('tol', params.tolerance);
    return window.location.origin + window.location.pathname + '?' + urlParams.toString();
  }

  function decodeFromURL() {
    var params = new URLSearchParams(window.location.search);
    if (!params.has('eq')) return null;
    return {
      equation: params.get('eq'),
      methodKey: params.get('method') || 'bisection',
      a: params.has('a') ? parseFloat(params.get('a')) : undefined,
      b: params.has('b') ? parseFloat(params.get('b')) : undefined,
      x0: params.has('x0') ? parseFloat(params.get('x0')) : undefined,
      x1: params.has('x1') ? parseFloat(params.get('x1')) : undefined,
      precision: parseInt(params.get('p')) || 6,
      maxIter: parseInt(params.get('max')) || 50,
      tolerance: parseFloat(params.get('tol')) || 0.00001
    };
  }

  // ===== THEME =====
  function initTheme() {
    var saved = localStorage.getItem('nc_theme') || 'light';
    document.documentElement.setAttribute('data-theme', saved);
    updateThemeIcon(saved);
  }

  function toggleTheme() {
    var current = document.documentElement.getAttribute('data-theme');
    var next = current === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    localStorage.setItem('nc_theme', next);
    updateThemeIcon(next);
  }

  function updateThemeIcon(theme) {
    var btn = document.getElementById('theme-toggle');
    if (btn) btn.textContent = theme === 'dark' ? '☀️' : '🌙';
  }

  // ===== UI HELPERS =====
  function updateMethodUI() {
    var method = document.getElementById('method-select').value;
    var intervalA = document.getElementById('interval-a-group');
    var intervalB = document.getElementById('interval-b-group');
    var x0Group = document.getElementById('x0-group');
    var x1Group = document.getElementById('x1-group');
    var note = document.getElementById('method-note');

    if (method === 'bisection' || method === 'regulafalsi') {
      intervalA.style.display = '';
      intervalB.style.display = '';
      x0Group.style.display = 'none';
      x1Group.style.display = 'none';
      if (method === 'bisection') {
        note.textContent = 'Bisection method requires f(a) and f(b) to have opposite signs. Converges linearly.';
      } else {
        note.textContent = 'Regula Falsi (False Position) uses linear interpolation. Requires f(a) and f(b) to have opposite signs.';
      }
      note.classList.add('visible');
    } else if (method === 'newton') {
      intervalA.style.display = 'none';
      intervalB.style.display = 'none';
      x0Group.style.display = '';
      x1Group.style.display = 'none';
      note.textContent = 'Newton-Raphson requires an initial guess x₀. Derivative is computed automatically. Converges quadratically near the root.';
      note.classList.add('visible');
    } else if (method === 'secant') {
      intervalA.style.display = 'none';
      intervalB.style.display = 'none';
      x0Group.style.display = '';
      x1Group.style.display = '';
      note.textContent = 'Secant method requires two initial guesses x₀ and x₁. Does not need the derivative. Convergence order ≈ 1.618.';
      note.classList.add('visible');
    }
  }

  function showError(msg) {
    document.getElementById('output').innerHTML = '<div class="error-message"><strong>Error:</strong> ' + msg + '</div>';
  }

  function escapeHTML(str) {
    return str.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  }

  function loadExample() {
    var method = document.getElementById('method-select').value;
    document.getElementById('equation-input').value = 'x^3 - x - 2';

    if (method === 'bisection' || method === 'regulafalsi') {
      document.getElementById('interval-a').value = '1';
      document.getElementById('interval-b').value = '2';
    } else if (method === 'newton') {
      document.getElementById('x0-input').value = '1.5';
    } else if (method === 'secant') {
      document.getElementById('x0-input').value = '1';
      document.getElementById('x1-input').value = '2';
    }

    document.getElementById('precision-input').value = '6';
    document.getElementById('max-iter-input').value = '50';
    document.getElementById('tolerance-input').value = '0.00001';
  }

  // ===== MAIN CALCULATION =====
  function calculate() {
    // Clear errors
    var errorEls = document.querySelectorAll('.input-error');
    for (var e = 0; e < errorEls.length; e++) errorEls[e].classList.remove('input-error');

    var exprStr = document.getElementById('equation-input').value.trim();
    if (!exprStr) {
      document.getElementById('equation-input').classList.add('input-error');
      showError('Please enter an equation.');
      return;
    }

    if (!compileExpression(exprStr)) {
      document.getElementById('equation-input').classList.add('input-error');
      showError('Invalid equation. Please check the syntax. Example: x^3 - x - 2');
      return;
    }

    var method = document.getElementById('method-select').value;
    var precision = parseInt(document.getElementById('precision-input').value) || 6;
    var maxIter = parseInt(document.getElementById('max-iter-input').value) || 50;
    var tolerance = parseFloat(document.getElementById('tolerance-input').value) || 0.00001;

    var result;
    var historyEntry = {
      date: Date.now(),
      equation: exprStr,
      methodKey: method,
      precision: precision,
      maxIter: maxIter,
      tolerance: tolerance
    };

    try {
      if (method === 'bisection' || method === 'regulafalsi') {
        var a = parseFloat(document.getElementById('interval-a').value);
        var b = parseFloat(document.getElementById('interval-b').value);
        if (isNaN(a)) { document.getElementById('interval-a').classList.add('input-error'); showError('Please enter a valid lower bound (a).'); return; }
        if (isNaN(b)) { document.getElementById('interval-b').classList.add('input-error'); showError('Please enter a valid upper bound (b).'); return; }
        if (a >= b) { showError('Lower bound (a) must be less than upper bound (b).'); return; }

        historyEntry.a = a;
        historyEntry.b = b;

        if (method === 'bisection') {
          result = bisection(a, b, tolerance, maxIter, precision);
        } else {
          result = regulaFalsi(a, b, tolerance, maxIter, precision);
        }
      } else if (method === 'newton') {
        var x0 = parseFloat(document.getElementById('x0-input').value);
        if (isNaN(x0)) { document.getElementById('x0-input').classList.add('input-error'); showError('Please enter a valid initial guess (x₀).'); return; }

        historyEntry.x0 = x0;
        result = newtonRaphson(x0, tolerance, maxIter, precision);
      } else if (method === 'secant') {
        var x0s = parseFloat(document.getElementById('x0-input').value);
        var x1s = parseFloat(document.getElementById('x1-input').value);
        if (isNaN(x0s)) { document.getElementById('x0-input').classList.add('input-error'); showError('Please enter a valid first guess (x₀).'); return; }
        if (isNaN(x1s)) { document.getElementById('x1-input').classList.add('input-error'); showError('Please enter a valid second guess (x₁).'); return; }

        historyEntry.x0 = x0s;
        historyEntry.x1 = x1s;
        result = secant(x0s, x1s, tolerance, maxIter, precision);
      }
    } catch (err) {
      showError(err.message);
      return;
    }

    // Build output
    var output = document.getElementById('output');
    var html = '';

    // Result panel
    html += '<div class="result-panel">';
    html += '<h2 class="result-title">' + result.method + ' Method</h2>';

    // Root summary
    html += '<div class="root-summary">';
    if (result.converged) {
      html += '<strong>Root found:</strong> x ≈ <span class="root-value">' + fmt(result.root, precision) + '</span>';
      html += '<div class="root-meta">Converged after ' + result.iterations.length + ' iteration' + (result.iterations.length > 1 ? 's' : '') + ' | f(root) ≈ ' + fmt(f(result.root), precision) + '</div>';
    } else {
      html += '<strong>Did not converge</strong> after ' + result.iterations.length + ' iterations.';
      html += '<div class="root-meta">Best approximation: x ≈ ' + fmt(result.root, precision) + ' | f(x) ≈ ' + fmt(f(result.root), precision) + '</div>';
    }
    html += '</div>';

    // Iteration table
    html += '<div class="result-section"><h3>Iteration Table</h3>';
    html += renderIterationTable(result, precision);
    html += '</div>';

    // Convergence graph
    html += '<div class="result-section"><h3>Convergence Graph</h3>';
    html += '<div class="convergence-graph-wrapper">';
    html += '<canvas id="convergence-canvas" class="convergence-canvas"></canvas>';
    html += '</div>';
    html += '<div class="graph-legend">';
    html += '<span class="legend-item"><span class="legend-dot" style="background:var(--danger)"></span> Error</span>';
    html += '<span class="legend-item"><span class="legend-dot" style="background:var(--accent)"></span> Root Approximation</span>';
    html += '</div>';
    html += '</div>';

    // Step-by-step
    html += '<details class="steps-details"><summary>Step-by-Step Solution</summary>';
    html += '<div class="steps-content">' + renderSteps(result, precision, exprStr) + '</div>';
    html += '</details>';

    // LaTeX
    var latex = generateLaTeX(result, precision, exprStr);
    html += '<details class="latex-details"><summary>LaTeX Export</summary>';
    html += '<div class="latex-content"><pre class="latex-code">' + escapeHTML(latex) + '</pre>';
    html += '<button class="btn btn-small btn-copy" id="copy-latex-btn">Copy LaTeX</button>';
    html += '</div></details>';

    html += '</div>'; // result-panel

    // Share
    var shareURL = encodeToURL(historyEntry);
    html += '<div class="share-section">';
    html += '<button class="btn btn-small btn-share" id="share-btn">Share Link</button>';
    html += '<button class="btn btn-small" id="print-btn">🖨 Print</button>';
    html += '<input type="text" class="share-url" id="share-url" value="' + escapeHTML(shareURL) + '" readonly>';
    html += '</div>';

    output.innerHTML = html;

    // Draw convergence graph
    var canvas = document.getElementById('convergence-canvas');
    if (canvas && result.iterations.length >= 2) {
      drawConvergenceGraph(canvas, result, precision);
    }

    // Event listeners
    document.getElementById('share-btn').addEventListener('click', function () {
      var urlInput = document.getElementById('share-url');
      navigator.clipboard.writeText(urlInput.value).then(function () {
        var btn = document.getElementById('share-btn');
        btn.textContent = 'Copied!';
        setTimeout(function () { btn.textContent = 'Share Link'; }, 1500);
      }).catch(function () {
        document.getElementById('share-url').select();
      });
    });

    document.getElementById('print-btn').addEventListener('click', function () {
      window.print();
    });

    document.getElementById('copy-latex-btn').addEventListener('click', function () {
      var pre = this.previousElementSibling;
      navigator.clipboard.writeText(pre.textContent).then(function () {
        var btn = document.getElementById('copy-latex-btn');
        btn.textContent = 'Copied!';
        setTimeout(function () { btn.textContent = 'Copy LaTeX'; }, 1500);
      });
    });

    // Save history
    historyEntry.method = result.method;
    historyEntry.iterations = result.iterations.length;
    historyEntry.root = result.root;
    saveHistory(historyEntry);
    renderHistoryPanel();
  }

  // ===== INITIALIZATION =====
  function init() {
    initTheme();

    document.getElementById('theme-toggle').addEventListener('click', toggleTheme);
    document.getElementById('method-select').addEventListener('change', updateMethodUI);
    document.getElementById('calculate-btn').addEventListener('click', calculate);
    document.getElementById('example-btn').addEventListener('click', function () {
      loadExample();
      setTimeout(calculate, 50);
    });
    document.getElementById('clear-btn').addEventListener('click', function () {
      document.getElementById('output').innerHTML = '';
      var errorEls = document.querySelectorAll('.input-error');
      for (var i = 0; i < errorEls.length; i++) errorEls[i].classList.remove('input-error');
    });
    document.getElementById('clear-history-btn').addEventListener('click', clearHistory);
    document.getElementById('history-toggle-btn').addEventListener('click', function () {
      document.getElementById('history-panel').classList.toggle('open');
    });

    updateMethodUI();
    renderHistoryPanel();

    // Load from URL
    var urlData = decodeFromURL();
    if (urlData) {
      document.getElementById('equation-input').value = urlData.equation;
      document.getElementById('method-select').value = urlData.methodKey;
      updateMethodUI();
      if (urlData.a !== undefined) document.getElementById('interval-a').value = urlData.a;
      if (urlData.b !== undefined) document.getElementById('interval-b').value = urlData.b;
      if (urlData.x0 !== undefined) document.getElementById('x0-input').value = urlData.x0;
      if (urlData.x1 !== undefined) document.getElementById('x1-input').value = urlData.x1;
      document.getElementById('precision-input').value = urlData.precision;
      document.getElementById('max-iter-input').value = urlData.maxIter;
      document.getElementById('tolerance-input').value = urlData.tolerance;
      setTimeout(calculate, 100);
    }

    // Keyboard shortcut
    document.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' && e.ctrlKey) calculate();
    });

    // Redraw graph on resize
    window.addEventListener('resize', function () {
      var canvas = document.getElementById('convergence-canvas');
      if (canvas && canvas._lastResult) {
        drawConvergenceGraph(canvas, canvas._lastResult, canvas._lastPrecision);
      }
    });
  }

  document.addEventListener('DOMContentLoaded', init);
})();
