# NC Root Finder

A modern root finding calculator with step-by-step solutions, convergence graphs, and LaTeX export.

## Methods

- **Bisection** — Interval halving with guaranteed convergence
- **Regula Falsi** — False position method with linear interpolation
- **Newton-Raphson** — Tangent line method with quadratic convergence
- **Secant** — Derivative-free variant with superlinear convergence

## Features

- Step-by-step iteration tables
- Convergence graph visualization
- Dark/light theme
- Configurable precision, tolerance, and max iterations
- History (localStorage)
- Share via URL parameters
- LaTeX export
- Print support
- Full SEO (meta tags, OG, JSON-LD, sitemap)
- PWA manifest

## Tech

Vanilla HTML/CSS/JS with [math.js](https://mathjs.org/) from CDN.

## Docker

```bash
docker compose up -d --build
```

## Live

[https://nc.arose-niazi.me/roots/](https://nc.arose-niazi.me/roots/)

## License

MIT
