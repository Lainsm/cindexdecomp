// Typesets the deck's three formula SVGs from TeX with MathJax. No LaTeX
// install is needed, and every glyph is written as a path, so the files render
// identically on a venue laptop, in a browser and when dropped into PowerPoint.
//
// Regenerate from the repository root:
//   npm install --prefix "$TMPDIR/mj" mathjax-full@3
//   NODE_PATH="$TMPDIR/mj/node_modules" node dev/talks/formulas/render.cjs

const fs = require("fs");
const path = require("path");
const { mathjax } = require("mathjax-full/js/mathjax.js");
const { TeX } = require("mathjax-full/js/input/tex.js");
const { SVG } = require("mathjax-full/js/output/svg.js");
const { liteAdaptor } = require("mathjax-full/js/adaptors/liteAdaptor.js");
const { RegisterHTMLHandler } = require("mathjax-full/js/handlers/html.js");
// Listing a TeX package by name only works once its configuration is loaded.
require("mathjax-full/js/input/tex/ams/AmsConfiguration.js");
require("mathjax-full/js/input/tex/color/ColorConfiguration.js");

// Mirrors symposium.scss. Colour is baked in rather than left as
// currentColor: an SVG shown through <img> cannot inherit the page's colour,
// and currentColor would fall back to black on the dark ground.
const FG = "#F8F8F2";
// MathJax 3 has no HTML colour model, so the hex values are given as RGB:
// #DC267F and #648FFF.
const colours = String.raw`
  \definecolor{ee}{RGB}{220,38,127}
  \definecolor{ec}{RGB}{100,143,255}
`;

// Pair sets follow the engine's tie rules (R/engine.R): i is the subject with
// the event; two events at the same time are not comparable, an event and a
// censoring at the same time are.
const formulas = {
  ee: String.raw`
    \begin{aligned}
      \textcolor{ee}{C_{ee}} &= \frac{\sum_{(i,j)\,\in\,\textcolor{ee}{\mathcal{P}_{ee}}} w_{ij}\,c_{ij}}{\textcolor{ee}{W_{ee}}} \\[6pt]
      \textcolor{ee}{W_{ee}} &= {\textstyle\sum}_{(i,j)\,\in\,\textcolor{ee}{\mathcal{P}_{ee}}} w_{ij} \\[6pt]
      \textcolor{ee}{\mathcal{P}_{ee}} &= \{(i,j) : T_i < T_j,\ \delta_i = 1,\ \delta_j = 1\}
    \end{aligned}`,
  ec: String.raw`
    \begin{aligned}
      \textcolor{ec}{C_{ec}} &= \frac{\sum_{(i,j)\,\in\,\textcolor{ec}{\mathcal{P}_{ec}}} w_{ij}\,c_{ij}}{\textcolor{ec}{W_{ec}}} \\[6pt]
      \textcolor{ec}{W_{ec}} &= {\textstyle\sum}_{(i,j)\,\in\,\textcolor{ec}{\mathcal{P}_{ec}}} w_{ij} \\[6pt]
      \textcolor{ec}{\mathcal{P}_{ec}} &= \{(i,j) : T_i \le T_j,\ \delta_i = 1,\ \delta_j = 0\}
    \end{aligned}`,
  // Shares the identity slide with `decomposition`; typeset here too so the
  // two formulas match in font and scale.
  pooled: String.raw`
    C_w = \frac{\sum_{ij} w_{ij}\,c_{ij}}{\sum_{ij} w_{ij}}`,
  decomposition: String.raw`
    C_w = \frac{\textcolor{ee}{W_{ee}\,C_{ee}} \;+\; \textcolor{ec}{W_{ec}\,C_{ec}}}{\textcolor{ee}{W_{ee}} + \textcolor{ec}{W_{ec}}}`,
};

const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const doc = mathjax.document("", {
  InputJax: new TeX({ packages: ["base", "ams", "color"] }),
  // fontCache "none" writes each glyph inline instead of via <use>, which
  // PowerPoint's SVG import does not always resolve.
  OutputJax: new SVG({ fontCache: "none" }),
});

// MathJax's boxes are tight to the ink, and calligraphic P overhangs its box.
const PAD = 80; // viewBox units; 1000 = 1em

for (const [name, tex] of Object.entries(formulas)) {
  const node = doc.convert(colours + tex, { display: true, em: 16, ex: 8 });
  let svg = adaptor.innerHTML(node);

  const [x, y, w, h] = svg.match(/viewBox="([^"]+)"/)[1].split(/\s+/).map(Number);
  const vb = [x - PAD, y - PAD, w + 2 * PAD, h + 2 * PAD];
  // Intrinsic size at 1em = 32px, so the file opens at a legible size.
  const px = (v) => ((v / 1000) * 32).toFixed(1);

  svg = svg
    .replace(/ style="[^"]*"/, "")
    .replace(/ role="img"/, "")
    .replace(/ focusable="false"/, "")
    .replace(/viewBox="[^"]+"/, `viewBox="${vb.join(" ")}"`)
    .replace(/width="[^"]+"/, `width="${px(vb[2])}"`)
    .replace(/height="[^"]+"/, `height="${px(vb[3])}"`)
    .replace(/currentColor/g, FG);

  const out = path.join(__dirname, `${name}.svg`);
  fs.writeFileSync(out, `<?xml version="1.0" encoding="UTF-8"?>\n${svg}\n`);
  console.log(`${name}.svg  ${(vb[2] / 1000).toFixed(2)}em x ${(vb[3] / 1000).toFixed(2)}em`);
}
