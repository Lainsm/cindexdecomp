# Symposium talk: "Check your C-index" (design)

**Date:** 2026-09-12
**Presented:** 2026-09-19
**Status:** Approved design, pending build
**Package version at time of writing:** 0.2.0 (development)

---

## 1. Purpose

A 15-minute talk (plus 5 minutes of questions) for a symposium of students
and early-career peers.

**Takeaway, in one sentence:** the pooled C-index averages two ranking tasks
of very different difficulty, and the mix between them is set by your
censoring rate rather than by your model, so check the split before you
trust the number.

The villain is the pooled C-index. `cindexdecomp` appears as the evidence
and the remedy, not as the subject. This is deliberately *not* a package
pitch: the room should leave suspicious of a statistic they report, not
merely aware that a package exists.

Relation to the publication path in
`2026-08-25-cindex-decomposition-design.md`: this talk is the argument of
the planned methods paper, delivered early to a friendly audience. Q&A is
therefore useful signal on which objections the paper must pre-empt.

## 2. Audience and constraints

| Constraint | Value |
|---|---|
| Audience | Students and early-career peers; mixed statistical depth |
| Slot | 15 min talk + 5 min Q&A |
| Slide budget | 18 counted slides, plus backups reachable only from Q&A |
| Toolchain | Quarto `revealjs` (`.qmd`) |
| Appearance | Dark ground (`#282A36`), package accent colours unchanged |
| Figures | Rendered live from packaged `survival` cohorts, seeded |
| House style | British English; no em dashes anywhere |

Because the audience is peers rather than methodologists, the estimator is
motivated before it is written down, and no slide opens with notation.

## 3. Narrative order and why

The order is **mechanism → evidence → machinery**:

1. Censoring, and what it does to the score.
2. Five cohorts.
3. The decomposition.

The alternative, opening with a hollow-looking result and revealing the
split as a twist, was rejected. It requires a cohort where the reveal
lands, and `survival::lung` is too lightly censored to supply one. Leading
with the mechanism inverts that weakness: once the room understands *why*
the pooled score drifts toward the easy task, `lung` stops being a weak
villain and becomes the **best case** in a ladder.

### 3.1 Evidence: five real cohorts, not a simulated sweep

The evidence is a ladder of five published cohorts ordered by censoring
rate, each with its own Cox model. All five ship with `survival`, so the
whole deck is reproducible with no data dependency beyond an existing
import.

The pattern that carries the talk, stated at exactly the strength the data
supports:

- **`C_ee` is flat and mediocre everywhere: 0.578 to 0.622.** Every one of
  these models is about equally poor at ordering real events. This is the
  load-bearing fact.
- **`C_global` rises monotonically with censoring**, 0.637 to 0.794. This is
  the only column that is strictly monotone across all five, and it is the
  only monotone claim to make out loud.
- The event-event share of pairs falls from 67.5% to 17.5% across the range,
  but **not monotonically** (`colon` 32.7% then `gbsg` 33.5%). Describe it as
  a collapse across the range, never as a trend row by row.
- The gap widens from 0.120 to 0.261 across the range, also **not
  monotonically** (`colon` 0.093 and `gbsg` 0.100 both sit below `lung`'s
  0.120).

So the pooled score tracks the censoring rate, not model quality. The
strongest single fact available: `flchain` reports a `C_global` of 0.794,
a number anyone would call a good model, on a `C_ee` of 0.578.

The two non-monotone columns are visible on the slide and a peer audience
will spot them. Slide 8 must therefore lead with the flat `C_ee` column and
the rising `C_global` column, and treat the gap as a derived consequence
rather than as the trend being asserted.

### 3.2 Recorded decision: the censoring curve is cut

`censoring_curve()` and its figure do not appear in the deck, not even as a
backup.

The known cost of this, recorded so it is not rediscovered under
questioning: the five cohorts differ in disease, model and predictor set,
so censoring is **confounded** with everything else. The ladder is
correlational, not a controlled comparison. The censoring curve was the
controlled version of the same claim: one cohort, one model, only the
cut-off moving. Cutting it removes the prepared answer to the most
likely objection ("different diseases and different covariates; of course
the numbers differ").

Mitigation, since the figure is gone: slide 9 must concede the confounding
explicitly rather than wait to be caught by it, and the presenter should be
ready to describe the within-cohort sweep verbally in Q&A. The honest
verbal answer is that `C_ee` being flat at 0.58 to 0.62 across five unrelated
diseases is itself evidence against a disease-specific explanation.

## 4. Slide-by-slide

18 counted slides. Build stages carry `visibility="uncounted"` so they do not
inflate the count or the numbering.

### Part 1. Censoring and what it does to the score (slides 1-6)

**1. Title.**

**2. The number we all report.** Every survival paper ends with a C-index.
Plant the question the talk answers: *which pairs?*

**3. The metric landscape.** Brier/IBS, time-dependent AUC, calibration slope
and CITL, Royston-Sauerbrei D. Scope-setting: this talk is about the C-index,
and the argument is not that C is wrong but that its pooled form hides its own
composition.

**4. Censoring and comparability.** A partial outcome, not a missing one.
Non-comparable pairs are dropped before any score is computed.

**5. Two pairs, two very different questions.** *Load-bearing slide; three
build stages as uncounted slides.* Ordering (both died) versus detection (one
died, one left alive), drawn as timelines in the package's own two accent
colours so the room associates magenta with event-event and blue with
event-censored long before `C_ee` and `C_ec` appear.

**6. The mix is set by your study, not your model.**

### Part 2. Five cohorts (section break, then slides 8-13)

**8. Five published cohorts, five Cox models.** Centred table: cohort,
setting, **years**, n, censoring. Years come from the data where the data
carries them (`rotterdam$year`, `flc$sample.yr`) and from the `survival`
documentation otherwise (lung 1994, colon 1989, gbsg 1984-1989). Caption:
"Unrelated diseases. Different covariates. Different decades."
The censoring ordering is deliberate but said **verbally**, not printed.

**9. How Harrell's C is computed.** Three steps: form every comparable pair;
score each 1 / 0.5 / 0; average with equal weight on every pair. Exists so
that `C_global` on the next slide is not a black box and so that `w = 1` in
the identity later lands as a *choice* rather than as notation. About 45
seconds.

**10. Read the global column first.** Table of cohort, censoring, `C_global`
only. Let the room read it as models improving. Followed by an uncounted
slide adding the `EE share`, `C_ee`, `C_ec` and gap columns.

**11. The pooled score tracks the censoring rate.** `C_global` alone against
censoring rate. Faint dotted connector in censoring order, not a fitted line:
`C_global` is strictly monotone across all five, so the connector traces an
ordering that is genuinely there. Climbs 0.637 to 0.794.

**12. The two halves, side by side.** Grouped bar chart, `C_ee` and `C_ec` per
cohort, cohorts in censoring order with their censoring rate in the axis
label. **Bars start at zero deliberately**; a baseline truncated at 0.5 would
exaggerate exactly the effect being claimed. "chance" rides a right-hand
secondary axis because every in-panel position collides with a bar.

**13. The same model, scored two ways.** `flchain`: 0.794 labelled
`C_global`, 0.578 labelled `C_ee`. Concede the confounding here.

### Part 3. The decomposition (section break, then slides 15-18)

**15. It is an identity, not an estimator.** Both display equations enlarged
(`.big-math`). Fits at 686px against a 700px slide, so the prose around it
cannot grow without pushing it over.

**16. You cannot escape this by...** Three moves: a better estimator, a better
model, and pooling being a convention rather than a law. The words **pools**
and **pooling** are set as filled chips, against italic *reweights*, because
that distinction is the point of the slide.

**17. What this does not do.** Limitations, volunteered before Q&A.

**18. Check your C-index.** A checklist of what to do to be right, each line
paired with the package call that does it: report the censoring rate
(printed in every result header), report `C_ee` (`fit$C_ee`), say whether the
gap is distinguishable from zero (`confint()`), show it is not an estimator
artefact (`weights_uno()`), compare models on the split
(`compare_decompositions()`).

### Backup slides (Q&A only)

- The ladder recomputed under Uno's weighting.
- Uno under tied event times: divergence of about 1e-04.
- The tied-times bug fixed in 0.2.0.
- The five model formulas.

## 5. Timing

| Part | Slides | Budget |
|---|---|---|
| 1. Censoring and the score | 1-6 | 5 min |
| 2. Five cohorts | 7-13 | 6 min |
| 3. The decomposition | 14-17 | 3 min |
| Close | 18 | 1 min |
| **Total** | **18** | **15 min** |

**This is now tight.** The deck grew from 13 counted slides to 18 without the
slot growing, which is roughly 50 seconds a slide against figure slides that
want 60 to 90.

Known redundancy, and the first place to cut: slides 10 and 11 both say
`C_global` rises with censoring, once as a table and once as a plot, and the
uncounted full table overlaps slide 12's bar chart. If a rehearsal runs long,
drop the two tables and keep the two plots (the plots were the point of
splitting this section), moving the full table to the backup pile. Slide 3 is
the next cut after that.

## 6. House style

- **British English throughout**, matching the package's `Language: en-GB`.
- **No em dashes.** Use a comma, colon or full stop instead. En dashes are
  retained only where they are not dashes-as-punctuation: joint-name pairs
  (Royston-Sauerbrei, Gonen-Heller, written with en dashes in the source).
  Numeric ranges are written as "0.578 to 0.622" rather than with a dash,
  because the deck is read aloud.

## 7. Build requirements

- Quarto `revealjs`, source at `dev/talks/2026-09-19-symposium.qmd`,
  following the existing `dev/lung-smoke-test.qmd` pattern.
- Render with `embed-resources: true`. This produces one self-contained
  `.html`, which `dev/.gitignore` (`*.html`) already ignores, avoids the
  `_files/` sidecar directory that `*.html` would *not* ignore, and
  satisfies the offline requirement below in the same stroke.
- `dev/` is already in `.Rbuildignore`, so nothing here reaches the package
  tarball.
- Every number and figure on a slide is computed at render time from
  `survival` cohorts. No transcribed values, matching the discipline
  already applied to `README.Rmd`.
- **Bootstrap chunks must be cached** (`cache: true`). Five cohorts at
  `n_boot = 1000` is minutes of compute, and rehearsal will re-render
  repeatedly. Seed every bootstrap so cached and fresh renders agree.
- Point estimates and pair counts are exact and reproduce deterministically;
  **interval bounds are seed-dependent** and move in the third decimal.
  Phrase gaps as "distinguishable from zero" rather than anything stronger,
  and let the rendered chunk supply the bounds rather than transcribing
  them.
- Slide 5 is built as a layered reveal: render the ggplot several times with
  successive layers added and stack the versions as reveal.js fragments, so
  the figure appears to draw itself. No animation dependency; degrades to a
  static final figure in the PDF export.
- Plots use `theme_cindex(dark = TRUE)` so the deck matches package output;
  check legibility at projector distance and raise base font size if needed.
- The deck's SCSS ground (`#282A36`) must stay identical to
  `cindex_palette(dark = TRUE)$bg`, so figures sit flush against the slide
  with no visible panel seam.
- **`theme_cindex()` never sets `strip.text`**, so a facetted plot inherits
  `theme_classic()`'s near-black `grey10` and the strip titles vanish against
  the dark ground. Slide 5 sets `strip.text` colour explicitly. (This is a
  latent gap in the package's dark variant, not a deck bug: the package's own
  `autoplot()` methods do not facet, so nothing there exposes it. Worth
  fixing upstream if `theme_cindex()` is ever used for facetted work.)
- **`auto-stretch: false` is required.** Quarto rewrites a slide's lone image
  as `.r-stretch` and sizes it to the leftover vertical space; on these
  slides that computed to height 0 and every figure silently vanished while
  the render still reported success. CSS caps image width instead.
- `code-line-numbers: false` and `code-overflow: wrap`: revealjs otherwise
  numbers a one-line snippet and draws a bright scrollbar across the dark
  slide.
- `navigation-mode: linear`: the `#` section breaks nest their slides as
  vertical stacks, so without this a right-arrow during the talk skips a
  whole section's contents.
- A footer carrying the canonical CRAN reference
  (`CRAN.R-project.org/package=cindexdecomp`) sits bottom-right on every
  slide, with the slide number moved to the bottom-left and clear of
  revealjs's menu button. **That URL is dead until CRAN accepts the package.**
  If presenting before then, swap the one `footer:` line for the live pkgdown
  vignette: `lainsm.github.io/cindexdecomp/articles/cindexdecomp.html`. Quarto ships its own `.reveal .footer` rule at
  equal specificity and later in the cascade, so the override needs
  `!important`.
- The deck must run offline from a local file, and a PDF export travels as a
  fallback in case the venue's display misbehaves. Verified working: append
  `?print-pdf` and print from the browser, which yields 23 pages with the
  dark theme, footer and every fragment stage preserved as its own page.
  Note the `?print-pdf` view looks blank on screen while initialising; that
  is expected and the printed output is correct.
- Slide 5's timeline figure is bespoke to the talk and is the one figure not
  produced by an exported package function.

## 8. Open questions

None blocking. To settle during rehearsal:

- Whether slide 6 needs a small figure showing pair composition shifting
  with censoring, or whether the sentence carries it. Build without; add
  only if rehearsal shows the point landing softly.
- Whether slide 8 reads better as a table or a slope chart. Build the table
  first; it is the safer object at projector distance.
- Whether slide 3 survives if rehearsal runs long. It is the first cut, but
  cutting it forfeits the scope-setting that keeps Q&A on topic.

## Appendix A: verified reference values

Computed 2026-09-12 against `survival` cohorts with `cindexdecomp` 0.2.0
(development), `set.seed(2026)`, `n_boot = 500`. Recorded so the built deck
can be checked against a known-good run, **not** for transcription onto
slides: the deck computes its own.

| cohort | n | cens | EE share | `C_ee` | `C_ec` | `C_global` | gap | gap 95% CI |
|---|---|---|---|---|---|---|---|---|
| lung | 167 | 28.1% | 67.5% | 0.598 | 0.717 | 0.637 | 0.120 | [0.015, 0.222] |
| colon | 888 | 51.6% | 32.7% | 0.599 | 0.692 | 0.661 | 0.093 | [0.050, 0.135] |
| gbsg | 686 | 56.4% | 33.5% | 0.622 | 0.721 | 0.688 | 0.100 | [0.054, 0.151] |
| rotterdam | 2982 | 57.3% | 31.0% | 0.588 | 0.740 | 0.693 | 0.151 | [0.129, 0.173] |
| flchain | 7871 | 72.5% | 17.5% | 0.578 | 0.840 | 0.794 | 0.261 | [0.246, 0.276] |

Models fitted:

- `lung`: `Surv(time, status - 1) ~ age + sex + ph.ecog`, complete cases.
- `colon`: `etype == 2` (death), complete cases on the covariates below;
  `Surv(time, status) ~ age + sex + obstruct + nodes + differ + extent + rx`.
- `gbsg`: `Surv(rfstime, status) ~ age + meno + size + grade + nodes + pgr +
  er + hormon`.
- `rotterdam`: `Surv(dtime, death) ~ age + meno + size + grade + nodes + pgr
  + er + hormon + chemo`.
- `flchain`: `futime > 0`, complete cases;
  `Surv(futime, death) ~ age + sex + kappa + lambda + mgus`.

Notes for the presenter:

- **Every gap excludes zero.** Interval width tracks sample size, as it
  should.
- **`lung` is the inferentially weakest row**: its lower bound is +0.015 on
  n = 167. That is fine for its role as the lightly censored anchor, but do
  not lean on `lung` alone if challenged: point at `rotterdam` and
  `flchain`, whose intervals are tight and nowhere near zero.
- **The gap is not monotone across all five.** `colon` (0.093) and `gbsg`
  (0.100) sit below `lung` (0.120). Do not claim a monotone trend in the
  gap. The monotone and defensible claims are that `C_global` rises with
  censoring and that the event-event share of pairs collapses; the
  load-bearing fact is that `C_ee` is flat at 0.578 to 0.622 across all five.
- Full `flchain` decomposes in about 0.4 s for a point estimate, so no
  subsampling is needed anywhere.
