# clinical-km-survival-analysis

A Claude skill for **Kaplan–Meier survival analysis** of clinical/medical time-to-event data. Produces NEJM-style survival curves and the descriptive statistics that accompany them.

This skill is **descriptive / unadjusted only**. For multivariable adjustment (adjusted HR, VIF, proportional-hazards diagnostics), use the companion skill [clinical-cox-regression](https://github.com/JeonKH81/clinical-cox-regression).

## What it does

- Single- and multi-group **Kaplan–Meier curves** (NEJM color style, number-at-risk table)
- **Log-rank test** for group comparison
- **Median survival** with 95% CI
- **Number-at-risk** table beneath the curve
- **Unadjusted HR** for the group variable via univariate Cox

## When it triggers

When an Excel file is provided together with terms like `생존분석`, `KM curve`, `Kaplan-Meier`, `생존곡선`, `log-rank`, `median survival`.

For `다변량 보정 / adjusted HR / VIF / Schoenfeld / PH 가정`, the request is routed to `clinical-cox-regression` instead.

## Usage

```bash
pip install lifelines openpyxl matplotlib pandas numpy --quiet --break-system-packages

python scripts/km_analysis.py \
  --file <filepath> --time <time_col> --event <event_col> \
  --group <group_col> --interval <30|365> --output result \
  --labels 0:Control 1:Treatment
```

**Arguments**

| Flag | Meaning |
|------|---------|
| `--file` | Excel file path |
| `--time` | Time-to-event column (days) |
| `--event` | Event column (1 = event, 0 = censored) |
| `--group` | Group column for comparison (optional) |
| `--interval` | X-axis tick interval in days (`30` or `365`) |
| `--output` | Output file prefix |
| `--labels` | Group labels as `value:label` (e.g. `0:non-DM 1:DM`, optional) |

## Output

| Group | N | Events | Median Survival (95% CI) |
|-------|---|--------|--------------------------|
| Control | 150 | 45 | 365 days (280–450) |
| Treatment | 148 | 32 | 520 days (410–NR) |

**Log-rank test**: χ² = 4.52, p = 0.033
**Unadjusted HR**: 0.68 (95% CI: 0.48–0.96), p = 0.028

## Installation as a Claude skill

Place this repository's contents (the `SKILL.md` and `scripts/`) into your Claude skills directory, or import it through the Claude Desktop / Claude Code skills interface.

## License

MIT © 2026 Kihyun Jeon
