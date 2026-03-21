#!/usr/bin/env python3
"""
plot_a247975.py — Logarithmic scatterplot of OEIS sequence A247975
a(n) = least positive integer m such that (m+n) | prime(m)^2 + prime(n)^2.

Generates: a247975_2.png
  Left y-axis:  A247975(n) — plain integers for 10, 1000, 100000;
                scientific notation (1e+07 etc.) for larger values.
                Labels rotated 90 deg CCW (parallel to axis).
  Right y-axis: log(A247975(n)), ticks at even integers 0,2,4,...,12.
  Resolved cases: blue dots. Unresolved (A247975(n)>2e11): red triangles.

Usage:
    python3 plot_a247975.py

Input files (same directory, or adjust paths below):
    phase1_results.csv   — Phase 1 (n=1..120000, bound 2e9)
    ext1_results.csv     — Extended search S1
    ext2_results.csv     — Extended search S2

Output:
    a247975_2.png        — PNG for OEIS upload (~512 KB, 1492x977 px)

Author: Carlo Corti, Mar 2026.
Code and data: https://github.com/carcorti/A247975
"""

import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE1_CSV = os.path.join(SCRIPT_DIR, 'phase1_results.csv')
EXT1_CSV   = os.path.join(SCRIPT_DIR, 'ext1_results.csv')
EXT2_CSV   = os.path.join(SCRIPT_DIR, 'ext2_results.csv')
OUTPUT_PNG = os.path.join(SCRIPT_DIR, 'a247975_2.png')

# ── Load data ──────────────────────────────────────────────────────────────
def load_phase1(path):
    rows = {}
    with open(path) as f:
        for _ in range(3): next(f)
        for row in csv.DictReader(f, skipinitialspace=True):
            m = int(row['m'])
            rows[m] = None if row['a_m'] == 'UNRESOLVED' else int(row['a_m'])
    return rows

def load_ext(path):
    rows = {}
    with open(path) as f:
        for _ in range(4): next(f)
        for row in csv.DictReader(f, skipinitialspace=True):
            m = int(row['m'])
            if row['a_m'] != 'UNRESOLVED':
                rows[m] = int(row['a_m'])
    return rows

data = load_phase1(PHASE1_CSV)
for path in [EXT1_CSV, EXT2_CSV]:
    for m, v in load_ext(path).items():
        data[m] = v

n_resolved   = sorted((m, v) for m, v in data.items() if v is not None)
n_unresolved = sorted(m for m, v in data.items() if v is None)
ns = np.array([x[0] for x in n_resolved])
vs = np.array([x[1] for x in n_resolved])

print(f"Resolved:   {len(ns)}")
print(f"Unresolved: {len(n_unresolved)} — {n_unresolved}")

# ── Figure ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6.5), dpi=150)

ax.scatter(ns, vs, s=0.5, c='#2166ac', alpha=0.45,
           linewidths=0, rasterized=True)

y_data_max = float(vs.max())
y_top_axis = y_data_max * 80
y_marker   = y_data_max * 18

for m in n_unresolved:
    ax.scatter(m, y_marker, marker='v', s=70,
               c='#d6604d', zorder=5, linewidths=0)

# ── Left y-axis ────────────────────────────────────────────────────────────
ax.set_yscale('log')
ax.set_xlim(0, 121000)
ax.set_ylim(0.8, y_top_axis)
ax.set_ylabel('A247975(n)', fontsize=12, rotation=90, labelpad=10)
ax.set_xlabel('n', fontsize=12)

# Automatic log ticks fall at: 1e+01, 1e+03, 1e+05, 1e+07, 1e+09, 1e+11
# First three shown as plain integers, rest in scientific notation
plain_vals = {1e1: '10', 1e3: '1000', 1e5: '100\u202f000'}

def left_fmt(y, pos):
    if y <= 0:
        return ''
    for val, label in plain_vals.items():
        if abs(y - val) / val < 1e-6:
            return label
    exp = int(np.floor(np.log10(y)))
    mant = y / 10**exp
    if abs(mant - 1.0) < 1e-9:
        return f'1e+{exp:02d}'
    return f'{mant:.0f}e+{exp:02d}'

ax.yaxis.set_major_formatter(ticker.FuncFormatter(left_fmt))
ax.tick_params(axis='y', which='both', labelrotation=90, labelsize=8)
ax.tick_params(axis='x', which='both', labelsize=8)

# ── Right y-axis ───────────────────────────────────────────────────────────
ax2 = ax.twinx()
ax2.set_yscale('log')
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks(10.0 ** np.arange(0, 13, 2))
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(
    lambda y, _: f'{int(round(np.log10(y)))}' if y >= 1 else ''))
ax2.set_ylabel('log(A247975(n))', fontsize=12, rotation=90, labelpad=10)
ax2.tick_params(axis='y', which='major', labelrotation=90, labelsize=8)
ax2.yaxis.set_minor_locator(ticker.NullLocator())

# ── x-axis ─────────────────────────────────────────────────────────────────
ax.xaxis.set_major_locator(ticker.MultipleLocator(20000))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(
    lambda x, _: f'{int(x):,}'.replace(',', '\u202f')))
ax.grid(True, which='major', linestyle='--', linewidth=0.4, alpha=0.5)
ax.grid(True, which='minor', linestyle=':', linewidth=0.2, alpha=0.3)

ax.set_title(
    'A247975: a(n) = least m such that (m+n) | prime(m)\u00b2 + prime(n)\u00b2\n'
    'Logarithmic scatterplot,  n = 1\u202f..\u202f120\u202f000',
    fontsize=11, pad=10)

# ── Annotate top-3 resolved Cramér records ─────────────────────────────────
highlights = [
    (23995, 47_572_503_750, 'n=23\u202f995',  42000, 8.0e10),
    (22802, 44_000_698_448, 'n=22\u202f802',   6000, 8.0e10),
    (82300, 81_961_045_941, 'n=82\u202f300',  96000, 4.0e10),
]
for m_h, v_h, label, tx, ty in highlights:
    ax.annotate(label, xy=(m_h, v_h), xytext=(tx, ty),
                fontsize=7.5, color='#222222', ha='center',
                arrowprops=dict(arrowstyle='->', color='#888888',
                                lw=0.7, shrinkB=3))

# ── Legend ─────────────────────────────────────────────────────────────────
h_blue = mlines.Line2D([], [], color='#2166ac', marker='o',
                       linestyle='None', markersize=5, alpha=0.8,
                       label=f'A247975(n) resolved  ({len(ns):,} values)')
h_red  = mlines.Line2D([], [], color='#d6604d', marker='v',
                       linestyle='None', markersize=8,
                       label=f'A247975(n) > 2\u00d710\u00b9\u00b9  ({len(n_unresolved)} unresolved)')
ax.legend(handles=[h_blue, h_red], loc='lower right', fontsize=9,
          framealpha=0.88, borderpad=0.8)

ax.annotate('A247975(n) > 2\u00d710\u00b9\u00b9',
            xy=(121000, y_marker), xytext=(121800, y_marker),
            fontsize=7, color='#d6604d', va='center', annotation_clip=False)

fig.text(0.99, 0.005, 'Carlo Corti, 2026 \u2014 github.com/carcorti/A247975',
         ha='right', va='bottom', fontsize=7, color='#999999')

plt.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=150, bbox_inches='tight',
            metadata={'Author': 'Carlo Corti',
                      'Title':  'Logarithmic scatterplot of A247975 (n=1..120000)'})
plt.close(fig)
print(f"Saved: {OUTPUT_PNG}")
