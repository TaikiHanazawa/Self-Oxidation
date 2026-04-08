#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path
from typing import Iterable

CI_MODEL_PPM = {
    'C': 35000.0,
    'N': 1500.0,
    'H': 6900.0,
}

CURRENT_BSE_RANGE_PPM = {
    'C': (42.0, 730.0),
    'N': (0.83, 2.5),
    'H': (44.0, 450.0),
}

CURRENT_BSE_MEAN_PPM = {
    key: 0.5 * (bounds[0] + bounds[1]) for key, bounds in CURRENT_BSE_RANGE_PPM.items()
}

SAKURABA_MAIN_TARGETS = [0.10, 0.30, 0.50, 0.70, 0.995]
SAKURABA_LATE_TARGETS = [0.995, 1.0]


def compile_makeearth(source: Path, binary: Path) -> None:
    subprocess.run([
        'c++', '-std=c++17', '-O2', str(source), '-o', str(binary)
    ], check=True)


def export_history(binary: Path, history_csv: Path) -> None:
    subprocess.run([
        str(binary), '--quiet', '--history-csv', str(history_csv)
    ], check=True)


def read_history(history_csv: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with history_csv.open(newline='') as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row: dict[str, object] = dict(raw)
            for key, value in raw.items():
                if key in {'label', 'phase', 'impactor_class'}:
                    continue
                row[key] = float(value)
            rows.append(row)
    return rows


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def nearest_row(rows: Iterable[dict[str, object]], target_mass: float) -> dict[str, object]:
    return min(
        rows,
        key=lambda row: (
            abs(float(row['planet_mass_after_earth']) - target_mass),
            int(float(row['step_index'])),
        ),
    )


def unique_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    seen: set[int] = set()
    unique: list[dict[str, object]] = []
    for row in rows:
        idx = int(float(row['step_index']))
        if idx in seen:
            continue
        seen.add(idx)
        unique.append(row)
    return unique


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>',
        'text { font-family: Georgia, "Times New Roman", serif; fill: #17202a; }',
        '.small { font-size: 12px; }',
        '.tick { font-size: 11px; }',
        '.label { font-size: 14px; }',
        '.title { font-size: 20px; font-weight: 600; }',
        '.panel { font-size: 18px; font-weight: 600; }',
        '</style>',
        '<rect width="100%" height="100%" fill="#fffdf8"/>',
    ]


def svg_footer(parts: list[str]) -> str:
    return '\n'.join(parts + ['</svg>'])


def map_linear(value: float, lo: float, hi: float, pixel_lo: float, pixel_hi: float) -> float:
    if hi == lo:
        return pixel_hi
    return pixel_lo + (value - lo) / (hi - lo) * (pixel_hi - pixel_lo)


def map_log10(value: float, lo: float, hi: float, pixel_lo: float, pixel_hi: float) -> float:
    if value <= 0.0:
        raise ValueError('log-scale values must be positive')
    log_lo = math.log10(lo)
    log_hi = math.log10(hi)
    log_value = math.log10(value)
    return pixel_lo + (log_value - log_lo) / (log_hi - log_lo) * (pixel_hi - pixel_lo)


def polyline(points: list[tuple[float, float]], stroke: str, width: float = 2.0,
             dasharray: str | None = None, fill: str = 'none') -> str:
    dash = f' stroke-dasharray="{dasharray}"' if dasharray else ''
    path = ' '.join(f'{x:.2f},{y:.2f}' for x, y in points)
    return (
        f'<polyline points="{path}" fill="{fill}" stroke="{stroke}" '
        f'stroke-width="{width}"{dash}/>'
    )


def circle(x: float, y: float, radius: float, fill: str, stroke: str = 'none') -> str:
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="{fill}" stroke="{stroke}"/>'


def text(x: float, y: float, content: str, klass: str = 'small', anchor: str = 'start') -> str:
    return f'<text x="{x:.2f}" y="{y:.2f}" class="{klass}" text-anchor="{anchor}">{content}</text>'


def rect(x: float, y: float, width: float, height: float, fill: str,
         stroke: str = 'none', stroke_width: float = 1.0, opacity: float = 1.0) -> str:
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" opacity="{opacity}"/>'
    )


def polygon(points: list[tuple[float, float]], fill: str, opacity: float = 1.0,
            stroke: str = 'none', stroke_width: float = 1.0) -> str:
    path = ' '.join(f'{x:.2f},{y:.2f}' for x, y in points)
    return (
        f'<polygon points="{path}" fill="{fill}" opacity="{opacity}" '
        f'stroke="{stroke}" stroke-width="{stroke_width}"/>'
    )


def format_mass_label(mass_earth: float) -> str:
    return f'{mass_earth * 100:.1f} wt%' if mass_earth < 1.0 else 'full accreted'


def render_sakuraba_fig2(rows: list[dict[str, object]], output_path: Path) -> None:
    width = 1080
    height = 620
    left = 90
    top = 90
    panel_w = 380
    panel_h = 360
    gap = 120
    ymin, ymax = 1.0e-4, 1.0e1
    categories = ['C', 'N', 'H']

    main_rows = [row for row in rows if float(row['planet_mass_after_earth']) <= 0.995 + 1.0e-12]
    late_rows = [row for row in rows if float(row['planet_mass_after_earth']) >= 0.995 - 1.0e-12]
    main_snapshots = unique_rows([nearest_row(main_rows, target) for target in SAKURABA_MAIN_TARGETS])
    late_snapshots = unique_rows([nearest_row(late_rows, target) for target in SAKURABA_LATE_TARGETS])

    main_palette = ['#8c8c8c', '#bdbdbd', '#f6b26b', '#f39c12', '#c0392b']
    late_palette = ['#f39c12', '#c0392b']

    parts = svg_header(width, height)
    parts.append(text(width / 2, 42, 'Sakuraba et al. (2021) Fig. 2 analogue from MakeEarth.cpp', 'title', 'middle'))
    parts.append(text(width / 2, 64, 'BSE is taken as mantle + atmosphere, normalized by total planetary mass and scaled by CI model abundances.', 'small', 'middle'))

    def draw_panel(panel_x: float, panel_y: float, panel_title: str,
                   snapshots: list[dict[str, object]], palette: list[str]) -> None:
        x_positions = [panel_x + 40 + i * (panel_w - 80) / (len(categories) - 1) for i in range(len(categories))]
        plot_left = panel_x
        plot_right = panel_x + panel_w
        plot_top = panel_y
        plot_bottom = panel_y + panel_h

        parts.append(text(panel_x - 24, panel_y - 12, panel_title, 'panel'))
        parts.append(rect(plot_left, plot_top, panel_w, panel_h, 'none', '#333333', 1.2))

        for tick_exp in range(-4, 2):
            tick_value = 10 ** tick_exp
            y = map_log10(tick_value, ymin, ymax, plot_bottom, plot_top)
            parts.append(f'<line x1="{plot_left:.2f}" y1="{y:.2f}" x2="{plot_right:.2f}" y2="{y:.2f}" stroke="#d0d0d0" stroke-width="1"/>')
            parts.append(text(plot_left - 10, y + 4, f'10^{tick_exp}', 'tick', 'end'))

        for x, label in zip(x_positions, categories):
            parts.append(f'<line x1="{x:.2f}" y1="{plot_top:.2f}" x2="{x:.2f}" y2="{plot_bottom:.2f}" stroke="#e1e1e1" stroke-width="1"/>')
            parts.append(text(x, plot_bottom + 24, label, 'label', 'middle'))

        upper_points = []
        lower_points = []
        mean_points = []
        for x, label in zip(x_positions, categories):
            lo_ppm, hi_ppm = CURRENT_BSE_RANGE_PPM[label]
            lo_ratio = lo_ppm / CI_MODEL_PPM[label]
            hi_ratio = hi_ppm / CI_MODEL_PPM[label]
            mean_ratio = CURRENT_BSE_MEAN_PPM[label] / CI_MODEL_PPM[label]
            upper_points.append((x, map_log10(hi_ratio, ymin, ymax, plot_bottom, plot_top)))
            lower_points.append((x, map_log10(lo_ratio, ymin, ymax, plot_bottom, plot_top)))
            mean_points.append((x, map_log10(mean_ratio, ymin, ymax, plot_bottom, plot_top)))
        band_polygon = upper_points + list(reversed(lower_points))
        parts.append(polygon(band_polygon, '#9fe2bf', opacity=0.75))
        parts.append(polyline(mean_points, '#2e8b57', 2.0, dasharray='6 4'))

        for snapshot, color in zip(snapshots, palette):
            points = []
            ratio_map = {
                'C': float(snapshot['bse_carbon_ppm']) / CI_MODEL_PPM['C'],
                'N': float(snapshot['bse_nitrogen_ppm']) / CI_MODEL_PPM['N'],
                'H': float(snapshot['bse_hydrogen_ppm']) / CI_MODEL_PPM['H'],
            }
            for x, label in zip(x_positions, categories):
                y = map_log10(ratio_map[label], ymin, ymax, plot_bottom, plot_top)
                points.append((x, y))
            parts.append(polyline(points, color, 2.4))
            for x, y in points:
                parts.append(circle(x, y, 4.0, color))

        legend_x = plot_right - 150
        legend_y = plot_bottom - 112
        parts.append(rect(legend_x - 12, legend_y - 22, 152, 20 + 18 * (len(snapshots) + 2), '#ffffff', '#cccccc', 1.0, 0.92))
        parts.append(polyline([(legend_x, legend_y), (legend_x + 20, legend_y)], '#2e8b57', 2.0, dasharray='6 4'))
        parts.append(text(legend_x + 28, legend_y + 4, 'BSE mean', 'small'))
        parts.append(rect(legend_x + 3, legend_y + 12, 14, 10, '#9fe2bf', opacity=0.75))
        parts.append(text(legend_x + 28, legend_y + 22, 'BSE range', 'small'))
        for idx, (snapshot, color) in enumerate(zip(snapshots, palette), start=2):
            y = legend_y + idx * 18
            parts.append(polyline([(legend_x, y), (legend_x + 20, y)], color, 2.4))
            parts.append(circle(legend_x + 10, y, 3.6, color))
            parts.append(text(legend_x + 28, y + 4, format_mass_label(float(snapshot['planet_mass_after_earth'])), 'small'))

        parts.append(text(plot_left - 60, plot_top + panel_h / 2, 'Abundance / CI chondrites', 'label'))

    draw_panel(left, top + 40, 'a', main_snapshots, main_palette)
    draw_panel(left + panel_w + gap, top + 40, 'b', late_snapshots, late_palette)

    if len(late_snapshots) < 3:
        parts.append(text(width / 2, height - 18, 'Late-accretion panel shows only start/end snapshots because the current scenario treats late veneer as a single step.', 'small', 'middle'))

    output_path.write_text(svg_footer(parts))


def render_shi_fig4b(rows: list[dict[str, object]], output_path: Path) -> None:
    width = 1120
    height = 620
    left = 110
    right = width - 80
    top = 90
    bottom = height - 90
    xmin, xmax = 0.01, 1.0
    ymin, ymax = -35.0, 5.0

    parts = svg_header(width, height)
    parts.append(text(width / 2, 42, 'Shi et al. (2022) Fig. 4b analogue from MakeEarth.cpp', 'title', 'middle'))
    parts.append(text(width / 2, 64, 'Atmosphere, mantle, and core nitrogen isotopic composition through accretion.', 'small', 'middle'))

    stage_bounds = [
        (0.01, 0.80, '#f8dfc5', 'EC-like growth'),
        (0.80, 0.90, '#d5e8ff', 'GI8 CC event'),
        (0.90, 1.00, '#dff3d6', 'Final assembly'),
    ]

    def x_map(value: float) -> float:
        return map_linear(value, xmin, xmax, left, right)

    def y_map(value: float) -> float:
        return map_linear(value, ymin, ymax, bottom, top)

    for x0, x1, color, label in stage_bounds:
        px0 = x_map(x0)
        px1 = x_map(x1)
        parts.append(rect(px0, top, px1 - px0, bottom - top, color, opacity=0.75))
        parts.append(text(0.5 * (px0 + px1), top - 10, label, 'small', 'middle'))

    parts.append(rect(left, top, right - left, bottom - top, 'none', '#333333', 1.2))

    for tick in [-35, -30, -25, -20, -15, -10, -5, 0, 5]:
        y = y_map(tick)
        parts.append(f'<line x1="{left:.2f}" y1="{y:.2f}" x2="{right:.2f}" y2="{y:.2f}" stroke="#d0d0d0" stroke-width="1"/>')
        parts.append(text(left - 10, y + 4, f'{tick:g}', 'tick', 'end'))

    for tick in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        x = x_map(tick)
        parts.append(f'<line x1="{x:.2f}" y1="{top:.2f}" x2="{x:.2f}" y2="{bottom:.2f}" stroke="#e1e1e1" stroke-width="1"/>')
        parts.append(text(x, bottom + 24, f'{tick:.1f}', 'tick', 'middle'))

    series = [
        ('Atmosphere', 'atmosphere_delta15n_permil', '#e67e22', '8 5'),
        ('Mantle', 'mantle_delta15n_permil', '#c0392b', None),
        ('Core', 'core_delta15n_permil', '#7b241c', '12 4 3 4'),
    ]

    for label, key, color, dash in series:
        points = [(x_map(float(row['planet_mass_after_earth'])), y_map(float(row[key]))) for row in rows]
        parts.append(polyline(points, color, 2.6, dasharray=dash))

    legend_x = left + 25
    legend_y = top + 25
    parts.append(rect(legend_x - 12, legend_y - 20, 164, 74, '#ffffff', '#cccccc', 1.0, 0.92))
    for idx, (label, _key, color, dash) in enumerate(series):
        y = legend_y + idx * 18
        parts.append(polyline([(legend_x, y), (legend_x + 26, y)], color, 2.6, dasharray=dash))
        parts.append(text(legend_x + 36, y + 4, label, 'small'))

    final_row = rows[-1]
    parts.append(text(right - 10, y_map(float(final_row['atmosphere_delta15n_permil'])) - 8, 'Atmosphere', 'small', 'end'))
    parts.append(text(right - 10, y_map(float(final_row['mantle_delta15n_permil'])) - 8, 'Mantle', 'small', 'end'))
    parts.append(text(right - 10, y_map(float(final_row['core_delta15n_permil'])) - 8, 'Core', 'small', 'end'))

    parts.append(text((left + right) / 2, bottom + 52, 'Mass of Earth accreted (M_E)', 'label', 'middle'))
    parts.append(text(left - 72, (top + bottom) / 2, 'δ15N (‰)', 'label'))
    output_path.write_text(svg_footer(parts))


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description='Generate Sakuraba Fig. 2 and Shi Fig. 4b analogue figures from MakeEarth.cpp output.')
    parser.add_argument('--source', type=Path, default=script_dir / 'MakeEarth.cpp')
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--history-csv', type=Path)
    parser.add_argument('--out-dir', type=Path, default=script_dir / 'benchmark_figures')
    parser.add_argument('--skip-build', action='store_true')
    parser.add_argument('--skip-run', action='store_true')
    args = parser.parse_args()

    ensure_dir(args.out_dir)
    if args.binary is None:
        args.binary = args.out_dir / 'MakeEarth_benchmark_plot'
    history_csv = args.history_csv or (args.out_dir / 'makeearth_history.csv')

    if not args.skip_build:
        compile_makeearth(args.source, args.binary)
    if not args.skip_run:
        export_history(args.binary, history_csv)

    rows = read_history(history_csv)
    render_sakuraba_fig2(rows, args.out_dir / 'sakuraba_fig2.svg')
    render_shi_fig4b(rows, args.out_dir / 'shi_fig4b.svg')

    print(f'history_csv={history_csv}')
    print(f'sakuraba_fig2_svg={args.out_dir / "sakuraba_fig2.svg"}')
    print(f'shi_fig4b_svg={args.out_dir / "shi_fig4b.svg"}')


if __name__ == '__main__':
    main()
