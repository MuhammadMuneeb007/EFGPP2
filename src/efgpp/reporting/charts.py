"""Plotly figure builders for the data report.

Colour follows the job: magnitude uses one sequential hue (blue); identity uses the
categorical slots in fixed order (never cycled); pass/fail uses the reserved status
colours and always carries a text label. Dark mode swaps to the dark steps client-side.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import plotly.graph_objects as go

CATEGORICAL_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
CATEGORICAL_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]
SEQUENTIAL_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
MAGNITUDE = CATEGORICAL_LIGHT[0]
SCATTER_MAX_SERIES = 3  # all-pairs CVD safety holds for the first three slots only

INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7"}
FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'


def _layout(fig: go.Figure, title: str, *, height: int = 320, x: str | None = None, y: str | None = None,
            showlegend: bool = False) -> go.Figure:
    fig.update_layout(
        title={"text": title, "x": 0, "xanchor": "left", "font": {"size": 14, "color": INK["primary"]}},
        height=height, margin={"l": 56, "r": 16, "t": 44, "b": 48},
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font={"family": FONT, "size": 12, "color": INK["secondary"]},
        showlegend=showlegend, bargap=0.25, hoverlabel={"font": {"family": FONT}},
        legend={"orientation": "h", "y": -0.2, "x": 0},
    )
    axis = {"gridcolor": INK["grid"], "linecolor": INK["axis"], "zerolinecolor": INK["axis"],
            "tickfont": {"color": INK["muted"]}, "title": {"font": {"color": INK["secondary"]}}}
    fig.update_xaxes(**axis, title_text=x, showgrid=False)
    fig.update_yaxes(**axis, title_text=y)
    return fig


def bar(labels: Sequence[str], values: Sequence[float], title: str, *, x: str | None = None,
        y: str | None = None, horizontal: bool = False, hover_unit: str = "") -> go.Figure:
    """Magnitude bars in one hue."""
    kw = {"x": list(values), "y": list(labels), "orientation": "h"} if horizontal else {"x": list(labels), "y": list(values)}
    fig = go.Figure(go.Bar(**kw, marker={"color": MAGNITUDE, "line": {"width": 0}, "cornerradius": 4},
                           hovertemplate=("%{y}: %{x:,}" if horizontal else "%{x}: %{y:,}") + hover_unit + "<extra></extra>"))
    height = max(260, 28 * len(labels) + 90) if horizontal else 320
    if horizontal:
        fig.update_yaxes(autorange="reversed")
    return _layout(fig, title, height=height, x=x, y=y)


def histogram(values: Sequence[float], title: str, *, x: str | None = None, bins: int = 40,
              threshold: float | None = None, threshold_label: str | None = None, log_y: bool = False) -> go.Figure:
    arr = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    fig = go.Figure(go.Histogram(x=arr, nbinsx=bins, marker={"color": MAGNITUDE, "line": {"width": 1, "color": "#fcfcfb"}},
                                 hovertemplate="%{x}: %{y:,} <extra></extra>"))
    if threshold is not None:
        fig.add_vline(x=threshold, line={"color": STATUS["critical"], "width": 2, "dash": "dot"},
                      annotation_text=threshold_label or f"threshold {threshold:g}",
                      annotation_font={"color": INK["secondary"]})
    if log_y:
        fig.update_yaxes(type="log")
    return _layout(fig, title, x=x, y="count")


def scatter(x: Sequence[float], y: Sequence[float], title: str, *, groups: Sequence[str] | None = None,
            xlabel: str = "", ylabel: str = "", hover: Sequence[str] | None = None) -> go.Figure:
    """Scatter; more than three groups fold the rest into 'Other' (all-pairs colour limit)."""
    fig = go.Figure()
    if groups is None:
        fig.add_trace(go.Scattergl(x=list(x), y=list(y), mode="markers", text=list(hover or []),
                                   marker={"size": 8, "color": MAGNITUDE, "opacity": 0.8,
                                           "line": {"width": 1, "color": "#fcfcfb"}},
                                   hovertemplate="%{text}<br>%{x:.3f}, %{y:.3f}<extra></extra>"))
        return _layout(fig, title, height=380, x=xlabel, y=ylabel)
    g = np.asarray(groups, dtype=object)
    counts = {k: int((g == k).sum()) for k in dict.fromkeys(g.tolist())}
    top = sorted(counts, key=lambda k: -counts[k])[:SCATTER_MAX_SERIES]
    labels = np.where(np.isin(g, top), g, "Other")
    order = [*top, *(["Other"] if (labels == "Other").any() else [])]
    xs, ys, hv = np.asarray(x), np.asarray(y), np.asarray(hover if hover is not None else [""] * len(g), dtype=object)
    for i, name in enumerate(order):
        m = labels == name
        color = INK["muted"] if name == "Other" else CATEGORICAL_LIGHT[i]
        fig.add_trace(go.Scattergl(x=xs[m], y=ys[m], mode="markers", name=f"{name} ({int(m.sum())})",
                                   text=hv[m], marker={"size": 8, "color": color, "opacity": 0.85,
                                                       "line": {"width": 1, "color": "#fcfcfb"}},
                                   hovertemplate="%{text}<br>%{x:.3f}, %{y:.3f}<extra>" + str(name) + "</extra>"))
    return _layout(fig, title, height=400, x=xlabel, y=ylabel, showlegend=True)


def status_bar(labels: Sequence[str], passed: Sequence[int], failed: Sequence[int], title: str) -> go.Figure:
    """Pass/fail counts: status colours with explicit labels in the legend and tooltips."""
    fig = go.Figure([
        go.Bar(name="✓ pass", x=list(labels), y=list(passed), marker={"color": STATUS["good"], "cornerradius": 4},
               hovertemplate="%{x}: %{y:,} pass<extra></extra>"),
        go.Bar(name="✗ fail", x=list(labels), y=list(failed), marker={"color": STATUS["critical"], "cornerradius": 4},
               hovertemplate="%{x}: %{y:,} fail<extra></extra>"),
    ])
    fig.update_layout(barmode="stack")
    fig.update_traces(marker_line_width=2, marker_line_color="#fcfcfb")
    return _layout(fig, title, showlegend=True)


def upset(patterns: list[tuple[tuple[str, ...], int]], sets: list[str], title: str, max_bars: int = 20) -> go.Figure:
    """UpSet-style chart: intersection sizes (bars) above a set-membership dot matrix."""
    from plotly.subplots import make_subplots

    patterns = sorted(patterns, key=lambda p: -p[1])[:max_bars]
    xs = list(range(len(patterns)))
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.62, 0.38], vertical_spacing=0.03)
    fig.add_trace(go.Bar(x=xs, y=[p[1] for p in patterns], marker={"color": MAGNITUDE, "cornerradius": 4},
                         customdata=[" + ".join(p[0]) or "(none)" for p in patterns],
                         hovertemplate="%{customdata}<br>%{y:,} participants<extra></extra>"), row=1, col=1)
    for j, s in enumerate(sets):
        members = [i for i, p in enumerate(patterns) if s in p[0]]
        fig.add_trace(go.Scatter(x=xs, y=[j] * len(xs), mode="markers", hoverinfo="skip",
                                 marker={"size": 10, "color": INK["grid"]}), row=2, col=1)
        fig.add_trace(go.Scatter(x=members, y=[j] * len(members), mode="markers", hoverinfo="skip",
                                 marker={"size": 10, "color": INK["primary"]}), row=2, col=1)
    for i, p in enumerate(patterns):
        ys = [sets.index(s) for s in p[0] if s in sets]
        if len(ys) > 1:
            fig.add_trace(go.Scatter(x=[i, i], y=[min(ys), max(ys)], mode="lines", hoverinfo="skip",
                                     line={"color": INK["primary"], "width": 2}), row=2, col=1)
    fig.update_yaxes(tickvals=list(range(len(sets))), ticktext=sets, row=2, col=1, showgrid=False)
    fig.update_xaxes(showticklabels=False)
    _layout(fig, title, height=300 + 22 * len(sets), y=None)
    fig.update_yaxes(title_text="participants", row=1, col=1)
    return fig


def heatmap(z: Sequence[Sequence[float]], labels: Sequence[str], title: str) -> go.Figure:
    fig = go.Figure(go.Heatmap(z=z, x=list(labels), y=list(labels), colorscale=[[i / 6, c] for i, c in enumerate(SEQUENTIAL_BLUE)],
                               xgap=2, ygap=2, hovertemplate="%{y} + %{x}: %{z:,}<extra></extra>",
                               colorbar={"thickness": 10, "outlinewidth": 0}))
    fig.update_yaxes(autorange="reversed")
    return _layout(fig, title, height=max(320, 30 * len(labels) + 120))


def to_html(fig: go.Figure, div_id: str) -> str:
    return fig.to_html(full_html=False, include_plotlyjs=False, div_id=div_id,
                       config={"displaylogo": False, "responsive": True})


def dark_mode_script() -> str:
    """Client-side swap of chart ink and categorical steps for the dark surface."""
    mapping = dict(zip(CATEGORICAL_LIGHT, CATEGORICAL_DARK, strict=True))
    mapping[INK["primary"]] = "#ffffff"
    mapping[INK["grid"]] = "#2c2c2a"
    mapping["#fcfcfb"] = "#1a1a19"
    import json

    return f"""
(function() {{
  const MAP = {json.dumps(mapping)};
  const swap = c => (typeof c === 'string' && MAP[c]) ? MAP[c] : c;
  function isDark() {{
    const t = document.documentElement.getAttribute('data-theme');
    if (t) return t === 'dark';
    return window.matchMedia('(prefers-color-scheme: dark)').matches;
  }}
  function apply() {{
    if (!isDark() || !window.Plotly) return;
    document.querySelectorAll('.plotly-graph-div').forEach(div => {{
      (div.data || []).forEach((tr, i) => {{
        const upd = {{}};
        if (tr.marker) {{
          if (tr.marker.color) upd['marker.color'] = [Array.isArray(tr.marker.color) ? tr.marker.color.map(swap) : swap(tr.marker.color)];
          if (tr.marker.line && tr.marker.line.color) upd['marker.line.color'] = [swap(tr.marker.line.color)];
        }}
        if (tr.line && tr.line.color) upd['line.color'] = [swap(tr.line.color)];
        if (Object.keys(upd).length) Plotly.restyle(div, upd, [i]);
      }});
      Plotly.relayout(div, {{'font.color': '#c3c2b7', 'title.font.color': '#ffffff',
        'xaxis.gridcolor': '#2c2c2a', 'yaxis.gridcolor': '#2c2c2a',
        'xaxis.linecolor': '#383835', 'yaxis.linecolor': '#383835'}});
    }});
  }}
  window.addEventListener('load', apply);
}})();
"""


def fmt(value: Any) -> str:
    if value is None:
        return "–"
    if isinstance(value, float):
        return f"{value:,.4g}"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)
