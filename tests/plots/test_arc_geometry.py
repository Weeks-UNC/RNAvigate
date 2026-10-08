"""Geometry tests for arc plot arcs (``max_arc_height`` option).

Nesting and crossing are checked on the exact half-pill curves that
``get_arc_shape`` samples, so polygon chords do not affect the comparisons.
"""

from __future__ import annotations

import matplotlib.patches as mp_patches
import matplotlib.pyplot as plt
import numpy as np
import pytest

from rnavigate import plots

MAX_ARC_HEIGHTS = [None, 10, 50, 150]
RNG_SEED = 0


def _edge(x, left, right, height):
    """Exact height of a half-pill outline at x (0 outside [left, right])."""
    x = np.asarray(x, dtype=float)
    t = np.clip(np.minimum(x - left, right - x), 0, None)
    curve = np.sqrt(np.clip(2 * height * t - t**2, 0, None))
    return np.where(t < height, curve, height)


def _outer_edge(x, i, j, max_arc_height):
    height = plots.get_arc_height(j - i + 1, max_arc_height)
    return _edge(x, i - 0.5, j + 0.5, height)


def _inner_edge(x, i, j, max_arc_height):
    height = plots.get_arc_height(j - i - 1, max_arc_height)
    return _edge(x, i + 0.5, j - 0.5, height)


def test_arc_height_without_max_arc_height_is_semicircle():
    for span in [1, 2, 7, 100, 1726]:
        assert plots.get_arc_height(span, None) == span / 2


def test_arc_height_approaches_max_arc_height():
    heights = [plots.get_arc_height(span, 150) for span in [10, 300, 3000, 1e6]]
    assert all(np.diff(heights) > 0)
    assert all(h < 150 for h in heights)
    assert heights[-1] == pytest.approx(150, rel=1e-3)


@pytest.mark.parametrize("max_arc_height", MAX_ARC_HEIGHTS)
def test_shape_vertices_lie_on_exact_edges(max_arc_height):
    for i, j in [(1, 2), (5, 30), (116, 1842)]:
        outer, inner = plots.get_arc_shape(i, j, max_arc_height)
        np.testing.assert_allclose(
            outer[:, 1], _outer_edge(outer[:, 0], i, j, max_arc_height), atol=1e-9
        )
        np.testing.assert_allclose(
            inner[:, 1], _inner_edge(inner[:, 0], i, j, max_arc_height), atol=1e-9
        )


@pytest.mark.parametrize("max_arc_height", MAX_ARC_HEIGHTS)
def test_stacked_pairs_tile_exactly(max_arc_height):
    for i, j in [(1, 4), (10, 60), (116, 1842), (850, 1720)]:
        _, inner = plots.get_arc_shape(i, j, max_arc_height)
        outer_next, _ = plots.get_arc_shape(i + 1, j - 1, max_arc_height)
        assert np.array_equal(inner, outer_next)


@pytest.mark.parametrize("max_arc_height", MAX_ARC_HEIGHTS)
def test_nested_arcs_never_overlap(max_arc_height):
    rng = np.random.default_rng(RNG_SEED)
    for _ in range(500):
        i, j = sorted(rng.choice(2000, size=2, replace=False))
        if j - i < 3:
            continue
        i2, j2 = sorted(rng.choice(np.arange(i + 1, j), size=2, replace=False))
        x = np.linspace(i2 - 0.5, j2 + 0.5, 2001)
        inner_arc_top = _outer_edge(x, i2, j2, max_arc_height)
        outer_arc_bottom = _inner_edge(x, i, j, max_arc_height)
        assert np.all(inner_arc_top <= outer_arc_bottom + 1e-9), (i, j, i2, j2)


@pytest.mark.parametrize("max_arc_height", MAX_ARC_HEIGHTS)
def test_crossing_arcs_cross(max_arc_height):
    """For i < i2 < j < j2, arc (i2, j2) passes from below to above (i, j)."""
    rng = np.random.default_rng(RNG_SEED)
    for _ in range(500):
        i, i2, j, j2 = sorted(rng.choice(2000, size=4, replace=False))
        if i2 - i < 2 or j2 - j < 2:
            continue  # legs are adjacent, so neither arc is fully above
        x_below, x_above = i2 - 0.5, j + 0.5
        assert _outer_edge(x_below, i2, j2, max_arc_height) < _inner_edge(
            x_below, i, j, max_arc_height
        )
        assert _inner_edge(x_above, i2, j2, max_arc_height) > _outer_edge(
            x_above, i, j, max_arc_height
        )


def test_shape_without_max_arc_height_matches_wedge():
    for i, j in [(1, 2), (5, 30), (116, 1842)]:
        wedge = mp_patches.Wedge(((i + j) / 2, 0), 0.5 + (j - i) / 2, 0, 180, width=1)
        outer, inner = plots.get_arc_shape(i, j, None)
        center = np.array(wedge.center)
        np.testing.assert_allclose(np.hypot(*(outer - center).T), wedge.r)
        np.testing.assert_allclose(
            np.hypot(*(inner - center).T), wedge.r - wedge.width, atol=1e-12
        )


@pytest.mark.parametrize("panel", ["top", "bottom"])
def test_plot_interactions_arcs_default_is_unchanged(tpp, panel):
    """Without max_arc_height, arcs are the same Wedges as before the option existed."""
    interactions = tpp.data["ss"].as_interactions()
    _, ax = plt.subplots()
    plots.plot_interactions_arcs(ax, interactions, panel=panel, yvalue=3)
    theta1, theta2 = (0, 180) if panel == "top" else (180, 360)
    expected = []
    for i, j, _ in zip(*interactions.get_ij_colors()):
        i, j = sorted((i, j))
        wedge = mp_patches.Wedge(
            ((i + j) / 2.0, 3), 0.5 + (j - i) / 2.0, theta1, theta2, width=1
        )
        expected.append(wedge.get_transform().transform_path(wedge.get_path()))
    paths = ax.collections[-1].get_paths()
    assert len(paths) == len(expected)
    for path, expected_path in zip(paths, expected):
        np.testing.assert_array_equal(path.vertices, expected_path.vertices)


@pytest.mark.parametrize("panel", ["top", "bottom"])
def test_plot_interactions_arcs_follow_panel(tpp, panel):
    interactions = tpp.data["ss"].as_interactions()
    _, ax = plt.subplots()
    plots.plot_interactions_arcs(
        ax, interactions, panel=panel, yvalue=3, max_arc_height=150
    )
    paths = ax.collections[-1].get_paths()
    assert len(paths) == len(interactions.get_ij_colors()[0])
    ys = np.concatenate([path.vertices[:, 1] for path in paths]) - 3
    assert np.all(ys >= 0) if panel == "top" else np.all(ys <= 0)


def test_ap_ylim_follows_max_arc_height():
    semicircle = plots.AP(1, nt_length=1845)
    pill = plots.AP(1, nt_length=1845, max_arc_height=150)
    for plot in (semicircle, pill):
        plot.set_axis(plot.axes[0, 0], sequence="A" * 1845, yticks=[], ylabels=[])
    assert semicircle.axes[0, 0].get_ylim() == (-301, 301)
    expected = plots.get_arc_height(1845, 150) + 1
    np.testing.assert_allclose(pill.axes[0, 0].get_ylim(), (-expected, expected))
