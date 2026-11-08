"""Shared panel styling for the paper figures.

Layout requested for this project (2026-10-22)
----------------------------------------------
* the **figure-level title is not drawn on the canvas** -- it lives in
  `figures/main/FIGURE_LEGENDS.md`, so every panel gets the full canvas width;
* each panel is labelled with a bold **UPPERCASE** letter (A, B, C ...);
* when a panel carries a subtitle, that subtitle is **centred** over the axes.

`panel(ax, "a where the gate puts every gene")` draws a bold "A" flush with the
left edge of the axes and centres "where the gate puts every gene" above it.

Collision handling
------------------
Centring the subtitle on the whole axes is only safe while it is short: a long
subtitle reaches left and runs into the bold letter.  That is what happened in
Fig 2 E, Fig 6 C and Fig 8 H after the 1.32x font enlargement.

Fixing it inside `panel()` alone is **not enough**: `panel()` runs before
`tight_layout()` / `subplots_adjust()`, so the axes it measures are wider than
the axes that are finally drawn, and the fit test passes on a width that never
materialises (Fig 6 C still collided that way).

Therefore:

1. `panel()` records the letter and the title and tags both with a gid;
2. `fit_panel_titles(fig)` re-runs the measurement on the *final* geometry and
   must be called immediately before `fig.savefig(...)`.

`fit_panel_titles` needs no extra draw, so it is cheap.
"""
from __future__ import annotations

#: point gap kept between the bold letter and the start of the subtitle
GAP_PT = 5.0
#: gids used to find the two artists again at fit time
LETTER_GID = "panel_letter"
TITLE_GID = "panel_title"


def panel(ax, text, pad=6.0, letter_size=None, letter_x=0.0, **kw):
    """Label one panel: bold uppercase letter at the left edge + centred subtitle.

    Parameters
    ----------
    ax
        The matplotlib axes to label.
    text
        A string whose first character is the panel letter, e.g. ``"a where ..."``.
        Everything after the letter (and any separator characters) is the subtitle.
    pad
        Points between the top of the axes and the baseline of the text.
    letter_size
        Font size for the letter; defaults to ``fontsize + 2``.
    letter_x
        Axes-fraction x for the letter. 0.0 puts it flush with the left spine.
    **kw
        Forwarded to `Axes.set_title` (notably ``fontsize``).
    """
    text = str(text).strip()
    letter, subtitle = text[:1], text[1:].strip(" \t:-—|")
    fs = kw.get("fontsize", 9.0)
    if letter_size is None:
        letter_size = fs + 2.0

    ann = ax.annotate(letter.upper(), xy=(letter_x, 1.0), xycoords="axes fraction",
                      xytext=(0.0, pad), textcoords="offset points",
                      ha="left", va="bottom",
                      fontsize=letter_size, fontweight="bold",
                      annotation_clip=False)
    ann.set_gid(LETTER_GID)
    if not subtitle:
        return

    title = ax.set_title(subtitle, loc="center", pad=pad, **kw)
    title.set_gid(TITLE_GID)
    # a first attempt now keeps the common case correct even if the caller
    # forgets fit_panel_titles(); the final pass below is authoritative.
    _resolve(ax, ann, title)


def _resolve(ax, ann, title):
    """Centre `title` in the space that is free to the right of `ann`."""
    try:
        fig = ax.figure
        renderer = fig.canvas.get_renderer()
    except Exception:                       # e.g. a non-Agg canvas
        return

    gap_px = GAP_PT * fig.dpi / 72.0
    ax_bb = ax.get_window_extent(renderer)

    for _ in range(10):
        bb_l = ann.get_window_extent(renderer)
        bb_t = title.get_window_extent(renderer)
        free_left = bb_l.x1 + gap_px
        free_w = ax_bb.x1 - free_left
        if bb_t.x0 >= free_left:
            return                          # already clear
        if free_w <= 4.0:
            return

        if bb_t.width > free_w:             # too wide -> shrink, then retry
            scale = max(0.6, 0.97 * free_w / max(bb_t.width, 1.0))
            new = max(title.get_fontsize() * scale, 5.0)
            if abs(new - title.get_fontsize()) < 0.05:
                return
            title.set_fontsize(new)
            continue

        centre_px = 0.5 * (free_left + ax_bb.x1)
        title.set_position(((centre_px - ax_bb.x0) / ax_bb.width,
                            title.get_position()[1]))
        return


def fit_panel_titles(fig):
    """Re-resolve every panel letter/title pair on the final geometry.

    Call this once, immediately before ``fig.savefig(...)``, after any
    ``tight_layout`` / ``subplots_adjust`` / ``constrained_layout`` has run.
    """
    for ax in fig.axes:
        ann = title = None
        for child in ax.texts:
            if child.get_gid() == LETTER_GID:
                ann = child
                break
        t = ax.get_title()
        if t:
            title = ax.title
        if ann is not None and title is not None:
            _resolve(ax, ann, title)
