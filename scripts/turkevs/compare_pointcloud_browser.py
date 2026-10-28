#!/usr/bin/env python3
"""
Browser-based viewer for one or two JSON point-cloud dataset files.

With one file, browse that dataset directly. With two files, records are
matched by:
    (additional_info["point_cloud_index"],
     additional_info["transformation"])
and shown side by side.

Navigation
----------
Left / Right : previous / next point_cloud_index
Up / Down    : previous / next transformation for the current index
Mouse/toolbar: pan, wheel zoom, box zoom, reset

In two-file mode the plots share x/y ranges, so zooming or panning either plot
updates both at once.

Direct jump
-----------
Type a point_cloud_index and transformation in the controls and press Enter
or click Go. Both fields provide autocomplete suggestions.

Run
---
    pip install bokeh numpy

    # Inspect one dataset
    python compare_pointcloud_browser.py file_a.json

    # Compare two datasets
    python compare_pointcloud_browser.py file_a.json file_b.json

Optional:
    python compare_pointcloud_browser.py file_a.json --port 5006
    python compare_pointcloud_browser.py file_a.json file_b.json --no-browser

This script does not use Tkinter.
"""

from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from bokeh.application import Application
from bokeh.application.handlers.function import FunctionHandler
from bokeh.document import Document
from bokeh.events import DocumentReady, ValueSubmit
from bokeh.layouts import column, row
from bokeh.models import (
    AutocompleteInput,
    Button,
    CheckboxGroup,
    ColumnDataSource,
    CustomJS,
    Div,
    Range1d,
    TextInput,
)
from bokeh.plotting import figure
from bokeh.server.server import Server


Key = tuple[int, str]


def extract_records(obj: Any) -> list[dict[str, Any]]:
    """Extract dataset records while deliberately ignoring declared counts."""
    if isinstance(obj, list):
        records = obj
    elif isinstance(obj, dict):
        for name in ("datasets", "data", "items", "records"):
            value = obj.get(name)
            if isinstance(value, list):
                records = value
                break
        else:
            candidates: list[list[Any]] = []
            for value in obj.values():
                if (
                    isinstance(value, list)
                    and value
                    and all(isinstance(x, dict) for x in value)
                ):
                    sample = value[0]
                    if "points" in sample and "additional_info" in sample:
                        candidates.append(value)

            if len(candidates) == 1:
                records = candidates[0]
            elif len(candidates) > 1:
                raise ValueError(
                    "The top-level JSON contains multiple dataset-like lists. "
                    "Rename the desired one to 'datasets'."
                )
            else:
                raise ValueError(
                    "Could not find a dataset list. Expected a top-level list or "
                    "an object containing a list such as 'datasets'."
                )
    else:
        raise ValueError("Top-level JSON must be a list or object.")

    out: list[dict[str, Any]] = []
    for i, record in enumerate(records):
        if not isinstance(record, dict):
            raise ValueError(f"Dataset #{i} is not a JSON object.")
        if "points" not in record:
            raise ValueError(f"Dataset #{i} has no 'points' field.")

        info = record.get("additional_info")
        if not isinstance(info, dict):
            raise ValueError(
                f"Dataset #{i} has no object-valued 'additional_info'."
            )
        if "point_cloud_index" not in info or "transformation" not in info:
            raise ValueError(
                f"Dataset #{i} is missing 'point_cloud_index' or "
                "'transformation' inside 'additional_info'."
            )
        out.append(record)

    return out


def load_json_records(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as f:
        return extract_records(json.load(f))


def make_index(
    records: Iterable[dict[str, Any]], source_name: str
) -> dict[Key, dict[str, Any]]:
    index: dict[Key, dict[str, Any]] = {}
    duplicates: list[Key] = []

    for record in records:
        info = record["additional_info"]
        try:
            pc_index = int(info["point_cloud_index"])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid point_cloud_index in {source_name}: "
                f"{info.get('point_cloud_index')!r}"
            ) from exc

        transformation = str(info["transformation"])
        key = (pc_index, transformation)
        if key in index:
            duplicates.append(key)
        else:
            index[key] = record

    if duplicates:
        preview = ", ".join(map(str, duplicates[:5]))
        extra = " ..." if len(duplicates) > 5 else ""
        raise ValueError(
            f"{source_name} contains duplicate "
            f"(point_cloud_index, transformation) keys: {preview}{extra}"
        )

    return index


def points_xy(record: dict[str, Any]) -> np.ndarray:
    points = np.asarray(record["points"], dtype=float)
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError("Each 'points' field must be an N x 2 (or wider) array.")
    return points[:, :2]


def pretty_info(record: dict[str, Any]) -> str:
    return json.dumps(
        record.get("additional_info", {}),
        indent=2,
        ensure_ascii=False,
        default=str,
    )


def make_app(path_a: Path, path_b: Path | None = None) -> Application:
    records_a = load_json_records(path_a)
    index_a = make_index(records_a, path_a.name)

    comparison_mode = path_b is not None

    if comparison_mode:
        assert path_b is not None
        records_b = load_json_records(path_b)
        index_b = make_index(records_b, path_b.name)

        available_keys = set(index_a) & set(index_b)
        if not available_keys:
            raise ValueError(
                "The two files have no shared "
                "(point_cloud_index, transformation) pairs."
            )

        only_a = len(set(index_a) - set(index_b))
        only_b = len(set(index_b) - set(index_a))
    else:
        records_b = None
        index_b = None
        available_keys = set(index_a)
        if not available_keys:
            raise ValueError("The dataset file contains no point-cloud records.")
        only_a = 0
        only_b = 0

    indices = sorted({idx for idx, _ in available_keys})

    # Preserve transformation order from file A instead of alphabetizing it.
    transform_order: list[str] = []
    seen_transforms: set[str] = set()
    for record in records_a:
        info = record["additional_info"]
        try:
            key = (int(info["point_cloud_index"]), str(info["transformation"]))
        except (TypeError, ValueError):
            continue
        transform = key[1]
        if key in available_keys and transform not in seen_transforms:
            seen_transforms.add(transform)
            transform_order.append(transform)

    transforms_by_index: dict[int, list[str]] = {}
    for idx in indices:
        present = {t for i, t in available_keys if i == idx}
        ordered = [t for t in transform_order if t in present]
        # Safety fallback for transformations not seen in the order above.
        ordered.extend(sorted(present - set(ordered)))
        transforms_by_index[idx] = ordered

    def modify_doc(doc: Document) -> None:
        state = {
            "index": indices[0],
            "transformation": transforms_by_index[indices[0]][0],
        }

        source_a = ColumnDataSource(data={"x": [], "y": []})
        source_b = (
            ColumnDataSource(data={"x": [], "y": []})
            if comparison_mode
            else None
        )

        shared_x = Range1d(0, 1)
        shared_y = Range1d(0, 1)

        tools = "pan,wheel_zoom,box_zoom,reset,save"
        plot_a = figure(
            title=path_a.name,
            x_range=shared_x,
            y_range=shared_y,
            tools=tools,
            active_scroll="wheel_zoom",
            match_aspect=True,
            sizing_mode="stretch_width",
            height=650,
        )
        plot_a.scatter("x", "y", source=source_a, size=5, alpha=0.75)

        plot_b = None
        if comparison_mode:
            assert path_b is not None
            assert source_b is not None
            plot_b = figure(
                title=path_b.name,
                x_range=shared_x,
                y_range=shared_y,
                tools=tools,
                active_scroll="wheel_zoom",
                match_aspect=True,
                sizing_mode="stretch_width",
                height=650,
            )
            plot_b.scatter("x", "y", source=source_b, size=5, alpha=0.75)

        for p in ([plot_a, plot_b] if plot_b is not None else [plot_a]):
            p.toolbar.autohide = True
            p.xaxis.axis_label = "x"
            p.yaxis.axis_label = "y"

        index_input = AutocompleteInput(
            title="Point cloud index",
            value=str(state["index"]),
            completions=[str(i) for i in indices],
            min_characters=0,
            width=190,
        )
        transform_input = AutocompleteInput(
            title="Transformation",
            value=state["transformation"],
            completions=transforms_by_index[state["index"]],
            min_characters=0,
            width=250,
        )

        go_button = Button(label="Go", button_type="primary", width=70)
        prev_index_button = Button(label="← index", width=85)
        next_index_button = Button(label="index →", width=85)
        prev_transform_button = Button(label="↑ transform", width=105)
        next_transform_button = Button(label="transform ↓", width=105)
        fit_button = Button(label="Fit view", width=85)
        auto_fit = CheckboxGroup(
            labels=["Auto-fit after navigation"], active=[0], width=210
        )

        status = Div(text="", sizing_mode="stretch_width")

        info_style = {
            "overflow-y": "auto",
            "border": "1px solid #ddd",
            "padding": "8px",
            "box-sizing": "border-box",
            "font-family": "monospace",
            "white-space": "pre-wrap",
        }
        info_a = Div(
            text="",
            sizing_mode="stretch_width",
            height=280,
            styles=info_style,
        )
        info_b = (
            Div(
                text="",
                sizing_mode="stretch_width",
                height=280,
                styles=info_style,
            )
            if comparison_mode
            else None
        )

        # Hidden bridge used by browser JavaScript for arrow-key commands.
        keyboard_command = TextInput(value="", visible=False)

        def current_key() -> Key:
            return state["index"], state["transformation"]

        def update_transform_completions() -> None:
            transform_input.completions = transforms_by_index[state["index"]]

        def fit_current_view() -> None:
            point_sets = [points_xy(index_a[current_key()])]
            if comparison_mode:
                assert index_b is not None
                point_sets.append(points_xy(index_b[current_key()]))

            all_points = np.vstack(point_sets)

            xmin = float(np.nanmin(all_points[:, 0]))
            xmax = float(np.nanmax(all_points[:, 0]))
            ymin = float(np.nanmin(all_points[:, 1]))
            ymax = float(np.nanmax(all_points[:, 1]))

            xspan = xmax - xmin
            yspan = ymax - ymin
            if not np.isfinite(xspan) or xspan <= 0:
                xspan = 1.0
            if not np.isfinite(yspan) or yspan <= 0:
                yspan = 1.0

            pad_x = 0.06 * xspan
            pad_y = 0.06 * yspan
            shared_x.start = xmin - pad_x
            shared_x.end = xmax + pad_x
            shared_y.start = ymin - pad_y
            shared_y.end = ymax + pad_y

        def refresh_view(*, fit: bool | None = None) -> None:
            key = current_key()
            record_a = index_a[key]
            pts_a = points_xy(record_a)
            source_a.data = {"x": pts_a[:, 0], "y": pts_a[:, 1]}

            index_input.value = str(state["index"])
            update_transform_completions()
            transform_input.value = state["transformation"]

            plot_a.title.text = f"{path_a.name}  —  {len(pts_a):,} points"
            info_a.text = (
                f"<b>{html.escape(path_a.name)} — additional_info</b>"
                f"<pre>{html.escape(pretty_info(record_a))}</pre>"
            )

            if comparison_mode:
                assert path_b is not None
                assert index_b is not None
                assert source_b is not None
                assert plot_b is not None
                assert info_b is not None

                record_b = index_b[key]
                pts_b = points_xy(record_b)
                source_b.data = {"x": pts_b[:, 0], "y": pts_b[:, 1]}
                plot_b.title.text = f"{path_b.name}  —  {len(pts_b):,} points"
                info_b.text = (
                    f"<b>{html.escape(path_b.name)} — additional_info</b>"
                    f"<pre>{html.escape(pretty_info(record_b))}</pre>"
                )

            transform_list = transforms_by_index[state["index"]]
            transform_position = transform_list.index(state["transformation"]) + 1
            index_position = indices.index(state["index"]) + 1

            unmatched = ""
            if comparison_mode and (only_a or only_b):
                unmatched = (
                    f" | unmatched pairs: {only_a} only in A, {only_b} only in B"
                )

            status.text = (
                f"<b>point_cloud_index:</b> {state['index']} "
                f"({index_position}/{len(indices)}) &nbsp;&nbsp; "
                f"<b>transformation:</b> {html.escape(state['transformation'])} "
                f"({transform_position}/{len(transform_list)})"
                f"{unmatched}"
            )

            should_fit = (0 in auto_fit.active) if fit is None else fit
            if should_fit:
                fit_current_view()

        def move_index(delta: int) -> None:
            pos = indices.index(state["index"])
            new_pos = (pos + delta) % len(indices)
            new_index = indices[new_pos]
            old_transform = state["transformation"]
            state["index"] = new_index

            available = transforms_by_index[new_index]
            state["transformation"] = (
                old_transform if old_transform in available else available[0]
            )
            refresh_view()

        def move_transformation(delta: int) -> None:
            available = transforms_by_index[state["index"]]
            pos = available.index(state["transformation"])
            state["transformation"] = available[(pos + delta) % len(available)]
            refresh_view()

        def apply_direct_jump() -> None:
            try:
                requested_index = int(index_input.value.strip())
            except ValueError:
                status.text = (
                    f"<b>Invalid point_cloud_index:</b> "
                    f"{html.escape(index_input.value)!s}"
                )
                return

            requested_transform = transform_input.value.strip()
            requested_key = (requested_index, requested_transform)
            if requested_key not in available_keys:
                if requested_index not in transforms_by_index:
                    noun = "shared point cloud" if comparison_mode else "point cloud"
                    status.text = (
                        f"<b>No {noun} with index {requested_index}.</b>"
                    )
                else:
                    available = ", ".join(
                        html.escape(x) for x in transforms_by_index[requested_index]
                    )
                    noun = "shared pair" if comparison_mode else "dataset entry"
                    status.text = (
                        f"<b>No {noun}</b> ({requested_index}, "
                        f"{html.escape(requested_transform)}). "
                        f"Available transformations for that index: {available}"
                    )
                return

            state["index"] = requested_index
            state["transformation"] = requested_transform
            refresh_view()

        def on_index_text_changed(attr: str, old: str, new: str) -> None:
            # Update transformation suggestions as soon as a valid index is typed,
            # without changing the displayed dataset until Go/Enter.
            try:
                idx = int(new.strip())
            except ValueError:
                return
            if idx in transforms_by_index:
                transform_input.completions = transforms_by_index[idx]

        def handle_keyboard_command(attr: str, old: str, new: str) -> None:
            command = new.split(":", 1)[0]
            if command == "left":
                move_index(-1)
            elif command == "right":
                move_index(+1)
            elif command == "up":
                move_transformation(-1)
            elif command == "down":
                move_transformation(+1)

        prev_index_button.on_click(lambda: move_index(-1))
        next_index_button.on_click(lambda: move_index(+1))
        prev_transform_button.on_click(lambda: move_transformation(-1))
        next_transform_button.on_click(lambda: move_transformation(+1))
        fit_button.on_click(fit_current_view)
        go_button.on_click(apply_direct_jump)
        index_input.on_change("value", on_index_text_changed)
        index_input.on_event(ValueSubmit, lambda event: apply_direct_jump())
        transform_input.on_event(ValueSubmit, lambda event: apply_direct_jump())
        keyboard_command.on_change("value", handle_keyboard_command)

        # Capture arrow keys anywhere on the page except while the user is
        # actively typing/selecting in a form control.
        doc.js_on_event(
            DocumentReady,
            CustomJS(
                args=dict(bridge=keyboard_command),
                code="""
                if (window.__pointcloud_compare_key_handler_installed) {
                    return;
                }
                window.__pointcloud_compare_key_handler_installed = true;

                document.addEventListener('keydown', (event) => {
                    const tag = (event.target && event.target.tagName)
                        ? event.target.tagName.toUpperCase()
                        : '';
                    if (['INPUT', 'TEXTAREA', 'SELECT'].includes(tag)) {
                        return;
                    }

                    const commands = {
                        ArrowLeft: 'left',
                        ArrowRight: 'right',
                        ArrowUp: 'up',
                        ArrowDown: 'down',
                    };
                    const command = commands[event.key];
                    if (!command) {
                        return;
                    }

                    event.preventDefault();
                    bridge.value = `${command}:${Date.now()}:${Math.random()}`;
                    bridge.change.emit();
                });
                """,
            ),
        )

        controls = row(
            index_input,
            transform_input,
            go_button,
            prev_index_button,
            next_index_button,
            prev_transform_button,
            next_transform_button,
            fit_button,
            auto_fit,
            sizing_mode="stretch_width",
        )

        plot_help = (
            "<b>Plots:</b> wheel zoom, pan, box zoom; ranges are linked"
            if comparison_mode
            else "<b>Plot:</b> wheel zoom, pan, box zoom"
        )
        instructions = Div(
            text=(
                "<b>Keyboard:</b> ←/→ point cloud &nbsp; "
                "↑/↓ transformation &nbsp; | &nbsp; "
                f"{plot_help}"
            ),
            sizing_mode="stretch_width",
        )

        plots = (
            row(plot_a, plot_b, sizing_mode="stretch_width")
            if plot_b is not None
            else plot_a
        )
        metadata = (
            row(info_a, info_b, sizing_mode="stretch_width")
            if info_b is not None
            else info_a
        )

        root = column(
            controls,
            instructions,
            status,
            plots,
            metadata,
            keyboard_command,
            sizing_mode="stretch_width",
        )

        doc.add_root(root)
        doc.title = (
            "Point-cloud JSON comparison"
            if comparison_mode
            else "Point-cloud JSON viewer"
        )
        refresh_view(fit=True)

    return Application(FunctionHandler(modify_doc))

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Browse one point-cloud JSON file, or compare two files side by side "
            "in a local browser UI."
        )
    )
    parser.add_argument("file_a", type=Path, help="JSON dataset file to inspect")
    parser.add_argument(
        "file_b",
        type=Path,
        nargs="?",
        default=None,
        help="Optional second JSON file for side-by-side comparison",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=5006,
        help="Local server port (default: 5006)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Start the server without automatically opening a browser tab",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    paths = [args.file_a]
    if args.file_b is not None:
        paths.append(args.file_b)

    for path in paths:
        if not path.is_file():
            print(f"error: file not found: {path}", file=sys.stderr)
            return 2

    try:
        app = make_app(args.file_a, args.file_b)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    server = Server({"/": app}, port=args.port)
    server.start()

    url = f"http://localhost:{server.port}/"
    mode = "comparison" if args.file_b is not None else "viewer"
    print(f"Point-cloud {mode} running at {url}")
    print("Press Ctrl+C in this terminal to stop it.")

    if not args.no_browser:
        server.io_loop.add_callback(server.show, "/")

    try:
        server.io_loop.start()
    except KeyboardInterrupt:
        pass
    finally:
        server.stop()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
