#!/usr/bin/env python3
"""
GPH Viewer - PyQt GUI to browse and edit GPH binary files (similar to HDFView).

Large meshes (e.g. laptop_simplified_voxel_v4/v6.gph, laptop_simplified_denser_v2_gph.gph)
are opened via mmap; LS_Links conn split status is shown in the status bar
(conn_split×N from parse_ls_links_summary).
"""

import struct
import sys
from pathlib import Path

try:
    from PyQt6.QtWidgets import (
        QApplication,
        QMainWindow,
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QTreeWidget,
        QTreeWidgetItem,
        QSplitter,
        QTextEdit,
        QTableWidget,
        QTableWidgetItem,
        QFileDialog,
        QLabel,
        QGroupBox,
        QLineEdit,
        QPushButton,
        QMessageBox,
        QHeaderView,
        QTabWidget,
        QAbstractItemView,
    )
    from PyQt6.QtCore import Qt, QPointF
    from PyQt6.QtGui import QAction, QColor, QPainter, QPen, QBrush, QPolygonF
    HAS_PYQT = "PyQt6"
except ImportError:
    try:
        from PyQt5.QtWidgets import (
            QApplication,
            QMainWindow,
            QWidget,
            QVBoxLayout,
            QHBoxLayout,
            QTreeWidget,
            QTreeWidgetItem,
            QSplitter,
            QTextEdit,
            QTableWidget,
            QTableWidgetItem,
                    QFileDialog,
                QLabel,
            QGroupBox,
            QLineEdit,
            QPushButton,
            QMessageBox,
            QHeaderView,
            QTabWidget,
            QAbstractItemView,
        )
        from PyQt5.QtCore import Qt, QPointF
        from PyQt5.QtGui import QColor, QPainter, QPen, QBrush, QPolygonF
        from PyQt5.QtWidgets import QAction
        HAS_PYQT = "PyQt5"
    except ImportError:
        HAS_PYQT = None

import numpy as np

from gph_model import (
    GphDocument,
    GphNode,
    build_mesh_preview,
    classify_volume_region_cells,
    parse_ls_cvol_ids,
    parse_ls_links_summary,
    parse_ls_nodes_vertices,
    parse_ls_parts,
    parse_ls_surface_regions,
    parse_ls_string_list,
    part_cvol_cell_mask,
    format_part_cvol_spec,
)


def read_f32_be(data: bytes, pos: int) -> float:
    return struct.unpack(">f", data[pos : pos + 4])[0]


def read_i32_be(data: bytes, pos: int) -> int:
    return int.from_bytes(data[pos : pos + 4], "big")


class MeshPreviewWidget(QWidget):
    """Simple Qt-painted 3D mesh preview with rotate / pan / zoom."""

    def __init__(self):
        super().__init__()
        self.preview = None
        self.title = "Open a GPH file to preview mesh regions"
        self.rot_x = -25.0
        self.rot_y = 35.0
        self.zoom = 1.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self._last_pos = None
        self.setMinimumHeight(320)

    def set_preview(self, preview: dict | None, title: str) -> None:
        self.preview = preview
        self.title = title
        self.update()

    def _event_xy(self, event):
        if hasattr(event, "position"):
            p = event.position()
            return p.x(), p.y()
        p = event.pos()
        return p.x(), p.y()

    def mousePressEvent(self, event):
        self._last_pos = self._event_xy(event)

    def mouseMoveEvent(self, event):
        if self._last_pos is None:
            return
        x, y = self._event_xy(event)
        dx = x - self._last_pos[0]
        dy = y - self._last_pos[1]
        buttons = event.buttons()
        if hasattr(Qt, "MouseButton"):
            middle = Qt.MouseButton.MiddleButton
            right = Qt.MouseButton.RightButton
        else:
            middle = Qt.MiddleButton
            right = Qt.RightButton
        if (buttons & middle) or (buttons & right):
            self.pan_x += dx
            self.pan_y += dy
        else:
            self.rot_y += dx * 0.5
            self.rot_x += dy * 0.5
        self._last_pos = (x, y)
        self.update()

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        self.zoom *= 1.15 if delta > 0 else 1 / 1.15
        self.zoom = max(0.05, min(self.zoom, 100.0))
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        if hasattr(QPainter, "RenderHint"):
            antialias = QPainter.RenderHint.Antialiasing
        else:
            antialias = QPainter.Antialiasing
        painter.setRenderHint(antialias)
        painter.fillRect(self.rect(), QColor(22, 24, 28))
        painter.setPen(QPen(QColor(230, 230, 230)))
        painter.drawText(12, 22, self.title)

        if not self.preview or not self.preview.get("faces"):
            summary = (
                self.preview.get("summary", "No mesh preview available")
                if self.preview else "No mesh preview available"
            )
            painter.setPen(QPen(QColor(180, 180, 180)))
            painter.drawText(12, 48, summary)
            painter.end()
            return

        faces = self.preview["faces"]
        all_pts = np.vstack(faces)
        center = all_pts.mean(axis=0)
        span = np.ptp(all_pts, axis=0)
        scale = max(float(span.max()), 1e-12)
        rx = np.deg2rad(self.rot_x)
        ry = np.deg2rad(self.rot_y)
        cx, sx = np.cos(rx), np.sin(rx)
        cy, sy = np.cos(ry), np.sin(ry)
        rot_x = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
        rot_y = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
        rot = rot_y @ rot_x
        size = min(self.width(), self.height()) * 0.78 * self.zoom
        origin = np.array([
            self.width() * 0.5 + self.pan_x,
            self.height() * 0.55 + self.pan_y,
        ])

        projected = []
        for poly in faces:
            pts = ((poly - center) / scale) @ rot.T
            screen = np.column_stack([pts[:, 0], -pts[:, 1]]) * size + origin
            projected.append((float(pts[:, 2].mean()), screen))
        projected.sort(key=lambda item: item[0])

        selected = self.preview.get("selection_active", False)
        fill = QColor(255, 155, 70, 150) if selected else QColor(90, 150, 255, 95)
        edge = QColor(255, 235, 190) if selected else QColor(175, 205, 255)
        painter.setPen(QPen(edge, 1.1))
        painter.setBrush(QBrush(fill))
        for _, screen in projected:
            polygon = QPolygonF([QPointF(float(x), float(y)) for x, y in screen])
            painter.drawPolygon(polygon)

        painter.setBrush(QBrush())
        painter.setPen(QPen(QColor(145, 145, 145)))
        painter.drawText(
            12,
            self.height() - 16,
            "Drag: rotate   Right/middle drag: pan   Wheel: zoom",
        )
        painter.end()


class GphViewerMain(QMainWindow):
    def __init__(self):
        super().__init__()
        self.doc = GphDocument()
        self.current_node: GphNode | None = None
        self._surface_regions = []
        self._cvol_ids = None
        self._parts = []
        self._volume_regions = []
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("GPH Viewer - Browse & Edit GPH Files")
        self.setMinimumSize(900, 600)
        self.resize(1100, 700)

        # Menu
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")
        open_act = QAction("Open...", self)
        open_act.setShortcut("Ctrl+O")
        open_act.triggered.connect(self.on_open)
        file_menu.addAction(open_act)

        save_act = QAction("Save", self)
        save_act.setShortcut("Ctrl+S")
        save_act.triggered.connect(self.on_save)
        file_menu.addAction(save_act)

        save_as_act = QAction("Save As...", self)
        save_as_act.triggered.connect(self.on_save_as)
        file_menu.addAction(save_as_act)

        file_menu.addSeparator()
        exit_act = QAction("Exit", self)
        exit_act.setShortcut("Ctrl+Q")
        exit_act.triggered.connect(self.close)
        file_menu.addAction(exit_act)

        self._mmap_readonly = False

        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        splitter = QSplitter()

        # Left: Tree
        tree_box = QGroupBox("Structure")
        tree_layout = QVBoxLayout(tree_box)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Name", "Type", "Address Range", "Size"])
        if hasattr(QAbstractItemView, "SelectionMode"):
            self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        else:
            self.tree.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tree.itemSelectionChanged.connect(self.on_selection_changed)
        tree_layout.addWidget(self.tree)
        splitter.addWidget(tree_box)

        # Right: Detail + Edit
        right = QWidget()
        right_layout = QVBoxLayout(right)

        self.tabs = QTabWidget()
        self.attrs_text = QTextEdit()
        self.attrs_text.setReadOnly(True)
        self.attrs_text.setPlaceholderText("Select a node to view details")
        self.tabs.addTab(self.attrs_text, "Attributes")

        self.data_table = QTableWidget()
        self.data_table.setAlternatingRowColors(True)
        self.data_table.itemChanged.connect(self.on_table_cell_changed)
        self._table_loading = False
        self.tabs.addTab(self.data_table, "Data")

        self.hex_text = QTextEdit()
        self.hex_text.setReadOnly(True)
        self.hex_text.setFontFamily("Consolas")
        self.tabs.addTab(self.hex_text, "Hex Dump")

        mesh_tab = QWidget()
        mesh_layout = QVBoxLayout(mesh_tab)
        self.mesh_info = QLabel("Select a surface or volume region in the tree")
        self.mesh_view = MeshPreviewWidget()
        mesh_layout.addWidget(self.mesh_info)
        mesh_layout.addWidget(self.mesh_view, 1)
        self.tabs.addTab(mesh_tab, "3D Regions")

        right_layout.addWidget(self.tabs)

        edit_box = QGroupBox("Edit")
        edit_layout = QHBoxLayout(edit_box)
        self.edit_label = QLabel("Value:")
        self.edit_input = QLineEdit()
        self.edit_input.setPlaceholderText("Enter new value and click Apply")
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.clicked.connect(self.on_apply_edit)
        edit_layout.addWidget(self.edit_label)
        edit_layout.addWidget(self.edit_input, 1)
        edit_layout.addWidget(self.apply_btn)
        right_layout.addWidget(edit_box)

        splitter.addWidget(right)
        splitter.setSizes([350, 650])
        layout.addWidget(splitter)

        self.statusBar().showMessage("Ready - Open a GPH file (Ctrl+O)")

    def on_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open GPH File",
            "",
            "GPH Files (*.gph);;All Files (*)",
        )
        if path:
            self.load_file(path)

    def load_file(self, path: str):
        if not self.doc.load(path):
            QMessageBox.critical(self, "Error", f"Failed to load: {path}")
            return
        self._mmap_readonly = getattr(self.doc, "_mmap_mode", False)
        self._load_region_index()
        self._build_tree()
        self._update_mesh_preview([])
        self.setWindowTitle(f"GPH Viewer - {Path(path).name}")
        summary = self._summarize_buffer(self.doc._raw_data)
        size_mb = len(self.doc._raw_data) / (1024 * 1024)
        mmap_note = " (mmap, read-only)" if self._mmap_readonly else ""
        self.statusBar().showMessage(
            f"Opened: {path} ({size_mb:.1f} MiB{mmap_note}) — {summary}"
        )

    def _load_region_index(self):
        data = self.doc._raw_data
        self._surface_regions = parse_ls_surface_regions(data)
        self._cvol_ids = parse_ls_cvol_ids(data)
        self._parts = parse_ls_parts(data, cvol_id=self._cvol_ids)
        self._volume_regions = parse_ls_string_list(data, "LS_VolumeRegions")

    def closeEvent(self, event):
        self.doc.close()
        super().closeEvent(event)

    def _summarize_buffer(self, data) -> str:
        parts = []
        _, _, nv = parse_ls_nodes_vertices(data)
        if nv:
            parts.append(f"{nv} verts")
        links = parse_ls_links_summary(data)
        if links:
            parts.append(f"{links['n_faces']} faces")
            parts.append(f"{links['n_cells']} cells")
            if links.get("polyhedral"):
                parts.append("polyhedral")
            if links.get("conn_split"):
                chunks = links.get("conn_chunks", 2)
                parts.append(f"conn_split×{chunks}")
            elif not links.get("conn_complete", True):
                parts.append("conn_INCOMPLETE")
        cvol = parse_ls_cvol_ids(data)
        pmeta = parse_ls_parts(data, cvol_id=cvol)
        n_cells = links["n_cells"] if links else (len(cvol) if cvol is not None else 0)
        if pmeta and cvol is not None and len(cvol) == n_cells:
            for pname, cv in pmeta:
                n = int((cvol == cv).sum())
                parts.append(f"{pname}={n}")
        elif pmeta:
            parts.append(f"{len(pmeta)} parts")
        regions = parse_ls_string_list(data, "LS_VolumeRegions")
        if regions:
            parts.append(f"vol_regions={len(regions)}")
        return ", ".join(parts) if parts else "structure only"

    def _build_tree(self):
        self.tree.clear()
        root = self.doc.root
        if not root:
            return
        item = QTreeWidgetItem([
            root.name,
            root.data_type,
            self._range_text(root),
            str(root.size),
        ])
        item.setData(0, 100, root)
        self._add_children(item, root)
        self._add_region_tree(item)
        self.tree.addTopLevelItem(item)
        self.tree.expandToDepth(1)
        for col in range(self.tree.columnCount()):
            self.tree.resizeColumnToContents(col)

    def _add_children(self, parent_item: QTreeWidgetItem, node: GphNode):
        for child in node.children:
            item = QTreeWidgetItem([
                child.name,
                child.data_type,
                self._range_text(child),
                str(child.size),
            ])
            item.setData(0, 100, child)
            parent_item.addChild(item)
            self._add_children(item, child)

    def _range_text(self, node: GphNode) -> str:
        if node.size <= 0:
            return "virtual"
        return f"0x{node.offset:04X}-0x{node.offset + node.size - 1:04X}"

    def _add_region_tree(self, root_item: QTreeWidgetItem):
        regions_root = GphNode(
            "3D Regions", 0, 0, "virtual", "Interactive mesh region selectors",
            metadata={"view_kind": "folder"},
        )
        regions_item = QTreeWidgetItem([regions_root.name, regions_root.data_type, "virtual", ""])
        regions_item.setData(0, 100, regions_root)
        root_item.addChild(regions_item)

        all_node = GphNode(
            "All mesh", 0, 0, "mesh preview", "Show the full mesh preview",
            metadata={"view_kind": "mesh"},
        )
        all_item = QTreeWidgetItem([all_node.name, all_node.data_type, "virtual", ""])
        all_item.setData(0, 100, all_node)
        regions_item.addChild(all_item)

        surf_root = GphNode("Surface Regions", 0, 0, "folder", metadata={"view_kind": "folder"})
        surf_item = QTreeWidgetItem([surf_root.name, surf_root.data_type, "virtual", str(len(self._surface_regions))])
        surf_item.setData(0, 100, surf_root)
        regions_item.addChild(surf_item)
        for name, face_ids in self._surface_regions:
            node = GphNode(
                name, 0, 0, "surface region", f"{len(face_ids)} faces",
                metadata={"view_kind": "surface", "face_ids": face_ids},
            )
            item = QTreeWidgetItem([node.name, node.data_type, "virtual", str(len(face_ids))])
            item.setData(0, 100, node)
            surf_item.addChild(item)

        vol_root = GphNode("Volume Regions", 0, 0, "folder", metadata={"view_kind": "folder"})
        vol_item = QTreeWidgetItem([vol_root.name, vol_root.data_type, "virtual", ""])
        vol_item.setData(0, 100, vol_root)
        regions_item.addChild(vol_item)
        for name in self._volume_regions:
            node = GphNode(
                name, 0, 0, "volume region", "Volume-region cell selection",
                metadata={"view_kind": "volume_region", "region_name": name},
            )
            item = QTreeWidgetItem([node.name, node.data_type, "virtual", ""])
            item.setData(0, 100, node)
            vol_item.addChild(item)
        for name, cvol_spec in self._parts:
            count = (int(part_cvol_cell_mask(self._cvol_ids, cvol_spec).sum())
                     if self._cvol_ids is not None else 0)
            node = GphNode(
                f"Part: {name}", 0, 0, "volume part",
                f"cvol={format_part_cvol_spec(cvol_spec)}, cells={count}",
                metadata={"view_kind": "volume_part", "cvol_spec": cvol_spec},
            )
            item = QTreeWidgetItem([node.name, node.data_type, "virtual", str(count)])
            item.setData(0, 100, node)
            vol_item.addChild(item)

    def on_selection_changed(self):
        items = self.tree.selectedItems()
        if not items:
            self.current_node = None
            self._show_empty()
            return
        node = items[0].data(0, 100)
        if isinstance(node, GphNode):
            self.current_node = node
            self._show_node(node)
            self._update_mesh_preview([
                item.data(0, 100) for item in items
                if isinstance(item.data(0, 100), GphNode)
            ])
        else:
            self.current_node = None

    def _show_empty(self):
        self.attrs_text.setPlainText("")
        self.data_table.setRowCount(0)
        self.hex_text.setPlainText("")
        self.edit_input.setEnabled(False)
        self.apply_btn.setEnabled(False)

    def _show_node(self, node: GphNode):
        # Attributes
        lines = [
            f"Name: {node.name}",
            f"Offset: 0x{node.offset:04X}",
            f"Address range: {self._range_text(node)}",
            f"Size: {node.size} bytes",
            f"Type: {node.data_type}",
            "",
        ]
        if node.value is not None:
            if isinstance(node.value, str) and "\n" in node.value:
                lines.append("Value:")
                lines.append(node.value)
            elif isinstance(node.value, (list, tuple)):
                if len(node.value) <= 20:
                    lines.append(f"Value: {node.value}")
                else:
                    lines.append(f"Value: [{len(node.value)} items] (see Data tab)")
            else:
                lines.append(f"Value: {node.value}")
        if node.modified:
            lines.append("\n(Modified)")
        extra = self._extra_node_lines(node)
        if extra:
            lines.append("")
            lines.extend(extra)
        self.attrs_text.setPlainText("\n".join(lines))

        # Data tab (table for arrays)
        self._table_loading = True
        if isinstance(node.value, (list, tuple)) and node.value:
            first = node.value[0]
            if isinstance(first, str):
                self.data_table.setColumnCount(1)
                self.data_table.setHorizontalHeaderLabels(["Entry"])
                self.data_table.setRowCount(len(node.value))
                for i, v in enumerate(node.value):
                    self.data_table.setItem(i, 0, QTableWidgetItem(str(v)))
            elif isinstance(first, (int, float)):
                self.data_table.setColumnCount(1)
                self.data_table.setHorizontalHeaderLabels(["Value"])
                self.data_table.setRowCount(len(node.value))
                for i, v in enumerate(node.value):
                    self.data_table.setItem(i, 0, QTableWidgetItem(str(v)))
            elif isinstance(first, (tuple, list)) and len(first) >= 2:
                ncols = len(first)
                self.data_table.setColumnCount(ncols)
                self.data_table.setHorizontalHeaderLabels([f"Col{i}" for i in range(ncols)])
                self.data_table.setRowCount(len(node.value))
                for i, row in enumerate(node.value):
                    for j, v in enumerate(row):
                        if j < ncols:
                            self.data_table.setItem(i, j, QTableWidgetItem(str(v)))
            else:
                self.data_table.setColumnCount(1)
                self.data_table.setRowCount(len(node.value))
                for i, v in enumerate(node.value):
                    self.data_table.setItem(i, 0, QTableWidgetItem(str(v)))
            hdr = self.data_table.horizontalHeader()
            if hasattr(QHeaderView.ResizeMode, "Stretch"):
                hdr.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            else:
                hdr.setSectionResizeMode(QHeaderView.Stretch)
        else:
            self.data_table.setRowCount(0)
        self._table_loading = False

        # Hex dump
        raw = node.get_raw(self.doc)
        hex_lines = []
        for i in range(0, min(len(raw), 4096), 16):
            chunk = raw[i : i + 16]
            hex_part = " ".join(f"{b:02X}" for b in chunk)
            ascii_part = "".join(chr(b) if 32 <= b < 127 else "." for b in chunk)
            hex_lines.append(f"{node.offset + i:04X}: {hex_part:<48} {ascii_part}")
        if len(raw) > 4096:
            hex_lines.append(f"... ({len(raw) - 4096} more bytes)")
        self.hex_text.setPlainText("\n".join(hex_lines))

        # Edit: enable for scalar I4, R4, C1
        can_edit = node.data_type in ("I4", "C1[8]", "C1") and node.value is not None
        if can_edit and not isinstance(node.value, (list, tuple)):
            self.edit_input.setEnabled(True)
            self.edit_input.setText(str(node.value))
            self.apply_btn.setEnabled(True)
        elif (node.data_type.startswith(("I4[", "R8["))
              or node.data_type in ("I4[]", "R8[]", "topology", "parts", "regions",
                                    "surf_regions", "assembly_xml")) and node.value:
            self.edit_input.setEnabled(True)
            self.edit_input.setPlaceholderText("Use Data tab; edit not yet implemented for arrays")
            self.edit_input.clear()
            self.apply_btn.setEnabled(False)
        else:
            self.edit_input.setEnabled(False)
            self.edit_input.clear()
            self.apply_btn.setEnabled(False)

    def _extra_node_lines(self, node: GphNode) -> list[str]:
        """Partition / topology details aligned with gph2cgns diagnostics."""
        data = self.doc._raw_data
        lines: list[str] = []

        if node.name == "LS_CvolIdOfElements":
            cvol = parse_ls_cvol_ids(data)
            if cvol is not None:
                lines.append(f"Per-cell array length: {len(cvol)}")
                for pname, cv in parse_ls_parts(data, cvol_id=cvol):
                    n = int(part_cvol_cell_mask(cvol, cv).sum())
                    lines.append(f"  {pname} (cvol={format_part_cvol_spec(cv)}): {n} cells")

        if node.name == "LS_Links":
            summary = parse_ls_links_summary(data)
            if summary:
                if summary.get("conn_split"):
                    lines.append(
                        f"conn split: {summary.get('conn_got', '?')}/"
                        f"{summary.get('conn_entries', '?')} entries "
                        f"({summary.get('conn_chunks', '?')} chunks)"
                    )
                elif not summary.get("conn_complete", True):
                    lines.append(
                        f"conn INCOMPLETE: {summary.get('conn_got', '?')}/"
                        f"{summary.get('conn_entries', '?')} entries"
                    )
                lines.append(
                    f"npe range: [{summary['npe_min']}..{summary['npe_max']}]"
                )

        if node.name == "LS_Parts":
            cvol = parse_ls_cvol_ids(data)
            if cvol is not None:
                for pname, cv in parse_ls_parts(data, cvol_id=cvol):
                    n = int(part_cvol_cell_mask(cvol, cv).sum())
                    lines.append(f"  {pname}: cvol={format_part_cvol_spec(cv)}, cells={n}")

        return lines

    def _update_mesh_preview(self, nodes: list[GphNode]):
        if not getattr(self.doc, "_raw_data", None):
            return
        selected_faces: list[np.ndarray] = []
        selected_cells: list[np.ndarray] = []
        labels: list[str] = []
        show_all = False

        for node in nodes:
            meta = node.metadata or {}
            kind = meta.get("view_kind")
            if kind == "mesh":
                show_all = True
                labels.append("All mesh")
            elif kind == "surface":
                selected_faces.append(np.asarray(meta.get("face_ids", []), dtype=np.int64))
                labels.append(node.name)
            elif kind == "volume_part" and self._cvol_ids is not None:
                cvol_spec = meta.get("cvol_spec", meta.get("cvol_id"))
                if cvol_spec is not None:
                    selected_cells.append(
                        np.flatnonzero(part_cvol_cell_mask(self._cvol_ids, cvol_spec))
                    )
                labels.append(node.name)
            elif kind == "volume_region":
                summary = parse_ls_links_summary(self.doc._raw_data)
                if summary:
                    mask = classify_volume_region_cells(
                        meta.get("region_name", node.name),
                        self._parts,
                        self._cvol_ids,
                        int(summary["n_cells"]),
                    )
                    selected_cells.append(np.flatnonzero(mask))
                    labels.append(node.name)

        face_ids = np.concatenate(selected_faces) if selected_faces and not show_all else None
        cell_ids = np.concatenate(selected_cells) if selected_cells and not show_all else None
        preview = build_mesh_preview(
            self.doc._raw_data,
            selected_face_ids=face_ids,
            selected_cell_ids=cell_ids,
        )
        if preview is None:
            self.mesh_info.setText("3D preview unavailable: LS_Nodes or LS_Links is incomplete")
            self.mesh_view.set_preview(None, "3D preview unavailable")
            return

        selected_label = " + ".join(labels) if labels else "All mesh"
        stats = [
            selected_label,
            preview.get("summary", ""),
            f"mesh={preview.get('n_vertices', 0)} vertices / {preview.get('n_faces', 0)} faces / {preview.get('n_cells', 0)} cells",
        ]
        if preview.get("dialect"):
            stats.append(preview["dialect"])
        self.mesh_info.setText(" | ".join(part for part in stats if part))
        self.mesh_view.set_preview(preview, selected_label)

    def on_table_cell_changed(self, item: QTableWidgetItem):
        # Vertex / topology / partition arrays are read-only in the viewer.
        # Coordinates: three separate axis blocks (float32 or float64 BE, optional
        # word-reversed f64); parsed via gph_model.parse_ls_nodes_xyz.  LS_Links
        # uses variable-length CSR connectivity for polyhedral meshes — inline
        # edits are not supported.
        if self._table_loading or not self.current_node:
            return
        # No-op: keep the table read-only for vertex arrays.
        self.statusBar().showMessage(
            "Vertex editing is disabled in the viewer (see DEV_SUMMARY.md)."
        )

    def on_apply_edit(self):
        node = self.current_node
        if not node or not self.edit_input.isEnabled():
            return
        text = self.edit_input.text().strip()
        raw = node.get_raw(self.doc)

        # Determine data location in section (heuristic per known fields)
        if node.name == "FileRevision" and node.data_type == "I4":
            local_offset = 64
            size = 4
            try:
                v = int(text)
                new_raw = v.to_bytes(4, "big")
            except ValueError:
                QMessageBox.warning(self, "Invalid", "Enter an integer")
                return
        elif node.name == "Dimension" and node.data_type == "I4":
            local_offset = 64
            size = 4
            try:
                v = int(text)
                new_raw = v.to_bytes(4, "big")
            except ValueError:
                QMessageBox.warning(self, "Invalid", "Enter an integer")
                return
        elif node.name == "Application" and node.data_type == "C1[8]":
            local_offset = 64  # value after descriptor
            size = 8
            s = text.encode("ascii", errors="replace")[:8].ljust(8)
            new_raw = s
        else:
            QMessageBox.information(self, "Info", f"Editing '{node.name}' not yet implemented")
            return

        file_offset = node.offset + local_offset
        if local_offset + size > len(raw):
            QMessageBox.warning(self, "Error", "Offset out of range")
            return

        self.doc.apply_patch(file_offset, new_raw)
        node.raw = raw[:local_offset] + new_raw + raw[local_offset + size :]
        node.value = int(text) if node.data_type == "I4" else text
        node.modified = True
        self._show_node(node)
        self.statusBar().showMessage(f"Modified {node.name} (unsaved)")

    def on_save(self):
        if getattr(self.doc, "_mmap_mode", False):
            QMessageBox.warning(
                self, "Read-only",
                "This file was opened via memory-map (large GPH). "
                "Saving in place is not supported.",
            )
            return
        if not self.doc.filepath:
            self.on_save_as()
            return
        if self.doc.save():
            self.statusBar().showMessage(f"Saved: {self.doc.filepath}")
            QMessageBox.information(self, "Saved", f"File saved: {self.doc.filepath}")
        else:
            QMessageBox.critical(self, "Error", "Failed to save")

    def on_save_as(self):
        if self._mmap_readonly:
            QMessageBox.warning(
                self, "Read-only",
                "This file was opened via memory-map (large GPH). "
                "Save As is not supported without loading the full file into RAM.",
            )
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save GPH As",
            self.doc.filepath or "",
            "GPH Files (*.gph);;All Files (*)",
        )
        if path:
            if self.doc.save(path):
                self.setWindowTitle(f"GPH Viewer - {Path(path).name}")
                self.statusBar().showMessage(f"Saved: {path}")
                QMessageBox.information(self, "Saved", f"File saved: {path}")
            else:
                QMessageBox.critical(self, "Error", "Failed to save")


def main():
    if HAS_PYQT is None:
        print("Error: PyQt5 or PyQt6 required. Install with: pip install PyQt6")
        sys.exit(1)

    app = QApplication(sys.argv)
    win = GphViewerMain()
    win.show()
    if len(sys.argv) > 1 and Path(sys.argv[1]).exists():
        win.load_file(sys.argv[1])
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
