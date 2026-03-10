#!/usr/bin/env python3
"""
GPH Viewer - PyQt GUI to browse and edit GPH binary files (similar to HDFView).
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
        QMenuBar,
        QMenu,
        QFileDialog,
        QStatusBar,
        QLabel,
        QGroupBox,
        QLineEdit,
        QPushButton,
        QMessageBox,
        QHeaderView,
        QTabWidget,
    )
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QAction
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
            QMenuBar,
            QMenu,
            QFileDialog,
            QStatusBar,
            QLabel,
            QGroupBox,
            QLineEdit,
            QPushButton,
            QMessageBox,
            QHeaderView,
            QTabWidget,
        )
        from PyQt5.QtCore import Qt
        from PyQt5.QtWidgets import QAction
        HAS_PYQT = "PyQt5"
    except ImportError:
        HAS_PYQT = None

from gph_model import GphDocument, GphNode


def read_f32_be(data: bytes, pos: int) -> float:
    return struct.unpack(">f", data[pos : pos + 4])[0]


def read_i32_be(data: bytes, pos: int) -> int:
    return int.from_bytes(data[pos : pos + 4], "big")


class GphViewerMain(QMainWindow):
    def __init__(self):
        super().__init__()
        self.doc = GphDocument()
        self.current_node: GphNode | None = None
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

        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        splitter = QSplitter()

        # Left: Tree
        tree_box = QGroupBox("Structure")
        tree_layout = QVBoxLayout(tree_box)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Name", "Type", "Size"])
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
        self._build_tree()
        self.setWindowTitle(f"GPH Viewer - {Path(path).name}")
        self.statusBar().showMessage(f"Opened: {path} ({len(self.doc._raw_data)} bytes)")

    def _build_tree(self):
        self.tree.clear()
        root = self.doc.root
        if not root:
            return
        item = QTreeWidgetItem([root.name, root.data_type, str(root.size)])
        item.setData(0, 100, root)
        self._add_children(item, root)
        self.tree.addTopLevelItem(item)
        self.tree.expandToDepth(1)

    def _add_children(self, parent_item: QTreeWidgetItem, node: GphNode):
        for child in node.children:
            item = QTreeWidgetItem([child.name, child.data_type, str(child.size)])
            item.setData(0, 100, child)
            parent_item.addChild(item)
            self._add_children(item, child)

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
            f"Size: {node.size} bytes",
            f"Type: {node.data_type}",
            "",
        ]
        if node.value is not None:
            if isinstance(node.value, (list, tuple)):
                if len(node.value) <= 20:
                    lines.append(f"Value: {node.value}")
                else:
                    lines.append(f"Value: [{len(node.value)} items] (see Data tab)")
            else:
                lines.append(f"Value: {node.value}")
        if node.modified:
            lines.append("\n(Modified)")
        self.attrs_text.setPlainText("\n".join(lines))

        # Data tab (table for arrays)
        self._table_loading = True
        if isinstance(node.value, (list, tuple)) and node.value:
            first = node.value[0]
            if isinstance(first, (int, float)):
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
        elif node.data_type in ("I4[135]", "R4[n,3]", "I4[]") and node.value:
            self.edit_input.setEnabled(True)
            self.edit_input.setPlaceholderText("Use Data tab; edit not yet implemented for arrays")
            self.edit_input.clear()
            self.apply_btn.setEnabled(False)
        else:
            self.edit_input.setEnabled(False)
            self.edit_input.clear()
            self.apply_btn.setEnabled(False)

    def on_table_cell_changed(self, item: QTableWidgetItem):
        if self._table_loading or not self.current_node:
            return
        node = self.current_node
        if node.data_type != "R4[n,3]" or not isinstance(node.value, (list, tuple)):
            return
        row, col = item.row(), item.column()
        if row >= len(node.value) or col >= 3:
            return
        try:
            val = float(item.text())
        except ValueError:
            return
        # LS_Nodes: data at 0x2750, vertex i has x,y,z at i*12, i*12+4, i*12+8
        file_offset = 0x2750 + row * 12 + col * 4
        import struct
        new_raw = struct.pack(">f", val)
        self.doc.apply_patch(file_offset, new_raw)
        old_row = list(node.value[row])
        old_row[col] = val
        node.value[row] = tuple(old_row)
        node.modified = True
        self.statusBar().showMessage(f"Modified {node.name}[{row},{col}] = {val} (unsaved)")

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
        if not self.doc.filepath:
            self.on_save_as()
            return
        if self.doc.save():
            self.statusBar().showMessage(f"Saved: {self.doc.filepath}")
            QMessageBox.information(self, "Saved", f"File saved: {self.doc.filepath}")
        else:
            QMessageBox.critical(self, "Error", "Failed to save")

    def on_save_as(self):
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
