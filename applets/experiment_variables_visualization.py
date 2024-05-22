from PyQt5 import QtWidgets, QtCore, QtGui
import sys, os

class KeyValueTable(QtWidgets.QTableWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.widge = parent
        self.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        self.is_selecting = False
        self.highlighted_cells = set()

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self.start_pos = self.indexAt(event.pos())
            self.end_pos = self.start_pos
            self.clearSelection()
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if event.buttons() & QtCore.Qt.LeftButton:
            self.end_pos = self.indexAt(event.pos())
            self.updateSelection()
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self.end_pos = self.indexAt(event.pos())
            self.updateSelection(final=True)
            self.start_pos = None
            self.end_pos = None
        super().mouseReleaseEvent(event)

    def updateSelection(self, final=False):
        try:
            if not self.start_pos or not self.end_pos:
                return

            top_left = self.start_pos
            bottom_right = self.end_pos

            if top_left.row() > bottom_right.row():
                top_left, bottom_right = bottom_right, top_left

            for row in range(top_left.row(), bottom_right.row() + 1):
                for col in range(top_left.column(), bottom_right.column() + 1):
                    index = self.model().index(row, col)
                    if index.isValid():
                        item = self.itemFromIndex(index)
                        if final:
                            if (row, col) in self.highlighted_cells:
                                self.highlighted_cells.remove((row, col))
                                item.setBackground(QtGui.QColor('white'))
                            else:
                                self.highlighted_cells.add((row, col))
                                item.setBackground(QtGui.QColor('yellow'))
                        else:
                            item.setBackground(QtGui.QColor('yellow'))

        except Exception as e:
            1+1



class KeyValueTableWidget(QtWidgets.QWidget):
    def __init__(self, key_value_pairs, pairs_per_row=3, path_to_db="dataset_db.pyon"):
        super().__init__()
        self.path = path_to_db
        self.key_value_pairs = key_value_pairs
        self.pairs_per_row = pairs_per_row
        self.table = KeyValueTable(parent=self)
        self.initUI()
        self.pressCounter = 0

    def initUI(self):
        layout = QtWidgets.QVBoxLayout()

        # Create a table
        self.table.setColumnCount(self.pairs_per_row * 2)
        headers = []
        for i in range(self.pairs_per_row):
            headers.extend([f"Expmt. Variable", f"Value"])
        self.table.setHorizontalHeaderLabels(headers)

        # Populate the table
        self.update_table(self.key_value_pairs)

        # Create a refresh button
        self.refresh_button = QtWidgets.QPushButton("Refresh")
        self.refresh_button.clicked.connect(self.refresh_data)
        self.clear_button = QtWidgets.QPushButton("Clear Highlights")
        self.clear_button.clicked.connect(self.clear_highlights)
        # Make the table scrollable
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidget(self.table)
        self.scroll.setWidgetResizable(True)

        # Add widgets to layout
        layout.addWidget(self.scroll)
        layout.addWidget(self.refresh_button)
        layout.addWidget(self.clear_button)
        self.setLayout(layout)

    def update_table(self, key_value_pairs):
        self.key_value_pairs = key_value_pairs
        keys = list(key_value_pairs.keys())
        values = list(key_value_pairs.values())
        row_count = (len(keys) + self.pairs_per_row - 1) // self.pairs_per_row
        self.table.setRowCount(row_count)

        for row in range(row_count):
            for col in range(self.pairs_per_row):
                key_index = row * self.pairs_per_row + col
                if key_index < len(keys):
                    key_item = QtWidgets.QTableWidgetItem(keys[key_index])
                    value_item = QtWidgets.QTableWidgetItem(str(values[key_index]))
                    key_item.setFlags(key_item.flags() & ~QtCore.Qt.ItemIsEditable)
                    value_item.setFlags(value_item.flags() & ~QtCore.Qt.ItemIsEditable)
                    self.table.setItem(row, col * 2, key_item)
                    self.table.setItem(row, col * 2 + 1, value_item)

    def clear_highlights(self):
        columns = self.pairs_per_row*2

        for i in range(len(self.key_value_pairs)*2):
            try:
                column = i%(self.pairs_per_row*2)
                row = int(i/(self.pairs_per_row*2))
                self.table.item(row, column).setBackground(QtGui.QColor('white'))
            except:
                1+1
    def refresh_data(self):
        # Store highlighted keys
        highlighted_keys = set()
        for row in range(self.table.rowCount()):
            for col in range(self.table.columnCount() // 2):
                key_item = self.table.item(row, col * 2)
                value_item = self.table.item(row, col * 2 + 1)
                
                if key_item and key_item.background() == QtGui.QColor('yellow') or \
                        (value_item and value_item.background() == QtGui.QColor('yellow')):
                    highlighted_keys.add(key_item.text())

        # Load new data and update the table
        simple_pairs_dict = load_and_extract(self.path)
        self.update_table(simple_pairs_dict)

        # Reapply highlights
        keys = list(simple_pairs_dict.keys())
        row_count = (len(keys) + self.pairs_per_row - 1) // self.pairs_per_row

        for row in range(row_count):
            for col in range(self.pairs_per_row):
                key_index = row * self.pairs_per_row + col
                if key_index < len(keys):
                    key = keys[key_index]
                    if key in highlighted_keys:
                        self.table.item(row, col * 2).setBackground(QtGui.QColor('yellow'))
                        self.table.item(row, col * 2 + 1).setBackground(QtGui.QColor('yellow'))
                        self.table.highlighted_cells.add((row, col))

        if self.pressCounter % 5 == 0:
            print(
                "Refreshing table: if not updated, wait before pressing the Refresh button again.\nEvery refresh updates the table, but this message only shows every 5 presses")
        self.pressCounter += 1

def extract_simple_pairs(input_dict):
    simple_pairs = {}
    for key, value in input_dict.items():
        if not isinstance(value, (dict, list)):
            simple_pairs[key] = value
    return simple_pairs

def load_and_extract(path_to_db):
    with open(path_to_db) as f:
        datasets_str = f.read()
        datasets_str = datasets_str.replace("true", "True")
        datasets_str = datasets_str.replace("false", "False")
        return extract_simple_pairs(eval(datasets_str))

def main():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    path_to_db = os.path.join(dir_path, "..", "..", "..", "dataset_db.pyon")
    app = QtWidgets.QApplication(sys.argv)

    simple_pairs_dict = load_and_extract(path_to_db)

    widget = KeyValueTableWidget(simple_pairs_dict, pairs_per_row=5, path_to_db=path_to_db)
    widget.setWindowTitle("Experiment Variable Tracker")
    widget.resize(800, 300)
    widget.show()

    def close_app():
        sys.exit(app.exec_())

    print("Refreshing table: if not updated, wait before pressing the Refresh button again.\n"
          "Every refresh updates the table, but this message only shows every 5 presses.")
    app.aboutToQuit.connect(close_app)
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
