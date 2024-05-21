from PyQt5 import QtWidgets, QtCore
import sys, os


#### add in applets with python "C:\\[your directory path]\\qn_artiq_routines\\applets\\experiment_variables_visualization.py"
class KeyValueTable(QtWidgets.QTableWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.widge = parent

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_R:  # Refresh on 'R' key press
            self.widge.refresh_data()
        else:
            super().keyPressEvent(event)

class KeyValueTableWidget(QtWidgets.QWidget):
    def __init__(self, key_value_pairs, pairs_per_row=3, path_to_db = "dataset_db.pyon"):
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

        # Make the table scrollable
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidget(self.table)
        self.scroll.setWidgetResizable(True)
          # Set fixed height for scrolling

        layout.addWidget(self.scroll)
        self.setLayout(layout)

    def update_table(self, key_value_pairs):
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

    def refresh_data(self):
        simple_pairs_dict = load_and_extract(self.path)
        if self.pressCounter % 5 == 0:
            print("Refreshing table: if not updated, wait before pressing [R] again.\nEvery [R] press refreshes but this"
              " message only shows ever 5 presses")
        self.pressCounter+=1
        self.update_table(simple_pairs_dict)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_R:  # Refresh on 'R' key press
            self.refresh_data()

def extract_simple_pairs(input_dict):
    simple_pairs = {}
    for key, value in input_dict.items():
        if not isinstance(value, (dict, list)):
            simple_pairs[key] = value
    return simple_pairs

def load_and_extract(path_to_db):
    with open(path_to_db) as f:
        datasets_str = f.read()
        # when the pyon file is saved python True and False are converted to lowercase...
        datasets_str = datasets_str.replace("true", "True")
        datasets_str = datasets_str.replace("false", "False")
        return extract_simple_pairs(eval(datasets_str))

def main():
    # File path-to
    dir_path = os.path.dirname(os.path.realpath(__file__))
    path_to_db = os.path.join(dir_path, "..", "..", "..", "dataset_db.pyon")
    app = QtWidgets.QApplication(sys.argv)  # Create the widget application

    # Initial load of the data
    #Extracts all simple key-value pairs for dataset_db.pyon
    simple_pairs_dict = load_and_extract(path_to_db)

    widget = KeyValueTableWidget(simple_pairs_dict, pairs_per_row = 5, path_to_db=path_to_db)
    widget.setWindowTitle("Key-Value Pairs")
    widget.resize(800, 300)
    widget.show()

    def close_app():
        sys.exit(app.exec_())

    print("Refreshing table: if not updated, wait before pressing [R] again.\nEvery [R] press refreshes but this"
          " message only shows ever 5 presses")
    app.aboutToQuit.connect(close_app)
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()

