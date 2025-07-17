import os, re
import sys
import numpy as np
import pandas as pd
import importlib.resources
from PIL import Image
import DashML.Basecall.run_basecall as run_basecall
import DashML.Landscape.Cluster.run_landscape as landscape
import DashML.Predict.run_predict as predict
import DashML.GUI
import DashML.Database_fx.Insert_DB as dbins
import DashML.Database_fx.Select_DB as dbsel
from PyQt6.QtWidgets import (
    QApplication, QWidget, QFrame, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QCheckBox, QLineEdit, QFileDialog, QSizePolicy,
    QScrollArea, QMessageBox, QComboBox
)
from PyQt6.QtCore import Qt,QTimer, QObject, pyqtSignal, QThread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class SelectSampleSection(QFrame):
    sample_changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.setFrameShape(QFrame.Shape.Box)
        self.setLineWidth(2)
        self.df = self.create_sample_df()
        self.lid = None
        self.contig = None
        self.init_ui()

    def create_sample_df(self):
        df = dbsel.select_library_full()
        return df

    def init_ui(self):
        layout = QVBoxLayout()
        title = QLabel("Add or Select Sample Sequence")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        layout.addWidget(title)

        id_contig_layout = QHBoxLayout()
        self.id_input = QLineEdit("")
        self.id_input.setFixedWidth(100)
        self.id_input.setReadOnly(True)
        self.id_input.setStyleSheet("background-color: #D3D3D3;")
        self.contig_combo = QComboBox()
        self.df['combo'] = self.df["ID"].astype(str) + " " + self.df["contig"] + " " + self.df["type1"] + " " + self.df["type2"]
        self.contig_combo.addItems(self.df['combo'])
        # Select first item properly
        if not self.df.empty:
            self.contig_combo.setCurrentText(self.df['combo'].iloc[0])
            self.lid = self.df["ID"].iloc[0]
            self.contig = self.df["contig"].iloc[0]

        self.contig_combo.currentTextChanged.connect(self.update_fields_from_contig)

        self.contig_input = QLineEdit()
        self.contig_input.hide()

        id_contig_layout.addWidget(QLabel("ID:"))
        id_contig_layout.addWidget(self.id_input)
        id_contig_layout.addSpacing(20)
        id_contig_layout.addWidget(QLabel("Contig (Select Existing Contig to Load Data):"))
        id_contig_layout.addWidget(self.contig_combo)
        id_contig_layout.addWidget(self.contig_input)
        layout.addLayout(id_contig_layout)

        # Sequence name comes first
        layout.addWidget(QLabel("Sequence Name (Reference Name in fasta file)"))
        self.sequence_name_input = QLineEdit(self.df["sequence_name"].iloc[0] if not self.df.empty else "")
        layout.addWidget(self.sequence_name_input)

        # Then nucleotide sequence
        layout.addWidget(QLabel("Sequence (Nucleotide Sequence)"))
        self.sequence_input = QLineEdit(self.df["sequence"].iloc[0] if not self.df.empty else "")
        self.sequence_input.setMinimumWidth(400)
        layout.addWidget(self.sequence_input)

        # structure and secondary
        layout.addWidget(QLabel("Secondary Structure (Optional Confirmed Secondary Structure w/o Psuedoknots)"))
        self.secondary_input = QLineEdit(self.df["secondary"].iloc[0] if not self.df.empty else "")
        self.secondary_input.setMinimumWidth(400)
        layout.addWidget(self.secondary_input)

        layout.addWidget(QLabel("Experiment (Optional Experimental Conf. eg Xray Crystallography)"))
        self.experiment_input = QLineEdit(self.df["experiment"].iloc[0] if not self.df.empty else "")
        self.experiment_input.setMinimumWidth(400)
        layout.addWidget(self.experiment_input)

        self.temperature_input = QLineEdit(str(self.df["temp"].iloc[0]) if not self.df.empty else "")
        self.condition_1_input = QLineEdit(self.df["type1"].iloc[0] if not self.df.empty else "")
        self.condition_2_input = QLineEdit(self.df["type2"].iloc[0] if not self.df.empty else "")
        self.experimental_run_input = QLineEdit(str(self.df["run"].iloc[0]) if not self.df.empty else "")

        self.fields = [
            ("Temperature (Default temperature is 37)", self.temperature_input),
            ("Condition 1 (Unmodified eg DMSO)", self.condition_1_input),
            ("Condition 2 (Probe Reagent eg 1M7, repeat Condition 1 if unmodified)", self.condition_2_input),
            ("Experimental Run (Multiple experiments on same sequence should have separate run numbers.)", self.experimental_run_input)
        ]
        for label, widget in self.fields:
            layout.addWidget(QLabel(label))
            layout.addWidget(widget)

        self.set_fields_enabled(False)

        self.new_button = QPushButton("New")
        self.new_button.clicked.connect(self.toggle_editable)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel_edit)
        self.cancel_button.hide()

        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.new_button)
        btn_layout.addWidget(self.cancel_button)
        layout.addLayout(btn_layout)

        self.setLayout(layout)
        # Initialize fields with first item
        if not self.df.empty:
            self.update_fields_from_contig(self.contig_combo.currentText())

    def validate_inputs_and_collect_errors(self):
        errors = []

        # Check required text fields
        sequence_txt = self.sequence_input.text().strip()
        if not re.fullmatch(r"^[ACGTU]+$", sequence_txt, re.IGNORECASE):
            errors.append("Sequence nucleotides must be ACGTU.")

        # Check required text fields
        if not self.validate_inputs():
            errors.append("All required fields must be filled.")

        # Check temperature is numeric
        temp_text = self.temperature_input.text().strip()
        if not re.fullmatch(r"-?\d+(\.\d+)?", temp_text):
            errors.append("Temperature must be a number.")

        # Check experimental run is numeric
        run_text = self.experimental_run_input.text().strip()
        if not re.fullmatch(r"-?\d+(\.\d+)?", run_text):
            errors.append("Experimental Run must be a number.")

        # Check secondary structure and experiment consistency
        secondary = self.secondary_input.text().strip()
        if not re.fullmatch(r"^[().]+$", secondary):
            errors.append("[red]Error: Secondary structure must contain only '.', '(', and ')' characters.[/red]")

        experiment = self.experiment_input.text().strip()

        if secondary:
            if not re.fullmatch(r"^[().]+$", secondary):
                errors.append("Secondary structure must be in valid dot-bracket notation (only '.', '(', ')').")

            if not experiment:
                errors.append("If specifying a control structure, an experiment must also be provided.")

            if len(self.sequence_input.text().strip()) != len(secondary):
                errors.append("Sequence and secondary structure must be the same length.")

        # Check if experiment is given without secondary (if needed)
        if experiment and not secondary:
            errors.append("If specifying an experiment, a secondary structure must also be provided.")

        return errors

    def toggle_editable(self):
        if self.new_button.text() == "New":
            self.clear_fields()
            self.set_fields_enabled(True)
            self.id_input.hide()
            self.contig_input.hide()
            self.contig_combo.hide()
            self.new_button.setText("Save")
            self.cancel_button.show()
        else:
            errors = self.validate_inputs_and_collect_errors()
            if errors:
                QMessageBox.warning(self, "Validation Error", "\n".join(errors))
                return
            self.add_seq()
            self.set_fields_enabled(False)
            self.id_input.setReadOnly(True)
            self.id_input.hide()
            self.contig_input.hide()
            self.contig_combo.show()
            self.new_button.setText("New")
            self.cancel_button.hide()

    def cancel_edit(self):
        self.set_fields_enabled(False)
        self.id_input.setReadOnly(True)
        self.id_input.show()
        self.contig_input.hide()
        self.contig_combo.show()
        self.new_button.setText("New")
        self.cancel_button.hide()
        self.update_fields_from_contig(self.contig_combo.currentText())

    def set_fields_enabled(self, enabled):
        for widget in [
            self.sequence_input, self.sequence_name_input, self.temperature_input,
            self.condition_1_input, self.condition_2_input,
            self.experimental_run_input, self.secondary_input, self.experiment_input
        ]:
            widget.setReadOnly(not enabled)

    def clear_fields(self):
        self.sequence_name_input.clear()
        self.sequence_input.clear()
        self.secondary_input.clear()
        self.experiment_input.clear()
        for _, widget in self.fields:
            widget.clear()

    def validate_inputs(self):
        widgets = [self.sequence_input, self.sequence_name_input, self.temperature_input,
                   self.condition_1_input, self.condition_2_input,
                   self.experimental_run_input]
        return all(w.text().strip() != "" for w in widgets)

    def add_seq(self):
        contig_value = self.contig_input.text().strip() if self.contig_input.isVisible() else self.contig_combo.currentText()
        try:
            dtr = pd.DataFrame.from_dict({
                "contig": [self.sequence_name_input.text()],
                "sequence": [self.sequence_input.text()],
                "secondary": [self.secondary_input.text()],
                "experiment": [self.experiment_input.text()],
                "sequence_name": [self.sequence_name_input.text()],
                "sequence_len": [len(self.sequence_input.text().strip())],
                "temp": [int(self.temperature_input.text())],
                "is_modified": [0 if self.condition_1_input.text() == self.condition_2_input.text() else 1],
                "type1": [self.condition_1_input.text()],
                "type2": [self.condition_2_input.text()],
                "complex": [0],  # TODO: Add complex input if needed
                "run": [int(self.experimental_run_input.text())]
            }, orient='columns')

            lid = dbins.insert_library(dtr)
            if lid==None:
                raise Exception("Error processing sequence. Please try again.")

            dtr["ID"] = lid
            #print(dtr.head())
            self.df = dbsel.select_library_full()

            # Compose new combo string
            contig_value = (str(dtr["ID"][0]) + " " + str(dtr["contig"][0]) + " " +
                            str(dtr["type1"][0]) + " " + str(dtr["type2"][0]))

            self.contig_combo.addItem(contig_value)
            self.contig_combo.setCurrentText(contig_value)
            QMessageBox.information(self, "Saved", "Sample added successfully.")
            self.sample_changed.emit()
        except Exception as error:
            QMessageBox.information(self, "Error", str(error))
            print("Failed to execute library: {}".format(error))
            self.cancel_edit()

    def update_fields_from_contig(self, contig):
        if not contig:
            return
        data = contig.split()
        if len(data) < 2:
            return
        try:
            lid = int(data[0].strip())
            contig_name = data[1].strip()
        except Exception:
            return
        if lid in self.df["ID"].values:
            row = self.df[self.df['ID'] == lid].iloc[0]
            self.lid = row["ID"]
            self.contig = row["contig"]
            self.id_input.setText(str(row["ID"]))
            self.sequence_input.setText(str(row["sequence"]))
            self.secondary_input.setText(str(row["secondary"]))
            self.experiment_input.setText(str(row["experiment"]))
            self.sequence_name_input.setText(str(row["sequence_name"]))
            self.temperature_input.setText(str(row["temp"]))
            self.condition_1_input.setText(str(row["type1"]))
            self.condition_2_input.setText(str(row["type2"]))
            self.experimental_run_input.setText(str(row["run"]))
            self.sample_changed.emit()

    def get_lid(self):
        return self.lid

    def get_contig(self):
        return self.contig

    def get_mod(self):
        if self.condition_1_input.text() == self.condition_2_input.text():
            return self.condition_1_input.text(), self.condition_2_input.text()
        else:
            return self.condition_2_input.text(), self.condition_1_input.text()


class SignalWorker(QObject):
    finished = pyqtSignal(str, object, object)  # source, gpath1, gpath2
    error = pyqtSignal(str)

    def __init__(self, source, lid, contig, signal_path, modification, modification2):
        super().__init__()
        self.lid = lid
        self.contig = contig
        self.signal_path = signal_path
        self.modification = modification
        self.modification2 = modification2

    def run(self):
        try:
            # Select only current contig for upload from tx
            cols = ['contig', 'position', 'reference_kmer', 'read_index',
                    'event_level_mean', 'event_length', 'event_stdv']
            df = pd.read_csv(self.signal_path, sep='\t', usecols=cols)
            df = df[df['contig'] == self.contig]
            df['LID'] = self.lid
            df['type1'] = self.modification
            df['type2'] = self.modification2

            dbins.insert_signal(df)

            # groups for display plotting
            df1 = df.groupby('position').agg(mean_val=('event_level_mean', 'mean')).reset_index()
            df2 = df.groupby('position').agg(mean_val=('event_length', 'mean')).reset_index()

            self.finished.emit("signal", df1, df2)
        except Exception as e:
            self.error.emit(str(e))


class BasecallWorker(QObject):
    finished = pyqtSignal(str, object, object)  # source, gpath1, gpath2
    error = pyqtSignal(str)

    def __init__(self, source, lid, contig, basecall_path, modification, modification2):
        super().__init__()
        self.lid = lid
        self.contig = contig
        self.basecall_path = basecall_path
        self.modification = modification
        self.modification2 = modification2

    def run(self):
        try:
            gpath1, gpath2 = run_basecall.get_modification(self.lid, self.contig,
                                                           self.basecall_path, self.modification, plot=True)
            self.finished.emit("basecall", gpath1, gpath2)
        except Exception as e:
            self.error.emit(str(e))


class LoadBasecall(QFrame):
    request_start_worker = pyqtSignal(dict)  # emit parameters as dict

    def __init__(self, sample_section):
        super().__init__()
        self.sample_section = sample_section
        self.sample_section.sample_changed.connect(self.update_sample_info)
        self.lid = self.sample_section.get_lid()
        self.contig = self.sample_section.get_contig()
        self.modification, self.modification2 = self.sample_section.get_mod()
        self.modification = self.modification.upper()
        self.modification2 = self.modification2.upper()
        self.path = None
        self.source = "basecall"

        self.basecall_thread = None
        self.basecall_worker = None

        self.init_ui()

    def init_ui(self):
        self.setFrameShape(QFrame.Shape.Box)
        self.setMinimumHeight(400)
        self.setLineWidth(2)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        layout = QVBoxLayout()

        # Title and ID row
        title_layout = QHBoxLayout()
        title = QLabel("Load Basecall Data")
        title.setStyleSheet("font-size: 20px; font-weight: bold;")
        title_layout.addWidget(title)

        title_layout.addStretch()  # Pushes the ID box to the right

        id_label = QLabel("Current ID:")
        self.id_display = QLineEdit()
        self.id_display.setText(str(self.lid))
        self.id_display.setReadOnly(True)
        self.id_display.setMaximumWidth(150)

        title_layout.addWidget(id_label)
        title_layout.addWidget(self.id_display)

        layout.addLayout(title_layout)

        # Basecall input group
        basecall_layout = QHBoxLayout()
        self.basecall_data_input = QLineEdit()
        basecall_button = QPushButton("Browse")
        basecall_button.clicked.connect(self.select_basecall_file)
        basecall_layout.addWidget(QLabel("Basecall Alignment Directory:"))
        basecall_layout.addWidget(self.basecall_data_input)
        basecall_layout.addWidget(basecall_button)

        layout.addLayout(basecall_layout)

        # Load buttons
        button_layout = QHBoxLayout()
        self.load_basecall_button = QPushButton("Load Basecall")
        self.load_basecall_button.clicked.connect(self.load_basecall)
        button_layout.addWidget(self.load_basecall_button)

        layout.addLayout(button_layout)

        # Basecall graph placeholders
        self.basecall_graphs_layout = QHBoxLayout()
        default_basecall1 = "default1.png"
        default_basecall2 = "default2.png"
        self.basecall_graph1 = self.create_sample_graph("Basecall: Avg Modification Rates", default_basecall1)
        self.basecall_graph2 = self.create_sample_graph("Basecall: Modification by Position", default_basecall2)
        self.basecall_graphs_layout.addWidget(self.basecall_graph1)
        self.basecall_graphs_layout.addWidget(self.basecall_graph2)

        layout.addLayout(self.basecall_graphs_layout)
        self.setLayout(layout)

    def select_basecall_file(self):
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Basecall Alignment Directory",
            options=QFileDialog.Option.ShowDirsOnly
        )
        if directory:
            self.basecall_data_input.setText(directory)

    def reset_ui(self):
        # Clear input field
        self.basecall_data_input.clear()

        # Replace graphs with placeholders
        for old_graph in [self.basecall_graph1, self.basecall_graph2]:
            self.basecall_graphs_layout.removeWidget(old_graph)
            old_graph.setParent(None)

        # Add default graphs
        default1 = "default1.png"
        default2 = "default2.png"
        self.basecall_graph1 = self.create_sample_graph("Basecall: Avg Modification Rates", default1)
        self.basecall_graph2 = self.create_sample_graph("Basecall: Modification by Position", default2)
        self.basecall_graphs_layout.addWidget(self.basecall_graph1)
        self.basecall_graphs_layout.addWidget(self.basecall_graph2)

        # Optionally disable load button until input is provided
        self.load_basecall_button.setEnabled(True)
    def update_sample_info(self):
        self.lid = self.sample_section.get_lid()
        self.contig = self.sample_section.get_contig()
        self.modification, self.modification2 = self.sample_section.get_mod()
        self.modification = self.modification.upper()
        self.modification2 = self.modification2.upper()
        self.id_display.setText(f"{self.lid} {self.contig}")
        self.id_display.setReadOnly(True)
        self.reset_ui()
    def load_basecall(self):
        self.load_basecall_button.setEnabled(False)
        self.lid = self.sample_section.get_lid()
        self.contig = self.sample_section.get_contig()
        self.modification, self.modification2 = self.sample_section.get_mod()
        self.modification = self.modification.upper()
        self.modification2 = self.modification2.upper()
        self.path = self.basecall_data_input.text()
        self.source = "basecall"

        self.basecall_thread = QThread()
        self.basecall_worker = BasecallWorker(
            self.source,
            self.sample_section.get_lid(),
            self.sample_section.get_contig(),
            self.path,
            self.modification,
            self.modification2
        )
        self.basecall_worker.moveToThread(self.basecall_thread)
        self.basecall_thread.started.connect(self.basecall_worker.run)
        self.basecall_worker.finished.connect(self.on_worker_finished)
        self.basecall_worker.error.connect(self.on_worker_error)
        self.basecall_worker.finished.connect(self.basecall_thread.quit)
        self.basecall_worker.finished.connect(self.basecall_worker.deleteLater)
        self.basecall_thread.finished.connect(self.basecall_thread.deleteLater)
        self.basecall_thread.start()

    def on_worker_finished(self, source, gpath1, gpath2):
        # Replace old graphs with new ones
        new_graph1 = self.create_sample_graph(f"{source.capitalize()}: Avg Modification Rates", gpath1)
        new_graph2 = self.create_sample_graph(f"{source.capitalize()}: Modification by Position", gpath2)

        # Remove old widgets
        for old_graph in [self.basecall_graph1, self.basecall_graph2]:
            self.basecall_graphs_layout.removeWidget(old_graph)
            old_graph.setParent(None)

        # Add new widgets
        self.basecall_graphs_layout.addWidget(new_graph1)
        self.basecall_graphs_layout.addWidget(new_graph2)

        # Save references
        self.basecall_graph1 = new_graph1
        self.basecall_graph2 = new_graph2

        QMessageBox.information(self, "Success", "Sample added successfully.")
        self.load_basecall_button.setEnabled(True)

    def on_worker_error(self, message):
        QMessageBox.critical(self, "Error", message)
        self.load_basecall_button.setEnabled(True)

    def create_sample_graph(self, label, gpath=None):
        fig = Figure(figsize=(4, 3))
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        ax.set_title(label, fontsize=10)
        ax.axis('off')

        image = None

        if isinstance(gpath, str):
            basename = os.path.basename(gpath)

            # Check if it's a default packaged image
            if basename.startswith("default"):
                try:
                    with importlib.resources.files("DashML.GUI").joinpath(basename).open("rb") as f:
                        image = mpimg.imread(f, format='png')
                except FileNotFoundError:
                    print(f"Default image '{basename}' not found in DashML.GUI package.")
                except Exception as e:
                    print(f"Error loading default image '{basename}': {e}")
            elif os.path.exists(gpath):
                try:
                    image = mpimg.imread(gpath)
                except Exception as e:
                    print(f"Error loading external image '{gpath}': {e}")

        if image is not None:
            ax.imshow(image, aspect='auto')
        else:
            ax.text(0.5, 0.5, "No Data", ha='center', va='center')

        fig.tight_layout()
        return canvas

class LoadSignal(QFrame):
    request_start_worker = pyqtSignal(dict)  # emit parameters as dict

    def __init__(self, sample_section):
        super().__init__()
        self.sample_section = sample_section
        self.sample_section.sample_changed.connect(self.update_sample_info)
        self.lid = self.sample_section.get_lid()
        self.contig = self.sample_section.get_contig()
        self.modification, self.modification2 = self.sample_section.get_mod()
        self.modification = self.modification.upper()
        self.modification2 = self.modification2.upper()
        self.path = None
        self.source = "signal"
        self.signal_thread = None
        self.signal_worker = None

        self.init_ui()

    def init_ui(self):
        self.setFrameShape(QFrame.Shape.Box)
        self.setMinimumHeight(400)
        self.setLineWidth(2)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        layout = QVBoxLayout()

        # Title and ID row
        title_layout = QHBoxLayout()
        title = QLabel("Load Signal Data")
        title.setStyleSheet("font-size: 20px; font-weight: bold;")
        title_layout.addWidget(title)

        title_layout.addStretch()  # Pushes the ID box to the right

        id_label = QLabel("Current ID:")
        self.id_display = QLineEdit()
        self.id_display.setText(str(self.lid))
        self.id_display.setReadOnly(True)
        self.id_display.setMaximumWidth(150)

        title_layout.addWidget(id_label)
        title_layout.addWidget(self.id_display)

        layout.addLayout(title_layout)

        # Label on one line
        layout.addWidget(QLabel("Signal Alignment File:"))

        # Input + Browse button on next line
        signal_input_layout = QHBoxLayout()
        self.signal_data_input = QLineEdit()
        signal_button = QPushButton("Browse")
        signal_button.clicked.connect(self.select_signal_file)
        signal_input_layout.addWidget(self.signal_data_input)
        signal_input_layout.addWidget(signal_button)
        layout.addLayout(signal_input_layout)

        # Load buttons
        button_layout = QHBoxLayout()
        self.load_signal_button = QPushButton("Load signal")
        self.load_signal_button.clicked.connect(self.load_signal)
        button_layout.addWidget(self.load_signal_button)

        layout.addLayout(button_layout)

        # signal graph placeholders
        self.signal_graphs_layout = QHBoxLayout()
        default_signal1 = "default3.png"
        default_signal2 = "default4.png"
        self.signal_graph1 = self.create_sample_graph("Signal: Avg Signal", default_signal1)
        self.signal_graph2 = self.create_sample_graph("Signal: Avg Dwell Time", default_signal2)
        self.signal_graphs_layout.addWidget(self.signal_graph1)
        self.signal_graphs_layout.addWidget(self.signal_graph2)

        layout.addLayout(self.signal_graphs_layout)
        self.setLayout(layout)

    def select_signal_file(self):
        file_dialog = QFileDialog(self, "Select Signal Alignment File")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setNameFilter("All Files (*)")

        if file_dialog.exec():
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                self.signal_data_input.setText(selected_files[0])

    def reset_ui(self):
        # Clear input field
        self.signal_data_input.clear()

        # Replace graphs with placeholders
        for old_graph in [self.signal_graph1, self.signal_graph2]:
            self.signal_graphs_layout.removeWidget(old_graph)
            old_graph.setParent(None)

        # Add default graphs
        default1 = "default3.png"
        default2 = "default4.png"
        self.signal_graph1 = self.create_sample_graph("Signal: Avg Signal", gpath=default1)
        self.signal_graph2 = self.create_sample_graph("Signal: Avg Dwell Time", gpath=default2)
        self.signal_graphs_layout.addWidget(self.signal_graph1)
        self.signal_graphs_layout.addWidget(self.signal_graph2)

        # Optionally disable load button until input is provided
        self.load_signal_button.setEnabled(True)

    def update_sample_info(self):
        self.lid = self.sample_section.get_lid()
        self.contig = self.sample_section.get_contig()
        self.modification, self.modification2 = self.sample_section.get_mod()
        self.modification = self.modification.upper()
        self.modification2 = self.modification2.upper()
        self.id_display.setText(f"{self.lid} {self.contig}")
        self.id_display.setReadOnly(True)
        self.reset_ui()
    def load_signal(self):
        self.load_signal_button.setEnabled(False)
        self.lid = self.sample_section.get_lid()
        self.contig = self.sample_section.get_contig()
        self.modification, self.modification2 = self.sample_section.get_mod()
        self.modification = self.modification.upper()
        self.modification2 = self.modification2.upper()
        self.path = self.signal_data_input.text()

        self.signal_thread = QThread()
        self.signal_worker = SignalWorker(self.source, self.lid, self.contig, self.path, self.modification, self.modification2)
        self.signal_worker.moveToThread(self.signal_thread)
        self.signal_thread.started.connect(self.signal_worker.run)
        self.signal_worker.finished.connect(self.on_worker_finished)
        self.signal_worker.error.connect(self.on_worker_error)
        self.signal_worker.finished.connect(self.signal_thread.quit)
        self.signal_worker.finished.connect(self.signal_worker.deleteLater)
        self.signal_thread.finished.connect(self.signal_thread.deleteLater)
        self.signal_thread.start()

    def on_worker_finished(self, source, gpath1, gpath2):
        # Replace old graphs with new ones
        new_graph1 = self.create_sample_graph(f"{source.capitalize()}: Avg Signal", gpath1)
        new_graph2 = self.create_sample_graph(f"{source.capitalize()}: Avg Dwell Time", gpath2)

        # Remove old widgets
        for old_graph in [self.signal_graph1, self.signal_graph2]:
            self.signal_graphs_layout.removeWidget(old_graph)
            old_graph.setParent(None)

        # Add new widgets
        self.signal_graphs_layout.addWidget(new_graph1)
        self.signal_graphs_layout.addWidget(new_graph2)

        # Save references
        self.signal_graph1 = new_graph1
        self.signal_graph2 = new_graph2

        QMessageBox.information(self, "Success", "Sample added successfully.")
        self.load_signal_button.setEnabled(True)

    def on_worker_error(self, message):
        QMessageBox.critical(self, "Error", message)
        self.load_signal_button.setEnabled(True)

    def create_sample_graph(self, label, gpath=None):
        fig = Figure(figsize=(4, 3))
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        ax.set_title(label, fontsize=10)

        image = None

        # Handle DataFrame case first
        if isinstance(gpath, pd.DataFrame):
            if not gpath.empty:
                if 'position' in gpath.columns and 'mean_val' in gpath.columns:
                    ax.scatter(gpath['position'], gpath['mean_val'], s=10, alpha=0.7)
                    ax.set_xlabel("Position", fontsize=10)
                    ax.set_ylabel("Mean Value", fontsize=10)
                    ax.tick_params(axis='both', which='major', labelsize=10)
                    ax.tick_params(axis='both', which='minor', labelsize=5)
                else:
                    ax.text(0.5, 0.5, "Invalid DataFrame format", ha='center', va='center')
                    ax.set_xticks([])
                    ax.set_yticks([])
            else:
                ax.text(0.5, 0.5, "Empty DataFrame", ha='center', va='center')
                ax.set_xticks([])
                ax.set_yticks([])
        elif isinstance(gpath, str):
            basename = os.path.basename(gpath)

            # Check if default image
            if basename.startswith("default"):
                try:
                    with importlib.resources.files("DashML.GUI").joinpath(basename).open("rb") as f:
                        image = mpimg.imread(f, format='png')
                except FileNotFoundError:
                    print(f"Default image '{basename}' not found in DashML.GUI package.")
                except Exception as e:
                    print(f"Error loading default image '{basename}': {e}")
            elif os.path.exists(gpath):
                try:
                    image = mpimg.imread(gpath)
                except Exception as e:
                    print(f"Error loading external image '{gpath}': {e}")

            if image is not None:
                ax.imshow(image, aspect='auto')
                ax.axis('off')
            else:
                ax.text(0.5, 0.5, "Image not found", ha='center', va='center')
                ax.set_xticks([])
                ax.set_yticks([])
        else:
            ax.text(0.5, 0.5, "No Data", ha='center', va='center')
            ax.set_xticks([])
            ax.set_yticks([])

        fig.tight_layout()
        return canvas


class PredictSection(QFrame):
    def __init__(self):
        super().__init__()
        self.setFrameShape(QFrame.Shape.Box)
        self.setMinimumHeight(400)
        self.setLineWidth(2)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.lid_unmod = None
        self.contig_unmod = None
        self.lid_mod = None
        self.contig_mod = None
        self.df= None

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        # Title
        title = QLabel("Predict")
        title.setStyleSheet("font-size: 20px; font-weight: bold;")
        layout.addWidget(title)

        # Dropdowns
        dropdown_layout = QHBoxLayout()
        self.unmod_combo = QComboBox()
        self.mod_combo= QComboBox()
        self.unmod_combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.mod_combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        self.df = dbsel.select_library_full()
        self.df['combo'] = (self.df["ID"].astype(str) + " " + self.df["contig"] + " "
                            + self.df["type1"] + " " + self.df["type2"])
        self.df_unmod = self.df.loc[self.df['is_modified'] == 0]
        self.df_mod = self.df.loc[self.df['is_modified'] == 1]
        self.unmod_combo.addItems(self.df_unmod['combo'])

        # Select first item properly
        if not self.df_unmod.empty:
            self.unmod_combo.setCurrentText(self.df_unmod['combo'].iloc[0])
            self.lid_unmod = self.df_unmod["ID"].iloc[0]
            self.contig_unmod = self.df_unmod["contig"].iloc[0]

        self.mod_combo.addItems(self.df_mod['combo'])
        # Select first item properly
        if not self.df_mod.empty:
            self.mod_combo.setCurrentText(self.df_mod['combo'].iloc[0])
            self.lid_mod = self.df_mod["ID"].iloc[0]
            self.contig_mod = self.df_mod["contig"].iloc[0]


        dropdown_layout.addWidget(QLabel("Select Unmodified:"))
        dropdown_layout.addWidget(self.unmod_combo)
        dropdown_layout.addWidget(QLabel("Select Modified:"))
        dropdown_layout.addWidget(self.mod_combo)
        layout.addLayout(dropdown_layout)

        # Buttons
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        self.btn_mod = QPushButton("Predict Modifications")
        self.btn_mod.clicked.connect(self.run_mod)

        # Add checkbox
        self.probability_checkbox = QCheckBox("Base Pairing Probability (Optional)")
        button_layout.addWidget(self.probability_checkbox)

        button_layout.addWidget(self.btn_mod)
        layout.addLayout(button_layout)

        # Graphs
        self.graph1_layout = QHBoxLayout()
        self.graph2_layout = QHBoxLayout()

        self.graph1 = self.create_sample_graph("default5.png")
        self.graph2 = self.create_sample_graph("default6.png")

        self.graph1_layout.addWidget(self.graph1)
        self.graph2_layout.addWidget(self.graph2)

        layout.addLayout(self.graph1_layout)
        layout.addLayout(self.graph2_layout)

        self.setLayout(layout)
        self.populate_dropdowns()

    def populate_dropdowns(self):
        if self.df_unmod is not None:
            items_unmod = self.df_unmod["combo"].astype(str).tolist()
            self.unmod_combo.clear()
            for item in items_unmod:
                self.unmod_combo.addItem(item)
            self.unmod_combo.setCurrentIndex(0)  # explicitly reset selection
            self.unmod_combo.update()  # force repaint

        if self.df_mod is not None:
            items_mod = self.df_mod["combo"].astype(str).tolist()
            self.mod_combo.clear()
            for item in items_mod:
                self.mod_combo.addItem(item)
            self.mod_combo.setCurrentIndex(0)  # explicitly reset selection
            self.mod_combo.update()  # force repaint

    def refresh_data(self):
        #print("PredictSection: refresh_data() called")
        self.df = dbsel.select_library_full()
        self.df['combo'] = (
                self.df["ID"].astype(str) + " " + self.df["contig"] + " "
                + self.df["type1"] + " " + self.df["type2"]
        )
        self.df_unmod = self.df.loc[self.df['is_modified'] == 0]
        self.df_mod = self.df.loc[self.df['is_modified'] == 1]

        # Clear and repopulate unmodified combo
        self.unmod_combo.clear()
        if not self.df_unmod.empty:
            self.unmod_combo.addItems(self.df_unmod["combo"])
            self.unmod_combo.setCurrentIndex(0)
            self.lid_unmod = self.df_unmod["ID"].iloc[0]
            self.contig_unmod = self.df_unmod["contig"].iloc[0]

        # Clear and repopulate modified combo
        self.mod_combo.clear()
        if not self.df_mod.empty:
            self.mod_combo.addItems(self.df_mod["combo"])
            self.mod_combo.setCurrentIndex(0)
            self.lid_mod = self.df_mod["ID"].iloc[0]
            self.contig_mod = self.df_mod["contig"].iloc[0]

    def create_sample_graph(self, image_path, label=None):
        image = None
        aspect_ratio = 1.5  # Default fallback

        # Check if image_path is a string
        if isinstance(image_path, str):
            basename = os.path.basename(image_path)

            # Check if this is a default image (starts with "default")
            if basename.startswith("default"):
                try:
                    # Use importlib.resources to get the image inside the package
                    with importlib.resources.files(DashML.GUI).joinpath(basename).open("rb") as f:
                        image = mpimg.imread(f, format='png')
                        height, width = image.shape[:2]
                        aspect_ratio = width / height if height != 0 else 1.5
                except FileNotFoundError:
                    print(f"Default image '{basename}' not found in DashML.GUI package.")
                except Exception as e:
                    print(f"Error loading default image '{basename}': {e}")
            else:
                # Load external file from filesystem
                if os.path.exists(image_path):
                    try:
                        image = mpimg.imread(image_path)
                        height, width = image.shape[:2]
                        aspect_ratio = width / height if height != 0 else 1.5
                    except Exception as e:
                        print(f"Error loading external image '{image_path}': {e}")
                else:
                    print(f"External image file '{image_path}' not found.")

        # Dynamically set figure size based on aspect ratio
        fig = Figure(figsize=(6 * aspect_ratio, 6))
        canvas = FigureCanvas(fig)
        canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        ax = fig.add_subplot(111)
        ax.axis('off')

        if label is not None:
            ax.set_title(label, fontsize=12, pad=10)

        if image is not None:
            ax.imshow(image, aspect='auto')
        else:
            ax.text(0.5, 0.5, "No image", ha='center', va='center')

        fig.tight_layout()

        return canvas
    def run_mod(self):
        self.run_prediction("Predict Modifications")

    def run_prediction(self, label):
        # Disable UI
        self.btn_mod.setEnabled(False)
        self.probability_checkbox.setEnabled(False)

        try:
            use_base_pairing = self.probability_checkbox.isChecked()
            self.lid_unmod = self.unmod_combo.currentText().strip().split()[0]
            self.contig_unmod = self.unmod_combo.currentText().strip().split()[1]
            self.lid_mod = self.mod_combo.currentText().strip().split()[0]
            self.contig_mod = self.mod_combo.currentText().strip().split()[1]

            seq_unmod = self.df_unmod.loc[self.df_unmod['ID'] == int(self.lid_unmod), 'sequence'].unique()[0]
            seq_mod = self.df_mod.loc[self.df_mod['ID'] == int(self.lid_mod), 'sequence'].unique()[0]

            if label == "Predict Modifications":
                if str(seq_mod) == str(seq_unmod):
                    image1, image2 = predict.run_predict(
                        unmod_lids=self.lid_unmod,
                        lids=self.lid_mod,
                        continue_reads=False,
                        vienna=use_base_pairing
                    )

                    if image1 is not None:
                        # Use singleShot if you want delayed UI update
                        QTimer.singleShot(3000, lambda: self.update_graphs(image1, image2, label))
                    else:
                        QMessageBox.critical(
                            self,
                            "Prediction Error",
                            "Prediction failed. Please check the input or try again."
                        )
                else:
                    QMessageBox.critical(
                        self,
                        "Prediction Error",
                        "Please select the unmod and modified versions of the same sequence and try again."
                    )
        except Exception as e:
            print(f"Prediction failed: {e}")
            QMessageBox.critical(self, "Prediction Error", f"An unexpected error occurred: {e}")
        finally:
            # Always re-enable UI controls
            self.btn_mod.setEnabled(True)
            self.probability_checkbox.setEnabled(True)
            self.probability_checkbox.setChecked(False)

    def update_graphs(self, image_path, image_path2, label):
        if label == "Predict Modifications":
            # Replace old graphs with new ones
            new_graph1 = self.create_sample_graph(image_path, label="Predicted Modifications")

            # Remove old graph1 from its layout
            self.graph1_layout.removeWidget(self.graph1)
            self.graph1.setParent(None)

            # Add the new graph1 to the same layout
            self.graph1_layout.addWidget(new_graph1)
            self.graph1 = new_graph1

            # Replace old graphs with new ones
            new_graph2 = self.create_sample_graph(image_path2, label="Predicted Modifications")

            # Remove old graph2 from its layout
            self.graph2_layout.removeWidget(self.graph2)
            self.graph2.setParent(None)

            # Add the new graph2 to the same layout
            self.graph2_layout.addWidget(new_graph2)
            self.graph2 = new_graph2

        # Show success message
        QMessageBox.information(self, "Prediction Complete", f"{label} prediction finished.")

        # Re-enable both buttons
        self.btn_mod.setEnabled(True)



class CreateLandscapeSection(QFrame):
    def __init__(self, parent=None):
        super().__init__()
        self.setFrameShape(QFrame.Shape.Box)
        self.setMinimumHeight(400)
        self.setLineWidth(2)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.lid_unmod = None
        self.contig_unmod = None
        self.lid_mod = None
        self.contig_mod = None
        self.df = dbsel.select_library_full()
        self.df['combo'] = (self.df["ID"].astype(str) + " " + self.df["contig"] + " "
                            + self.df["type1"] + " " + self.df["type2"])
        self.df_unmod = self.df.loc[self.df['is_modified'] == 0]
        self.df_mod = self.df.loc[self.df['is_modified'] == 1]

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        # Title
        title = QLabel("Create Landscape")
        title.setStyleSheet("font-size: 20px; font-weight: bold;")
        layout.addWidget(title)

        # Dropdowns
        dropdown_layout = QHBoxLayout()
        self.unmod_combo = QComboBox()
        self.mod_combo = QComboBox()
        self.unmod_combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.mod_combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.unmod_combo.addItems(self.df_unmod['combo'])
        # Select first item properly
        if not self.df_unmod.empty:
            self.unmod_combo.setCurrentText(self.df_unmod['combo'].iloc[0])
            self.lid_unmod = self.df_unmod["ID"].iloc[0]
            self.contig_unmod = self.df_unmod["contig"].iloc[0]

        self.mod_combo.addItems(self.df_mod['combo'])
        # Select first item properly
        if not self.df_mod.empty:
            self.mod_combo.setCurrentText(self.df_mod['combo'].iloc[0])
            self.lid_mod = self.df_mod["ID"].iloc[0]
            self.contig_mod = self.df_mod["contig"].iloc[0]

        dropdown_layout.addWidget(QLabel("Select Unmodified:"))
        dropdown_layout.addWidget(self.unmod_combo)
        dropdown_layout.addWidget(QLabel("Select Modified:"))
        dropdown_layout.addWidget(self.mod_combo)
        layout.addLayout(dropdown_layout)

        # Create button
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        # Add checkbox
        self.optimize_checkbox = QCheckBox("Optimize Cluster (Optional Default 20)")
        button_layout.addWidget(self.optimize_checkbox)

        self.create_button = QPushButton("Create Landscape")
        self.create_button.setFixedSize(200, 28)
        self.create_button.clicked.connect(self.create_landscape)
        button_layout.addWidget(self.create_button)

        layout.addLayout(button_layout)

        # Graphs
        self.graphs_layout = QHBoxLayout()
        self.graphs_layout.setSpacing(6)

        self.graph1 = self.create_sample_graph("default7.png")
        self.graph2 = self.create_sample_graph("default7.png")
        self.graph3 = self.create_sample_graph("default7.png")
        self.graph4 = self.create_sample_graph("default8.png")

        self.graphs_layout.addWidget(self.graph1)
        self.graphs_layout.addWidget(self.graph2)
        self.graphs_layout.addWidget(self.graph3)
        self.graphs_layout.addWidget(self.graph4)

        layout.addLayout(self.graphs_layout)
        self.setLayout(layout)

    def populate_dropdowns(self):
        if self.df_unmod is not None:
            self.unmod_combo.clear()
            self.unmod_combo.addItems(self.df_unmod["combo"].astype(str).tolist())
            if not self.df_unmod.empty:
                self.unmod_combo.setCurrentIndex(0)
        if self.df_mod is not None:
            self.mod_combo.clear()
            self.mod_combo.addItems(self.df_mod["combo"].astype(str).tolist())
            if not self.df_mod.empty:
                self.mod_combo.setCurrentIndex(0)

    def refresh_data(self):
        # print("PredictSection: refresh_data() called")
        self.df = dbsel.select_library_full()
        self.df['combo'] = (
                self.df["ID"].astype(str) + " " + self.df["contig"] + " "
                + self.df["type1"] + " " + self.df["type2"]
        )
        self.df_unmod = self.df.loc[self.df['is_modified'] == 0]
        self.df_mod = self.df.loc[self.df['is_modified'] == 1]

        # Clear and repopulate unmodified combo
        self.unmod_combo.clear()
        if not self.df_unmod.empty:
            self.unmod_combo.addItems(self.df_unmod["combo"])
            self.unmod_combo.setCurrentIndex(0)
            self.lid_unmod = self.df_unmod["ID"].iloc[0]
            self.contig_unmod = self.df_unmod["contig"].iloc[0]

        # Clear and repopulate modified combo
        self.mod_combo.clear()
        if not self.df_mod.empty:
            self.mod_combo.addItems(self.df_mod["combo"])
            self.mod_combo.setCurrentIndex(0)
            self.lid_mod = self.df_mod["ID"].iloc[0]
            self.contig_mod = self.df_mod["contig"].iloc[0]

    def create_landscape(self):
        self.create_button.setEnabled(False)
        try:
            optimize_clusters = self.optimize_checkbox.isChecked()
            self.lid_unmod = self.unmod_combo.currentText().strip().split()[0]
            self.contig_unmod = self.unmod_combo.currentText().strip().split()[1]
            self.lid_mod = self.mod_combo.currentText().strip().split()[0]
            self.contig_mod = self.mod_combo.currentText().strip().split()[1]
            seq_unmod = self.df_unmod.loc[self.df_unmod['ID'] == int(self.lid_unmod), 'sequence'].unique()[0]
            seq_mod = self.df_mod.loc[self.df_mod['ID'] == int(self.lid_mod), 'sequence'].unique()[0]

            if str(seq_mod) == str(seq_unmod):
                images = landscape.run_landscape(
                    unmod_lid=self.lid_unmod,
                    lid=self.lid_mod,
                    optimize_clusters=optimize_clusters
                )

                if images is not None:
                    QTimer.singleShot(1000, lambda: self.update_graphs(images))
                else:
                    QMessageBox.critical(
                        self,
                        "Landscape Error",
                        "Landscape failed. Please check the input or try again."
                    )
            else:
                QMessageBox.critical(
                    self,
                    "Landscape Error",
                    "Please select the unmod and modified versions of the same sequence and try again."
                )
        except Exception as e:
            print(f"Landscape creation failed: {e}")
            QMessageBox.critical(self, "Landscape Error", f"An unexpected error occurred: {e}")
        finally:
            self.create_button.setEnabled(True)

    def update_graphs(self, images):
        # Remove old graphs from layout and UI
        for graph in [self.graph1, self.graph2, self.graph3, self.graph4]:
            self.graphs_layout.removeWidget(graph)
            graph.setParent(None)

        # Create new graphs using updated images
        self.graph1 = self.create_sample_graph(images[0], label="HeatMap")
        self.graph2 = self.create_sample_graph(images[1], label="Dendrogram")
        self.graph3 = self.create_sample_graph(images[2], label="Read Corr")
        self.graph4 = self.create_sample_graph(images[3], label="Cluster Opt.")

        # Add them back to the layout
        self.graphs_layout.addWidget(self.graph1)
        self.graphs_layout.addWidget(self.graph2)
        self.graphs_layout.addWidget(self.graph3)
        self.graphs_layout.addWidget(self.graph4)

        # Done
        QMessageBox.information(self, "Landscape Complete", "Landscape analysis finished.")
        self.create_button.setEnabled(True)

    def create_sample_graph(self, image_path, label=None):
        # Check if this is a default image
        if os.path.basename(image_path).startswith("default"):
            try:
                # Load from DashML.GUI package
                with importlib.resources.files("DashML.GUI").joinpath(os.path.basename(image_path)).open("rb") as f:
                    image = mpimg.imread(f, format='png')
            except FileNotFoundError:
                print(f"Default image {image_path} not found in DashML.GUI package.")
                image = None
        else:
            # Load from local filesystem
            if os.path.exists(image_path):
                image = mpimg.imread(image_path)
            else:
                print(f"Image file {image_path} not found.")
                image = None

        # Fixed figure size to ensure all graphs have the same size
        fig = Figure(figsize=(6, 6))  # 6x6 inches fixed size
        canvas = FigureCanvas(fig)
        canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        ax = fig.add_subplot(111)
        ax.axis('off')  # Hide axes lines and ticks

        if label is not None:
            ax.set_title(label, fontsize=12, fontweight='bold', pad=10)

        if image is not None:
            ax.imshow(image, aspect='auto')  # Scale naturally
        else:
            ax.text(0.5, 0.5, "No image", ha='center', va='center')

        # Adjust layout to make room for title
        fig.tight_layout()
        fig.subplots_adjust(top=0.85)  # Adjust top margin so title doesn't get cut off

        return canvas


class MainApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Dashing Turtle RNA Landscape Analysis")
        self.setGeometry(100, 100, 1200, 1000)

        # Create shared sample section
        self.sample_section = SelectSampleSection()

        # Pass shared sample section to LoadBasecall
        self.load_basecall_section = LoadBasecall(self.sample_section)
        self.load_basecall_section.request_start_worker.connect(self.handle_start_worker)

        # Pass shared sample section to LoadSignal
        self.load_signal_section = LoadSignal(self.sample_section)
        self.load_signal_section.request_start_worker.connect(self.handle_start_worker)

        #predict
        self.predict_section = PredictSection()
        self.sample_section.sample_changed.connect(self.predict_section.refresh_data)

        #landscape
        self.landscape_section = CreateLandscapeSection()
        self.sample_section.sample_changed.connect(self.landscape_section.refresh_data)

        self.init_ui()

        # Keep references to threads and workers to avoid premature GC
        # self.basecall_thread = None
        # self.basecall_worker = None
        # self.signal_thread = None
        # self.signal_worker = None

    def init_ui(self):
        scroll = QScrollArea()
        container = QWidget()
        layout = QVBoxLayout()

        # Add sections
        layout.addWidget(self.sample_section)
        layout.addWidget(self.load_basecall_section)
        layout.addWidget(self.load_signal_section)
        layout.addWidget(self.predict_section)
        layout.addWidget(self.landscape_section)

        container.setLayout(layout)

        # Apply responsive sizing
        container.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        container.adjustSize()

        scroll.setWidget(container)
        scroll.setWidgetResizable(True)

        main_layout = QVBoxLayout()
        main_layout.addWidget(scroll)
        self.setLayout(main_layout)

    def handle_start_worker(self, params):
        thread = QThread()

        source = params.get('source')
        if source == 'basecall':
            worker = BasecallWorker(
                params['lid'], params['contig'], params['basecall_path'],
                params['modification'], params['modification2']
            )
        elif source == 'signal':
            worker = SignalWorker(
                params['lid'], params['contig'], params['signal_path'],
                params['modification'], params['modification2']
            )
        else:
            print(f"Unknown source: {source}")
            return

        self._start_worker(thread, worker)

    def _start_worker(self, thread, worker):
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self.on_worker_finished)
        worker.error.connect(self.on_worker_error)
        worker.finished.connect(thread.quit)
        thread.finished.connect(thread.deleteLater)
        worker.finished.connect(worker.deleteLater)
        thread.start()

    def on_worker_finished(self, source, gpath1, gpath2):
        self.load_basecall_section.load_basecall_button.setEnabled(True)
        self.load_signal_section.load_signal_button.setEnabled(True)

        if source == "basecall":
            layout = self.load_basecall_section.basecall_graphs_layout
            old1, old2 = self.load_basecall_section.basecall_graph1, self.load_basecall_section.basecall_graph2
            new_graph1 = self.load_basecall_section.create_sample_graph(f"{source.capitalize()}: Avg Modification Rates",
                                                                    gpath1)
            new_graph2 = self.load_basecall_section.create_sample_graph(f"{source.capitalize()}: Modification by Position",
                                                                    gpath2)
            self.load_basecall_section.basecall_graph1 = new_graph1
            self.load_basecall_section.basecall_graph2 = new_graph2
        else:
            layout = self.load_signal_section.signal_graphs_layout
            old1, old2 = self.load_signal_section.signal_graph1, self.load_basecall_section.signal_graph2
            new_graph1 = self.load_signal_section.create_sample_graph(
                f"{source.capitalize()}: Average Signal Rates by Position", gpath1)
            new_graph2 = self.load_signal_section.create_sample_graph(
                f"{source.capitalize()}: Average Dwell Time by Position", gpath2)
            self.load_signal_section.signal_graph1 = new_graph1
            self.load_signal_section.signal_graph2 = new_graph2

        layout.replaceWidget(old1, new_graph1)
        layout.replaceWidget(old2, new_graph2)

        old1.setParent(None)
        old2.setParent(None)

        # Clear references to allow garbage collection
        if source == "basecall":
            self.basecall_worker = None
            self.basecall_thread = None
        elif source == "signal":
            self.signal_worker = None
            self.signal_thread = None

    def on_worker_error(self, message):
        self.load_basecall_section.load_basecall_button.setEnabled(True)
        self.load_signal_section.load_signal_button.setEnabled(True)

def main():
     print("Launching GUI...")
     app = QApplication(sys.argv)
     window = MainApp()
     window.show()
     sys.exit(app.exec())

if __name__ == "__main__":
    main()
