# Python 3

import os
import re
import sys
import _thread
from PyQt5 import QtCore, QtGui, QtWidgets, uic

import structs.sequence as sequence
import structs.spectra as spectra

import structs.mass_histogram as histogram

ui_path = os.path.join('ui', 'mac', 'main.ui')
if sys.platform == 'win32':
    ui_path = os.path.join('ui', 'win', 'main.ui')

form_class, base_class = uic.loadUiType(ui_path)

class MyTableWidgetItem(QtWidgets.QTableWidgetItem):
    def __init__(self, *args):
        super(MyTableWidgetItem, self).__init__(*args)
        self.id = ''

    def get_id(self):
        return self.id

class SortTableWidgetItem(QtWidgets.QTableWidgetItem):
    def __init__(self, *args):
        super(SortTableWidgetItem, self).__init__(*args)

    def __lt__(self, other):
        if isinstance(other, SortTableWidgetItem):
            self_data_value = float(str(self.data(QtCore.Qt.EditRole)))
            other_data_value = float(str(other.data(QtCore.Qt.EditRole)))
            return self_data_value < other_data_value
        else:
            return QtGui.QTableWidgetItem.__lt__(self, other)
                

class InterLink(QtWidgets.QMainWindow, form_class):
    def __init__(self, *args):
        super(InterLink, self).__init__(*args)
        self.setupUi(self)
        self.protein = None
        self.spectra = None
        self.window = self.window()
        self.window.setGeometry(5, 60, self.width(), self.height())
        self.setup_initial_toolbar()
        self.setup_protease_list()
        self.setup_signals()
        self.search_masses = False
        # self.peptides_list.setColumnWidth(0, 1060)
        # self.peptides_list.setColumnWidth(1, 150)
        # self.peptides_list.setColumnWidth(2, 150)
        self.histogram = histogram.View(self, width=5, height=4, dpi=100)
        self.histogram_layout.addWidget(self.histogram)
        # self.load_test_files()

    def setup_signals(self):
        self.protein_sequence.textChanged.connect(self.protein_sequence_edited)

        self.protease_list.currentIndexChanged.connect(self.cleave_protein)
        self.max_missed_cleavage.valueChanged.connect(self.cleave_protein)
        self.mass_list.itemSelectionChanged.connect(self.mass_spectra_selected)
        self.peptides_list.itemSelectionChanged.connect(self.peptide_selected)

        for checkbox in self.ion_type_group.findChildren(QtWidgets.QCheckBox):
            checkbox.clicked.connect(self.fragment_peptides)
        for radiobutton in self.charge_state_group.findChildren(QtWidgets.QRadioButton):
            radiobutton.clicked.connect(self.fragment_peptides)
        for radiobutton in self.crosslinker_group.findChildren(QtWidgets.QRadioButton):
            radiobutton.clicked.connect(self.cleave_protein)

        self.actionShow_Fragments_List.triggered.connect(self.show_hide_fragments_list)
        self.actionImport_Peptide.triggered.connect(self.import_peptide)
        self.actionImport_Spectra.triggered.connect(self.import_mass_spectra)


    def show_hide_fragments_list(self):
        if self.side_dock.isHidden():
            self.side_dock.setHidden(False)
        else:
            self.side_dock.setHidden(True)

    def protein_sequence_edited(self):
        self.peptides_list.clearContents()
        self.fragments_list.clearContents()
        self.fragments_list.setRowCount(0)
        if not self.protein:
            self.protein = sequence.protein()
        text = str(self.protein_sequence.toPlainText())
        self.protein.sequence = re.sub('\n| ', '', text)
        self.protease_selected()

    def setup_protease_list(self):
        self.protease_list.addItems(['Trypsin', 'Chymotrypsin high specificity', 'Chymotrypsin low specificity'])
        
    def cleave_protein(self):
        selection = str(self.protease_list.currentText())
        self.protein.cleave_sequence(selection.lower(), self.max_missed_cleavage.value(),
                                     crosslinker='bpa' if self.bpa.isChecked() else 'bpa_alk')
        self.update_peptides_list()
    
    def update_peptides_list(self):
        self.peptides_list.clearContents()
        self.fragments_list.clearContents()
        self.peptides_list.setRowCount(self.protein.x_count)
        row = 0
        for peptide_id in self.protein.peptides.keys():
            temp_pep = self.protein.peptides[peptide_id]
            if re.match('.*X.*', temp_pep.sequence):
                pep = MyTableWidgetItem("%s" % temp_pep.sequence)
                pep.id = peptide_id
                self.peptides_list.setItem(row, 0, pep)
                self.peptides_list.setItem(row, 1, QtWidgets.QTableWidgetItem("+%d" % temp_pep.charge_state))
                self.peptides_list.setItem(row, 2, QtWidgets.QTableWidgetItem("%.3f" % temp_pep.mass))
                row += 1
        self.peptides_list.resizeColumnsToContents()
        if row > 0:
            self.peptides_list.selectRow(0)
        
    def protease_selected(self):
        self.cleave_protein()
        

    def peptide_selected(self):
        self.update_fragment_table()

    def mass_spectra_selected(self):
        try:
            row = self.mass_list.selectedItems()[0]
            _thread.start_new_thread(self.histogram.plot_figure, (self.spectra.reads[row.id],))
            self.update_binding_list()
        except IndexError:
            pass
    
    def update_consolidated_list(self):
        consolidated = {}
        for mass_id in self.spectra.reads.keys():
            if self.spectra.reads[mass_id].score >= self.pass_window.value():
                read = self.spectra.reads[mass_id]
                for row, peptide_id in enumerate(read.peptide_matches):
                    peptide = self.protein.peptides[peptide_id]
                    try:
                        consolidated[peptide.sequence] += read.peptide_scores[row]
                    except KeyError:
                        consolidated[peptide.sequence] = read.peptide_scores[row]
        
        self.consolidated_list.clearContents()
        self.consolidated_list.setRowCount(len(consolidated.keys()))
        i = 0
        for key in consolidated.keys():
            normalized_score = float(consolidated[key])/len(key)
            self.consolidated_list.setItem(i, 0, QtWidgets.QTableWidgetItem(key))
            self.consolidated_list.setItem(i, 1, SortTableWidgetItem("%.2f" % normalized_score))
            i += 1
        self.consolidated_list.resizeColumnsToContents()
        self.consolidated_list.sortItems(1, order=1)
        
        
    def update_binding_list(self):
        mass_id = self.mass_list.selectedItems()[0].get_id()
        self.binding_list.clearContents()
        self.binding_list.setRowCount(0)
        rows = len(self.spectra.reads[mass_id].peptide_matches)
        if rows != 0:
            self.binding_list.setRowCount(rows)
            i = 0
            for peptide_id in self.spectra.reads[mass_id].peptide_matches:
                self.binding_list.setItem(i, 0, QtWidgets.QTableWidgetItem(self.protein.peptides[peptide_id].sequence))
                self.binding_list.setItem(i, 1, SortTableWidgetItem("%d" % self.spectra.reads[mass_id].peptide_scores[i]))
                i += 1
            self.binding_list.resizeColumnsToContents()
            self.binding_list.sortItems(1, order=1)

    def update_mass_list(self):
        self.mass_list.clearContents()
        self.mass_list.setRowCount(self.spectra.hits)
        i = 0
        for id in self.spectra.reads.keys():
            read = self.spectra.reads[id]
            if read.score >= self.pass_window.value():
                temp_item = MyTableWidgetItem("%f" % read.charge_adjusted_mass)
                temp_item.id = id
                self.mass_list.setItem(i, 0, temp_item)
                self.mass_list.setItem(i, 1, QtWidgets.QTableWidgetItem("+%d" % read.charge_state))
                self.mass_list.setItem(i, 2, QtWidgets.QTableWidgetItem("%f" % read.retention_time))
                score_item = SortTableWidgetItem("%d" % read.score)
                if read.score >= self.pass_window.value():
                    score_item.setBackground(QtGui.QColor(110, 235, 131))
                self.mass_list.setItem(i, 3, score_item)
                i += 1
        self.mass_list.resizeColumnsToContents()
        self.mass_list.sortItems(3, order=1)
        if i > 0:
            self.mass_list.selectRow(0)

    def setup_initial_toolbar(self):
        self.toolBar.setContentsMargins(0, 5, 5, 0)
        self.toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)

        #Import peptide sequence action
        importProteinAction = QtWidgets.QAction(QtGui.QIcon('icons/48/Files/add_file-48.png'),
                                       'Import Protein',
                                       self)

        importProteinAction.setShortcut('Ctrl+O')
        importProteinAction.triggered.connect(self.import_peptide)
        self.toolBar.addAction(importProteinAction)

        #Import mass spec data action
        importSpectraAction = QtWidgets.QAction(QtGui.QIcon('icons/48/Files/add_file-48.png'),
                                     'Import Spectra',
                                     self)

        importSpectraAction.setShortcut('Ctrl+I')
        importSpectraAction.triggered.connect(self.import_mass_spectra)
        self.toolBar.addAction(importSpectraAction)

        saveAction = QtWidgets.QAction(QtGui.QIcon('icons/48/User_Interface/save-48.png'),
                                   'Save',
                                   self)

        saveAction.setShortcut('Ctrl+S')
        self.toolBar.addAction(saveAction)

        right_spacer = QtWidgets.QWidget()
        right_spacer.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.toolBar.addWidget(right_spacer)
        settingsAction = QtWidgets.QAction(QtGui.QIcon('icons/48/User_Interface/vertical_settings_mixer-48.png'),
                                       'Settings',
                                   self)
        settingsAction.setShortcut('Ctrl+X')
        settingsAction.triggered.connect(self.open_settings)
        self.toolBar.addAction(settingsAction)

    def score_matches(self, peptide_id):
        def _get_matches(frags, read):
            score = 0.0
            for frag in frags:
                for mi in read.mass_int:
                    read_mass = mi[0]
                    ppm_difference = abs((read_mass-frag.mass)*1000000/frag.mass)
                    if ppm_difference <= self.ppm_cutoff.value():
                        score += 1
            read.score = score
        
        fragments = [x for x in self.protein.peptides[peptide_id].fragments
                     if x.ion_type in self.get_selected_ion_types()]
        
        for mass_id in self.spectra.reads.keys():
            _get_matches(fragments, self.spectra.reads[mass_id])
        
        self.spectra.hits = 0
        for mass_id in self.spectra.reads.keys():
            if self.spectra.reads[mass_id].score >= self.pass_window.value():
                self.spectra.hits += 1
                self.score_peptides_for_mass_id(mass_id)
        
        self.update_mass_list()
        self.update_binding_list()
        self.update_consolidated_list()
        self.show_ui_elements(True)

    def score_peptides_for_mass_id(self, mass_id):
        def _get_matches(frags, read, pep_id):
            score = 0.0
            for frag in frags:
                for mi in read.mass_int:
                    read_mass = mi[0]
                    ppm_difference = abs((read_mass - frag.mass) * 1000000 / frag.mass)
                    if ppm_difference <= self.ppm_cutoff.value():
                        score += 1
            read.peptide_scores.append(score)
            read.peptide_matches.append(pep_id)

        self.spectra.reads[mass_id].peptide_scores = []
        self.spectra.reads[mass_id].peptide_matches = []
        for peptide_id in self.protein.peptides.keys():
            fragments = [x for x in self.protein.peptides[peptide_id].fragments
                         if x.ion_type in self.get_selected_ion_types()]
            _get_matches(fragments, self.spectra.reads[mass_id], peptide_id)

    @QtCore.pyqtSlot()
    def on_search_reads_clicked(self):
        self.statusbar.showMessage("Searching all peptides against all MS2 read")
        self.show_ui_elements(False)
        peptide_id = self.peptides_list.selectedItems()[0].get_id()
        _thread.start_new_thread(self.score_matches, (peptide_id,))
        
    def show_ui_elements(self, state=True):
        for widget in self.centralwidget.children():
            try:
                widget.setEnabled(state)
            except AttributeError:
                pass
        for widget in self.bottom_dock.children():
            try:
                widget.setEnabled(state)
            except AttributeError:
                pass
                
    
    def update_fragment_table(self):
        try:
            peptide_id = self.peptides_list.selectedItems()[0].get_id()
            self.fragments_list.setRowCount(len(self.protein.peptides[peptide_id].fragments))
            self.fragments_list.setColumnCount(4)
            self.fragments_list.clearContents()
            i = 0
            for fragment in self.protein.peptides[peptide_id].fragments:
                self.fragments_list.setItem(i, 0, QtWidgets.QTableWidgetItem(fragment.sequence))
                self.fragments_list.setItem(i, 1, QtWidgets.QTableWidgetItem(fragment.ion_type))
                self.fragments_list.setItem(i, 2, QtWidgets.QTableWidgetItem("+%d" % fragment.charge))
                self.fragments_list.setItem(i, 3, QtWidgets.QTableWidgetItem(str(fragment.mass)))
                i += 1
            self.fragments_list.resizeColumnsToContents()
        except IndexError:
            print("Bypassed")
            pass
    
    def fragment_peptides(self):
        self.protein.fragment_all_peptides(types=tuple(self.get_selected_ion_types()),
                                           maxcharge=self.get_max_fragment_charge(),
                                           crosslinker='bpa' if self.bpa.isChecked() else 'bpa_alk')
        self.update_fragment_table()

    def import_peptide(self, file_name=None):
        if not file_name:
            dlg = QtWidgets.QFileDialog()
            dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
            dlg.setWindowTitle("Import Peptide Sequence")
            dlg.setNameFilters(["Peptide Sequence Files (*.fasta)"])
            if dlg.exec_():
                self.protein = sequence.protein(dlg.selectedFiles()[0])
        else:
            self.protein = sequence.protein("/Users/Venky/Work/Softwares/Interlink/test_files/HDH_wt.fasta")
        self.protein.parse_fasta()
        self.protein_sequence_title.setText(self.protein.name)
        self.protein_sequence.setText(self.protein.sequence)
        self.fragment_peptides()

    def import_mass_spectra(self, file_name=None):
        if not file_name:
            dlg = QtWidgets.QFileDialog()
            dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
            dlg.setWindowTitle("Import Mass Spectra")
            dlg.setNameFilters(["Mass Spectra Files (*.ms2)"])
            if dlg.exec_():
                self.spectra = spectra.spectra(dlg.selectedFiles()[0])
        else:
            self.spectra = spectra.spectra("/Users/Venky/Work/Softwares/Interlink/test_files/CM_151128_Elite_32_N94HDH_gel_band_UV_HCD_rerun.ms2")
        self.spectra.parse_spectra()

    def open_settings(self):
        print("settings")
        
    def get_selected_ion_types(self):
        ion_types = []
        for checkbox in self.ion_type_group.findChildren(QtWidgets.QCheckBox):
            if checkbox.isChecked():
                ion_types.append(str(checkbox.text()))
        return ion_types
    
    def get_max_fragment_charge(self):
        charge = 1
        for radiobutton in self.charge_state_group.findChildren(QtWidgets.QRadioButton):
            if radiobutton.isChecked():
                txt = str(radiobutton.text())
                charge_string = txt.split(',')
                charge = int(charge_string[0].replace('+', ''))
        return charge
    
    def load_test_files(self):
        self.import_peptide(file_name="default")
        self.import_mass_spectra(file_name="default")

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    form = InterLink()
    form.show()
    app.exec_()
    app.deleteLater()
    sys.exit()
