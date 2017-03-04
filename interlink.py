# Python 2.7

import os
import re
import sys
import thread
from PyQt4 import QtCore, QtGui, uic

import structs.sequence as sequence
import structs.spectra as spectra

import structs.mass_histogram as histogram

# ui_path = os.path.join('ui', 'mac', 'main.ui')
# if sys.platform == 'win32':
#     ui_path = os.path.join('ui', 'win', 'main.ui')

ui_path = os.path.join('ui', 'main.ui')
main_class, _ = uic.loadUiType(ui_path)

ui_path = os.path.join('ui', 'workflow.ui')
workflow_class, _ = uic.loadUiType(ui_path)

ui_path = os.path.join('ui', 'spectrum.ui')
spectrum_class, _ = uic.loadUiType(ui_path)


class MyTableWidgetItem(QtGui.QTableWidgetItem):
    def __init__(self, *args):
        super(MyTableWidgetItem, self).__init__(*args)
        self.id = ''

    def get_id(self):
        return self.id

class SortTableWidgetItem(QtGui.QTableWidgetItem):
    def __init__(self, *args):
        super(SortTableWidgetItem, self).__init__(*args)

    def __lt__(self, other):
        if isinstance(other, SortTableWidgetItem):
            self_data_value = float(str(self.data(QtCore.Qt.EditRole).toPyObject()))
            other_data_value = float(str(other.data(QtCore.Qt.EditRole).toPyObject()))
            return self_data_value < other_data_value
        else:
            return QtGui.QTableWidgetItem.__lt__(self, other)
                
class Workflow(QtGui.QWidget, workflow_class):
    def __init__(self, parent, *args):
        super(Workflow, self).__init__(*args)
        self.setupUi(self)
        self.parent = parent

class Spectrum(QtGui.QWidget, spectrum_class):
    def __init__(self, parent, *args):
        super(Spectrum, self).__init__(*args)
        self.setupUi(self)
        self.parent = parent

class InterLink(QtGui.QMainWindow, main_class):
    def __init__(self, *args):
        super(InterLink, self).__init__(*args)
        self.setupUi(self)
        self.spectrum_window = Spectrum(self)
        self.workflow_window = Workflow(self)
        self.protein = None
        self.spectra = None
        self.window = self.window()
        self.window.setGeometry(5, 60, self.width(), self.height())
        self.setup_initial_toolbar()
        self.setup_protease_list()
        self.setup_signals()
        self.search_masses = False
        self.histogram = histogram.View(self, width=5, height=4, dpi=100)
        self.window.setGeometry(10, 30, self.width(), self.height())
        self.spectrum_window.histogram_layout_2.addWidget(self.histogram)
        self.spectrum_window.setWindowFlags(QtCore.Qt.CustomizeWindowHint |
                                            QtCore.Qt.WindowTitleHint |
                                            QtCore.Qt.WindowMinMaxButtonsHint)
        self.spectrum_window.setGeometry(10, self.height() - self.spectrum_window.height(),
                                         self.spectrum_window.width(), self.spectrum_window.height())
        self.spectrum_window.show()
        self.workflow_window.setWindowFlags(QtCore.Qt.CustomizeWindowHint |
                                            QtCore.Qt.WindowTitleHint |
                                            QtCore.Qt.WindowMinMaxButtonsHint)
        self.workflow_window.setGeometry(10, self.height() + 50,
                                         self.workflow_window.width(), self.workflow_window.height())
        self.workflow_window.show()
        self.load_test_files()

    def setup_signals(self):
        self.protein_sequence.textChanged.connect(self.protein_sequence_edited)

        self.protease_list.currentIndexChanged.connect(self.cleave_protein)
        self.max_missed_cleavage.valueChanged.connect(self.cleave_protein)
        self.workflow_window.mass_list.itemSelectionChanged.connect(self.mass_spectra_selected)
        self.peptides_list.itemSelectionChanged.connect(self.peptide_selected)
        self.workflow_window.binding_list.itemSelectionChanged.connect(self.binding_peptide_selected)
        self.workflow_window.search_reads.clicked.connect(self.on_search_reads_clicked)
        self.connect(self, QtCore.SIGNAL("FINISHED_LOADING"), self.show_ui_elements)

        for checkbox in self.ion_type_group.findChildren(QtGui.QCheckBox):
            checkbox.clicked.connect(self.fragment_peptides)
        for radiobutton in self.charge_state_group.findChildren(QtGui.QRadioButton):
            radiobutton.clicked.connect(self.fragment_peptides)
        for radiobutton in self.crosslinker_group.findChildren(QtGui.QRadioButton):
            radiobutton.clicked.connect(self.cleave_protein)

        self.actionImport_Protein.triggered.connect(self.import_peptide)
        self.actionImport_UV_Spectra_1.triggered.connect(self.import_mass_spectra_uv)
        self.actionImport_UV_Spectra_2.triggered.connect(self.import_mass_spectra_no_uv)

    def protein_sequence_edited(self):
        self.peptides_list.clearContents()
        self.workflow_window.fragments_list.clearContents()
        self.workflow_window.fragments_list.setRowCount(0)
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
        self.workflow_window.fragments_list.clearContents()
        self.peptides_list.setRowCount(self.protein.x_count)
        row = 0
        for peptide_id in self.protein.peptides.keys():
            temp_pep = self.protein.peptides[peptide_id]
            if re.match('.*X.*', temp_pep.sequence):
                pep = MyTableWidgetItem("%s" % temp_pep.sequence)
                pep.id = peptide_id
                self.peptides_list.setItem(row, 0, pep)
                self.peptides_list.setItem(row, 1, QtGui.QTableWidgetItem("+%d" % temp_pep.charge_state))
                self.peptides_list.setItem(row, 2, QtGui.QTableWidgetItem("%.3f" % temp_pep.mass))
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
            row = self.workflow_window.mass_list.selectedItems()[0]
            read = self.spectra.uv_reads[row.id]
            thread.start_new_thread(self.histogram.plot_figure, (read, 'b',))
            # no_uv_read1 = self.spectra.no_uv_reads[self.spectra.no_uv_tree.floor_item(read.charge_adjusted_mass)[1]]
            # no_uv_read2 = self.spectra.no_uv_reads[self.spectra.no_uv_tree.ceiling_item(read.charge_adjusted_mass)[1]]
            # no_uv_read = no_uv_read1
            # if (read.charge_adjusted_mass - no_uv_read1.charge_adjusted_mass) > \
            #         (read.charge_adjusted_mass - no_uv_read2.charge_adjusted_mass):
            #     no_uv_read = no_uv_read2
            # thread.start_new_thread(self.histogram.plot_figure, (read, 'b',))
            self.update_binding_list()
        except IndexError:
            pass

    def binding_peptide_selected(self):
        peptide_id = self.workflow_window.binding_list.selectedItems()[0].get_id()
        thread.start_new_thread(self.histogram.plot_fragment, (self.protein.peptides[peptide_id].fragments,))
    
    # def update_consolidated_list(self):
    #     consolidated = {}
    #     for mass_id in self.spectra.reads.keys():
    #         if self.spectra.reads[mass_id].score >= self.pass_window.value():
    #             read = self.spectra.reads[mass_id]
    #             for row, peptide_id in enumerate(read.peptide_matches):
    #                 peptide = self.protein.peptides[peptide_id]
    #                 try:
    #                     consolidated[peptide.sequence] += read.peptide_scores[row]
    #                 except KeyError:
    #                     consolidated[peptide.sequence] = read.peptide_scores[row]
    #
    #     self.consolidated_list.clearContents()
    #     self.consolidated_list.setRowCount(len(consolidated.keys()))
    #     i = 0
    #     for key in consolidated.keys():
    #         normalized_score = float(consolidated[key])/len(key)
    #         self.consolidated_list.setItem(i, 0, QtGui.QTableWidgetItem(key))
    #         self.consolidated_list.setItem(i, 1, SortTableWidgetItem("%.2f" % normalized_score))
    #         i += 1
    #     self.consolidated_list.resizeColumnsToContents()
    #     self.consolidated_list.sortItems(1, order=1)
        
        
    def update_binding_list(self):
        mass_id = self.workflow_window.mass_list.selectedItems()[0].get_id()
        self.workflow_window.binding_list.clearContents()
        self.workflow_window.binding_list.setRowCount(0)
        rows = len(self.spectra.uv_reads[mass_id].peptide_matches)
        if rows != 0:
            self.workflow_window.binding_list.setRowCount(rows)
            i = 0
            for peptide_id in self.spectra.uv_reads[mass_id].peptide_matches:
                pep = MyTableWidgetItem("%s" % self.protein.peptides[peptide_id].sequence)
                pep.id = peptide_id
                self.workflow_window.binding_list.setItem(i, 0, pep)
                score = self.spectra.uv_reads[mass_id].peptide_scores[i]
                self.workflow_window.binding_list.setItem(i, 1, SortTableWidgetItem("%d" % score[0]))
                self.workflow_window.binding_list.setItem(i, 2, SortTableWidgetItem("%d" % score[1]))
                i += 1
            self.workflow_window.binding_list.resizeColumnsToContents()
            self.workflow_window.binding_list.sortItems(1, order=1)

    def update_mass_list(self):
        self.workflow_window.mass_list.clearContents()
        self.workflow_window.mass_list.setRowCount(self.spectra.hits)
        i = 0
        for id in self.spectra.uv_reads.keys():
            read = self.spectra.uv_reads[id]
            if read.score >= self.workflow_window.pass_window.value() and read.charge_state >= 3:
                temp_item = MyTableWidgetItem("%f" % read.charge_adjusted_mass)
                temp_item.id = id
                self.workflow_window.mass_list.setItem(i, 0, temp_item)
                self.workflow_window.mass_list.setItem(i, 1, QtGui.QTableWidgetItem("+%d" % read.charge_state))
                self.workflow_window.mass_list.setItem(i, 2, QtGui.QTableWidgetItem("%f" % read.retention_time))
                score_item = SortTableWidgetItem("%d" % read.score)
                if read.score >= self.workflow_window.pass_window.value():
                    score_item.setBackground(QtGui.QColor(110, 235, 131))
                self.workflow_window.mass_list.setItem(i, 3, score_item)
                i += 1
        self.workflow_window.mass_list.resizeColumnsToContents()
        self.workflow_window.mass_list.sortItems(3, order=1)
        if i > 0:
            self.workflow_window.mass_list.selectRow(0)

    def setup_initial_toolbar(self):
        self.toolBar.setContentsMargins(0, 5, 5, 0)
        self.toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)

        #Import peptide sequence action
        importProteinAction = QtGui.QAction(QtGui.QIcon('icons/48/Files/add_file-48.png'),
                                       'Import Protein',
                                       self)

        importProteinAction.setShortcut('Ctrl+O')
        importProteinAction.triggered.connect(self.import_peptide)
        self.toolBar.addAction(importProteinAction)

        #Import mass spec data action
        importSpectraAction = QtGui.QAction(QtGui.QIcon('icons/48/Files/add_file-48.png'),
                                     'UV+ Spectra',
                                     self)

        importSpectraAction.setShortcut('Ctrl+I')
        importSpectraAction.triggered.connect(self.import_mass_spectra_uv)
        self.toolBar.addAction(importSpectraAction)

        # Import mass spec data action
        importSpectraAction2 = QtGui.QAction(QtGui.QIcon('icons/48/Files/add_file-48.png'),
                                            'UV- Spectra',
                                            self)

        importSpectraAction2.setShortcut('Ctrl+Y')
        importSpectraAction2.triggered.connect(self.import_mass_spectra_no_uv)
        self.toolBar.addAction(importSpectraAction2)

        right_spacer = QtGui.QWidget()
        right_spacer.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.toolBar.addWidget(right_spacer)
        saveAction = QtGui.QAction(QtGui.QIcon('icons/48/User_Interface/save-48.png'),
                                   'Save',
                                   self)

        saveAction.setShortcut('Ctrl+S')
        self.toolBar.addAction(saveAction)

        # right_spacer = QtGui.QWidget()
        # right_spacer.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        # self.toolBar.addWidget(right_spacer)
        # settingsAction = QtGui.QAction(QtGui.QIcon('icons/48/User_Interface/vertical_settings_mixer-48.png'),
        #                                'Settings',
        #                            self)
        # settingsAction.setShortcut('Ctrl+X')
        # settingsAction.triggered.connect(self.open_settings)
        # self.toolBar.addAction(settingsAction)

    def score_matches(self, peptide_id):
        def _get_crosslinked_fragment_matches(frags):
            for frag in frags:
                for key in self.spectra.uv_reads.keys():
                    for mass_int in self.spectra.uv_reads[key].mass_int:
                        ppm_difference = abs(float(mass_int[0] - frag.mass) * 1000000.0 / frag.mass)
                        if ppm_difference <= self.ppm_cutoff.value():
                            self.spectra.uv_reads[key].score += 1
        
        fragments = [x for x in self.protein.peptides[peptide_id].fragments
                     if x.ion_type in self.get_selected_ion_types()]
        
        # for mass_id in self.spectra.uv_reads.keys():
        _get_crosslinked_fragment_matches(fragments)

        self.spectra.hits = 0
        for mass_id in self.spectra.uv_reads.keys():
            if self.spectra.uv_reads[mass_id].score >= self.workflow_window.pass_window.value() and \
                    self.spectra.uv_reads[mass_id].charge_state >= 3:
                self.spectra.hits += 1
                self.score_peptides_for_mass_id(self.spectra.uv_reads, mass_id)
        self.update_mass_list()
        self.update_binding_list()
        # self.update_consolidated_list()
        self.show_ui_elements(True)

    def score_peptides_for_mass_id(self, uv_reads, mass_id):
        def _get_matches(frags, read, pep_id):
            uv_score = 0.0
            no_uv_score = 0.0
            for frag in frags:
                for mi in read.mass_int:
                    read_mass = mi[0]
                    ppm_difference = abs((read_mass - frag.mass) * 1000000 / frag.mass)
                    if ppm_difference <= self.ppm_cutoff.value():
                        uv_score += 1
                        no_uv_read1 = self.spectra.no_uv_reads[self.spectra.no_uv_tree.floor_item(read.charge_adjusted_mass)[1]]
                        no_uv_read2 = self.spectra.no_uv_reads[self.spectra.no_uv_tree.ceiling_item(read.charge_adjusted_mass)[1]]
                        no_uv_read = no_uv_read1
                        if (read.charge_adjusted_mass - no_uv_read1.charge_adjusted_mass) > \
                                (read.charge_adjusted_mass - no_uv_read2.charge_adjusted_mass):
                            no_uv_read = no_uv_read2
                        for nmi in no_uv_read.mass_int:
                            nread_mass = nmi[0]
                            nppm_diff = abs((nread_mass - frag.mass) * 1000000 / frag.mass)
                            if nppm_diff <= self.ppm_cutoff.value():
                                no_uv_score += 1
            read.peptide_scores.append((uv_score, no_uv_score))
            read.peptide_matches.append(pep_id)

        uv_reads[mass_id].peptide_scores = []
        uv_reads[mass_id].peptide_matches = []
        for peptide_id in self.protein.peptides.keys():
            fragments = [x for x in self.protein.peptides[peptide_id].fragments
                         if x.ion_type in self.get_selected_ion_types()]
            _get_matches(fragments, uv_reads[mass_id], peptide_id)

    def on_search_reads_clicked(self):
        self.statusbar.showMessage("Searching all peptides against all MS2 read")
        self.show_ui_elements(False)
        peptide_id = self.peptides_list.selectedItems()[0].get_id()
        thread.start_new_thread(self.score_matches, (peptide_id,))
        
    def show_ui_elements(self, state=True):
        for widget in self.centralwidget.children():
            try:
                widget.setEnabled(state)
            except AttributeError:
                pass
    
    def update_fragment_table(self):
        try:
            peptide_id = self.peptides_list.selectedItems()[0].get_id()
            self.workflow_window.fragments_list.setRowCount(len(self.protein.peptides[peptide_id].fragments))
            self.workflow_window.fragments_list.setColumnCount(4)
            self.workflow_window.fragments_list.clearContents()
            i = 0
            for fragment in self.protein.peptides[peptide_id].fragments:
                self.workflow_window.fragments_list.setItem(i, 0, QtGui.QTableWidgetItem(fragment.sequence))
                self.workflow_window.fragments_list.setItem(i, 1, QtGui.QTableWidgetItem(fragment.ion_type))
                self.workflow_window.fragments_list.setItem(i, 2, QtGui.QTableWidgetItem("+%d" % fragment.charge))
                self.workflow_window.fragments_list.setItem(i, 3, QtGui.QTableWidgetItem(str(fragment.mass)))
                i += 1
            self.workflow_window.fragments_list.resizeColumnsToContents()
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
            dlg = QtGui.QFileDialog()
            dlg.setFileMode(QtGui.QFileDialog.AnyFile)
            dlg.setWindowTitle("Import Peptide Sequence")
            dlg.setNameFilters(["Peptide Sequence Files (*.fasta)"])
            if dlg.exec_():
                self.protein = sequence.protein(dlg.selectedFiles()[0])
        else:
            self.protein = sequence.protein("./test_files/HDH_wt.fasta")
        self.protein.parse_fasta()
        self.protein_sequence.setText(self.protein.sequence)
        self.fragment_peptides()

    def import_mass_spectra_uv(self, file_name=None):
        if not file_name:
            dlg = QtGui.QFileDialog()
            dlg.setFileMode(QtGui.QFileDialog.AnyFile)
            dlg.setWindowTitle("UV+ Mass Spectra")
            dlg.setNameFilters(["Mass Spectra Files (*.ms2)"])
            if dlg.exec_():
                if not self.spectra:
                    self.spectra = spectra.spectra(self)
                self.spectra.uv_filename = dlg.selectedFiles()[0]
        else:
            self.spectra = spectra.spectra(self)
            self.spectra.uv_filename = "./test_files/CM_151128_Elite_32_N94HDH_gel_band_UV_CID.ms2"
        self.show_ui_elements(False)
        thread.start_new_thread(self.spectra.parse_spectra, (0,))

    def import_mass_spectra_no_uv(self, file_name=None):
        if not file_name:
            dlg = QtGui.QFileDialog()
            dlg.setFileMode(QtGui.QFileDialog.AnyFile)
            dlg.setWindowTitle("UV- Mass Spectra")
            dlg.setNameFilters(["Mass Spectra Files (*.ms2)"])
            if dlg.exec_():
                if not self.spectra:
                    self.spectra = spectra.spectra(self)
                self.spectra.no_uv_filename = dlg.selectedFiles()[0]
        else:
            self.spectra = spectra.spectra(self)
            self.spectra.no_uv_filename = "./test_files/CM_151128_Elite_32_N94HDH_gel_band_UV_CID.ms2"
        self.show_ui_elements(False)
        thread.start_new_thread(self.spectra.parse_spectra, (1,))

    def open_settings(self):
        print("settings")
        
    def get_selected_ion_types(self):
        ion_types = []
        for checkbox in self.ion_type_group.findChildren(QtGui.QCheckBox):
            if checkbox.isChecked():
                ion_types.append(str(checkbox.text()))
        return ion_types
    
    def get_max_fragment_charge(self):
        charge = 1
        for radiobutton in self.charge_state_group.findChildren(QtGui.QRadioButton):
            if radiobutton.isChecked():
                txt = str(radiobutton.text())
                charge_string = txt.split(',')
                charge = int(charge_string[0].replace('+', ''))
        return charge
    
    def load_test_files(self):
        self.import_peptide(file_name="default")
        # self.import_mass_spectra_uv(file_name="default")
        # self.import_mass_spectra_no_uv(file_name="default")

    def closeEvent(self, event):
        event.ignore()
        self.spectrum_window.close()
        self.workflow_window.close()
        super(InterLink, self).closeEvent(event)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    form = InterLink()
    form.show()
    app.exec_()
    app.deleteLater()
    sys.exit()
