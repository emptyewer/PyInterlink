import re
from libs.bintrees import AVLTree

class ms_read():
    def __init__(self):
        self.scan_number = 0
        self.observed_mass = 0.0
        self.retention_time = 0.0
        self.intensity_ms1 = 0.0
        self.charge_state = 0
        self.charge_adjusted_mass = 0.0
        self.mass_int = []
        self.score = 0.0
        self.peptide_matches = []
        self.peptide_scores = []
        self.no_uv_peptide_matches = []
        self.no_uv_peptide_scores = []

    def __repr__(self):
        return "Intensity: %f, Score: %f" % (self.intensity_ms1, self.score)

class spectra():
    def __init__(self, parent):
        self.parent = parent
        self.uv_filename = None
        self.no_uv_filename = None
        self.uv_reads = {}
        self.no_uv_reads = {}
        self.no_uv_ids = []
        self.no_uv_tree = None
        self.hits = 0
        self.parse_options = {
            0: self.s_matched,
            1: self.r_matched,
            2: self.i_matched,
            3: self.z_matched,
            4: self.m_matched,
        }
        
    def s_matched(self, read, line):
        split = line.split()
        read.scan_number = int(split[1])
        read.observed_mass = float(split[3])

    def r_matched(self, read, line):
        split = line.split()
        read.retention_time = float(split[2])

    def i_matched(self, read, line):
        split = line.split()
        read.intensity_ms1 = float(split[2])

    def z_matched(self, read, line):
        split = line.split()
        read.charge_state = int(split[1])
        read.charge_adjusted_mass = float(split[2])

    def m_matched(self, read, line):
        split = tuple(line.split())
        read.mass_int.append((float(split[0]), float(split[1])))

    def parse_spectra(self, tag):
        s_match = re.compile(r'^S\s+')
        r_match = re.compile(r'^I\s+RTime')
        i_match = re.compile(r'^I\s+ms1int')
        z_match = re.compile(r'^Z\s+')
        m_match = re.compile(r'^\d+.\d+\s+\d+')
        read = None
        handle = None
        if tag == 0:
            handle = open(self.uv_filename)
        elif tag == 1:
            handle = open(self.no_uv_filename)
        id = 0
        kv_pairs = []
        for line in handle.readlines():
            line = line.rstrip()
            if s_match.match(line):
                if read != None:
                    if tag == 0:
                        self.uv_reads[id] = read
                    elif tag == 1:
                        self.no_uv_reads[id] = read
                    id += 1
                read = ms_read()
                self.parse_options[0](read, line)
            elif r_match.match(line):
                self.parse_options[1](read, line)
            elif i_match.match(line):
                self.parse_options[2](read, line)
            elif z_match.match(line):
                self.parse_options[3](read, line)
                kv_pairs.append((read.charge_adjusted_mass, id))
            elif m_match.match(line):
                self.parse_options[4](read, line)
        if tag == 1:
            self.no_uv_tree = AVLTree(kv_pairs)
        self.parent.show_ui_elements()

