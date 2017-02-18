import re

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

    def __repr__(self):
        return "Intensity: %f, Score: %f" % (self.intensity_ms1, self.score)

class spectra():
    def __init__(self, filename=None):
        self.filename = filename
        self.reads = {}
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

    def parse_spectra(self):
        s_match = re.compile(r'^S\s+')
        r_match = re.compile(r'^I\s+RTime')
        i_match = re.compile(r'^I\s+ms1int')
        z_match = re.compile(r'^Z\s+')
        m_match = re.compile(r'^\d+.\d+\s+\d+')
        read = None
        handle = open(self.filename)
        id = 0
        for line in handle.readlines():
            line = line.rstrip()
            if s_match.match(line):
                if read != None:
                    self.reads[id] = read
                    id += 1
                read = ms_read()
                self.parse_options[0](read, line)
            elif r_match.match(line):
                self.parse_options[1](read, line)
            elif i_match.match(line):
                self.parse_options[2](read, line)
            elif z_match.match(line):
                self.parse_options[3](read, line)
            elif m_match.match(line):
                self.parse_options[4](read, line)


