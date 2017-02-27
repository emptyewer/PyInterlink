import re
from structs.pyteomics import mass, electrochem
from structs.pyteomics.parser import cleave, expasy_rules


class peptide():
    def __init__(self):
        self.sequence = 0
        self.pI = 0.0
        self.mass = 0.0
        self.charge_state = 0.0
        self.fragments = []


class fragment():
    def __init__(self, seq, type, charge, mass):
        self.sequence = seq
        self.ion_type = type
        self.charge = charge
        self.mass = mass
        
    def __repr__(self):
        return "%s (%f)" % (self.sequence, self.mass)


class protein():
    def __init__(self, filename=None):
        self.filename = filename
        self.sequence = ''
        self.name = ''
        self.peptides = {}
        self.x_count = 0
        self.charge = 0
        self.pI = 0

    def parse_fasta(self):
        filehandle = open(self.filename, 'r')
        seq = ''
        for line in filehandle.readlines():
            if re.match(r'^>', line):
                self.name = line.rstrip().replace('>', '')
            else:
                seq += line.rstrip()
        self.sequence = seq

    def properties(self):
        self.charge = electrochem.charge(self.sequence, 7.0)
        self.pI = electrochem.pI(self.sequence)

    def cleave_sequence(self, enzyme, max_missed, crosslinker='bpa'):
        if crosslinker == 'bpa':
            mass.std_aa_comp['X'] = mass.Composition({'H': 13, 'C': 16, 'O': 2, 'N': 1})
            mass.std_aa_mass['X'] = 251.0946
        elif crosslinker == 'bpa_alk':
            mass.std_aa_comp['X'] = mass.Composition({'H': 23, 'C': 26, 'O': 3, 'N': 5})
            mass.std_aa_mass['X'] = 453.18009

        self.peptides = {}
        self.x_count = 0
        id = 0
        for pep_seq in list(cleave(self.sequence, expasy_rules[enzyme], missed_cleavages=max_missed, min_length=None)):
            if re.match(r'.*X.*', pep_seq):
                self.x_count += 1
            temp_pep = peptide()
            temp_pep.sequence = pep_seq
            temp_pep.charge_state = 1
            temp_pep.pI = electrochem.pI(pep_seq)
            temp_pep.mass = mass.calculate_mass(sequence=pep_seq, charge=1)
            self.peptides[id] = temp_pep
            id += 1
        
    def fragment_all_peptides(self, types=('b', 'y'), maxcharge=1, crosslinker='bpa'):
        for pep_id in self.peptides.keys():
            self._fragment_peptide(pep_id, types=types, maxcharge=maxcharge, crosslinker=crosslinker)

    def _fragment_peptide(self, pep_id, types=('b', 'y'), maxcharge=1, crosslinker='bpa'):
        """
        The function generates all possible m/z for fragments of types
        `types` and of charges from 1 to `maxcharge`.
        """
        peptide = self.peptides[pep_id].sequence
        if crosslinker == 'bpa':
            mass.std_aa_comp['X'] = mass.Composition({'H': 13, 'C': 16, 'O': 2, 'N': 1})
            mass.std_aa_mass['X'] = 251.0946
        elif crosslinker == 'bpa_alk':
            mass.std_aa_comp['X'] = mass.Composition({'H': 23, 'C': 26, 'O': 3, 'N': 5})
            mass.std_aa_mass['X'] = 453.18009
        self.peptides[pep_id].fragments = []
        for i in range(1, len(peptide) - 1):
            for ion_type in types:
                for charge in range(1, maxcharge + 1):
                    if ion_type[0] in 'abc':
                        frag = fragment(peptide[:i], ion_type, charge, mass.calculate_mass(sequence=peptide[:i],
                                                                                ion_type=ion_type,
                                                                                charge=charge))
                        self.peptides[pep_id].fragments.append(frag)
                    else:
                        frag = fragment(peptide[i:], ion_type, charge, mass.calculate_mass(sequence=peptide[i:],
                                                                                ion_type=ion_type,
                                                                                charge=charge))
                        self.peptides[pep_id].fragments.append(frag)

