import Bio.PDB
from rnavigate import data
import numpy as np


class PDB(data.Data):
    def __init__(self, filepath, chain, sequence=None):
        self.chain = chain
        if sequence is None:
            sequence = self.get_sequence(filepath)
        super().__init__(sequence=sequence)
        self.read_pdb(filepath)
        self.path = filepath
        self.distance_matrix = {}

    def get_sequence(self, pdb):
        seqres = []
        with open(pdb) as file:
            for line in file.readlines():
                line = [field.strip() for field in line.split()]
                if line[0] == "SEQRES" and line[2] == self.chain:
                    seqres.extend(line[4:])
        return self.get_sequence_from_seqres(seqres)

    def get_sequence_from_seqres(self, seqres):
        sequence = []
        for nt in seqres:
            valid = nt[0].upper() in "GUACT"
            if valid and (len(nt) == 3) and nt.endswith("TP"):
                sequence.append(nt[0])
            elif valid and (len(nt) == 1):
                sequence.append(nt[0])
            else:
                raise ValueError("invalid nt in SEQRES: " + nt)
        return ''.join(sequence)

    def read_pdb(self, pdb):
        if pdb.split('.')[-1] == "pdb":
            parser = Bio.PDB.PDBParser(QUIET=True)
        elif pdb.split('.')[-1] == 'cif':
            parser = Bio.PDB.MMCIFParser(QUIET=True)
        self.pdb = parser.get_structure('RNA', pdb)
        self.pdb_idx = []
        self.pdb_seq = []
        for res in self.pdb[0][self.chain].get_residues():
            res_id = res.get_id()
            res_seq = res.get_resname()
            if res_id[0] == " ":
                self.pdb_idx.append(res_id[1])
                self.pdb_seq.append(res_seq.strip())
        self.pdb_idx = np.array(self.pdb_idx)
        self.set_indices()

    def set_indices(self):
        for i in range(len(self.sequence)):
            correct_offset = True
            for pdb_nt, pdb_idx in zip(self.pdb_seq, self.pdb_idx):
                if self.sequence[pdb_idx-i-1] != pdb_nt:
                    correct_offset = False
                    break
            if correct_offset:
                self.offset = i
                break
        if not correct_offset:
            print('PDB entries could not be matched to sequence.')

    def get_pdb_idx(self, seq_idx):
        return seq_idx + self.offset

    def get_seq_idx(self, pdb_idx):
        return pdb_idx - self.offset

    def is_valid_idx(self, pdb_idx=None, seq_idx=None):
        if pdb_idx is not None and pdb_idx in self.pdb_idx:
            return True
        elif seq_idx is not None and (seq_idx in (self.pdb_idx - self.offset)):
            return True
        else:
            return False

    def get_xyz_coord(self, nt, atom):
        pdb_idx = self.get_pdb_idx(nt)
        if atom == "DMS":
            if self.sequence.upper()[nt-1] in "AG":
                atom = "N1"
            elif self.sequence.upper()[nt-1] in "UC":
                atom = "N3"
        xyz = [float(c) for c in self.pdb[0][self.chain]
               [int(pdb_idx)][atom].get_coord()]
        return xyz

    def get_distance(self, i, j, atom="O2'"):
        if atom in self.distance_matrix.keys():
            return self.distance_matrix[atom][i-1, j-1]
        try:
            xi, yi, zi = self.get_xyz_coord(i, atom)
            xj, yj, zj = self.get_xyz_coord(j, atom)
            distance = ((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)**0.5
        except KeyError:
            distance = np.nan
        return distance

    def get_distance_matrix(self, atom="O2'"):
        if atom in self.distance_matrix.keys():
            return self.distance_matrix[atom]
        x = np.full(self.length, np.nan)
        y = np.full(self.length, np.nan)
        z = np.full(self.length, np.nan)
        for i in (self.pdb_idx - self.offset):
            try:
                x[i-1], y[i-1], z[i-1] = self.get_xyz_coord(i, atom)
            except KeyError:
                pass
        a = x - x[:, np.newaxis]
        b = y - y[:, np.newaxis]
        c = z - z[:, np.newaxis]
        matrix = np.sqrt(a*a + b*b + c*c)
        self.distance_matrix[atom] = np.nan_to_num(matrix, nan=1000)
        return matrix
