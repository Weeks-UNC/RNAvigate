import Bio.PDB
from .data import Data
import numpy as np


class PDB(Data):

    def __init__(self, filepath, chain=None, datatype="pdb", offset=None, fasta=None):
        self.datatype = datatype
        get_offset = offset is None
        get_seq = fasta is None
        self.sequence = ""
        if chain is None:
            self.chain = " "
        else:
            self.chain = chain
        if get_offset or get_seq:
            self.get_sequence_offset(filepath, get_seq, get_offset)
        if not get_offset:
            self.offset = offset
        if not get_seq:
            super().__init__(filepath=fasta)
        self.read_pdb(filepath)
        self.path = filepath
        self.distance_matrix = {}

    def get_sequence_offset(self, pdb, get_seq=True, get_offset=True):
        with open(pdb) as file:
            for line in file.readlines():
                line = [field.strip() for field in line.split()]
                if line[0] == "DBREF" and line[2] == self.chain and get_offset:
                    self.offset = int(line[3])-1
                if line[0] == "SEQRES" and line[2] == self.chain and get_seq:
                    self.sequence += ''.join(line[4:])

    def read_pdb(self, pdb):
        parser = Bio.PDB.PDBParser()
        self.pdb = parser.get_structure('RNA', pdb)
        self.validres = []
        for res in self.pdb[0][self.chain].get_residues():
            res_id = res.get_id()
            if res_id[0] == " ":
                self.validres.append(res_id[1]-self.offset)
        self.validres = np.array(self.validres)

    def get_xyz_coord(self, nt, atom):
        if atom == "DMS":
            if self.sequence.upper()[nt-1] in "AG":
                atom = "N1"
            elif self.sequence.upper()[nt-1] in "UC":
                atom = "N3"
        xyz = [float(c) for c in self.pdb[0][self.chain]
               [int(nt)][atom].get_coord()]
        return xyz

    def get_distance(self, i, j, atom="O2'"):
        if atom in self.distance_matrix.keys():
            return self.distance_matrix[atom][i-1, j-1]
        valid = [v - self.offset for v in self.validres]
        if i in valid and j in valid:
            xi, yi, zi = self.get_xyz_coord(i, atom)
            xj, yj, zj = self.get_xyz_coord(j, atom)
            distance = ((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)**0.5
        else:
            distance = 1000
        return distance

    def get_distance_matrix(self, atom="O2'"):
        if atom in self.distance_matrix.keys():
            return self.distance_matrix[atom]
        x = np.full(self.length, np.nan)
        y = np.full(self.length, np.nan)
        z = np.full(self.length, np.nan)
        for i in self.validres:
            x[i-1], y[i-1], z[i-1] = self.get_xyz_coord(i, atom)
        a = x - x[:, np.newaxis]
        b = y - y[:, np.newaxis]
        c = z - z[:, np.newaxis]
        matrix = np.sqrt(a*a + b*b + c*c)
        self.distance_matrix[atom] = np.nan_to_num(matrix, nan=1000)
        return matrix
