import Bio.PDB
from data.data import Data


class PDB(Data):

    def __init__(self, filepath, datatype="pdb"):
        self.datatype = datatype
        self.read_pdb(filepath)

    def read_pdb(self, pdb):
        parser = Bio.PDB.PDBParser()
        self.pdb = parser.get_structure('RNA', pdb)
        self.sequence = ''
        with open(pdb) as file:
            for line in file.readlines():
                line = [field.strip() for field in line.split()]
                if line[0] == "SEQRES":
                    self.sequence += ''.join(line[4:])
        self.validres = []
        for res in self.pdb[0]["A"].get_residues():
            res_id = res.get_id()
            if res_id[0] == " ":
                self.validres.append(res_id[1])

    # def get_xyz_coord(self, nt, atom="O2'"):
    #     xyz = [float(c) for c in self.pdb[0]["A"][int(nt)][atom].get_coord()]
    #     return xyz

    # def get_3d_distance(self, i, j):
    #     valid = [nt.get_id()[1] for nt in self.pdb[0]["A"].get_residues()]
    #     if i in valid and j in valid:
    #         xi, yi, zi = self.get_xyz_coord(i)
    #         xj, yj, zj = self.get_xyz_coord(j)
    #         distance = ((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)**0.5
    #     else:
    #         distance = 1000
    #     return distance
