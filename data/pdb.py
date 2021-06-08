import Bio.PDB


class PDB():

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
        self.length = len(self.sequence)
        self.validres = []
        for res in self.pdb[0]["A"].get_residues():
            res_id = res.get_id()
            if res_id[0] == " ":
                self.validres.append(res_id[1])
