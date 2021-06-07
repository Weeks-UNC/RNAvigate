import pandas as pd


class Log():
    def __init__(self, filepath, datatype="log"):
        self.datatype = datatype
        self.read_log(filepath)

    def read_log(self, log):
       with open(log, 'r') as f:
            flist = list(f)
            log_format_test = 0
            for i, line in enumerate(flist):
                if line.startswith("  |MutationCounter_Modified"):
                    log_format_test += 1
                    modlength = []
                    for x in flist[i+6:i+27]:
                        modlength.append(float(x.strip().split('\t')[1]))
                    modmuts = []
                    for x in flist[i+32:i+53]:
                        modmuts.append(float(x.strip().split('\t')[1]))
                if line.startswith("  |MutationCounter_Untreated"):
                    log_format_test += 1
                    untlength = []
                    for x in flist[i+6:i+27]:
                        untlength.append(float(x.strip().split('\t')[1]))
                    untmuts = []
                    for x in flist[i+32:i+53]:
                        untmuts.append(float(x.strip().split('\t')[1]))
        message = ("Histogram data missing from log file. Requires" +
                   " --per-read-histogram flag when running ShapeMapper.")
        assert log_format_test >= 2, message
        data = {'Read_length': ['0-49', '50-99', '100-149', '150-199',
                                '200-249', '250-299', '300-349', '350-399',
                                '400-449', '450-499', '500-549', '550-599',
                                '600-649', '650-699', '700-749', '750-799',
                                '800-849', '850-899', '900-949', '950-999',
                                '>1000'],
                'Mutation_count': list(range(21)),
                'Modified_read_length': modlength,
                'Modified_mutations_per_molecule': modmuts,
                'Untreated_read_length': untlength,
                'Untreated_mutations_per_molecule': untmuts}
        self.data = pd.DataFrame(data)
