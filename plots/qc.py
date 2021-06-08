
class qc():
    def __init__(self):
        self.logs = []
        self.profiles = []
        self.labels

##############################################################################
# Mutations per molecule plotting
#
#
##############################################################################

   def plot_MutsPerMol(self, sample="Modified", upper_limit=10):
        x = self.log['Mutation_count'][:upper_limit]
        y = self.log[sample+'_mutations_per_molecule'][:upper_limit]
        ax.plot(x, y, label=self.sample+": "+sample)

    def set_MutsPerMol(self, ax):
        ax.legend(title="Samples")
        ax.set(xlabel="Mutations per molecule",
               ylabel="Percentage of Reads",
               title='Mutations per molecule distribution')

    def make_MutsPerMol(self, ax):
        self.plot_log_MutsPerMol(ax, sample="Modified")
        self.plot_log_MutsPerMol(ax, sample="Untreated")
        self.set_log_MutsPerMol(ax)

    def plot_log_ReadLength(self, ax, sample="Modified", upper_limit=10,
                            n=1, of=1):
        width = 0.8/of
        x = np.arange(upper_limit) - 0.4 - (width/2) + (width*n)
        y = self.log[sample+'_read_length'][:upper_limit]
        ax.bar(x, y, width, label=self.sample+": "+sample)

    def set_log_ReadLength(self, ax, upper_limit=10):
        ax.legend(title="Samples")
        ax.set(xticks=range(upper_limit),
               xlabel='Read Length',
               ylabel='Percentage of Reads',
               title='Read length distribution')
        ax.set_xticklabels(self.log["Read_length"][:upper_limit],
                           rotation=45, ha='right', label=self.sample)

    def make_log_ReadLength(self, ax):
        self.plot_log_ReadLength(ax, sample="Modified", n=1, of=2)
        self.plot_log_ReadLength(ax, sample="Untreated", n=2, of=2)
        self.set_log_ReadLength(ax)

    def make_log_qc(self):
        fig, ax = plt.subplots(1, 3, figsize=(21, 7))
        self.make_log_MutsPerMol(ax[0])
        self.make_log_ReadLength(ax[1])
        self.plot_boxplot(ax[2])
        plt.tight_layout()

    def get_boxplot_data(self, sample=1,
                         cols=["Modified_rate", "Untreated_rate"]):
        data = self.profile[cols].copy()
        data = data.assign(Sample=sample)
        return data

    def plot_boxplot(self, ax, other_samples=[]):
        data = [self.get_boxplot_data(sample=1)]
        xticklabels = [self.sample]
        for i, sample in enumerate(other_samples):
            data.append(sample.get_boxplot_data(sample=i+2))
            xticklabels.append(sample.sample)
        data = pd.concat(data)
        data = pd.melt(data, id_vars=['Sample'], var_name=['Rate'])
        data.head()
        ax = sns.boxplot(x="Sample", y="value",
                         hue='Rate', data=data, orient='v')
        ax.set(yscale='log',
               ylim=(0.0005, 0.5),
               ylabel="Mutation Rate",
               xticklabels=xticklabels)
