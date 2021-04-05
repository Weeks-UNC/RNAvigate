import matplotlib as mp
import seaborn as sns
import pandas as pd
from scipy import stats

sns.set_style("ticks")
sns.set_context("talk")

# profile.txt header names and indices
# 0 Nucleotide
# 1 Sequence
# 2 Modified_mutations
# 3 Modified_read_depth
# 4 Modified_effective_depth
# 5 Modified_rate
# 6 Modified_off_target_mapped_depth
# 7 Modified_low_mapq_mapped_depth
# 8 Modified_primer_pair_1_mapped_depth
# 9 Untreated_mutations
# 10 Untreated_read_depth
# 11 Untreated_effective_depth
# 12 Untreated_rate
# 13 Untreated_off_target_mapped_depth
# 14 Untreated_low_mapq_mapped_depth
# 15 Untreated_primer_pair_1_mapped_depth
# 16 Denatured_mutations
# 17 Denatured_read_depth
# 18 Denatured_effective_depth
# 19 Denatured_rate
# 20 Denatured_off_target_mapped_depth
# 21 Denatured_low_mapq_mapped_depth
# 22 Denatured_primer_pair_1_mapped_depth
# 23 Reactivity_profile
# 24 Std_err
# 25 HQ_profile
# 26 HQ_stderr
# 27 Norm_profile
# 28 Norm_stderr


def plotRegression(ax, p1, p2, ctfile='None'):
    ax.plot([0, 1], [0, 1], color='black')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(p1, p2)
    ax.text(0.1, 0.8,
            'R^2 = {:.2f}\nslope = {:.2f}'.format(r_value**2, gradient),
            transform=ax.transAxes)
    if ctfile != 'None':
        ct = pd.read_csv(ctfile, sep='\s+', usecols=[4], names=['j'], header=1)
        paired = ct.j != 0
        unpaired = ct.j == 0
        ax.scatter(p1[paired], p2[paired], label="Paired")
        ax.scatter(p1[unpaired], p2[unpaired], label="Unpaired")
    else:
        ax.scatter(p1, p2)


def readHistograms(logfile):
    with open(logfile, 'r') as f:
        flist = list(f)
        for i, line in enumerate(flist):
            if line.startswith("  |MutationCounter_Modified"):
                modlength = []
                for x in flist[i+6:i+27]:
                    modlength.append(float(x.strip().split('\t')[1]))
                modmuts = []
                for x in flist[i+32:i+53]:
                    modmuts.append(float(x.strip().split('\t')[1]))
            if line.startswith("  |MutationCounter_Untreated"):
                untlength = []
                for x in flist[i+6:i+27]:
                    untlength.append(float(x.strip().split('\t')[1]))
                untmuts = []
                for x in flist[i+32:i+53]:
                    untmuts.append(float(x.strip().split('\t')[1]))
    data = {'Read_length': ['(0,49)', '(50,99)', '(100,149)', '(150,199)',
                            '(200,249)', '(250,299)', '(300,349)', '(350,399)',
                            '(400,449)', '(450,499)', '(500,549)', '(550,599)',
                            '(600,649)', '(650,699)', '(700,749)', '(750,799)',
                            '(800,849)', '(850,899)', '(900,949)', '(950,999)',
                            '>1000'],
            'Mutation_count': list(range(21)),
            'Modified_read_length': modlength,
            'Modified_mutations_per_molecule': modmuts,
            'Untreated_read_length': untlength,
            'Untreated_mutations_per_molecule': untmuts}
    return pd.DataFrame(data)


def plotSkyline(axis, profile, label=None, column='Reactivity_profile'):
    x = []
    y = []
    for n, r in zip(profile['Nucleotide'], profile[column]):
        x.extend([n - 0.5, n + 0.5])
        y.extend([r, r])
    axis.plot(x, y, label=label)


def addSeqBar(axis, profile, yvalue=-0.017):
    font_prop = mp.font_manager.FontProperties(
        family="monospace", style="normal", weight="bold", size="12")
    color_dict = {"A": "#f20000", "U": "#f28f00",
                  "G": "#00509d", "C": "#00c200"}
    for i, row in profile.iterrows():
        col = color_dict[row['Sequence'].upper()]
        axis.annotate(row['Sequence'], xy=(i + 1, yvalue),
                      fontproperties=font_prop,
                      color=col, annotation_clip=False,
                      horizontalalignment="center")


def getWidth(sample):
    left_inches = 0.9
    right_inches = 0.4
    sp_width = len(sample['Nucleotide']) * 0.1
    fig_width = max(7, sp_width + left_inches + right_inches)
    return fig_width
