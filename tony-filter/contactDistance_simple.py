import RNAtools2 as RNAtools

# To adapt for plotmapper dance groups
# sig and neg only require self.dance[i].ij_data['rings']
# ss only requires self.dance[i].ct
# cd filter on all requires [dance.ct for dance in self.dance]

# Load files
ringfile = 'deep/Mod100_newnull-0-rings.txt'
ctlist = [RNAtools.CT('Deep/Mod100-0.f.ct'), RNAtools.CT('Deep/Mod100-1.f.ct')]

# Set filters
filterneg = True  # require correlations to be positive
cdfilter = 15    # contact filter distance
sigfilter = 20   # G cutoff
ssfilter = True  # require both nts to be SS in the 1st file in ctlist

# Loop through lines of ring file and print if all filters pass
with open(ringfile) as inp:

    print(inp.readline(), end='')
    print(inp.readline(), end='')

    for line in inp:
        spl = line.split()
        i = int(spl[0])
        j = int(spl[1])

        # Skip line if any filters don't pass
        if ssfilter and (ctlist[0].ct[i-1] != 0 or ctlist[0].ct[j-1] != 0):
            continue
        if filterneg and int(spl[3]) < 0:
            continue
        if float(spl[2]) < sigfilter:
            continue
        cd = False
        for ct in ctlist:
            if ct.contactDistance(i, j) <= cdfilter:
                cd = True
                break
        if cd:
            continue

        print(line, end='')
