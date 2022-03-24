#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy as np
import scipy.stats
import argparse 
 
# I'm using this command to import ReactivityProfile from where it lives on my computer
import sys
sys.path.append('/Users/anthonymustoe/Dropbox/Code/RingMapper/')
from ReactivityProfile import ReactivityProfile

    

def filterProfilesByNt(profile1, profile2, nts=None, name=None, exclhighbg=None):
    """Return matched reactivities from ReactivityProfile objects prof1 and prof2,
    filtering out any NaNs
    If nts argument is None, return all nts. Otherwise, return only these nts   
    name argument is passed to ReactivityProfile to get desired profile    
    """
    
    # lets check to make sure the sequences are the same!
    if not np.array_equal(profile1.sequence, profile2.sequence):
        raise ValueError('Sequences are not the same!')
        
    
    # initiliaze my mask (all false)
    mask = np.zeros(len(profile1.sequence), dtype=bool)
    
    if nts is None:
        mask[:] = True #reassign all values to True
    
    else:
        for n in nts:
            mask = mask | (profile1.sequence == n)
    
    if exclhighbg is not None:
        with np.errstate(invalid='ignore'):
            mask = mask & (profile1.backprofile < exclhighbg) & (profile2.backprofile < exclhighbg)
    
    
    # get the desired profiles
    r1 = profile1.profile(name)     
    r2 = profile2.profile(name)
    
    mask = mask & np.isfinite(r1) & np.isfinite(r2)
    
    return r1[mask], r2[mask]
    
    

def plotCorrelation(var1, var2, ax, title='', xlabel='', ylabel=''):
    
    regress = scipy.stats.linregress(var1, var2)
    
    # get min and max values for plotting
    xmin, ymin = min(var1), min(var2)
    gmin = min(xmin, ymin)
    xmax, ymax = max(var1), max(var2)
    gmax = max(xmax, ymax)
    
    ax.scatter(var1,var2)
    
    # plot the diagonal
    ax.plot([gmin, gmax], [gmin, gmax], 'k--', label='diagonal')

    # plot the fit
    x = np.linspace(xmin, xmax, num=100)
    y = x*regress.slope + regress.intercept
    ax.plot(x, y, 'r', label='fit')
    
    # add labels and stuff
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('{} : R={:.3f} ; m={:.1f}'.format(title, regress.rvalue, regress.slope))
    
    
#############################################################################

    


if __name__=='__main__':
    
    parser=argparse.ArgumentParser()
    parser.add_argument('profile1', help='Path to first profile file')
    parser.add_argument('profile2', help='Path to second profile file')
    parser.add_argument('output', help='Path to output file')
    parser.add_argument('--exclhighbg', type=float, help='Exclude nts with bg values greater than this cutoff')
    parser.add_argument('--comparison', default='all', help='Type of comparison to perform. Options are raw/sub/back/norm/all (default=all)')
    
    args=parser.parse_args()

    if args.comparison not in ('raw','sub', 'back', 'norm','all'):
        exit('Invalid comparison selected. Options are raw/sub/back/norm/all')



    # read in the desired files
    p1 = ReactivityProfile(args.profile1)
    p1.normalize(DMS=True)
    p2 = ReactivityProfile(args.profile2)
    p2.normalize(DMS=True)

    # create the matplotlib figure object
	
    if args.comparison == 'all':
	
        fig, ax = plot.subplots(4, 5, figsize=(20,16))

        for i,name in enumerate(('raw','sub','back','norm')):
	
    	    # this loop will compute correlations for all nts combined, and then individually
    	    for j, n in enumerate( ('All', 'A', 'C', 'U', 'G')):

                if n=='All':
                    x,y = filterProfilesByNt(p1, p2, name=name, exclhighbg=args.exclhighbg)
                else:
                    x,y = filterProfilesByNt(p1, p2, nts=n, name=name, exclhighbg=args.exclhighbg)
                
                #if name =='back':
                #    print("{} {} {}".format(n, np.median(x), np.median(y)))
        
                plotCorrelation(x,y, ax[i,j], title=name+' '+n)
	
	# label the axes

        ax[3,2].set_xlabel(args.profile1.split('/')[-1])
        ax[1,0].set_ylabel(args.profile2.split('/')[-1])


    else:
	
        fig, ax = plot.subplots(1, 5, figsize=(20,4))

        # this loop will compute correlations for all nts combined, and then individually
        for i, n in enumerate( ('All', 'A', 'C', 'U', 'G')):

       	   if n=='All':
               x,y = filterProfilesByNt(p1, p2, name=args.comparison, exclhighbg=args.exclhighbg)
           else:
               x,y = filterProfilesByNt(p1, p2, nts=n, name=args.comparison, exclhighbg=args.exclhighbg)
        
           plotCorrelation(x,y, ax[i], title=n)
	

        # label the axes
        ax[2].set_xlabel(args.profile1.split('/')[-1])
        ax[0].set_ylabel(args.profile2.split('/')[-1])

    #make a bit prettier
    fig.tight_layout()

    # save the data
    fig.savefig(args.output)



