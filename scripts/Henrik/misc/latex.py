import matplotlib.pyplot as plt
import warnings
import matplotlib.colors as mcolors

from S2S.local_configuration import config

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def diverge_map(high=(0.565, 0.392, 0.173), low=(0.094, 0.310, 0.635)):
    '''
    low and high are colors that will be used for the two
    ends of the spectrum. they can be either color strings
    or rgb color tuples
    '''
    c = mcolors.ColorConverter().to_rgb
    if isinstance(low, str): low = c(low)
    if isinstance(high, str): high = c(high)
    return make_colormap([low, c('white'), 0.5, c('white'), high])

def diverge_map_grey(high=(0.565, 0.392, 0.173), low=(0.094, 0.310, 0.635)):
    '''
    low and high are colors that will be used for the two
    ends of the spectrum. they can be either color strings
    or rgb color tuples
    '''
    c = mcolors.ColorConverter().to_rgb
    if isinstance(low, str): low = c(low)
    if isinstance(high, str): high = c(high)
    return make_colormap([low, (226/255, 226/255, 226/255), 0.5, (226/255, 226/255, 226/255), high])

def cm_ETH_rb():
    return diverge_map(high=(168/255, 50/255, 45/255), low=(0/255, 122/255, 150/255))

def cm_ETH_yb():
    return diverge_map(high=(149/255, 96/255, 19/255), low=(0/255, 122/255, 150/255))

def cm_ETH_rgb():
    return diverge_map_grey(high=(168/255, 50/255, 45/255), low=(0/255, 122/255, 150/255))

def set_style(style='white'):

    if style=='white':
        basestyle = 'seaborn-whitegrid'
    else:
        basestyle = 'seaborn'

    plt.style.use(basestyle)
    plt.style.use(config['MPL_STYLE'])


def get_width():

    textwidth = 'textwidth'
    spec = 'No ' + textwidth + ' found'
    with open('/home/hauestad/thesis/main.log') as file:
        for line in file:
            if line[0]=='*':
                if line[3:12]==textwidth:
                    return float(line[13:-3])

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 441.01775
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

def save_figure(fig,path):

    with warnings.catch_warnings(record=True) as w:

        fig.savefig(path+'.pdf', format='pdf', bbox_inches='tight')

        if len(w)==0:

            return path+'.pdf',True
        else:
            print('constrained layout did not converge')
            return path+'.pdf',False

if __name__=='__main__':

    import os.path as path
    import matplotlib as mpl
    print("Your style sheets are located at: {}".format(path.join(mpl.__path__[0], 'mpl-data', 'stylelib')))
