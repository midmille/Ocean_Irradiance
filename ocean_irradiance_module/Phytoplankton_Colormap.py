
from ocean_irradiance_module.PARAMS import Param_Init
import matplotlib as mpl

def Get_Phy_Cmap_Dict(): 
    """
    Get the color map for the plots

    Parameters
    ----------
    phy_species: List
        The number of different colors in the list.

    Returns 
    -------
    colors: List,  [N]
        The list of colors from the color map.
    """

    phy_species = Param_Init().phy_species
    
    N = len(phy_species)

    ## [The number of colors in the cmap.]
    Ncolors = 8

    ## [Can only do up to 10 colors.]
    if N > Ncolors:
        raise ValueError(f'Limit to N .LE. {Ncolors} colors, or increase number of colors in cmap')

#    cmap = mpl.cm.get_cmap('nipy_spectral')
    cmap = mpl.cm.get_cmap('Dark2')

    color_dict = {}
    colors = ['g', 'tab:purple', 'b', 'm', 'c', 'tab:brown','r', 'y']
    for k,phy in enumerate(phy_species): 
#        color_dict[phy] = (cmap(((k+ 0.5*(1/10))/ (10))))
        color_dict[phy] = colors[k]
        

    return color_dict


