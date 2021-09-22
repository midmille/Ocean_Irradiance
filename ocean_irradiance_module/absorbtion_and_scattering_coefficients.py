                                # -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 08:56:50 2020

@author: Miles Miller
"""

"""
Here I am making a function that  contains all phytoplankton absorbtion and 
scattering by wavelength
--> All values will be eyeballed form the 
--> going to create first dictionaries for each phy tht has each wavelength as key and 
each corresponding value of a tuple (a,b)
-->Then nest all those dictionaries into a big dictionary with the name of the 
plankton as the key, and its wavelength dictionary as the value. 
"""

"""

"""
## dictionaries for each plankton type 
## {wavelength:(a,b)}
## HLPro 

def absorbtion_scattering(wavelength, constituent):
    """
    
    

    Parameters
    ----------
    wavelength : int
        --> Must be one of the following wavelengths:
            410,443,486,551,638,671
    constituent : string
        --> Must be one of the following words: 
            water, HLPro, Cocco, Diat, Generic, Syn 

    Returns
    -------
    None.

    """
    
    ## Necessary unit change 
    ## ----------------------
    ## unit change for back scattering coefficients
    ## from units of Carbon to units of chl-a
    ## Large phytoplankton ie diatoms 
    LC2chla = 50
    # LC2chla = 1
    ## Small phytoplankton ie nanophytoplankton
    SC2chla = 100
    # SC2chla = 1
    
    
    water = {}

    water[410] = (0.001,.007)
    water[443] = (0.01,.004)
    water[486] = (.02,.003)
    water[551] = (.07,.002)
    water[638] = (.3,.0015)
    water[671] = (.45,.001)
    
    
    HLPro = {} 
    
    HLPro[410] = (.03,.0038)
    HLPro[443] = (.05,.004)
    HLPro[486] = (.025,.004)
    HLPro[551] = (.003,.0032)
    HLPro[638] = (.005,.0017)
    HLPro[671] = (.02,.0014)
    
    Cocco = {}
    
    Cocco[410] = (.023,.007)
    Cocco[443] = (.028,.006)
    Cocco[486] = (.026,.0079)
    Cocco[551] = (.009,.009)
    Cocco[638] = (.005,.009)
    Cocco[671] = (.02,.008)
    
    Diat = {}
    
    Diat[410] = (.015,.0038*LC2chla)
    Diat[443] = (.014,.0038*LC2chla)
    Diat[486] = (.01,.0039*LC2chla)
    Diat[551] = (.005,.004*LC2chla)
    Diat[638] = (.005,.0039*LC2chla)
    Diat[671] = (.012,.0038*LC2chla)
    
    Generic = {}
    
    Generic[410] = (.03,.0062)
    Generic[443] = (.036,.006)
    Generic[486] = (.025,.006)
    Generic[551] = (.014,.0059)
    Generic[638] = (.008,.0046)
    Generic[671] = (.017,.004)
    
    Syn = {}
    
    Syn[410] = (.036,.0078*SC2chla)
    Syn[443] = (.043,.007*SC2chla)
    Syn[486] = (.031,.0065*SC2chla)
    Syn[551] = (.031,.0051*SC2chla)
    Syn[638] = (.008,.004*SC2chla)
    Syn[671] = (.017,.0035*SC2chla)
    

    if constituent == 'water':
        return water[wavelength]
    if constituent == 'HLPro': 
        return HLPro[wavelength]
    if constituent == 'Cocco': 
        return Cocco[wavelength]
    if constituent == 'Diat':
        return Diat[wavelength]
    if constituent == 'Generic': 
        return Generic[wavelength]
    if constituent == 'Syn':
        return Syn[wavelength]
    
 



