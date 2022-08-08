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

import numpy as np
from ocean_irradiance_module import PARAMS
import pandas

def Load_Coefficient(wavelength, constituent, opt_wat, opt_det, opt_phy, C2chla, LC2chla, SC2chla):
    """
    This loads the coefficient for the desired wavelength and constituent
    """
    if constituent == 'water': 

        a = opt_wat[opt_wat["wavelength"] == wavelength].values[0,1]
        b = opt_wat[opt_wat["wavelength"] == wavelength].values[0,2]

        return a,b

    if constituent == 'detritus': 
        
        ## [Retrieve the absorption, scattering, and back scattering from the data frame.]
        a = opt_det[opt_det['wavelength'] == wavelength].values[0,1]
        b = opt_det[opt_det['wavelength'] == wavelength].values[0,2]
        b_b = opt_det[opt_det['wavelength'] == wavelength].values[0,3]

        return a, b, b_b

    else: 

        ## [The if constituent == species_name syntax is used because we name some of 
        ##  the species differently then the text file does.]
        if constituent == 'Micromonas':
            spec_name = "Micromonas"
        if constituent == 'Syn':
            spec_name = "Syn"
            if C2chla == 'default': 
                C2chla = SC2chla
        if constituent == 'HLPro':
            spec_name = "HLPro"
        if constituent == 'LLPro':
            spec_name = "LLPro"
        if constituent == 'Diat':
            spec_name = "Diatom"
            if C2chla == 'default': 
                C2chla = LC2chla
        if constituent == 'Cocco':
            spec_name = "Coccoli"
        if constituent == 'Tricho':
            spec_name = "Tricho"
        if constituent == 'n2unicell':
            spec_name = "n2unicell"
        if constituent == 'Lgeuk':
            spec_name = "Lgeuk"
        if constituent == 'Generic':
            spec_name = "MEAN"
        if constituent == 'zooplankton':
            spec_name = "zooplankton"
        if constituent == 'heterobacteria':
            spec_name = "heterobacteria"
        ## [Retrieveing the species specific absorption and scattering for desired wavelength.]
        ## [Species specific portion of the data frame.]
        spec_df = opt_phy[opt_phy["species"] == spec_name]
        ## [Getting the wavelength specific values.]
        a = spec_df[spec_df["wavelength"] == wavelength].values[0,2]
        ## [Total Scattering is still in units of Carbon]
        b = spec_df[spec_df["wavelength"] == wavelength].values[0,4] * C2chla
            
        return a, b


    return

def absorbtion_scattering(wavelength, constituent, C2chla = 'default'):
    """
    Parameters
    ----------
    wavelength : int
        --> Must be one of the following wavelengths:
            410,443,486,551,638,671
    constituent : string
        --> Must be one of the following words: 
            water, HLPro, Cocco, Diat, Generic, Syn 
    C2chla: Float or 'default'
        The carbon to chla value. 
    dut_txt: Boolean 
        The flag for using the Dutkiwiecz values sent via an email on April 8, 2022 from Stephanie. 

    Returns
    -------
    None.

    """
    
    if C2chla == 'default':
        ## [Import the Carbon to Chla from the PARAMS.py file.] 
        PI = PARAMS.Param_Init()
        ## [Large phytoplankton C2Chla ratio.] 
        LC2chla = PI.LC2chla
        ## [Small phytoplankton C2Chla ration.]
        SC2chla = PI.SC2chla
        ## [Medium phytoplankton C2Chla ration.]
        C2chla = PI.C2chla
    else: 
        LC2chla = C2chla
        SC2chla = C2chla

    path = PI.path

    ## [The wavelengths of which the coefficients are explicitly defined.]
    wavelengths = np.arange(400, 725, 25)

    ## [The water data frame.]
    ## [Water data fram column names.]
    wat_cols = ['wavelength', 'a', 'b']
    opt_wat = pandas.read_csv(path + "optics_water.txt", header=None, sep=' ', skiprows=7, names=wat_cols)

    ## [Get the detritus data frame.]
    ## [The column labels for the data frame.]
    det_cols = ['wavelength', 
                'absorption (m2 particle-1)', 
                'scattering (m2 particle-1)', 
                'backscattering (m2 particle-1)']
    opt_det = pandas.read_csv(path+'optics_detritus.txt', header=None, sep=' ', skiprows=6, names=det_cols)

    ## [The phytoplankton data frame.]
    ## [Phytoplankton column names.]
    phy_cols = ["species", 
                "wavelength",
                "a (m2 mgChl-1)",
                "a_photosynthetic (m2 mgChl-1)",
                "b_tot (m2 mgC-1)",
                "non-spec bb coef (m2 mgC-1)",
                "a (m2 mgC-1)" ]
    opt_phy = pandas.read_csv(path + "optics_plankton.txt", 
                                header=None, 
                                sep=" ", 
                                skiprows=8, 
                                comment="*", 
                                names=phy_cols)


    ## [if the requested wavelenth is in the waveband of definde coefficients.]
    if np.any(wavelengths == wavelength): 
        if constituent == 'detritus':
            ## [Detritus returns thebackscattering coefficient as well.]
            a, b, b_b = Load_Coefficient(wavelength, constituent, opt_wat, opt_det, opt_phy, C2chla, LC2chla, SC2chla)
            return a, b, b_b
        else: 
            a, b = Load_Coefficient(wavelength, constituent, opt_wat, opt_det, opt_phy, C2chla, LC2chla, SC2chla)
            return a, b

    ## [If the requested wavelength is NOT in the waveband of defined coefficients
    ##  then an interpolation is performed.]
    else: 
        ## [The number of defined wavlengths.]
        Nlam = len(wavelengths)

        if constituent == 'detritus': 
            ## [The 1-D Arrays of the absorption, scattering, and back scattering for detritus.]
            a_arr = np.zeros(Nlam)
            b_arr = np.zeros(Nlam)
            b_b_arr = np.zeros(Nlam)

            ## [Loop over the wavelengths.]
            for k, lam in enumerate(wavelengths): 
                a, b, b_b = Load_Coefficient(lam, constituent, opt_wat, opt_det, opt_phy, C2chla, LC2chla, SC2chla)
                a_arr[k] = a
                b_arr[k] = b
                b_b_arr[k] = b_b

            ## [interpolate to the desired wavelength.]
            a = np.interp(wavelength, wavelengths, a_arr)
            b = np.interp(wavelength, wavelengths, b_arr)
            b_b = np.interp(wavelength, wavelengths, b_b_arr)

            return a, b, b_b

        ## [If not detritus.]
        else:
            ## [The 1-D Arrays of the absorption, scattering, and back scattering for detritus.]
            a_arr = np.zeros(Nlam)
            b_arr = np.zeros(Nlam)

            ## [Loop over the wavelengths.]
            for k, lam in enumerate(wavelengths): 
                a, b = Load_Coefficient(lam, constituent, opt_wat, opt_det, opt_phy, C2chla, LC2chla, SC2chla)
                a_arr[k] = a
                b_arr[k] = b

            ## [interpolate to the desired wavelength.]
            a = np.interp(wavelength, wavelengths, a_arr)
            b = np.interp(wavelength, wavelengths, b_arr)

            return a, b


def equivalent_spherical_diameter(species):
    """
    Retrieves the equivalent spherical diameter for each species of phytoplankton
    
    The ESD is returned in units of micro meters.

    See Notes for 12/7
    """ 

    if species == 'HLPro': 
        return 0.6 
    if species == 'LLPro': 
        return 0.6 
    if species == 'Cocco': 
        return 4
    if species == 'Diat':
        return 17
    if species == 'Generic': 
        return 10 
    if species == 'Syn':
        return 0.98 
    if species == 'Lgeuk':
        return 27.64
    if species == 'Tricho':
       return 6.00
    






"""

The old eyeballed coefficients
------------------------------

    if dut_txt == False: 
        ## [The start of the different coefficient dictionaries.]
        water = {}
    
        ## These values were edited on 1/28/2022 with Chris and Jonathan. 
        ## Absorption and Scattering values edited to be values taken from optical properties of water
        ## Light and Water Chapter 3: Optical Properties of Water Table 3.5 page 90. 
        ## http://misclab.umeoce.maine.edu/boss/classes/RT_Weizmann/Chapter3.pdf
        water[410] = np.array((0.0162,.007))  
        water[443] = np.array((0.0145,.0049))
        water[486] = np.array((.0196,.0031))
        water[551] = np.array((.0638,.0019))
        water[638] = np.array((.329,.0010))
        water[671] = np.array((.430,.0008))
     
     
        HLPro = {} 
     
        HLPro[410] = np.array((.03,.0038*C2chla))
        HLPro[443] = np.array((.05,.004*C2chla))
        HLPro[486] = np.array((.025,.004*C2chla))
        HLPro[551] = np.array((.003,.0032*C2chla))
        HLPro[638] = np.array((.005,.0017*C2chla))
        HLPro[671] = np.array((.02,.0014*C2chla))
     
        Cocco = {}
     
        Cocco[410] = np.array((.023,.007*C2chla))
        Cocco[443] = np.array((.028,.006*C2chla))
        Cocco[486] = np.array((.026,.0079*C2chla))
        Cocco[551] = np.array((.009,.009*C2chla))
        Cocco[638] = np.array((.005,.009*C2chla))
        Cocco[671] = np.array((.02,.008*C2chla))
     
        Diat = {}
     
        Diat[410] = np.array((.015,.0038*LC2chla))
        Diat[443] = np.array((.014,.0038*LC2chla))
        Diat[486] = np.array((.01,.0039*LC2chla))
        Diat[551] = np.array((.005,.004*LC2chla))
        Diat[638] = np.array((.005,.0039*LC2chla))
        Diat[671] = np.array((.012,.0038*LC2chla))
     
        Generic = {}
     
        Generic[410] = np.array((.03,.0062*C2chla))
        Generic[443] = np.array((.036,.006*C2chla))
        Generic[486] = np.array((.025,.006*C2chla))
        Generic[551] = np.array((.014,.0059*C2chla))
        Generic[638] = np.array((.008,.0046*C2chla))
        Generic[671] = np.array((.017,.004*C2chla))
     
        Syn = {}
     
        Syn[410] = np.array((.036,.0078*SC2chla))
        Syn[443] = np.array((.043,.007*SC2chla))
        Syn[486] = np.array((.031,.0065*SC2chla))
        Syn[551] = np.array((.031,.0051*SC2chla))
        Syn[638] = np.array((.008,.004*SC2chla))
        Syn[671] = np.array((.017,.0035*SC2chla))


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
"""

