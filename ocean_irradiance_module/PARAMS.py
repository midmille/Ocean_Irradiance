# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 17:10:18 2021

@author: Miles Miller

This file defines all variables and main constants. 
"""






class Param_Init:
    """
    This class defines the default constants and ROMS parameters.
    """
    def __init__(self):
        #####################################INITIAL PARAMS###########################
        ##############################################################################

        ## The default boundary conditions of problem following John Wilkin.
        self.Ed0 = .7
        self.Es0 = 1 - self.Ed0
        self.Euh = 0 
        
        ## Wavelengths taken from John Wilkin email regarding satelite wavelengths.
        self.wavelengths = [410,443,486,551,638,671] 
        ## The Chlorophyll to Nitrogen unit change. 
        """
        Unit Change Information:
        -----------------------
        To change units from 'mmol Nitrogen m^-3' to 'mg chl-a m^-3'
        
        --> Convert from Nitrogen to Carbon => multiply by 106[Carbon]/16[Nitrogen]
            -- This is called the Redfield Ratio.

        --> Convert from mmols to mgs => multiply by 12[mg]/1[mmol]

        --> Convert from Carbon to Chl-a => multiply by:
            -- For Large Phytoplankton such as diatoms:
               - 1[gram chl-a]/50[gram carbon]
            -- For Small Phytoplankton such as nanophytoplankton: 
               - 1[gram chl-a]/( 100 to 200 [gram carbon])  
        
        Overall unit change: 
        NL2Chl = (106/16) * (12/1) * (1/50) = 1.59
        NS2chl = (106/16) * (12/1) * (1/100) = 0.795
        """
        self.NL2Chl = 1.59
        self.NS2Chl = .7950
        
        ## Average cosines and such taken from Dutkiewicz 2015. 
        self.r_s = 1.5 
        self.r_u = 3.0 
    
        self.v_d = .9 
        self.v_s = .83 
        self.v_u = .4 
        
        self.coefficients = (self.r_s, self.r_u, self.v_d, self.v_s, self.v_u)

        ##############################################################################
        ##############################################################################
        
    
