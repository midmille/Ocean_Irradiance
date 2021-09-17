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
        self.Chl2NL = 1.59
        self.Chl2NS = .7950
        
        ## Average cosines and such taken from Dutkiewicz 2015. 
        self.r_s = 1.5 
        self.r_u = 3.0 
    
        self.v_d = .9 
        self.v_s = .83 
        self.v_u = .4 
        
        self.coefficients = (self.r_s, self.r_u, self.v_d, self.v_s, self.v_u)

        ##############################################################################
        ##############################################################################
        
    
