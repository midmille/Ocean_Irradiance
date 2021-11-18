"""
Created: November 1 2021 9:24am 

Author: Miles Miller 

"""


## External Mods
import numpy as np
## User made mod
from Read_CSV_Dat import Make_Data_Dict 
import ocean_irradiance_module.Ocean_Irradiance as OI
import ocean_irradiance_module.Ocean_Irradiance_ROMS as OIR
from ocean_irradiance_module.PARAMS import Param_Init
from ocean_irradiance_module.absorbtion_and_scattering_coefficients import absorbtion_scattering as abscat

def Argo_Data_to_Dict(file): 

    """
    This function reads in the data from the csv data file and returns the data 
    in the form of numpy arrays.
    """

     
    field_names = ('PLATFORM_CODE,DATE (YYYY-MM-DDTHH:MI:SSZ),DATE_QC,LATITUDE (degree_north),LONGITUDE (degree_east),'
                   'POSITION_QC,PRES (decibar),PRES_QC,PSAL (psu),PSAL_QC,TEMP (degree_Celsius),TEMP_QC,PRES_ADJUSTED (decibar),'
                   'PRES_ADJUSTED_QC,TEMP_ADJUSTED (degree_Celsius),TEMP_ADJUSTED_QC,PSAL_ADJUSTED (psu),PSAL_ADJUSTED_QC,'
                   'BBP700 (m-1),BBP700_QC,CPHL_ADJUSTED (milligram/m3),CPHL_ADJUSTED_QC,BBP700_ADJUSTED (m-1),BBP700_ADJUSTED_QC')
    
    data_dict = Make_Data_Dict(file, field_names, 1)

    return data_dict

    
def Calculate_Irradiance_Field(file, wavelength): 
    """
    This is to calculate the irradiance filed from the argo chla data profile
    """
    
    data_dict = Argo_Data_to_Dict(file) 
 
    ## The pressure is used to get the depth
    ## must change data type from object to float
    press = data_dict['PRES_ADJUSTED (decibar)'].astype(np.float64)
    ## The chla profile
    chla = data_dict['CPHL_ADJUSTED (milligram/m3)'].astype(np.float64)

    ## Changing pressure to depth
    ## pressure unit change to pascals
    press_pa = press * 10000
    ## The density of seewater [kg m^-3]
    rho_sea = 1023.6
    ## Equation for depth
    z = press_pa / (rho_sea * (-9.81))

    ## getting rid of the nans in the arrays
    z = np.flip(z[~np.isnan(chla)])
    chla = np.flip(chla[~np.isnan(chla)])

    ## Making the negative chla values zero
    chla[chla<0] = 0

    ## Params class
    PI = Param_Init()

    ## Using diatom coefficients
    ab_diat = abscat(wavelength, 'Diat')

    phy = OI.Phy(z, chla, ab_diat[0], ab_diat[1])

    irr_out = OI.ocean_irradiance_shoot_up(z[0],
                                           PI.Ed0, 
                                           PI.Es0,
                                           PI.Euh, 
                                           abscat(wavelength, 'water'),
                                           PI.coefficients, 
                                           phy = phy, 
                                           N = 100, 
                                           pt1_perc_zbot = True) 

    return irr_out 


def Calculate_RRS_Ratio(file):
    """
    Calculates the ratio of rrs for the two wavelengths required of the OCx algorithiom.
    """


    wavelengths = [443, 551]
    PI = Param_Init()

    RRS = np.zeros(len(wavelengths))

    for k, lam in enumerate(wavelengths): 
        irr_out = Calculate_Irradiance_Field(file, lam)

        RRS[k] = OIR.R_RS(PI.Ed0, PI.Es0, irr_out[2][-1])
        
    RRS_ratio = RRS[0] / RRS[1]

    chl_a = OIR.OCx_alg(RRS[0], RRS[1])

    return RRS_ratio, chl_a, RRS


















