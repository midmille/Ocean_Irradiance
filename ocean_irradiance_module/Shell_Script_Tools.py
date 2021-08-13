# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:39:08 2021

@author: miles

This file will contain all the tools that help to run experiments on with ROMS 
"""
import re 
import os
import sys
import subprocess


def Edit_ROMS_In_File(file, old_inst, new_inst ):
    """
    

    Parameters
    ----------
    file : String
        Absolute path to file to be edited. 
    old_inst : String
        The string instance to be replaced. 
    new_inst : String
        The string instance the will replace the old instance.

    Returns
    -------
    None.

    """
    
    with open(file, 'r+') as f:
        txt = f.read()
        txt = re.sub(old_inst, new_inst, txt)
        f.seek(0)
        f.write(txt)
        f.truncate()
        
    return 
    
    
def Make_Out_Dir(out_dir):
        """
        This File makes the out directory for all the roms_his nc files for each 
        run case to go into. There is some user input required. 

        Parameters
        ----------
        out_dir : String
            Absoilute path and name of directory to be created. 

        Returns
        -------
        None.

        """
        if os.path.exists(out_dir) == True: 
            yon = input(f'Output directory [{out_dir}] already exists... Would you like to destory and overwrite? [yes/no] ')
            ## If no, exit
            if yon=='no': 
                sys.exit(f'You said {yon}. Exiting Run. Please Try again later.')
            if yon=='yes': 
                print("FABULOUS! Starting Experiment.... and destroying old files")
                ## Remove direcetory and contents.
                out = subprocess.run(['rm', '-rfv', out_dir])
                ## Check for completion.
                if out.returncode != 0:
                    sys.exit('Non-Zero Retun Code... Exiting')
                ## Make new empty directory.
                os.mkdir(out_dir)
        
        elif os.path.exists(out_dir) == False:
            os.mkdir(out_dir)
        
        return 
    