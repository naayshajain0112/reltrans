import numpy as np
from astropy.io import fits
import os
import sys
########################################################################
# Utility function to handle yes/no user input
def query_yes_no(question, default="no"):
     """
    Ask a yes/no question and return the user's response.

    Parameters:
        question (str): The question displayed to the user.
        default (str): Default answer if the user presses Enter.
                       Options are "yes", "no", or None.

    Returns:
        bool: True for "yes", False for "no"."""
    # Mapping of user inputs to boolean values
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        # Raise error if default value is invalid
        raise ValueError("invalid default answer: '%s'" % default)
    # Loop until user provides valid input
    while True:
        sys.stdout.write(question + prompt)
        # Get user input and convert to lowercase
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
########################################################################
xillver_table_name = ['xillverCp_v3.4.fits', 'xillver-a-Ec5.fits']

if (query_yes_no('Have you set the environment variable RELTRANS_TABLES?')):
    reltrans_table = os.environ['RELTRANS_TABLES']
    print(f'The RELTRANS_TABLES env variable is {reltrans_table}')
    print
else:
    print('Indicate the path of the folder where either you keep the xillver tables or you want to download them.')
    print('(the script creates a new folder if the given path point to a non-existing one)')
    reltrans_table = input('input the path: ')
    print(reltrans_table)
    print
    
print()
print(f'The first part of the script gives the option to download the xillver tables')
print(f'(reltrans needs the xillver tables: {xillver_table_name[0]} and {xillver_table_name[1]})')
print()
if (query_yes_no('Do you need to download the xillver tables?')):
    command = 'mkdir -p '+ str(reltrans_table)
    os.system(command)
    command = 'wget -c --retry-connrefused --waitretry=5 --timeout=30 -t 10 -P ' + str(reltrans_table) + ' https://sites.srl.caltech.edu/~javier/xillver/tables/'

    for name_table in xillver_table_name:
        print
        command_temp = command + name_table
        print(command_temp)
        os.system(command_temp)
        print


print()
print('This script is going to create new version of the xillver tables')
print('in the same folder xillver tables (keep in mind the disk space)')
print()
if (query_yes_no('Do you agree?')):

    for table_name in xillver_table_name:    
        fits_image_filename = str(reltrans_table) + '/' + str(table_name)
        try:
            with fits.open(fits_image_filename) as hdul:
                spectra_extension = hdul['SPECTRA'].data
                length = len(spectra_extension.field(0))
                if (table_name == 'xillverCp_v3.4.fits'):
                    logxi = spectra_extension.field(0)[:,2].reshape(length,1) # logxi ionisation values of the spectra
                    logne = spectra_extension.field(0)[:,4].reshape(length,1) #  logne density values of the spectra
                    spectra_extension['INTPSPEC'] = spectra_extension['INTPSPEC'] / 10**(logxi + logne - 15)
                if (table_name == 'xillver-a-Ec5.fits'):
                    logxi = spectra_extension.field(0)[:,2].reshape(length,1) # logxi ionisation values of the spectra
                    spectra_extension['INTPSPEC'] = spectra_extension['INTPSPEC'] / 10**(logxi)
                if (table_name == 'xillverD-5.fits'):
                    logxi = spectra_extension.field(0)[:,2].reshape(length,1) # logxi ionisation values of the spectra
                    logne = spectra_extension.field(0)[:,3].reshape(length,1) #  logne density values of the spectra
                    spectra_extension['INTPSPEC'] = spectra_extension['INTPSPEC'] / 10**(logxi + logne - 15)
                norm_table_name = fits_image_filename[:-5] + '_normalised.fits'
                try: 
                    hdul.writeto(norm_table_name)
                    print(f'The script created a new table called {norm_table_name}')
                    print()
                except:
                    print(f'Table {norm_table_name} has not been created since it already exists')
        except:
            print(f'There is no table called {fits_image_filename}')
else:
    print('No table has been created')        
