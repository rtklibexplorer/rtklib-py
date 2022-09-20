"""
Run multiple GSDC 2022 PPK solutions with rtklib_py
"""

import sys
if 'rtklib_py/src' not in sys.path:
    sys.path.append('rtklib_py/src')
if 'android_rinex/src' not in sys.path:
    sys.path.append('android_rinex/src')

import os, shutil
from os.path import join, isdir, isfile
from glob import glob
from multiprocessing import Pool
import subprocess
import gnsslogger_to_rnx as rnx
from time import time

# set run parameters
maxepoch = None # max number of epochs, used for debug, None = no limit

# Set solution choices
ENABLE_PY = True        # Use rtklib_py to generate solutions 
OVERWRITE_RINEX = False  # overwrite existing rinex filex
OVERWRITE_SOL = False   # overwrite existing solution files

# specify location of input folder and files
datadir = r'C:\gps\GSDC_2022_py\data\test'
basefiles = '../base.obs' # rinex3, use this for python
navfiles = '../BRDM*MN.rnx' # navigation files with wild cards

PHONES = ['GooglePixel4', 'GooglePixel4XL', 'GooglePixel4Modded', 'GooglePixel5', 'GooglePixel6Pro', 'XiaomiMi8', 'SamsungGalaxyS20Ultra']

# Setup for rtklib_py
cfgfile_py = r'd:\gps\GSDC_2022_orig\config\ppk_phone_0625.py'
soltag_py = '_py0625'  # postfix for solution file names

# updated 9/18/22 to match bases.sta
BASE_POS = {'slac' : [-2703116.3161, -4291766.7898, 3854248.0816],  # WGS84 XYZ coordinates
            'vdcy' : [-2497836.8418, -4654543.0025, 3563029.0634],
            'p222' : [-2689640.5799, -4290437.1653, 3865051.0923]}


# input structure for rinex conversion
class args:
    def __init__(self):
        # Input parameters for conversion to rinex
        self.slip_mask = 0 # overwritten below
        self.fix_bias = True
        self.timeadj = 1e-7
        self.pseudorange_bias = 0
        self.filter_mode = 'sync'
        # Optional hader values for rinex files
        self.marker_name = ''
        self.observer = ''
        self.agency = ''
        self.receiver_number = ''
        self.receiver_type = ''
        self.receiver_version = ''
        self.antenna_number = ''
        self.antenna_type = ''

# Copy and read config file
if ENABLE_PY:
    shutil.copyfile(cfgfile_py, '__ppk_config.py')
    import __ppk_config as cfg
    import rinex as rn
    import rtkcmn as gn
    from rtkpos import rtkinit
    from postpos import procpos, savesol

# function to convert single rinex file
def convert_rnx(folder, rawFile, rovFile, slipMask):
    os.chdir(folder)
    argsIn = args()
    argsIn.input_log = rawFile
    argsIn.output = os.path.basename(rovFile)
    argsIn.slip_mask = slipMask
    rnx.convert2rnx(argsIn)

# function to run single rtklib_py solution
def run_ppk(folder, rovfile, basefile, navfile, solfile):
    # init solution
    os.chdir(folder)
    gn.tracelevel(0)
    nav = rtkinit(cfg)
    nav.maxepoch = maxepoch
    print(folder)

    # load rover obs
    rov = rn.rnx_decode(cfg)
    print('    Reading rover obs...')
    if nav.filtertype == 'backward':
        maxobs = None   # load all obs for backwards
    else:
        maxobs = maxepoch
    rov.decode_obsfile(nav, rovfile, maxobs)

    # load base obs and location
    base = rn.rnx_decode(cfg)
    print('   Reading base obs...')
    base.decode_obsfile(nav, basefile, None)
    
    # determine base location from original base obs file name
    if len(BASE_POS) > 1:
        baseName = glob('../*.2*o')[0][-12:-8]
        nav.rb[0:3]  = BASE_POS[baseName]
    elif nav.rb[0] == 0:
        nav.rb = base.pos # from obs file
        
    # load nav data from rover obs
    print('   Reading nav data...')
    rov.decode_nav(navfile, nav)

    # calculate solution
    print('    Calculating solution...')
    sol = procpos(nav, rov, base)

    # save solution to file
    savesol(sol, solfile)
    return rovfile


####### Start of main code ##########################

def main():

    # get list of data sets in data path
    datasets = os.listdir(datadir)

    # loop through data set folders
    rinexIn = []
    ppkIn = []
    rtklibIn = []
    for dataset in datasets:
        for phone in PHONES:
            # skip if no folder for this phone
            folder = join(datadir, dataset, phone)
            if not isdir(folder):  
                continue
            os.chdir(folder)
            rawFile = join('supplemental', 'gnss_log.txt')
            rovFile = join('supplemental', 'gnss_log.obs')

            rinex = False
            # check if need rinex conversion
            if OVERWRITE_RINEX or not isfile(rovFile):
                # generate list of input parameters for each rinex conversion
                if phone == 'SamsungS20Ultra': # Use cycle slip flags for Samsung phones
                    slipMask = 0 # 1 to unmask recevier cycle slips
                else:
                    slipMask = 0 
                rinexIn.append((folder, rawFile, rovFile, slipMask))
                print(rawFile, '->', rovFile) 
                rinex = True
            
            # check if need to create PPK solution
            try:
                baseFile = glob(basefiles)[0]
                navFile = glob(navfiles)[0]
                solFile = rovFile[:-4] + soltag_py + '.pos'
            except:
                print(folder,'  Error: Missing file')
                continue
            if ENABLE_PY and (OVERWRITE_SOL == True or len(glob(solFile)) == 0 
                              or rinex == True):
                # generate list of input/output files for each python ppk solution
                print('PY: ', join(dataset, phone))
                ppkIn.append((folder, rovFile, baseFile, navFile, solFile))

    if len(rinexIn) > 0:
        print('\nConvert rinex files...')
        # generate rinx obs files in parallel, does not give error messages
        # with Pool() as pool: # defaults to using cpu_count for number of procceses
        #     res = pool.starmap(convert_rnx, rinexIn)
        # run sequentially, use for debug
        for input in rinexIn:
            convert_rnx(input[0],input[1],input[2],input[3])

    if ENABLE_PY and len(ppkIn) > 0:
        print('Calculate PPK solutions...')
        # run PPK solutions in parallel, does not give error messages
        # with Pool() as pool: # defaults to using cpu_count for number of procceses
        #     res = pool.starmap(run_ppk, ppkIn)
        # run sequentially, use for debug
        for input in ppkIn:
            run_ppk(input[0],input[1],input[2],input[3],input[4])



if __name__ == '__main__':
    t0 = time()
    main()
    print('Runtime=%.1f' % (time() - t0))







