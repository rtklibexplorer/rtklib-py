"""
RTKLIB-Py: Top level code for post-processing solution from rinex data

Copyright (c) 2022 Tim Everett
"""

import sys, os, shutil

# set run parameters
maxepoch = None # max # of epochs, used for debug, None = no limit
trace_level = 3  # debug trace level
basepos = []  # default to not specified here

######## specify input files ######################################

# cellphone example
# datadir = r'C:\gps\python\rtklib-py\data\phone'
# navfile = 'nav_1350.nav'
# rovfile = 'Pixel4_GnssLog.obs'
# basefile = 'slac1350.obs'
# cfgfile = 'config_phone.py' # must be in src folder or absolute path

# Google Smartphone Decimeter Challenge example
#datadir = r'C:\gps\data\cellphone\Google\2022\train\2021-12-08-US-LAX-1'
#navfile = 'AC0300USA_R_20213420000_01D_MN.rnx'
#rovfile = r'C:\gps\data\cellphone\Google\2022\train\2021-12-08-US-LAX-1\GooglePixel6Pro\supplemental\gnss_log.obs'
#basefile = 'base.obs'
#cfgfile = r'C:\gps\data\cellphone\Google\2022\train\config_phone_0401.py' # must be in src folder or absolute path
#basepos = [-2497836.5139, -4654543.2609, 3563028.9379] # LAX

# u-blox example
datadir = r'C:\gps\python\rtklib-py\data\u-blox'
navfile = 'rover.nav'
rovfile = 'rover.obs'
basefile = 'tmg23590.obs'
cfgfile = 'config_f9p.py'  # must be in src folder or absolute path

###################################################################

# Copy config file
shutil.copyfile(cfgfile, '__ppk_config.py')

# import rtklib files
import __ppk_config as cfg
import rinex as rn
import rtkcmn as gn
from rtkpos import rtkinit
from postpos import procpos, savesol

# generate output file names
solfile = rovfile[:-4] + '.pos'
if trace_level > 0:
    trcfile = os.path.join(datadir, rovfile[:-4] + '.trace')
    sys.stderr = open(trcfile, "w")
    
# Read config file
shutil.copyfile(cfgfile, '__ppk_config.py')


# init solution
os.chdir(datadir)
gn.tracelevel(trace_level)
nav = rtkinit(cfg)
nav.maxepoch = maxepoch

# load rover obs
rov = rn.rnx_decode(cfg)
print('Reading rover obs...')
if nav.filtertype == 'backward':
    maxepoch = None   # load all obs for 
rov.decode_obsfile(nav, rovfile, maxepoch)

# load base obs
base = rn.rnx_decode(cfg)
print('Reading base obs...')
base.decode_obsfile(nav, basefile, None)
if basepos != []:
    nav.rb = basepos
elif nav.rb[0] == 0:
    nav.rb = base.pos
    
# load nav data from rover obs
print('Reading nav data...')
rov.decode_nav(navfile, nav)

# calculate solution
print('Calculating solution ...\n')
sol = procpos(nav, rov, base)

# save solution to file
savesol(sol, solfile)





