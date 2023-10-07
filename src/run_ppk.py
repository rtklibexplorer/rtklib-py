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

# cell phone example from 2021 Google Smartphone Decimeter Challenge
# datadir = r'C:\gps\python\rtklib-py\data\phone'
# navfile = 'nav_1350.nav'
# rovfile = 'Pixel4_GnssLog.obs'
# basefile = 'slac1350.obs'
# cfgfile = 'config_phone.py' # must be in src folder or absolute path

# u-blox example
datadir = '../data/u-blox'
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
statfile = os.path.join(datadir, rovfile[:-4] + '.pos.stat')
fp_stat = open(statfile, 'w')
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
sol = procpos(nav, rov, base, fp_stat)

# save solution to file
savesol(sol, solfile)
fp_stat.close()
