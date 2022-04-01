"""
 post-processing solution from rinex data
"""
import numpy as np
from copy import copy, deepcopy
from rtkpos import rtkpos, rtkinit
import __ppk_config as cfg
import rinex as rn
from pntpos import pntpos
import rtkcmn as gn
#import rtkcmn as gn

def combres(solf, solb):
    # combine forward/backward solutions 
    i, j, solc = 0, len(solb) - 1, []
    while i < len(solf) and j >= 0:
        tt = gn.timediff(solf[i].t, solb[j].t)
        if tt < -gn.DTTOL:
            sol = deepcopy(solf[i])
            j += 1
        elif tt > gn.DTTOL:
            sol = deepcopy(solb[j])
            i -= 1
        elif solf[i].stat < solb[j].stat:
            sol = deepcopy(solf[i])
        elif solf[i].stat > solb[j].stat:
            sol = deepcopy(solb[i]) 
        else:
            sol = deepcopy(solf[i])
            sol.t = gn.timeadd(sol.t, -tt / 2)
            sol.rr[0:3], sol.qr[0:3,0:3] = gn.smoother(solf[i].rr[0:3], 
                solb[j].rr[0:3], solf[i].qr, solb[j].qr)
        solc.append(sol)
        i, j = i + 1, j - 1
    return solc


def firstpos(nav, rov, base, dir):
    # find rover position from first obs, 
    obsr, obsb = rn.first_obs(nav, rov, base, dir)
    sol = pntpos(obsr, nav)
    # repeat until get solution
    while sol.stat == gn.SOLQ_NONE:
        obsr, obsb = rn.next_obs(nav, rov, base, dir)
        sol = pntpos(obsr, nav)
    nav.x[0:6] = copy(sol.rr[0:6])
    nav.rr[0:3] = copy(sol.rr[0:3])
    
def sqrtvar(cov):
    " sqrt of covariance "
    return np.sqrt(abs(cov)) * np.sign(cov)
    
def savesol(sol, solfile):
    D2R = gn.rCST.D2R
    solhdr = '%  GPST          latitude(deg) longitude(deg)  height(m)   Q  ' \
        'ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio\n'
    with open(solfile, 'w') as outfile:
        outfile.write(solhdr)
        for s in sol:
                wk, sec = gn.time2gpst(s.t)
                llh = gn.ecef2pos(s.rr[0:3])
                std = sqrtvar(gn.covenu(llh, s.qr))
                fmt = '%4d %10.3f %14.9f %14.9f %10.4f %3d %3d %8.4f' + \
                    '  %8.4f %8.4f %8.4f %8.4f %8.4f %6.2f %6.1f\n'
                outfile.write(fmt % (wk, sec, llh[0]/D2R, llh[1]/D2R, llh[2], 
                    s.stat, s.ns, std[1,1], std[0,0], std[2,2], std[0,1],
                    std[2,0], std[1,2], s.age, s.ratio))

def procpos(nav, rov, base):
 
    try:
        if nav.filtertype != 'backward':
            firstpos(nav, rov, base, dir=1)
            rtkpos(nav, rov, base, dir=1) # run forward solution
            sol0 = deepcopy(nav.sol)
            savesol(sol0,'forward.pos')
        if nav.filtertype != 'forward':
            if nav.filtertype == 'combined':
                # reset filter states
                rb = nav.rb.copy()
                eph = nav.eph.copy()
                maxepoch = nav.maxepoch
                nav = rtkinit(cfg)
                nav.rb = rb
                nav.eph = eph
                nav.maxepoch = maxepoch
            elif nav.filtertype == 'combined_noreset':
                nav.sol = []
            firstpos(nav, rov, base, dir=-1)
            rtkpos(nav, rov, base, dir=-1)  # run backward solution
            savesol(nav.sol,'backward.pos')
        if nav.filtertype == 'combined':
            sol = combres(sol0, nav.sol)
            savesol(sol,'combined.pos')
            return sol
    except KeyboardInterrupt:
        pass
    return nav.sol
