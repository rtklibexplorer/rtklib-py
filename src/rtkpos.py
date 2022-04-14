"""
module for PPK positioning

Copyright (c) 2021 Rui Hirokawa (from CSSRLIB)
Copyright (c) 2022 Tim Everett
"""

import numpy as np
from numpy.linalg import inv, norm
from sys import stdout
from copy import copy, deepcopy
import rtkcmn as gn
from rtkcmn import rCST, DTTOL, sat2prn, sat2freq, timediff, xyz2enu
import rinex as rn
from pntpos import pntpos
from ephemeris import satposs
from mlambda import mlambda
from rtkcmn import trace, tracemat, uGNSS
import __ppk_config as cfg
#rom ppp import tidedisp

MAX_VAR_EPH = 300**2

def rtkinit(cfg):
    nav = gn.Nav(cfg)
    """ initalize RTK-GNSS parameters from config file """
    nav.gnss_t = cfg.gnss_t
    nav.pmode = cfg.pmode
    nav.filtertype = cfg.filtertype
    # add rover vel and accel states for kinematic solution
    nav.na = nav.nq = 3 if nav.pmode == 'static' else 9
    nav.nx = nav.na + uGNSS.MAXSAT * nav.nf
    nav.x = np.zeros(nav.nx)
    nav.P = np.zeros((nav.nx, nav.nx))
    nav.xa = np.zeros(nav.na)
    nav.Pa = np.zeros((nav.na, nav.na))
    nav.el = np.zeros(uGNSS.MAXSAT)
    nav.gf = np.zeros(uGNSS.MAXSAT)
    nav.ph = np.zeros((2, uGNSS.MAXSAT, nav.nf))
    nav.pt = np.empty((2, uGNSS.MAXSAT, nav.nf), dtype=object)
    nav.nfix = nav.neb = nav.tt = 0
    
    nav.rb = cfg.rb

    # parameter for RTK/PPK    
    nav.use_sing_pos = cfg.use_sing_pos
    nav.cnr_min = cfg.cnr_min
    nav.maxout = cfg.maxout  # maximum outage [epoch]
    nav.elmin = np.deg2rad(cfg.elmin)
    nav.nf = cfg.nf
    nav.excsats = cfg.excsats
    nav.freq = cfg.freq
    nav.dfreq_glo = cfg.dfreq_glo
    nav.interp_base = cfg.interp_base
    nav.gnss_t = cfg.gnss_t
    nav.maxinno = cfg.maxinno
    nav.thresdop = cfg.thresdop
    nav.thresslip = cfg.thresslip
    nav.maxage = cfg.maxage
    nav.accelh = cfg.accelh
    nav.accelv = cfg.accelv
    nav.prnbias = cfg.prnbias
    
    # ambiguity resolution
    nav.armode = cfg.armode
    nav.ratio = 0
    nav.thresar = cfg.thresar
    nav.thresar1 = cfg.thresar1
    nav.var_holdamb = cfg.var_holdamb
    nav.elmaskar = np.deg2rad(cfg.elmaskar)
    nav.mindropsats = cfg.mindropsats
    nav.excsat_ix = 0
    
    # statistics
    nav.efact = cfg.efact
    nav.eratio = cfg.eratio
    nav.err = np.array(cfg.err) 
    nav.sig_p0 = cfg.sig_p0
    nav.sig_v0 = cfg.sig_v0
    nav.sig_n0 = cfg.sig_n0
    
    # solution parameters
    nav.sol = []

    dP = np.diag(nav.P)
    dP.flags['WRITEABLE'] = True
    dP[0:3] = nav.sig_p0**2
    if nav.pmode == 'kinematic':
        dP[3:9] = nav.sig_v0**2

    # obs index
    ix0, ix1 = cfg.freq_ix0, cfg.freq_ix1
    freq0 = {k: cfg.freq[ix0[k]] for k in ix0.keys()}
    freq1 = {k: cfg.freq[ix1[k]] for k in ix1.keys()}
    nav.obs_idx = [ix0, ix1]
    nav.obs_freq = [freq0, freq1]

    # sat index
    nav.sysprn = {i: gn.sat2prn(i) for i in range(1, uGNSS.MAXSAT+1)}
        
    return nav

def zdres_sat(nav, obs, r, rtype, dant, ix):
    _c = rCST.CLIGHT
    nf = nav.nf
    y = np.zeros(nf * 2)
    for f in range(nf):
        freq = sat2freq(obs.sat[ix], f, nav)
        if obs.S[ix,f] < nav.cnr_min:
            continue
        # residuals = observable - estimated range (phase and code)
        y[f] = obs.L[ix,f] * _c / freq - r - dant[f] if obs.L[ix,f] else 0
        y[f+nf] = obs.P[ix,f] - r - dant[f] if obs.P[ix,f] else 0
#        if obs.L[ix,f] != 0 or obs.P[ix,f] != 0:
#            trace(3, 'zdres_sat: %d: L=%.6f P=%.6f r=%.6f f=%.0f\n' % 
#                  (obs.sat[ix],obs.L[ix,f], obs.P[ix,f],r,freq))
    return y
        
def zdres(nav, obs, rs, dts, svh, var, rr, rtype):
    """ undifferenced phase/code residuals ----------------------------------------
    calculate zero diff residuals [observed pseudorange - range] 
        output is in y[0:nu-1], only shared input with base is nav 
 args:  I   obs  = sat observations
        I   n    = # of sats
        I   rs = sat position {x,y,z} (m)
        I   dts = sat clock {bias,drift} (s|s/s)
        I   var  = variance of ephemeris
        I   svh  = sat health flags
        I   nav  = sat nav data
        I   rr   = rcvr pos (x,y,z)
        I   rtype:  0=base,1=rover 
        O   y[] = zero diff residuals {phase,code} (m)
        O   e    = line of sight unit vectors to sats
        O   azel = [az, el] to sats  """
    if obs == []:
        return [], [], []
    _c = rCST.CLIGHT
    nf = nav.nf
    n = len(obs.P)
    y = np.zeros((n, nf * 2))
    el = np.zeros(n)
    e = np.zeros((n, 3))
    rr_ = rr.copy()
    trace(3, 'zdres: n=%d rr=%.2f %.2f %.2f\n' % (n, rr[0], rr[1], rr[2]))
    pos = gn.ecef2pos(rr_)
    # loop through satellites
    ix = np.argsort(obs.sat)
    for i in ix:
        if gn.satexclude(obs.sat[i], var[i], svh[i], nav):
            continue
        # compute geometric-range and azimuth/elevation angle
        r, e[i,:] = gn.geodist(rs[i,0:3], rr_)
        _, el[i] = gn.satazel(pos, e[i,:])
        if el[i] < nav.elmin:
            continue
        # adjust range for satellite clock-bias
        r += -_c * dts[i]
        # adjust range for troposphere delay model (hydrostatic)
        trophs, tropw, _ = gn.tropmodel(obs.t, pos, np.deg2rad(90.0), 0.0)
        zhd = trophs + tropw
        mapfh, _ = gn.tropmapf(obs.t, pos, el[i])
        r += mapfh * zhd
        # calc receiver antenna phase center correction
        dant = gn.antmodel(nav, el[i], nav.nf, rtype)
        trace(4,'sat=%d r=%.6f c*dts=%.6f zhd=%.6f map=%.6f\n' % (obs.sat[i],r,_c*dts[i],zhd,mapfh))
        # calc undifferenced phase/code residual for satellite
        y[i] = zdres_sat(nav, obs, r, rtype, dant, i)
    
    for i in ix:
        if obs.L[i,0] != 0 and rtype == 0:
            trace(4, 'sat=%2d %13.3f %13.3f %13.3f %13.10f %5.1f\n' %
                  (obs.sat[i], rs[i,0], rs[i,1], rs[i,2], dts[i], 
                   np.rad2deg(el[i])))
    
    tracemat(3, 'y=', y[ix,:].T, '13.3f')
    return y, e, el


def ddcov(nb, n, Ri, Rj, nv):
    """ double-differenced measurement error covariance ---------------------------
*
*   nb[n]:  # of sat pairs in group
*   n:      # of groups (2 for each system, phase and code)
*   Ri[nv]: variances of first sats in double diff pairs
*   Rj[nv]: variances of 2nd sats in double diff pairs
*   nv:     total # of sat pairs 
*   R[nv][nv]:  double diff measurement err covariance matrix """
    R = np.zeros((nv, nv))
    k = 0
    for b in range(n):
        block = R[k:nb[b]+k, k:nb[b]+k] # define subarray
        block += Ri[k:nb[b]+k]
        block[range(nb[b]), range(nb[b])] += Rj[k:nb[b]+k]
        k += nb[b]
    return R


def sysidx(satlist, sys_ref):
    """ return index of satellites with sys=sys_ref """
    idx = []
    for k, sat in enumerate(satlist):
        sys, _ = sat2prn(sat)
        if sys == sys_ref:
            idx.append(k)
    return idx


def IB(s, f, na=3):
    """ return index of phase ambguity """
    return na + uGNSS.MAXSAT * f + s - 1


def varerr(nav, sys, el, f, dt, rcvstd):
    """ variation of measurement """
    code = 1 * (f >= nav.nf) # 0 = phase, 1 = code
    freq = f % nav.nf
    s_el = np.sin(el)
    if s_el <= 0.0: 
        return 0.0
    fact = nav.eratio[freq] if code else 1
    fact *= nav.efact[sys]
    a, b = fact * nav.err[1:3]
    c = fact * 0  # nav.err[4]*bl/1E4  # TODO: add baseline term
    d = rCST.CLIGHT * nav.err[5] * dt # clock term
    var = 2.0 * (a**2 + (b / s_el)**2 + c**2) + d**2
    # TODO: add SNR term
    # add scaled stdevs from receiver
    if nav.err[3] > 0:
        var += (nav.err[3] * rcvstd)**2
    return var


def ddres(nav, x, P, yr, er, yu, eu, sat, el, dt, obsr):
    """ /* double-differenced residuals and partial derivatives  -----------------------------------
        I nav  = sat nav data
        I dt = time diff between base and rover observations
        I x = rover pos & vel and sat phase biases (float solution)
        I P = error covariance matrix of float states
        I sat = list of common sats
        I y = zero diff residuals (code and phase, base and rover)
        I e = line of sight unit vectors to sats
        I el = el to sats
        O v = double diff innovations (measurement-model) (phase and code)
        O H = linearized translation from innovations to states (az/el to sats)
        O R = measurement error covariances """
    _c = rCST.CLIGHT
    nf = nav.nf
    ns = len(el)
    ny = ns * nf * 2  # phase and code
    nb = np.zeros(2 * len(nav.gnss_t) * nf, dtype=int)
    Ri = np.zeros(ny)
    Rj = np.zeros(ny)
    H = np.zeros((nav.nx, ny))
    P_init = nav.sig_n0**2 # value used to initialize P states
    trace(3,"ddres   : dt=%.4f ns=%d\n" % (dt, ns))

    nv = b = 0
    v = np.zeros(ny)
    # step through sat systems
    for sys in nav.gnss_t:
        # step through phases/codes
        for f in range(0, nf*2):
            frq = f % nf
            code = 1 * (f >= nf)
            idx = sysidx(sat, sys) # find sats in sys
            # remove sats with missing base or rover residuals
            nozero = np.where((yr[:,f] != 0) & (yu[:,f] != 0))[0]
            idx = np.intersect1d(idx, nozero)
            if len(idx) == 0: 
                continue  # no common sats
            # find sat with max el and not just reset for reference
            i_el = idx[np.argsort(el[idx])]
            for i in i_el[::-1]:
                ii = IB(sat[i], frq, nav.na)
                # check if sat just reset
                if  P[ii,ii] != nav.sig_n0**2: 
                    break
            else:
                i = i_el[0] # use highest sat if none without reset
            # calculate double differences of residuals (code/phase) for each sat
            lami = _c / sat2freq(sat[i], frq, nav)
            for j in idx: # loop through sats
                if i == j: continue  # skip ref sat
                #  double-differenced measurements from 2 receivers and 2 sats in meters 
                v[nv] = (yu[i,f] - yr[i,f]) - (yu[j,f] - yr[j,f])
                # partial derivatives by rover position, combine unit vectors from two sats
                H[0:3, nv] = -eu[i,:] + er[j,:]
                
                jj = IB(sat[j], frq, nav.na)
                if not code:  # carrier phase
                    # adjust phase residual by double-differenced phase-bias term
                    lamj = _c / sat2freq(sat[j], frq, nav)
                    v[nv] -= lami * x[ii] - lamj * x[jj]
                    H[ii, nv], H[jj, nv] = lami, -lamj
                
                # if residual too large, flag as outlier
                thres = nav.maxinno
                # use larger thresh for code or just initialized phase
                if code or P[ii,ii] == P_init or P[jj,jj] == P_init:
                    thres *= nav.eratio[frq] 
                if abs(v[nv]) > thres:
                    nav.vsat[sat[j]-1,frq] = 0
                    nav.rejc[sat[j]-1,frq] += 1
                    trace(3,"outlier rejected: (sat=%3d-%3d %s%d v=%13.3f x=%13.3f %13.3f P=%.6f %.6f)\n" 
                          % (sat[i], sat[j], 'LP'[code],frq+1, v[nv], x[ii], x[jj], P[ii,ii],P[jj,jj]))
                    H[ii, nv], H[jj, nv] = 0, 0
                    continue 
                # single-differenced measurement error variances (m)
                si, sj = sat[i] - 1, sat[j] - 1
                Ri[nv] = varerr(nav, sys, el[i], f, dt, nav.rcvstd[si,f])
                Rj[nv] = varerr(nav, sys, el[j], f, dt, nav.rcvstd[sj,f])
                if not code:
                    nav.vsat[si,frq] = 1
                    nav.vsat[sj,frq] = 1
                #trace(3,'sys=%d f=%d,i=%d,j=%d\n' % (sys,f,i,j))
                trace(3,"sat=%3d-%3d %s%d v=%13.3f R=%9.6f %9.6f lock=%2d x=%13.3f\n" %
                      (sat[i], sat[j], 'LP'[code], frq+1, v[nv], Ri[nv], Rj[nv], nav.lock[sat[j]-1,frq],x[jj]))
                nv += 1
                nb[b] += 1
            b += 1
    R = ddcov(nb, b, Ri[:nv], Rj[:nv], nv)

    return v[:nv], H[:,:nv], R


def valpos(nav, v, R, thres=4.0):
    """ post-file residual test """
    trace(3, 'valpos  : nv=%d thres=%.1f\n' % (len(v), thres))
    nv = len(v)
    fact = thres**2
    for i in range(nv):
        if v[i]**2 > fact * R[i, i]:
            trace(3, 'large residual (ix_sat=%d  v=%.3f sig=%.3f)\n' % 
                  (i, v[i], np.sqrt(R[i, i])))
    return True

def intpres(time, nav, y0, y1, obs0, obs1):
    """ time-interpolation of residuals """
    tt, ttb = timediff(time, obs1.t), timediff(time, obs0.t)
    if len(y0) == 0 or abs(ttb) > nav.maxage or abs(tt) < DTTOL:
        return y1, tt
    # find common sats
    _, ix0, ix1 = np.intersect1d(obs0.sat, obs1.sat, return_indices=True)
    for i in range(len(ix0)):
        for j in range(4):
            i0, i1 = ix0[i], ix1[i]
            if y1[i1,j] == 0:
                y1[i1,j] = y0[i0,j]
            elif y0[i0,j] != 0:
                y1[i1,j] = (ttb * y1[i1,j] - tt * y0[i0,j]) / (ttb - tt)
    dt = min(abs(tt), abs(ttb)) / np.sqrt(2)
    return y1, dt

    

def ddidx(nav, sats):
    """ index for single to double-difference transformation matrix (D') """
    nb, fix, ref = 0, [], []
    ns = uGNSS.MAXSAT
    #na = nav.na
    ix = np.zeros((ns, 2), dtype=int)
    # clear fix flag for all sats (1=float, 2=fix)
    nav.fix[:,:] = 0
    # step through constellations
    for m in range(uGNSS.GNSSMAX):
        k = nav.na  # state index for first sat
        # step through freqs
        for f in range(nav.nf):
            # look for first valid sat (i=state index, i-k=sat index)
            for i in range(k, k + ns):
                sati = i - k + 1
                # if sati not in sats:
                #     xxx=1
                sys = nav.sysprn[sati][0]
                # skip if sat not active
                if nav.x[i] == 0.0  or sys != m or nav.vsat[sati-1,f] == 0:
                    continue
                if nav.lock[sati-1,f] >= 0 and nav.slip[sati-1,f] & 2 == 0 and \
                        nav.el[sati-1] >= nav.elmaskar:
                    # set sat to use for fixing ambiguity if meets criteria
                    nav.fix[sati-1,f] = 2 # fix
                    break # break out of loop if find good sat
                else: # don't use this sat for fixing ambiguity
                    nav.fix[sati-1,f] = 1 # float
            if nav.fix[sati-1,f] != 2: # no good sat found
                continue
            n = 0  # count of sat pairs for this freq/constellation
            # step through all sats (j=state index, j-k=sat index, i-k=first good sat)
            for j in range(k, k + ns):
                satj = j - k + 1
                sys = nav.sysprn[satj][0]
                if i == j or nav.x[j] == 0.0 or sys != m or nav.vsat[satj-1,f] <= 0:
                    continue
                if nav.lock[satj-1,f] >= 0 and nav.slip[satj-1,f] & 2 == 0 and \
                        nav.el[satj-1] >= nav.elmaskar:
                    # set D coeffs to subtract sat j from sat i
                    ix[nb, :] = [i,j] # state indices of ref bias and target bias
                    ref.append(sati)
                    fix.append(satj)
                    nav.fix[satj-1,f] = 2 # fix
                    nb += 1 # increment total count 
                    n += 1 # inc count in freq/constellation
                else: # don't use this sat for fixing ambiguity
                    nav.fix[satj-1,f] = 1 # float
            if n == 0: # don't use ref sat if no sat pairs
                nav.fix[sati-1,f] = 1 
            k += ns
    ix = np.resize(ix, (nb, 2))
    if nb > 0:
        tracemat(3,'refSats= ', np.array(ref), '7d')
        tracemat(3,'fixSats= ', np.array(fix), '7d')
    return ix


def restamb(nav, bias, nb):
    """ restore SD ambiguity """
    nv = 0
    xa = nav.x.copy()
    xa[0:nav.na] = nav.xa[0:nav.na]

    for m in range(uGNSS.GNSSMAX):
        for f in range(nav.nf):
            n = 0
            index = []
            for i in range(uGNSS.MAXSAT):
                sys = nav.sysprn[i+1][0]
                if sys != m or (sys not in nav.gnss_t) or nav.fix[i, f] != 2:
                    continue
                index.append(IB(i+1, f, nav.na))
                n += 1
            if n < 2:
                continue
            xa[index[0]] = nav.x[index[0]]
            for i in range(1, n):
                xa[index[i]] = xa[index[0]] - bias[nv]
                nv += 1
    return xa


def resamb_lambda(nav, sats):
    """ resolve integer ambiguity using LAMBDA method """
    nx = nav.nx
    na = nav.na
    xa = np.zeros(na)
    ix = ddidx(nav, sats)
    nav.nb_ar = nb = len(ix)
    if nb <= 0:
        trace(3, 'resamb_lambda: no valid DD\n')
        return -1, -1

    # y=D*xc, Qb=D*Qc*D', Qab=Qac*D'
    y = nav.x[ix[:, 0]] - nav.x[ix[:, 1]]
    DP = nav.P[ix[:, 0], na:nx] - nav.P[ix[:, 1], na:nx]
    Qb = DP[:, ix[:, 0] - na] - DP[:, ix[:, 1] - na]
    Qab = nav.P[0:na, ix[:, 0]] - nav.P[0:na, ix[:, 1]]
    tracemat(3,'N(0)=      ', y, '7.3f')
    tracemat(3, 'Qb*1000=   ', 1000 * np.diag(Qb[0:nb]), '7.4f')

    # MLAMBDA ILS
    b, s = mlambda(y, Qb)
    tracemat(3,'N(1)=      ', b[:,0], '7.3f')
    tracemat(3,'N(2)=      ', b[:,1], '7.3f')
    nav.ratio = s[1] / s[0]
    if s[0] <= 0.0 or nav.ratio >= nav.thresar:
        trace(3,'resamb : validation OK (nb=%d ratio=%.2f\n s=%.2f/%.2f\n' 
              % (nb, nav.ratio, s[1], s[0]))
        nav.xa = nav.x[0:na].copy()
        nav.Pa = nav.P[0:na, 0:na].copy()
        bias = b[:, 0]
        y -= b[:, 0]
        K = Qab @ inv(Qb)
        nav.xa -= K @ y
        nav.Pa -= K @ Qab.T

        # restore single diff ambiguity
        xa = restamb(nav, bias, nb)
    else:
        trace(3,'ambiguity validation failed (nb=%d ratio=%.2f\n s=%.2f/%.2f'
              % (nb, nav.ratio, s[1], s[0]))
        nb = 0
    return nb, xa


def manage_amb_LAMBDA(nav, sats, stat, posvar):
    """ resolve integer ambiguity by LAMBDA using partial fix techniques and 
    multiple attempts """
    
    # skip AR if don't meet criteria 
    if stat != gn.SOLQ_FLOAT or posvar > nav.thresar1:
        nav.ratio, nav.prev_ratio1, nav.prev_ratio2, nav.nb_ar = 0, 0, 0, 0
        return 0, []

    # if no fix on previous sample and enough sats, exclude next sat in list
    excflag = False
    if nav.prev_ratio2 < nav.thresar and nav.nb_ar >= nav.mindropsats:
        # find and count sats used last time for AR
        arsats = np.where(nav.prev_fix == 2)[0]
        excflag = 0
        if  nav.excsat_ix < len(arsats):
            excsat = arsats[nav.excsat_ix] + 1
            lockc = copy(nav.lock[excsat-1]) # save lock count
            # remove sat from AR long enough to enable hold if stays fixed
            nav.lock[excsat-1] = -nav.nb_ar
            trace(3, 'AR: exclude sat %d\n' % excsat);
            excflag = True
            nav.excsat_ix += 1
        else:
            nav.excsat_ix = 0 # exclude none and reset to beginning of list 
    
    # initial ambiguity resolution attempt, include all enabled sats
    nb, xa = resamb_lambda(nav, sats)
    ratio1 = nav.ratio
    rerun = False
    # if results are much poorer than previous epoch or dropped below AR ratio 
    # thresh, remove new sats
    trace(3, 'lambda: nb=%d r1= %.3f r2=%.3f r=%.3f\n' % ((nb, nav.prev_ratio1, nav.prev_ratio2, nav.ratio)))
    if nb >= 0 and nav.prev_ratio2 >= nav.thresar and (nav.ratio < nav.thresar 
           or (nav.ratio < nav.thresar * 1.1 and nav.ratio < nav.prev_ratio1 / 2.0)):
        trace(3, 'low ratio: check for new sat\n')
        dly = 2
        ix = np.where((nav.fix == 2) & (nav.lock == 0))
        for i,f in zip(ix[0],ix[1]):
            nav.lock[i,f] = -dly
            dly +=2
            trace(3, 'remove sat %d:%d lock=%d\n' % (i+1, f, nav.lock[i,f]))

            rerun = True
    
    # rerun if filter removed any sats
    if rerun:
        trace(3, 'rerun AR with new sat removed\n')
        nb, xa = resamb_lambda(nav, sats)
        
    # restore excluded sat if still no fix or significant increase in ar ratio 
    if excflag and nav.ratio < nav.thresar and nav.ratio < 1.5* nav.prev_ratio2:
        nav.lock[excsat-1] = lockc
        trace(3, 'AR: restore sat %d\n' % excsat)
       
    nav.prev_ratio1, nav.prev_ratio2  = ratio1, nav.ratio
    return nb, xa


def initx(nav, x0, v0, i):
    """ initialize x and P for index i """
    nav.x[i] = x0
    nav.P[i,:] = 0
    nav.P[:,i] = 0
    nav.P[i,i] = v0


def detslp_dop(rcv, nav, obs, ix):
    """ detect cycle slip with doppler measurement """
    if nav.thresdop <= 0:
        return
    # calculate doppler differences for all sats and freqs
    ns  = len(ix)
    mean_dop = ndop = 0
    dopdif = np.zeros((ns, nav.nf))
    tt = np.zeros((ns, nav.nf))
    for i, ii in enumerate(ix):
        sat = obs.sat[ii] - 1
        for f in range(nav.nf):
            if obs.L[ii,f] == 0.0 or obs.D[ii,f] == 0.0 or nav.ph[rcv,sat,f] == 0.0 \
                or nav.pt[rcv,sat,f] == None:
                continue
            tt[i,f] = timediff(obs.t, nav.pt[rcv,sat,f])
            if abs(tt[i,f]) < DTTOL:
                continue
            # calc phase difference and doppler x time (cycle)
            dph = (obs.L[ii,f] - nav.ph[rcv,sat,f]) / tt[i,f]
            dpt = -obs.D[ii,f]
            dopdif[i,f] = dph - dpt

            # if not outlier, use this to calculate mean
            if abs(dopdif[i,f]) < 3 * nav.thresdop:
                mean_dop += dopdif[i,f]
                ndop += 1
    # calc mean doppler diff, most likely due to clock error
    if ndop == 0:
        trace(4, 'detslp_dop rcv=%d: no valid doppler diffs\n' % (rcv+1))
        return # unable to calc mean doppler, usually very large clock err
    mean_dop /= ndop

    # set slip if doppler difference with mean removed exceeds threshold
    for i, ii in enumerate(ix):
        sat = obs.sat[ii] - 1
        for f in range(nav.nf):
            if dopdif[i,f] == 0.0:
                continue
            if abs(dopdif[i,f] - mean_dop) > nav.thresdop:
                nav.slip[sat,f] |= 1
                trace(3, "slip detected doppler (sat=%2d rcv=%d dL%d=%.3f off=%.3f tt=%.2f)\n"
                      % (sat+1, rcv+1, f+1, dopdif[i,f] - mean_dop, mean_dop, tt[i,f]))


def detslp_gf(nav, obsb, obsr, iu, ir):
    """ detect cycle slip with geometry-free LC """
    ns = len(iu)
    _c = rCST.CLIGHT
    for i in range(ns):
        sat = obsr.sat[iu[i]] - 1
        # skip check if slip already detected
        if (nav.slip[sat,0] & 1) or (nav.slip[sat,1] & 1):
            continue
        # calc SD geomotry free LC of phase between freq0 and freq1
        L1R = obsr.L[iu[i],0]
        L2R = obsr.L[iu[i],1]
        L1B = obsb.L[ir[i],0]
        L2B = obsb.L[ir[i],1]
        if L1R == 0.0 or L1B == 0.0 or L2R == 0 or L2B == 0:
            #trace(3, 'gf: skip sat %d, L=0\n' % sat)
            continue
        freq0 = sat2freq(sat + 1, 0, nav)
        freq1 = sat2freq(sat + 1, 1, nav)
        gf1 = ((L1R - L1B) * _c / freq0 - (L2R - L2B) * _c / freq1)
        gf0 = nav.gf[sat]    #retrieve previous gf
        nav.gf[sat] = gf1    # save current gf for next epoch
        if gf0 !=0.0 and abs(gf1 - gf0) > nav.thresslip:
            nav.slip[sat,0] |= 1
            nav.slip[sat,1] |= 1
            trace(3, "slip detected GF jump (sat=%2d L1-L2 dGF=%.3f)\n" %
                (sat + 1, gf0 - gf1))
            
def detslp_ll(nav, obs, ix, rcv):
    """ detect cycle slip from rinex file flags """
    
    # retrieve previous LLI
    LLI = nav.prev_lli[:,:,rcv]
    
    ixsat = obs.sat[ix] - 1 
    for f in range(nav.nf):
        ixL = np.where(obs.L[ix,f] != 0)[0]
        if nav.tt >= 0: # forward
            nav.slip[ixsat[ixL],f] |= (obs.lli[ix[ixL],f] & 3)
        else: # backward
            nav.slip[ixsat[ixL],f] |= (LLI[ixsat[ixL],f] & 3)
        
        # detect slip by parity unknown flag transition in LLI 
        hc_slip = np.where((obs.lli[ix[ixL],f] & 2) != 
                  (LLI[ixsat[ixL],f] & 2))[0]
        if len(hc_slip) > 0:
            nav.slip[ixsat[ixL[hc_slip]],f] |= 1
        
        # output results to trace
        ixslip = np.where((nav.slip[ixsat[ixL],f] & 1) != 0)[0]
        slipsats = ixsat[ixL[ixslip]] + 1
        if len(slipsats) > 0:
            trace(3, 'slip detected from LLI flags: f=%d, sats=%s\n' % (f, str(slipsats)))
            trace(3, '   slip=%s\n' % str(nav.slip[ixsat[ixL[ixslip]], f]))


def udpos(nav):
    """ states propagation for kalman filter """
    tt = nav.tt
    trace(3, 'udpos : tt=%.3f\n' % tt)
    
    if nav.pmode == 'static':
        return
    
    # check variance of estimated position
    posvar = np.sum(np.diag(nav.P[0:3])) / 3
    if posvar > nav.sig_p0**2:
        #reset position with large variance
        for i in range(3):
            initx(nav, nav.rr[i], nav.sig_p0**2, i)
            initx(nav, 0, nav.sig_v0**2, i + 3)
            initx(nav, 1e-6, nav.sig_v0**2, i + 6)
        trace(2, 'reset rtk position due to large variance: var=%.3f\n' % posvar)
        return

    # state transition of position/velocity/acceleration
    F = np.eye(nav.nx)
    F[0:6, 3:9] += np.eye(6) * tt
    # include accel terms if filter is converged
    if posvar < nav.thresar1:
        F[3:6, 6:9] += np.eye(3) * tt**2 / 2
    else:
        trace(3, 'pos var too high for accel term: %.4f\n' % posvar)
    # x=F*x, P=F*P*F
    nav.x = F @ nav.x
    nav.P = F @ nav.P @ F.T


    # process noise added to accel
    Q = np.zeros((3,3))
    Q[0,0] = Q[1,1] = nav.accelh**2 * abs(tt)
    Q[2,2] = nav.accelv**2 * abs(tt)    
    E = gn.xyz2enu(gn.ecef2pos(nav.x[0:3]))
    Qv = E.T @ Q @ E
    nav.P[6:9,6:9] += Qv
    
def udbias(nav, obsb, obsr, iu, ir):
    
    trace(3, 'udbias  : tt=%.3f ns=%d\n' % (nav.tt, len(iu)))
    
    # cycle slip detection from receiver flags
    detslp_ll(nav, obsb, ir, 0)
    detslp_ll(nav, obsr, iu, 1)
    # cycle slip detection by doppler and geom-free
    detslp_dop(0, nav, obsb, ir) # base
    detslp_dop(1, nav, obsr, iu) # rover
    detslp_gf(nav, obsb, obsr, iu, ir)

    # init sat and sys arrays
    ns = len(iu)
    sat = obsr.sat[iu]

    # update outage counters and reset phase-biases for sats with outage
    nav.outc += 1
    for f in range(nav.nf):
        for i in range(uGNSS.MAXSAT):
            j = IB(i+1, f, nav.na)
            if nav.outc[i,f] > nav.maxout and nav.x[j] != 0.0:
                trace(3, '  obs outage counter overflow ( sat=%d L%d: n=%d\n' 
                      % (i+1, f+1, nav.outc[i,f]))
                initx(nav, 0, 0, j)
            # TODO: set AR minlock
        # update phase bias noise and check for cycle slips and outliers
        for i in range(ns):
            j = IB(sat[i], f, nav.na)
            nav.P[j,j] += nav.prnbias**2 * abs(nav.tt)
            if (nav.slip[sat[i]-1,f] & 1) or nav.rejc[sat[i]-1,f] >= 2:
                trace(4, 'flag phase for reset: sat=%d f=%d slip=%d rejc=%d\n' % 
                      (sat[i], f, nav.slip[sat[i]-1,f], nav.rejc[sat[i]-1,f]))
                initx(nav, 0, 0, j)
            
        # estimate approximate phase-bias by delta phase - delta code
        bias = np.zeros(ns)
        offset = namb = 0
        for i in range(ns):
            freq = sat2freq(sat[i], f, nav)
            if obsr.L[iu[i], f] == 0 or obsb.L[ir[i], f] == 0 or \
                obsr.P[iu[i], f] == 0 or obsb.P[ir[i], f] == 0:
                    continue
            # calc single differences
            cp = obsr.L[iu[i], f] - obsb.L[ir[i], f]
            pr = obsr.P[iu[i], f] - obsb.P[ir[i], f]
            
            if cp == 0 or pr == 0 or freq == 0:
                continue
            # estimate bias in cycles
            bias[i] = cp - pr * freq / rCST.CLIGHT
            # offset = sum of (bias - phase-bias) for all valid sats in meters
            x = nav.x[IB(sat[i], f, nav.na)]
            if x != 0.0:
                offset += bias[i] - x
                namb += 1
                
        # correct phase-bias offset to ensure phase-code coherency
        offset = offset / namb if namb > 0 else 0
        trace(4, 'phase-code coherency adjust=%.2f, n=%d\n' % (offset, namb))
        ib1 = IB(1, f, nav.na)  # state offset for first sat
        # add offsest to non-zero states
        ix = np.where(nav.x[ib1:] != 0)[0]
        nav.x[ix+ib1] += offset
        
        # set initial states of phase-bias
        for i in range(ns):
            j = IB(sat[i], f, nav.na)
            if bias[i] == 0.0 or nav.x[j] != 0.0:
                continue
            freq = sat2freq(sat[i], f, nav)
            initx(nav, bias[i], nav.sig_n0**2, j)
            nav.outc[sat[i]-1,f] = 0
            nav.rejc[sat[i]-1,f] = 0
            nav.lock[sat[i]-1,f] = 0
            trace(3,"     sat=%3d, F=%d: init phase=%.3f\n" % (sat[i],f+1, bias[i]))


def udstate(nav, obsr, obsb, iu, ir):
    """ temporal update of states """
    trace(3, 'udstate : ns=%d\n' % len(iu))
    # temporal update of position/velocity/acceleration
    tracemat(4, 'before udstate x=', nav.x[0:9])
    udpos(nav) # updates nav.x and nav.P
    tracemat(4, 'after udstate x=', nav.x[0:9])
    # temporal update of phase-bias
    udbias(nav, obsb, obsr, iu, ir) # updates outxnav.x and nav.P

def selsat(nav, obsr, obsb, elb):
    """ select common satellite between rover and base station """
    
    trace(3, 'selsat  : nu=%d nr=%d\n' % (len(obsr.sat), len(obsb.sat)));
    # exclude satellite with missing pseudorange or low elevation
    idx_u = np.unique(np.where(obsr.P!=0)[0])
    idx_r = np.unique(np.where(obsb.P!=0)[0])
    idx_r = list(set(idx_r).intersection(np.where(elb >= nav.elmin)[0]))

    idx = np.intersect1d(obsr.sat[idx_u], obsb.sat[idx_r], return_indices=True)
    k = len(idx[0])
    iu = np.array(idx_u)[idx[1]]
    ir = np.array(idx_r)[idx[2]]
    return k, iu, ir


def holdamb(nav, xa):
    """ hold integer ambiguity """
    nb = nav.nx-nav.na
    v = np.zeros(nb)
    H = np.zeros((nb, nav.nx))
    nv = 0
    for m in range(uGNSS.GNSSMAX):
        for f in range(nav.nf):
            n = 0
            index = []
            for i in range(uGNSS.MAXSAT):
                sys = nav.sysprn[i+1][0]
                if sys != m or nav.fix[i, f] != 2:
                    continue
                index.append(IB(i+1, f, nav.na))
                n += 1
                nav.fix[i, f] = 3  # hold
            # constraint to fixed ambiguity
            for i in range(1, n):
                v[nv] = (xa[index[0]] - xa[index[i]]) - \
                    (nav.x[index[0]] - nav.x[index[i]])
                H[nv, index[0]] = 1.0
                H[nv, index[i]] = -1.0
                nv += 1
    if nv > 0:
        R = np.eye(nv) * nav.var_holdamb
        # update states with constraints
        nav.x, nav.P, _ = gn.filter(nav.x, nav.P, H[0:nv, :], v[0:nv], R)
    return 0
        
        
def relpos(nav, obsr, obsb, sol):
    """ relative positioning for PPK """
    
    # time diff between rover and base
    nav.dt = timediff(obsr.t, obsb.t)
    ep = gn.time2epoch(obsr.t)
    trace(1,"\n---------------------------------------------------------\n")
    trace(1, "relpos: nu=%d nr=%d\n" % (len(obsr.sat), len(obsb.sat)))
    trace(1, '        teph= %04d/%02d/%02d %02d:%02d:%06.3f\n' %           
          (ep[0], ep[1], ep[2], ep[3], ep[4], ep[5])) 
    trace(1,"---------------------------------------------------------\n")
    if abs(nav.dt) > nav.maxage:
        trace(3, 'Age of differential too large: %.2f\n' % nav.dt)
        return
    
    # clear valid sat status
    nav.vsat[:,:] = 0

    # compute satellite positions, velocities and clocks
    rs, var, dts, svh = satposs(obsr, nav)
    rsb, varb, dtsb, svhb = satposs(obsb, nav)

    # undifferenced residuals for base
    trace(3, 'base station:\n')
    yr, er, elr = zdres(nav, obsb, rsb, dtsb, svhb, varb, nav.rb, 0)
    if nav.interp_base:
        # get residuals for previous base station obs
        yr0, _, _ = zdres(nav, nav.obsb, nav.rsb, nav.dtsb, nav.svhb, nav.varb,
                          nav.rb, 0)
        # time-interpolation of base residuals
        yr, nav.dt = intpres(obsr.t, nav, yr0, yr, nav.obsb, obsb)
    
    # find common sats between base and rover
    ns, iu, ir = selsat(nav, obsr, obsb, elr)
    if ns <= 0:
        trace(3, 'no common sats: %d\n' % ns)
        return
    
    # kalman filter time propagation
    udstate(nav, obsr, obsb, iu, ir)

    # undifferenced residuals for rover
    trace(3, 'rover: dt=%.3f\n' % nav.dt)
    yu, eu, el = zdres(nav, obsr, rs, dts, svh, var, nav.x[0:3], 1)
    # decode stdevs from receiver
    rn.rcvstds(nav, obsr)

    # remove non-common residuals
    yr, er = yr[ir,:], er[ir,:]
    yu, eu = yu[iu,:], eu[iu,:]
    sats = obsr.sat[iu]
    els = nav.el[sats-1] = el[iu]
    
    # double differenced residuals
    v, H, R = ddres(nav, nav.x, nav.P, yr, er, yu, eu, sats, els, nav.dt, obsr)
    
    if len(v) < 4:
        trace(3, 'not enough double-differenced residual\n')
        stat = gn.SOLQ_NONE
    else:
        stat = gn.SOLQ_FLOAT
    
    if stat != gn.SOLQ_NONE:
        # kalman filter measurement update
        tracemat(3, 'before filter x=', nav.x[0:9])
        xp, Pp = gn.filter(nav.x, nav.P, H, v, R)
        tracemat(3, 'after filter x=', xp[0:9])
        posvar = np.sum(np.diag(Pp[0:3])) / 3
        trace(3,"posvar=%.6f \n" % posvar)

        # check validity of solution 
        # non-differencial residual for rover after measurement update
        yu, eu, _ = zdres(nav, obsr, rs, dts, svh, var, xp[0:3], 1)
        yu, eu = yu[iu,:], eu[iu,:]
        # residual for float solution
        v, H, R = ddres(nav, xp, Pp, yr, er, yu, eu, sats, els, nav.dt, obsr)
        # validation of float solution, always returns 1, msg to trace file if large residual
        valpos(nav, v, R)
        
        # save results of kalman filter update
        nav.x = xp.copy()
        nav.P = Pp.copy()
        
    # update sat status
    for f in range(nav.nf):
        ix = np.where(nav.vsat[:,f] > 0)[0]
        nav.outc[ix,f] = 0
        if f == 0:
            nav.ns = len(ix)
            
    # resolve integer ambiguities
    if nav.armode > 0:
        nb, xa = manage_amb_LAMBDA(nav, sats, stat, posvar)
        if nb > 0:
            yu, eu, el = zdres(nav, obsr, rs, dts, svh, var, xa[0:3], 1)
            yu, eu, el = yu[iu, :], eu[iu, :], el[iu]
            v, H, R = ddres(nav, xa, nav.P, yr, er, yu, eu, sats, el, nav.dt, obsr)
            if valpos(nav, v, R):
                if nav.armode == 3:
                    holdamb(nav, xa)
                stat = gn.SOLQ_FIX
        
    # save solution unless none
    if stat != gn.SOLQ_NONE:
        if stat == gn.SOLQ_FIX: 
            sol.rr = nav.xa[0:6]
            sol.qr = nav.Pa[0:3,0:3]
            sol.qv = nav.Pa[3:6,3:6]
        else: # SOLQ_FLOAT
            sol.rr = nav.x[0:6]
            sol.qr = nav.P[0:3,0:3]
            sol.qv = nav.P[3:6,3:6]
        sol.stat = stat
        sol.ratio = nav.ratio
        sol.age = nav.dt
        nav.sol.append(sol)
        nav.rr = sol.rr[0:3]
        tracemat(3, 'sol_rr= ', sol.rr, '15.3f')
                
    # save phases and times for cycle slip detection
    for i, sat in enumerate(sats):
        for f in range(nav.nf):
            if obsb.L[ir[i],f] != 0:
                nav.pt[0,sat-1,f] = obsb.t
                nav.ph[0,sat-1,f] = obsb.L[ir[i],f]
            if obsr.L[iu[i],f] != 0:
                nav.pt[1,sat-1,f] = obsr.t
                nav.ph[1,sat-1,f] = obsr.L[iu[i],f]
    # save current LLI and fix status
    nav.prev_lli[:,:], nav.slip[:,:] = 0, 0
    for f in range(nav.nf):
        nav.prev_lli[obsb.sat-1,f,0] = obsb.lli[:,f]
        nav.prev_lli[obsr.sat-1,f,1] = obsr.lli[:,f]
    if nav.armode > 0:
        nav.prev_fix = copy(nav.fix)
        # update lock counts for sats used in fix and disabled sats (lock < 0)
        for f in range(nav.nf):
            ix = np.where(((nav.nb_ar > 0) & (nav.fix[:,f] == 2)) | (nav.lock[:,f] < 0))[0]
            nav.lock[ix,f] += 1

            
def rtkpos(nav, rov, base, dir):
    """ relative positioning for PPK """
    trace(3, 'rtkpos: start solution, dir=%d\n' % dir)
    n = 0
    # loop through all epochs
    while True:
        if n== 0:
            # first epoch
            obsr, obsb = rn.first_obs(nav, rov, base, dir)
            t = obsr.t
            # force initial rover position (for align to RTKLIB)
            if dir == 1:
                if cfg.rr_f[0] != 0:
                    nav.x[0:6] = cfg.rr_f
            else:
                if cfg.rr_b[0] != 0:
                    nav.x[0:6] = cfg.rr_b
            nav.x[6:9] = 1E-6  # match RTKLIB
        else:
            # get next rover obs and next base obs if required
            if len(nav.sol) > 0:
                t = nav.sol[-1].t # previous epoch
            obsr, obsb = rn.next_obs(nav, rov, base, dir)
        if obsr == []:
            break
        # single precision solution, used to update solution time
        if nav.use_sing_pos:
            sol = pntpos(obsr, nav)
        else:
            sol = gn.Sol()
        if sol.t.time == 0:
            sol.t = obsr.t
        nav.tt = timediff(sol.t, t) # timediff from previous epoch
        # relative solution
        relpos(nav, obsr, obsb, sol)
        ep = gn.time2epoch(sol.t)
        stdout.write('\r   %2d/%2d/%4d %02d:%02d:%05.2f: %d' % (ep[1], ep[2], ep[0],
                ep[3], ep[4], ep[5], sol.stat))
        n += 1
        if nav.maxepoch != None and n > nav.maxepoch:
            break
    trace(3, 'rtkpos: end solution\n')
                


