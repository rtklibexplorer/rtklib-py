"""
module for PPK positioning

Copyright (c) 2021 Rui Hirokawa (from CSSRLIB)
Copyright (c) 2022 Tim Everett
"""

import numpy as np
from numpy.linalg import inv, norm
import rtkcmn as gn
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
    nav.nf = cfg.nf
    nav.pmode = cfg.pmode
    nav.filtertype = cfg.filtertype
    # add rover vel and accel states for kinematic solution
    nav.na = nav.nq = 3 if nav.pmode == 'static' else 9
    nav.ratio = 0
    nav.thresar = [cfg.thresar]
    nav.thresar1 = cfg.thresar1
    nav.nx = nav.na + uGNSS.MAXSAT * nav.nf
    nav.x = np.zeros(nav.nx)
    nav.P = np.zeros((nav.nx, nav.nx))
    nav.xa = np.zeros(nav.na)
    nav.Pa = np.zeros((nav.na, nav.na))
    nav.el = np.zeros(uGNSS.MAXSAT)
    nav.gf = np.zeros(uGNSS.MAXSAT)
    nav.ph = np.zeros((2, uGNSS.MAXSAT, nav.nf))
    nav.pt = np.empty((2, uGNSS.MAXSAT, nav.nf), dtype=object)
    nav.nfix = nav.neb = nav.db = nav.tt = 0
    
    nav.rb = cfg.rb

    # parameter for RTK/PPK    
    nav.eratio = cfg.eratio
    nav.err = np.array(cfg.err) 
    nav.use_sing_pos = cfg.use_sing_pos
   
    nav.armode = cfg.armode
    nav.elmaskar = np.deg2rad(cfg.elmaskar)
    nav.gnss_t = cfg.gnss_t
    nav.maxinno = cfg.maxinno
    nav.var_holdamb = cfg.var_holdamb
    nav.thresdop = cfg.thresdop
    nav.thresslip = cfg.thresslip
    nav.maxage = cfg.maxage
    nav.accelh = cfg.accelh
    nav.accelv = cfg.accelv
    nav.prnbias = cfg.prnbias
    
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
    ix0 = cfg.freq_ix0
    ix1 = cfg.freq_ix1
    freq0 = {}; freq1 = {}
    for k in ix0.keys():
        freq0[k] = cfg.freq[ix0[k]]
    for k in ix1.keys():
        freq1[k] = cfg.freq[ix1[k]]
    nav.obs_idx = [ix0, ix1]
    nav.obs_freq = [freq0, freq1]

    # sat index
    nav.sysprn = {}
    for i in range(uGNSS.MAXSAT):
        nav.sysprn[i+1] = gn.sat2prn(i+1)
        
    return nav

def zdres_sat(nav, obs, r, rtype, dant, ix):
    sys = nav.sysprn[obs.sat[ix]][0]
    _c = gn.rCST.CLIGHT
    nf = nav.nf
    y = np.zeros(nf * 2)
    for f in range(nf):
        j = nav.obs_idx[f][sys]
        # TODO: check SNR mask
        # residuals = observable - estimated range (phase and code)
        y[f] = obs.L[ix,f] * _c / nav.freq[j] - r - dant[f] if obs.L[ix,f] else 0
        y[f+nf] = obs.P[ix,f] - r - dant[f] if obs.P[ix,f] else 0
        #if obs.L[ix,f] == 0 or obs.P[ix,f] == 0:
        #    continue
        #trace(3, 'zdres_sat: %d: %.6f %.6f %.6f %.6f\n' % (obs.sat[ix],obs.L[ix,f],
        #    obs.P[ix,f],r,dant[f]))
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
    _c = gn.rCST.CLIGHT
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
        if svh[i] > 0 or var[i] > MAX_VAR_EPH or obs.sat[i] in nav.excl_sat:
            trace(3, 'exclude sat %d: svh=%d ura=%.2f\n' % (obs.sat[i], svh[i],
                        np.sqrt(var[i])))
            continue
        # compute geometric-range and azimuth/elevation angle
        r, e[i,:] = gn.geodist(rs[i,0:3], rr_)
        _, el[i] = gn.satazel(pos, e[i,:])
        if el[i] < nav.elmin:
            continue
        # TODO: check for excluded sat
        
        # adjust range for satellite clock-bias
        r += -_c * dts[i]
        # adjust range for troposphere delay model (hydrostatic)
        trophs, tropw, _ = gn.tropmodel(obs.t, pos, np.deg2rad(90.0), 0.0)
        zhd = trophs + tropw
        mapfh, _ = gn.tropmapf(obs.t, pos, el[i])
        r += mapfh * zhd
        # calc receiver antenna phase center correction
        dant = gn.antmodel(nav, el[i], nav.nf, rtype)
        trace(3,'sat=%d %.6f %.6f %.6f %.6f\n' % (obs.sat[i],r,_c*dts[i],zhd,mapfh))
        # calc undifferenced phase/code residual for satellite
        y[i] = zdres_sat(nav, obs, r, rtype, dant, i)
        trace(4, 'sat=%d: y=%.3f %.3f %.3f %.3f trop=%.3f\n' % (obs.sat[i], 
            y[i,0], y[i,2],y[i,1], y[i,3], mapfh * zhd))
    
    for i in ix:
        if obs.L[i,0] != 0 and rtype == 0:
            trace(3, 'sat=%2d %13.3f %13.3f %13.3f %13.10f %5.1f\n' %
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
        for i in range(nb[b]):
            for j in range(nb[b]):
                R[k+i, k+j] = Ri[k+i]
                if i == j:
                    R[k+i, k+j] += Rj[k+i]
        k += nb[b]
    return R


def sysidx(satlist, sys_ref):
    """ return index of satellites with sys=sys_ref """
    idx = []
    for k, sat in enumerate(satlist):
        sys, _ = gn.sat2prn(sat)
        if sys == sys_ref:
            idx.append(k)
    return idx


def IB(s, f, na=3):
    """ return index of phase ambguity """
    idx = na + uGNSS.MAXSAT * f + s - 1
    return idx


def varerr(nav, el, f, dt, rcvstd):
    """ variation of measurement """
    code = 1 * (f >= nav.nf) # 0 = phase, 1 = code
    freq = f % nav.nf
    s_el = np.sin(el)
    #if s_el <= 0.0: return 0.0
    fact = nav.eratio[freq] if code else 1
    a, b = fact * nav.err[1:3]
    c = fact * 0  # nav.err[4]*bl/1E4  # TODO: add baseline term
    d = gn.rCST.CLIGHT * nav.err[5] * dt # clock term
    var = 2.0 * (a**2 + (b / s_el)**2 + c**2) + d**2
    # TODO: add SNR term
    # add scaled stdevs from receiver
    if nav.err[3] > 0:
        var += (nav.err[3] * rcvstd)**2
    return var


def ddres(nav, x, yr, er, yu, eu, sat, el, dt, obsr):
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
    _c = gn.rCST.CLIGHT
    nf = nav.nf
    ns = len(el)
    ny = ns * nf * 2  # phase and code
    nb = np.zeros(2 * len(nav.gnss_t) * nf, dtype=int)
    Ri = np.zeros(ny)
    Rj = np.zeros(ny)
    H = np.zeros((nav.nx, ny))
    trace(3,"\n\nddres   : dt=%.2f nx=%d ns=%d\n" % (dt, nav.nx, ns))

    nv = b = 0
    v = np.zeros(ny)
    # step through sat systems
    for sys in nav.gnss_t:
        # step through phases/codes
        for f in range(0, nf*2):
            frq = f % nf
            code = 1 * (f >= nf)
            k = nav.obs_idx[frq][sys]  # index of frq
            lami = _c / nav.freq[k]
            idx = sysidx(sat, sys) # find sats in sys
            # remove sats with missing base or rover residuals
            nozero = np.where((yr[:,f] != 0) & (yu[:,f] != 0))[0]
            idx = np.intersect1d(idx, nozero)
            if len(idx) == 0: 
                continue  # no common sats
            # find reference satellite with highest elevation, set to i

            # find sat with max el and no slip for reference
            noslip=idx[np.where(nav.slip[sat[idx]-1,frq]==0)[0]]
            idx_noslip = np.intersect1d(idx, noslip)
            if len(idx_noslip) > 0:
                i = idx[np.argmax(el[idx_noslip])]
            else: # use sat with slip if no sats without slip
                i = idx[np.argmax(el[idx])]
            # calculate double differences of residuals (code/phase) for each sat
            for j in idx: # loop through sats
                if i == j: continue  # skip ref sat
                if yu[i,f] == 0 or yr[i,f] == 0 or yu[j,f] == 0 or yr[j,f] == 0:
                    continue
                #  double-differenced measurements from 2 receivers and 2 sats in meters 
                v[nv] = (yu[i,f] - yr[i,f]) - (yu[j,f] - yr[j,f])
                # partial derivatives by rover position, combine unit vectors from two sats
                H[0:3, nv] = -eu[i,:] + er[j,:]
                
                ii = IB(sat[i], frq, nav.na)
                jj = IB(sat[j], frq, nav.na)
                if not code:  # carrier phase
                    # adjust phase residual by double-differenced phase-bias term
                    v[nv] -= lami * (x[ii] - x[jj])
                    H[ii, nv] = lami
                    H[jj, nv] = -lami
                
                # if residual too large, flag as outlier
                thres = nav.maxinno
                # use larger thresh for code or just initialized phase
                if code or nav.P[ii,ii] == nav.sig_n0**2 or \
                        nav.P[jj,jj] == nav.sig_n0**2:
                    thres *= nav.eratio[frq] 
                if abs(v[nv]) > thres:
                    nav.vsat[sat[j]-1,frq] = 0
                    nav.rejc[sat[j]-1,frq] += 1
                    trace(3,"outlier rejected: (sat=%3d-%3d %s%d v=%13.3f)\n" 
                          % (sat[i], sat[j], 'LP'[code],frq+1, v[nv]))
                    continue 
                # single-differenced measurement error variances (m)
                si = sat[i] - 1; sj = sat[j] - 1
                Ri[nv] = varerr(nav, el[i], f, dt, nav.rcvstd[si,f])
                Rj[nv] = varerr(nav, el[j], f, dt, nav.rcvstd[sj,f])
                nav.vsat[si,frq] = 1
                nav.vsat[sj,frq] = 1
                #trace(3,'sys=%d f=%d,i=%d,j=%d\n' % (sys,f,i,j))
                trace(3,"sat=%3d-%3d %s%d v=%13.3f R=%9.6f %9.6f x=%13.3f\n" %
                      (sat[i], sat[j], 'LP'[code], frq+1, v[nv], Ri[nv], Rj[nv], x[jj]))
                nv += 1
                nb[b] += 1
            b += 1
    R = ddcov(nb, b, Ri, Rj, nv)

    return v[:nv], H[:,:nv], R


def valpos(nav, v, R, thres=4.0):
    """ post-file residual test """
    nv = len(v)
    fact = thres**2
    for i in range(nv):
        if v[i]**2 > fact * R[i, i]:
            trace(3, 'large residual (i=%d  v=%.3f)' % (i, v[i]))
    return True


def ddidx(nav, sats):
    """ index for single to double-difference transformation matrix (D') """
    nb, fix, ref = 0, [], []
    ns = uGNSS.MAXSAT
    #na = nav.na
    ix = np.zeros((ns, 2), dtype=int)
    # clear fix flag for all sats (1=float, 2=fix)
    nav.fix = np.zeros((ns, nav.nf), dtype=int)
    for m in range(uGNSS.GNSSMAX):
        k = nav.na
        for f in range(nav.nf): # step through freqs
            for i in range(k, k + ns):
                sat_i = i - k + 1
                sys = nav.sysprn[sat_i][0]
                if (sys != m): # or sys not in nav.gnss_t:
                    continue
                if sat_i not in sats or nav.x[i] == 0.0 \
                   or nav.vsat[sat_i-1, f] <= 0:
                    continue
                if nav.el[sat_i-1] >= nav.elmaskar:
                    # set sat to use for fixing ambiguity if meets criteria
                    nav.fix[sat_i-1, f] = 2 # fix
                    break
                else:
                    nav.fix[sat_i-1, f] = 1 # float
            for j in range(k, k + ns):
                sat_j = j - k + 1
                sys = nav.sysprn[sat_j][0]
                if (sys != m): # or sys not in nav.gnss_t:
                    continue
                if i == j or sat_j not in sats or nav.x[j] == 0.0 \
                   or nav.vsat[sat_j-1, f] <= 0:
                    continue
                if nav.el[sat_j-1] >= nav.elmaskar:
                    # set D coeffs to subtract sat j from sat i
                    ix[nb, :] = [i, j]
                    ref.append(sat_i)
                    fix.append(sat_j)
                    # inc # of sats used for fix
                    nb += 1
                    nav.fix[sat_j-1, f] = 2 # fix
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
    nb = len(ix)
    if nb <= 0:
        print("no valid DD")
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
    ratio = s[1] / s[0]
    if s[0] <= 0.0 or ratio >= nav.thresar[0]:
        trace(3,'resamb : validation OK (nb=%d ratio=%.2f\n s=%.2f/%.2f' 
              % (nb, ratio, s[1], s[0]))
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
              % (nb, ratio, s[1], s[0]))
        nb = 0
    return nb, xa


def initx(nav, x0, v0, i):
    """ initialize x and P for index i """
    nav.x[i] = x0
    for j in range(nav.nx):
        nav.P[j, i] = nav.P[i, j] = v0 if i == j else 0

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
            tt[i,f] = abs(gn.timediff(obs.t, nav.pt[rcv,sat,f]))
            if tt[i,f] < gn.DTTOL:
                continue
            # calc phase difference and doppler x time (cycle)
            dph = obs.L[ii,f] - nav.ph[rcv,sat,f]
            dpt = -obs.D[ii,f] * tt[i,f]
            dopdif[i,f] = (dph-dpt) / tt[i,f]

            # if not outlier, use this to calculate mean
            if abs(dopdif[i,f]) < 3 * nav.thresdop:
                mean_dop += dopdif[i,f]
                ndop += 1
    # calc mean doppler diff, most likely due to clock error
    if ndop == 0:
        trace(3, 'detslp_dop rcv=%d: no valid doppler diffs\n' % (rcv+1))
        return # unable to calc mean doppler, usually very large clock err
    mean_dop /= ndop;

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
    _c = gn.rCST.CLIGHT
    for i in range(ns):
        sat = obsr.sat[iu[i]] - 1
        sys = nav.sysprn[sat][0]
        # skip check if slip already detected
        if (nav.slip[sat,0] & 1) or (nav.slip[sat,1] & 1):
            #trace(3, 'gf: skip sat %d, LLI=1\n' % sat)
            continue
        # calc SD geomotry free LC of phase between freq0 and freq1
        j0 = nav.obs_idx[0][sys]
        j1 = nav.obs_idx[1][sys]
        L1R = obsr.L[iu[i],0]
        L2R = obsr.L[iu[i],1]
        L1B = obsb.L[ir[i],0]
        L2B = obsb.L[ir[i],1]
        if L1R == 0.0 or L1B == 0.0 or L2R == 0 or L2B == 0:
            #trace(3, 'gf: skip sat %d, L=0\n' % sat)
            continue
        gf1 = ((L1R - L1B) * _c / nav.freq[j0] - (L2R - L2B) * _c / nav.freq[j1])
        gf0 = nav.gf[sat]    #retrieve previous gf
        nav.gf[sat] = gf1    # save current gf for next epoch
        if gf0 !=0.0 and abs(gf1 - gf0) > nav.thresslip:
            nav.slip[sat,0] |= 1
            nav.slip[sat,1] |= 1
            trace(3, "slip detected GF jump (sat=%2d L1=%.3f L2=%.3f dGF=%.3f)\n" %
                (sat + 1, gf0, gf1, gf0 - gf1))
            
def detslp_ll(nav, obs, ix):
    for i in ix:
        sat = obs.sat[i] - 1
        for f in range(nav.nf):
            if obs.L[i,f] != 0:
                nav.slip[sat,f] |= (obs.lli[i,f] & 1)
            # TODO: add half-cycle ambiguity handling


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
        trace(3, 'ignore high accel: %.4f\n' % posvar)
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
    detslp_ll(nav, obsb, ir)
    detslp_ll(nav, obsr, iu)
    # cycle slip detection by doppler and geom-free
    detslp_dop(0, nav, obsb, ir) # base
    detslp_dop(1, nav, obsr, iu) # rover
    detslp_gf(nav, obsb, obsr, iu, ir)

    # init sat and sys arrays
    ns = len(iu)
    sat = obsr.sat[iu]

    # reset phase-biases for sats with outage
    for f in range(nav.nf):
        for i in range(uGNSS.MAXSAT):
            nav.outc[i, f] += 1
            j = IB(i+1, f, nav.na)
            if nav.outc[i,f] > nav.maxout and nav.x[j] != 0.0:
                trace(3, '  obs outage counter overflow ( sat=%d L%d: n=%d\n' 
                      % (i+1, f+1, nav.outc[i,f]))
                initx(nav, 0, 0, j)
                
        # update phase bias noise and check for cycle slips and outliers
        for i in range(ns):
            j = IB(sat[i], f, nav.na)
            nav.P[j,j] += nav.prnbias**2 * abs(nav.tt)
            if (nav.slip[sat[i]-1,f] & 1) or nav.rejc[sat[i]-1,f] > 1:
                initx(nav, 0, 0, j)
                nav.rejc[sat[i]-1,f] = 0
                nav.slip[sat[i]-1,f] = 0
            
        # estimate approximate phase-bias by delta phase - delta code
        bias = np.zeros(ns)
        offset = na = 0
        for i in range(ns):
            sys = nav.sysprn[sat[i]-1][0]
            freq = nav.obs_freq[f][sys]
            if obsr.L[iu[i], f] == 0 or obsb.L[ir[i], f] == 0 or \
                obsr.P[iu[i], f] == 0 or obsb.P[ir[i], f] == 0:
                    continue
            # calc single differences
            cp = obsr.L[iu[i], f] - obsb.L[ir[i], f]
            pr = obsr.P[iu[i], f] - obsb.P[ir[i], f]
            
            if cp == 0 or pr == 0 or freq == 0:
                continue
			# translate cycles diff to meters and subtract pseudorange diff
            bias[i] = cp * (gn.rCST.CLIGHT / freq) - pr
            # offset = sum of (bias - phase-bias) for all valid sats in meters
            x = nav.x[IB(sat[i], f, nav.na)]
            if x != 0.0:
                dbias = bias[i] - x * (gn.rCST.CLIGHT / freq)
                offset += dbias
                na += 1
                trace(4,'     sat:%d: bias=%.2f\n' % (sat[i], dbias))
                
        # adjust phase-code coherency
        nav.db = offset / na if na > 0 else 0
        trace(4, 'phase-code coherency adjust=%.2f, n=%d\n' % (nav.db, na))
        #for i in range(uGNSS.MAXSAT):
        #    if nav.x[IB(i+1, f, nav.na)] != 0.0:
        #        nav.x[IB(i+1, f, nav.na)] += db
        
        # initialize ambiguities
        for i in range(ns):
            j = IB(sat[i], f, nav.na)
            if bias[i] == 0.0 or nav.x[j] != 0.0:
                continue
            sys = nav.sysprn[sat[i]-1][0]
            freq = nav.obs_freq[f][sys]
            adjbias = (bias[i] - nav.db) * freq / gn.rCST.CLIGHT
            initx(nav, adjbias, nav.sig_n0**2, j)
            trace(3,"     sat=%3d, F=%d: init phase=%.3f\n" % (sat[i],f+1, adjbias))


def udstate(nav, obsr, obsb, iu, ir):
    """ temporal update of states """
    trace(3, 'udstate : ns=%d\n' % len(iu))
    # temporal update of position/velocity/acceleration
    tracemat(3, 'before udstate x=', nav.x[0:9])
    udpos(nav) # updates nav.x and nav.P
    tracemat(3, 'after udstate x=', nav.x[0:9])
    # temporal update of phase-bias
    udbias(nav, obsb, obsr, iu, ir) # updates outxnav.x and nav.P

def selsat(nav, obsr, obsb, elb):
    """ select common satellite between rover and base station """
    # exclude satellite with missing observation and cycle slip for rover
    idx_u = []
    for k, sat in enumerate(obsr.sat):
        if (obsr.P[k, 0] == 0.0 and obsr.P[k, 1] == 0.0):
            continue
        idx_u.append(k)

    # exclude satellite with missing observation and cycle slip for base
    idx_r = []
    for k, sat in enumerate(obsb.sat):
        if (obsb.P[k, 0] == 0.0 and obsb.P[k, 1] == 0.0) or elb[k] < nav.elmin:
            continue
        idx_r.append(k)

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
                v[nv] = (xa[index[0]]-xa[index[i]]) - \
                    (nav.x[index[0]]-nav.x[index[i]])
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
    nav.dt = gn.timediff(obsr.t, obsb.t) 
    trace(3,"\nrelpos  : nx=%d, dt=%.3f, nu=%d nr=%d\n" % (nav.nx, nav.dt,
            len(obsr.sat), len(obsb.sat)))
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
    
    # find common sats between base and rover
    ns, iu, ir = selsat(nav, obsr, obsb, elr)
    if ns <= 0:
        trace(3, 'no common sats: %d\n' % ns)
        return
    
    # kalman filter time propagation
    udstate(nav, obsr, obsb, iu, ir)

    # undifferenced residuals for rover
    trace(3, 'rover:\n')
    yu, eu, el = zdres(nav, obsr, rs, dts, svh, var, nav.x[0:3], 1)
    # decode stdevs from receiver
    rn.rcvstds(nav, obsr)

    # remove non-common residuals
    yr, er = yr[ir,:], er[ir,:]
    yu, eu = yu[iu,:], eu[iu,:]
    sats = obsr.sat[iu]
    els = nav.el[sats-1] = el[iu]
    
    # double differenced residuals
    v, H, R = ddres(nav, nav.x, yr, er, yu, eu, sats, els, nav.dt, obsr)
    
    if len(v) < 4:
        trace(3, 'not enough double-differenced residual\n')
        stat = gn.SOLQ_NONE
    else:
        stat = gn.SOLQ_FLOAT
    
    if stat != gn.SOLQ_NONE:
        # kalman filter measurement update
        tracemat(3, 'before filter x=', nav.x[0:9])
        #tracemat(3, 'before filter v=', v)
        xp, Pp = gn.filter(nav.x, nav.P, H, v, R)
        tracemat(3, 'after filter x=', xp[0:9])
        posvar = np.sum(np.diag(Pp[0:3])) / 3
        trace(3,"posvar=%.6f \n" % posvar)

        # check validity of solution (optional, doesn't affect result)
        # non-differencial residual for rover after measurement update
        # yu, eu, _ = zdres(nav, obsr, rs, dts, svh, xp[0:3], 1)
        # yu, eu = yu[iu,:], eu[iu,:]
        # # residual for float solution
        # v, H, R = ddres(nav, xp, yr, er, yu, eu, sat, el, dt, obsr)
        #valpos(nav, v, R):
        
        # save results of kalman filter update
        nav.x = xp
        nav.P = Pp
        
    # update sat status
    for f in range(nav.nf):
        ix = np.where(nav.vsat[:,f] > 0)[0]
        nav.outc[ix,f] = 0
        if f == 0:
            nav.ns = len(ix)

    # ambiguity resolution
    if stat == gn.SOLQ_FLOAT and nav.armode > 0 and posvar < nav.thresar1:
        nb, xa = resamb_lambda(nav, sats)
        if nb > 0:
            yu, eu, _ = zdres(nav, obsr, rs, dts, var, svh, xa[0:3], 1)
            yu, eu = yu[iu, :], eu[iu, :]
            v, H, R = ddres(nav, xa, yr, er, yu, eu, sats, el, nav.dt, obsr)
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
                
    # save phases and times for cycle slip detection
    for i, sat in enumerate(sats):
        for f in range(nav.nf):
            if obsb.L[ir[i],f] != 0:
                nav.pt[0,sat-1,f] = obsb.t
                nav.ph[0,sat-1,f] = obsb.L[ir[i],f]
            if obsr.L[iu[i],f] != 0:
                nav.pt[1,sat-1,f] = obsr.t
                nav.ph[1,sat-1,f] = obsr.L[iu[i],f]
           
def rtkpos(nav, rov, base, dir):
    """ relative positioning for PPK """

    n = 0
    while True:
        print(n)
        if n== 0:
            obsr, obsb = rn.first_obs(nav, rov, base, dir)
            t = obsr.t
            # force initial rover position (for align to RTKLIB)
            if dir == 1:
                if cfg.rr_f[0] != 0:
                    nav.x[0:6] = cfg.rr_f
            else:
                if cfg.rr_b[0] != 0:
                    nav.x[0:6] = cfg.rr_b
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
        nav.tt = gn.timediff(sol.t, t)
        # relative solution
        relpos(nav, obsr, obsb, sol)
        n += 1
        if nav.maxepoch != None and n > nav.maxepoch:
            break
                


