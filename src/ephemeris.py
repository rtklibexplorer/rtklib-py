"""
module for ephemeris processing

Copyright (c) 2021 Rui Hirokawa (from CSSRLIB)
Copyright (c) 2022 Tim Everett
"""

import numpy as np
from rtkcmn import uGNSS, rCST, timediff, timeadd, vnorm, time2epoch
from rtkcmn import sat2prn, trace

# ephemeris parameters
MAX_ITER_KEPLER = 30
RTOL_KEPLER = 1e-13
TSTEP = 120 #60.0  # time step for Glonass orbital calcs
ERREPH_GLO = 5.0


def seleph(nav, t, sat):
    """ select ephemeric for sat, assumes ephemeris is sorted by sat, then time """
    dt_p = 1e10 # timediff(t, nav.eph[nav.eph_index[sat]].toe)
    eph = None
    sys = sat2prn(sat)[0]
    i_p = 0
    if sys != uGNSS.GLO:
        # start with previous index for this sat
        for i, eph_ in enumerate(nav.eph[nav.eph_index[sat]:]):
            if eph_.sat != sat:
                continue
            # bit 8 set=E5a, bit 9 set=E5b
            if sys == uGNSS.GAL:
                # TODO: abstract hard coded freq
                if nav.obs_idx[1][uGNSS.GAL] == 2 and (eph_.code >> 8) & 1 == 0:
                    continue
                elif nav.obs_idx[1][uGNSS.GAL] == 3 and (eph_.code >> 9) & 1 == 0:
                    continue
            dt = timediff(t, eph_.toe)
            if abs(dt) <= dt_p:
                dt_p = abs(dt)
                i_p = i
                eph = eph_
            else:
                break
    else: # GLONASS
        # start with previous index for this sat
        for i, eph_ in enumerate(nav.geph[nav.eph_index[sat]:]):
            if eph_.sat != sat:
                continue
            dt = timediff(t, eph_.toe)
            if abs(dt) <= dt_p:
                dt_p = abs(dt)
                i_p = i
                eph = eph_
            else:
                break
    trace(4, 'seleph: sat=%d dt=%.0f\n' % (sat,dt_p))
    nav.eph_index[sat] = max(nav.eph_index[sat] + i_p - 1, 0) # save index for next time
    return eph


def dtadjust(t1, t2, tw=604800):
    """ calculate delta time considering week-rollover """
    dt = timediff(t1, t2)
    if dt > tw:
        dt -= tw
    elif dt < -tw:
        dt += tw
    return dt

def sva2ura(sys, sva):
    """ variance by ura ephemeris """
    ura_nominal = [2.0, 2.8, 4.0, 5.76, 8.0, 11.3, 16.0, 32.0, 64.0, 128.0,
                   256.0, 512.0, 1024.0, 2048.0,4096.0, 8192.0]
    if sys == uGNSS.GAL:  #galileo sisa 
       if sva < 0 or sva > 6: return 500**2
       return sva**2
    else: # gps ura
        if sva < 0 or sva > 15: return 500**2
        return ura_nominal[int(sva) + 1]

def eph2pos(t, eph):
    """ broadcast ephemeris to satellite position and clock bias -------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd) """
    tk = dtadjust(t, eph.toe)
    sys, _ = sat2prn(eph.sat)
    if sys == uGNSS.GAL:
        mu = rCST.MU_GAL
        omge = rCST.OMGE_GAL
    else:  # GPS,QZS
        mu = rCST.MU_GPS
        omge = rCST.OMGE

    M = eph.M0 + (np.sqrt(mu / eph.A**3) + eph.deln) * tk
    E, Ek = M, 0
    for _ in range(MAX_ITER_KEPLER):
        if abs(E - Ek) < RTOL_KEPLER:
            break
        Ek = E
        E -= (E - eph.e * np.sin(E) - M) / (1.0 - eph.e * np.cos(E))

    sinE, cosE = np.sin(E), np.cos(E)
    nus = np.sqrt(1.0 - eph.e**2) * sinE
    nuc = cosE - eph.e
    nue = 1.0 - eph.e * cosE
    u = np.arctan2(nus, nuc) + eph.omg
    r = eph.A * nue 
    i = eph.i0 + eph.idot * tk
    sin2u, cos2u = np.sin(2*u), np.cos(2*u)
    u += eph.cus * sin2u + eph.cuc * cos2u
    r += eph.crs * sin2u +eph.crc * cos2u
    i += eph.cis * sin2u + eph.cic * cos2u
    x = r * np.cos(u)
    y = r * np.sin(u)
    cosi = np.cos(i)
    O = eph.OMG0 + (eph.OMGd - omge) * tk - omge * eph.toes
    sinO, cosO = np.sin(O), np.cos(O)
    rs = [x * cosO - y * cosi * sinO, x * sinO + y * cosi * cosO, y * np.sin(i)]
    tk = dtadjust(t, eph.toc)
    dts = eph.f0 + eph.f1 * tk + eph.f2 * tk**2
    # relativity correction
    dts -= 2 *np.sqrt(mu * eph.A) * eph.e * sinE / rCST.CLIGHT**2
    var = sva2ura(sys, eph.sva)
    trace(4, 'eph2pos: sat=%d, dts=%.10f rs=%.4f %.4f %.4f var=%.3f\n' % 
          (eph.sat, dts,rs[0],rs[1],rs[2], var))
    return rs, var, dts

def deq(x, acc):
    """glonass orbit differential equations """
    xdot = np.zeros(6)
    r2 = np.dot(x[0:3], x[0:3])
    if r2 <= 0.0:
        return xdot
    r3 = r2 * np.sqrt(r2)
    omg2 = rCST.OMGE_GLO**2

    a = 1.5 * rCST.J2_GLO * rCST.MU_GLO * rCST.RE_GLO**2 / r2 / r3 
    b = 5.0 * x[2]**2 / r2 
    c = -rCST.MU_GLO / r3 - a * (1.0 - b)
    xdot[0:3] = x[3:6]
    xdot[3] = (c + omg2) * x[0] + 2.0 * rCST.OMGE_GLO * x[4] + acc[0]
    xdot[4] = (c + omg2) * x[1] - 2.0 * rCST.OMGE_GLO * x[3] + acc[1]
    xdot[5] = (c - 2.0 * a) * x[2] + acc[2]
    return xdot

def glorbit(t, x, acc):
    """ glonass position and velocity by numerical integration """
    k1 = deq(x, acc)
    w =x + k1 * t / 2
    k2 = deq(w, acc)
    w = x + k2 * t / 2
    k3 = deq(w, acc)
    w = x + k3 * t
    k4 = deq(w, acc)
    x += (k1 + 2 * k2 + 2 * k3 + k4) * t / 6
    return x
    

def geph2pos(time, geph):
    """ GLONASS ephemeris to satellite position and clock bias """
    t = timediff(time, geph.toe)
    dts = -geph.taun + geph.gamn * t
    x = np.array((*geph.pos, *geph.vel))
    
    trace(4, 'geph2pos: sat=%d\n' % geph.sat)
    tt = -TSTEP if t < 0 else TSTEP
    while abs(t) > 1E-5:  #1E-9
        if abs(t) < TSTEP:
            tt = t
        x = glorbit(tt, x, geph.acc)
        t -= tt

    var = ERREPH_GLO**2
    return x[0:3], var, dts
    
def ephpos(time, eph):
    tt = 1e-3  # delta t to calculate velocity
    rs = np.zeros(6)
    
    if sat2prn(eph.sat)[0] != uGNSS.GLO:
        rs[0:3], var, dts = eph2pos(time, eph)
        # use delta t to determine velocity
        t = timeadd(time, tt)
        rs[3:6], _, dtst = eph2pos(t, eph)
    else: # GLONASS
        rs[0:3], var, dts = geph2pos(time, eph)
        # use delta t to determine velocity
        t = timeadd(time, tt)
        rs[3:6], _, dtst = geph2pos(t, eph)
    rs[3:6] = (rs[3:6] - rs[0:3]) / tt
    return rs, var, dts

def satpos(t, eph):
    return ephpos(t, eph)

def eph2clk(time, eph):
    """ calculate clock offset based on ephemeris """
    t = ts = timediff(time, eph.toc)
    for _ in range(2):
        t = ts - (eph.f0 + eph.f1 * t + eph.f2 * t**2)
    dts = eph.f0 + eph.f1*t + eph.f2 * t**2
    trace(4, 'ephclk: t=%.12f ts=%.12f dts=%.12f f0=%.12f f1=%.9f f2=%.9f\n' % (t,ts,dts,eph.f0,eph.f1,eph.f2))
    return dts

def geph2clk(time, geph):
    """ calculate GLONASS clock offset based on ephemeris """
    t = ts = timediff(time, geph.toe)
    for _ in range(2):
        t = ts - (-geph.taun + geph.gamn * t)
    trace(4, 'geph2clk: t=%.12f ts=%.12f taun=%.12f gamn=%.12f\n' % (t, ts, 
        geph.taun, geph.gamn))
    return -geph.taun + geph.gamn * t

def ephclk(time, eph):
    if sat2prn(eph.sat)[0] != uGNSS.GLO:
        return eph2clk(time, eph)
    else:
        return geph2clk(time, eph)

def satposs(obs, nav):
    """ satellite positions and clocks ----------------------------------------------
    * compute satellite positions, velocities and clocks
    * args     obs_t obs       I   observation data
    *          nav_t  nav      I   navigation data
    *          double rs       O   satellite positions and velocities (ecef)
    *          double dts      O   satellite clocks
    *          double var      O   sat position and clock error variances (m^2)
    *          int    svh      O   sat health flag (-1:correction not available)
    * return : none
    * notes  : rs [0:2] = obs[i] sat position {x,y,z} (m)
    *          rs [3:5] = obs[i] sat velocity {vx,vy,vz} (m/s)
    *          dts[0:1] = obs[i] sat clock {bias,drift} (s|s/s)
    *          var[i]   = obs[i] sat position and clock error variance (m^2)
    *          svh[i]    = obs[i] sat health flag
    *          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
    *          satellite position and clock are values at signal transmission time
    *          satellite position is referenced to antenna phase center
    *          satellite clock does not include code bias correction (tgd or bgd)
    *          any pseudorange and broadcast ephemeris are always needed to get
    *          signal transmission time """
    n = obs.sat.shape[0]
    rs = np.zeros((n, 6))
    dts = np.zeros(n)
    var = np.zeros(n)
    svh = np.zeros(n, dtype=int)
    
    ep = time2epoch(obs.t)
    trace(3, 'satposs  : teph= %04d/%02d/%02d %02d:%02d:%06.3f n=%d ephopt=%d\n' %           
          (ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], n, 0)) 
    
    for i in np.argsort(obs.sat):
        sat = obs.sat[i]
        # search any pseudorange
        pr = obs.P[i,0] if obs.P[i,0] != 0 else obs.P[i,1]
        # transmission time by satellite clock
        t = timeadd(obs.t, -pr / rCST.CLIGHT)

        eph = seleph(nav, t, sat)
        if eph is None:
            svh[i] = 1
            trace(2, 'No broadcast ephemeris: sat=%d\n' % sat)
            continue
        svh[i] = eph.svh
        # satellite clock bias by broadcast ephemeris
        dt = ephclk(t, eph)
        t = timeadd(t, -dt)
        # satellite position and clock at transmission time 
        rs[i], var[i], dts[i] = satpos(t, eph)
        trace(4,'satposs: %d,time=%.9f dt=%.9f pr=%.3f rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f\n' 
              % (obs.sat[i],t.sec,dt,pr,*rs[i,0:3],dts[i]*1e9,var[i]))


    return rs, var, dts, svh
