"""
module for RINEX 3.0x processing

Copyright (c) 2021 Rui Hirokawa (from CSSRLIB)
Copyright (c) 2022 Tim Everett
"""

import numpy as np
from copy import deepcopy
from rtkcmn import uGNSS, rSIG, Eph, Geph, prn2sat, gpst2time, time2gpst, Obs, \
                    epoch2time, timediff, timeadd, utc2gpst
import rtkcmn as gn
from ephemeris import satposs

class rnx_decode:
    """ class for RINEX decoder """
    MAXSAT = uGNSS.GPSMAX+uGNSS.GLOMAX+uGNSS.GALMAX+uGNSS.BDSMAX+uGNSS.QZSMAX

    def __init__(self, cfg):
        self.ver = -1.0
        self.fobs = None
        self.gnss_tbl = {'G': uGNSS.GPS, 'E': uGNSS.GAL, 'R': uGNSS.GLO, 'J': uGNSS.QZS}
        self.sig_tbl = cfg.sig_tbl
        self.skip_sig_tbl = cfg.skip_sig_tbl
        self.nf = 4
        self.sigid = np.ones((uGNSS.GNSSMAX, rSIG.SIGMAX*3), dtype=int) * rSIG.NONE
        self.typeid = np.ones((uGNSS.GNSSMAX, rSIG.SIGMAX*3), dtype=int) * rSIG.NONE
        self.nsig = np.zeros((uGNSS.GNSSMAX), dtype=int)
        self.nband = np.zeros((uGNSS.GNSSMAX), dtype=int)
        self.pos = np.array([0, 0, 0])

    def flt(self, u, c=-1):
        if c >= 0:
            u = u[19*c+4:19*(c+1)+4]
        try:
            return float(u.replace("D", "E"))
        except:
            return 0
        
    
    def adjday(self, t, t0):
        """" adjust time considering week handover  """
        tt = timediff(t, t0)
        if tt < -43200.0:
            return timeadd(t, 86400.0)
        if tt > 43200.0:
            return timeadd(t,-86400.0)
        return t

    def decode_nav(self, navfile, nav):
        """decode RINEX Navigation message from file """
        nav.eph = []
        nav.geph = []
        with open(navfile, 'rt') as fnav:
            for line in fnav:
                if line[60:73] == 'END OF HEADER':
                    break
                elif line[60:80] == 'RINEX VERSION / TYPE':
                    self.ver = float(line[4:10])
                    if self.ver < 3.02:
                        return -1
                elif line[60:76] == 'IONOSPHERIC CORR':
                    if line[0:4] == 'GPSA' or line[0:4] == 'QZSA':
                        for k in range(4):
                            nav.ion[0, k] = self.flt(line[5+k*12:5+(k+1)*12])
                    if line[0:4] == 'GPSB' or line[0:4] == 'QZSB':
                        for k in range(4):
                            nav.ion[1, k] = self.flt(line[5+k*12:5+(k+1)*12])

            for line in fnav:
                if line[0] not in self.gnss_tbl:
                    continue
                sys = self.gnss_tbl[line[0]]
                prn = int(line[1:3])
                if sys == uGNSS.QZS:
                    prn += 192
                sat = prn2sat(sys, prn)
                year = int(line[4:8])
                month = int(line[9:11])
                day = int(line[12:14])
                hour = int(line[15:17])
                minute = int(line[18:20])
                sec = int(line[21:23])
                toc = epoch2time([year, month, day, hour, minute, sec])
                if sys != uGNSS.GLO:
                    eph = Eph(sat)
                    eph.toc = toc
                    eph.f0 = self.flt(line, 1)
                    eph.f1 = self.flt(line, 2)
                    eph.f2 = self.flt(line, 3)
    
                    line = fnav.readline() #3:6
                    eph.iode = int(self.flt(line, 0)) 
                    eph.crs = self.flt(line, 1)
                    eph.deln = self.flt(line, 2)
                    eph.M0 = self.flt(line, 3)
    
                    line = fnav.readline() #7:10
                    eph.cuc = self.flt(line, 0)
                    eph.e = self.flt(line, 1)
                    eph.cus = self.flt(line, 2)
                    sqrtA = self.flt(line, 3)
                    eph.A = sqrtA**2
    
                    line = fnav.readline() #11:14
                    eph.toes = int(self.flt(line, 0))
                    eph.cic = self.flt(line, 1)
                    eph.OMG0 = self.flt(line, 2)
                    eph.cis = self.flt(line, 3)
    
                    line = fnav.readline() #15:18
                    eph.i0 = self.flt(line, 0)
                    eph.crc = self.flt(line, 1)
                    eph.omg = self.flt(line, 2)
                    eph.OMGd = self.flt(line, 3)
    
                    line = fnav.readline() #19:22
                    eph.idot = self.flt(line, 0)
                    eph.code = int(self.flt(line, 1))  # source for GAL NAV type
                    eph.week = int(self.flt(line, 2))
    
                    line = fnav.readline() #23:26
                    eph.sva = self.flt(line, 0)
                    eph.svh = int(self.flt(line, 1))
                    tgd = np.zeros(2)
                    tgd[0] = float(self.flt(line, 2))
                    if sys == uGNSS.GAL:
                        tgd[1] = float(self.flt(line, 3))
                    else:
                        eph.iodc = int(self.flt(line, 3))
                    eph.tgd = tgd
    
                    line = fnav.readline() #27:30
                    tot = int(self.flt(line, 0))
                    if len(line) >= 42:
                        eph.fit = int(self.flt(line, 1))
    
                    eph.toe = gpst2time(eph.week, eph.toes)
                    eph.tot = gpst2time(eph.week, tot)
                    nav.eph.append(eph)
                else:  # GLONASS
                    if prn > uGNSS.GLOMAX:
                        print('Reject nav entry: %s' % line[:3])
                        break
                    geph = Geph(sat)
                    # Toc rounded by 15 min in utc 
                    week, tow = time2gpst(toc)
                    toc = gpst2time(week,np.floor((tow + 450.0) / 900.0) * 900)
                    dow = int(np.floor(tow  / 86400.0))
                    # time of frame in UTC 
                    tod = self.flt(line, 2) % 86400
                    tof = gpst2time(week ,tod + dow * 86400.0)
                    tof = self.adjday(tof, toc)
                    geph.toe = utc2gpst(toc)
                    geph.tof = utc2gpst(tof)
                    # IODE = Tb (7bit), Tb =index of UTC+3H within current day
                    geph.iode = int(((tow + 10800.0) % 86400) / 900.0 + 0.5)
                    geph.taun = -self.flt(line, 1)
                    geph.gamn = self.flt(line, 2)
                    
                    line = fnav.readline() #3:6
                    pos =np.zeros(3)
                    vel = np.zeros(3)
                    acc = np.zeros(3)
                    pos[0] = self.flt(line, 0)
                    vel[0] = self.flt(line, 1)
                    acc[0] = self.flt(line, 2)
                    geph.svh = self.flt(line, 3)
                    
                    line = fnav.readline() #7:10
                    pos[1] = self.flt(line, 0)
                    vel[1] = self.flt(line, 1)
                    acc[1] = self.flt(line, 2)
                    geph.frq = self.flt(line, 3)
                    nav.glofrq[sat - uGNSS.GPSMAX - 1] = int(geph.frq)

                    line = fnav.readline() #11:14
                    pos[2] = self.flt(line, 0)
                    vel[2] = self.flt(line, 1)
                    acc[2] = self.flt(line, 2)                                      
                    geph.age = self.flt(line, 2)
                    
                    geph.pos = pos * 1000
                    geph.vel = vel * 1000
                    geph.acc = acc * 1000
                    
                    nav.geph.append(geph)
    
        #nav.eph.sort(key=lambda x: (x.sat, x.toe.time))
        nav.eph.sort(key=lambda x: x.toe.time)
        nav.geph.sort(key=lambda x: x.toe.time)
        return nav

    def decode_obsh(self, obsfile):
        self.fobs = open(obsfile, 'rt')
        for line in self.fobs:
            if line[60:73] == 'END OF HEADER':
                break
            if line[60:80] == 'RINEX VERSION / TYPE':
                self.ver = float(line[4:10])
                if self.ver < 3.02:
                    return -1
            elif line[60:79] == 'APPROX POSITION XYZ':
                self.pos = np.array([float(line[0:14]),
                                     float(line[14:28]),
                                     float(line[28:42])])
            elif line[60:79] == 'SYS / # / OBS TYPES':
                if line[0] in self.gnss_tbl:
                    sys = self.gnss_tbl[line[0]]
                else:
                    continue
                self.nsig[sys] = int(line[3:6])
                s = line[7:7+4*13]
                if self.nsig[sys] >= 14:
                    line2 = self.fobs.readline()
                    s += line2[7:7+4*13]

                for k in range(self.nsig[sys]):
                    sig = s[4*k:3+4*k]
                    if sig[1:3] not in self.sig_tbl:
                        continue
                    if self.sig_tbl[sig[1:3]] in self.skip_sig_tbl[sys]:
                        continue
                    if sig[0] == 'C':
                        self.typeid[sys][k] = 0
                    elif sig[0] == 'L':
                        self.typeid[sys][k] = 1
                    elif sig[0] == 'S':
                        self.typeid[sys][k] = 2
                    elif sig[0] == 'D':
                        self.typeid[sys][k] = 3
                    else:
                        continue
                    self.sigid[sys][k] = self.sig_tbl[sig[1:3]]
                self.nband[sys] = len(np.where(self.typeid[sys]==1)[0])
        return 0

    def decode_obs(self, nav, maxepoch):
        """decode RINEX Observation message from file """

        self.obslist = []
        nepoch = 0
        for line in self.fobs:
            if line == '':
                break
            if line[0] != '>':
                continue
            obs = Obs()
            nsat = int(line[32:35])
            year = int(line[2:6])
            month = int(line[7:9])
            day = int(line[10:12])
            hour = int(line[13:15])
            minute = int(line[16:18])
            sec = float(line[19:29])
            obs.t = epoch2time([year, month, day, hour, minute, sec])
            obs.P = np.zeros((nsat, gn.MAX_NFREQ))
            obs.L = np.zeros((nsat, gn.MAX_NFREQ))
            obs.D = np.zeros((nsat, gn.MAX_NFREQ))
            obs.S = np.zeros((nsat, gn.MAX_NFREQ))
            obs.lli = np.zeros((nsat, gn.MAX_NFREQ), dtype=int)
            obs.Pstd = np.zeros((nsat, gn.MAX_NFREQ), dtype=int)
            obs.Lstd = np.zeros((nsat, gn.MAX_NFREQ), dtype=int)
            obs.mag = np.zeros((nsat, gn.MAX_NFREQ))
            obs.sat = np.zeros(nsat, dtype=int)
            n = 0
            for k in range(nsat):
                line = self.fobs.readline()
                if line[0] not in self.gnss_tbl:
                    continue
                sys = self.gnss_tbl[line[0]]
                if sys not in nav.gnss_t:
                    continue
                prn = int(line[1:3])
                if sys == uGNSS.QZS:
                    prn += 192
                obs.sat[n] = prn2sat(sys, prn)
                if obs.sat[n] == 0:
                    continue
                nsig_max = (len(line) - 4 + 2) // 16
                for i in range(self.nsig[sys]):
                    if i >= nsig_max:
                        break
                    obs_ = line[16*i+4:16*i+17].strip()
                    if obs_ == '' or self.sigid[sys][i] == 0:
                        continue
                    try:
                        obsval = float(obs_)
                    except:
                        obsval = 0
                    f = i // (self.nsig[sys] // self.nband[sys])
                    if f >= gn.MAX_NFREQ:
                        print('Obs file too complex, please use RTKCONV to remove unused signals')
                        raise SystemExit
                    if self.typeid[sys][i] == 0:  # code
                        obs.P[n, f] = obsval
                        Pstd = line[16*i+18]
                        obs.Pstd[n, f] = int(Pstd) if Pstd != " " else 0
                    elif self.typeid[sys][i] == 1:  # carrier
                        obs.L[n, f] = float(obs_)
                        lli = line[16*i+17]
                        obs.lli[n, f] = int(lli) if lli != " " else 0
                        Lstd = line[16*i+18]
                        obs.Lstd[n, f] = int(Lstd) if Lstd != " " else 0
                    elif self.typeid[sys][i] == 2:  # C/No
                        obs.S[n, f] = obsval
                    elif self.typeid[sys][i] == 3:  # Doppler
                            obs.D[n, f] = obsval
                n += 1
            obs.P = obs.P[:n, :]
            obs.L = obs.L[:n, :]
            obs.Pstd = obs.Pstd[:n, :]
            obs.Lstd = obs.Lstd[:n, :]
            obs.D = obs.D[:n, :]
            obs.S = obs.S[:n, :]
            obs.lli = obs.lli[:n, :]
            obs.mag = obs.mag[:n, :]
            obs.sat = obs.sat[:n]
            self.obslist.append(obs)
            nepoch += 1
            if maxepoch != None and nepoch >= maxepoch:
                break
        self.index = 0
        self.fobs.close()

    
    def decode_obsfile(self, nav, obsfile, maxepoch):
        self.decode_obsh(obsfile)
        self.decode_obs(nav, maxepoch)

def first_obs(nav, rov, base, dir):
    if dir == 1: # forward solution
        rov.index = base.index = 0
    else: # backward solution
        rov.index = len(rov.obslist) - 1
        base.index = len(base.obslist) - 1
    # sync base and rover, step one obs to sync
    _, _ =next_obs(nav, rov, base, dir)
    # step back to first obs
    obsr, obsb = next_obs(nav, rov, base, -dir)
    return obsr, obsb

def next_obs(nav, rov, base, dir):
    """ sync observations between rover and base """
    rov.index += dir   # 1=forward, -1=backward
    if abs(dir) != 1 or rov.index < 0 or rov.index >= len(rov.obslist):
        return [], []
    obsr, obsb = rov.obslist[rov.index], base.obslist[base.index]
    dt = timediff(obsr.t, obsb.t)
    baseChange = False
    ixb = base.index + dir
    while True:
        if ixb < 0 or ixb >= len(base.obslist):
            ixb -= dir
            dt_next = dt
            break # hit end of obs list
        dt_next = timediff(obsr.t, base.obslist[ixb].t)
        if abs(dt_next) >= abs(dt):
            break # next base obs is not closer
        else:
            base.index = ixb
            ixb += dir
            baseChange = True
            dt = dt_next
        
    if baseChange and nav.interp_base and len(nav.sol) > 0:
        # save base residuals for next epoch
        nav.obsb = deepcopy(obsb)
        nav.rsb, nav.varb, nav.dtsb, nav.svhb = satposs(obsb, nav)
    obsb = base.obslist[base.index]
    return obsr, obsb

def rcvstds(nav, obs):
    """ decode receiver stdevs from rinex fields """
    # skip if weighting factor is zero
    if nav.err[5] == 0:
        return
    for i in np.argsort(obs.sat):
        for f in range(nav.nf):
            s = obs.sat[i] - 1
            # decode receiver stdevs, 
            # Lstd: 0.004 cycles -> m
            nav.rcvstd[s,f] = obs.Lstd[i,f] * 0.004 * 0.2
            # Pstd: 0.01*2^(n+5)
            nav.rcvstd[s,f+nav.nf] = 0.01 * (1 << (obs.Pstd[i,f] + 5))

