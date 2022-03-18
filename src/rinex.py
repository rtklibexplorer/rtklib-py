"""
module for RINEX 3.0x processing

Copyright (c) 2021 Rui Hirokawa (from CSSRLIB)
Copyright (c) 2022 Tim Everett
"""

import numpy as np
from rtkcmn import uGNSS, rSIG, Eph, prn2sat, gpst2time, Obs, epoch2time, timediff
import time

class rnx_decode:
    """ class for RINEX decoder """
    MAXSAT = uGNSS.GPSMAX+uGNSS.GLOMAX+uGNSS.GALMAX+uGNSS.BDSMAX+uGNSS.QZSMAX

    def __init__(self, cfg):
        self.ver = -1.0
        self.fobs = None
        self.freq_tbl = {uGNSS.GPS: {rSIG.L1C: 0, rSIG.L1X: 0, rSIG.L2W: 1, rSIG.L2L: 1,  
                                     rSIG.L2X: 1, rSIG.L5Q: 2, rSIG.L5X: 2},
                         uGNSS.GAL: {rSIG.L1C: 0, rSIG.L1X: 0, rSIG.L5Q: 2, rSIG.L5X: 2, 
                                     rSIG.L7Q: 3, rSIG.L7X: 3},
                         uGNSS.QZS: {rSIG.L1C: 0, rSIG.L1X: 0, rSIG.L2W: 1, rSIG.L2L: 1,  
                                     rSIG.L2X: 1, rSIG.L5Q: 2, rSIG.L5X: 2}}
        self.gnss_tbl = {'G': uGNSS.GPS, 'E': uGNSS.GAL, 'J': uGNSS.QZS}
        self.sig_tbl = {'1C': rSIG.L1C, '1X': rSIG.L1X, '1W': rSIG.L1W,
                        '2W': rSIG.L2W, '2L': rSIG.L2L, '2X': rSIG.L2X,
                        '5Q': rSIG.L5Q, '5X': rSIG.L5X, '7Q': rSIG.L7Q,
                        '7X': rSIG.L7X}
        self.skip_sig_tbl = cfg.skip_sig_tbl
        self.nf = 4
        self.sigid = np.ones((uGNSS.GNSSMAX, rSIG.SIGMAX*3), dtype=int)*rSIG.NONE
        self.typeid = np.ones((uGNSS.GNSSMAX, rSIG.SIGMAX*3), dtype=int)*rSIG.NONE
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

    def decode_nav(self, navfile, nav):
        """decode RINEX Navigation message from file """
        nav.eph = []
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
                eph = Eph(sat)

                year = int(line[4:8])
                month = int(line[9:11])
                day = int(line[12:14])
                hour = int(line[15:17])
                minute = int(line[18:20])
                sec = int(line[21:23])
                eph.toc = epoch2time([year, month, day, hour, minute, sec])
                eph.af0 = self.flt(line, 1)
                eph.af1 = self.flt(line, 2)
                eph.af2 = self.flt(line, 3)

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
                eph.code = int(self.flt(line, 1))  # source for GAL
                eph.week = int(self.flt(line, 2))

                line = fnav.readline() #23:26
                eph.sva = self.flt(line, 0)
                eph.svh = int(self.flt(line, 1))
                eph.tgd = float(self.flt(line, 2))
                if sys == uGNSS.GAL:
                    tgd_b = float(self.flt(line, 3))
                    if (eph.code >> 9) & 1:
                        eph.tgd = tgd_b
                else:
                    eph.iodc = int(self.flt(line, 3))

                line = fnav.readline() #27:30
                tot = int(self.flt(line, 0))
                if len(line) >= 42:
                    eph.fit = int(self.flt(line, 1))

                eph.toe = gpst2time(eph.week, eph.toes)
                eph.tot = gpst2time(eph.week, tot)
                nav.eph.append(eph)
        nav.eph.sort(key=lambda x: x.toe.time)
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
        nf = nav.nf
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
            obs.P = np.zeros((nsat, nf))
            obs.L = np.zeros((nsat, nf))
            obs.D = np.zeros((nsat, nf))
            obs.S = np.zeros((nsat, nf))
            obs.lli = np.zeros((nsat, nf), dtype=int)
            obs.Pstd = np.zeros((nsat, nf), dtype=int)
            obs.Lstd = np.zeros((nsat, nf), dtype=int)
            obs.mag = np.zeros((nsat, nf))
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
                    f = i // (self.nsig[sys] // self.nband[sys])
                    if self.typeid[sys][i] == 0:  # code
                        obs.P[n, f] = float(obs_)
                        Pstd = line[16*i+18]
                        obs.Pstd[n, f] = int(Pstd) if Pstd != " " else 0
                    elif self.typeid[sys][i] == 1:  # carrier
                        obs.L[n, f] = float(obs_)
                        lli = line[16*i+17]
                        obs.lli[n, f] = int(lli) if lli != " " else 0
                        Lstd = line[16*i+18]
                        obs.Lstd[n, f] = int(Lstd) if Lstd != " " else 0
                    elif self.typeid[sys][i] == 2:  # C/No
                        obs.S[n, f] = float(obs_)
                    elif self.typeid[sys][i] == 3:  # Doppler
                            obs.D[n, f] = float(obs_)
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
    # sync base and rover, step in one obs to sync
    _, _ =next_obs(nav, rov, base, dir)
    # step back to first obs
    obsr, obsb = next_obs(nav, rov, base, -dir)
    return obsr, obsb

def next_obs(nav, rov, base, dir):
    """ sync observation beteen rover and base """
    if abs(dir) != 1:
        return [], []
    rov.index += dir   # 1=forward, -1=backward
    while True:
        if rov.index < 0 or rov.index >= len(rov.obslist) or \
                base.index < 0 or base.index >= len(base.obslist):
            return [], []
        obsr = rov.obslist[rov.index]
        obsb = base.obslist[base.index]
        dt = timediff(obsr.t, obsb.t)
        if dt * dir > nav.maxage:
            base.index += dir
            if base.index < 0 or base.index > len(base.obslist) - 1:
                base.index -= dir
                break
        if np.abs(dt) <= nav.maxage:
            break
        elif dt * dir < 0:
            rov.index += dir
            if rov.index < 0 or rov.index > len(rov.obslist) - 1:
                rov.index -= dir
                break
    return obsr, obsb

def rcvstds(nav, obs):
    """ decode receiver stdevs from rinex fields """
    # skip if weighting factor is zero
    if nav.err[3] == 0:
        return
    for i in np.argsort(obs.sat):
        for f in range(nav.nf):
            s = obs.sat[i] - 1
            # decode receiver stdevs, 
            # Lstd: 0.004 cycles -> m
            nav.rcvstd[s,f] = obs.Lstd[i,f] * 0.004 * 0.2
            # Pstd: 0.01*2^(n+5)
            nav.rcvstd[s,f+nav.nf] = 0.01 * (1 << (obs.Pstd[i,f] + 5))

