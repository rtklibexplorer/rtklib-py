# configuration settings

from rtkcmn import uGNSS

# ----------- PPK options -------------------------------
nf = 2                   # num frequencies ( 1 or 2)
pmode = 'kinematic'      # static, kinematic
filtertype = 'forward'  # forward, backward, combined, combined_noreset
use_sing_pos = False     # run initial single precision sol each epoch, not
                         # necessary unless receiever clock errors are large
elmin = 15
cnr_min = 20             # min signal strength, not currently supported
excl_sat = []            # excluded sats

maxinno = 1              # outlier threshold for phase
maxage = 0.6              # mag age of differential, set at half base sample rate for now
maxout = 10               # maximum outage [epoch]
thresdop = 0             # cycle slip detection by doppler method
thresslip = 0.10         # cycle slip detection by geom-free LC

# ------------  Kalman Filter Statistics ------------------------
eratio = [300, 300]
err = [0, 0.003, 0.003, 0, 0, 5e-12]  # err sigmas [-, base, el, rcvstd, bl, satclk]
#err = [0, 0.00, 0.00, 0.5, 0, 5e-12] 
accelh = 3
accelv = 1
prnbias = 1e-4
sig_p0 = 30.0            # initial pos sigma
sig_v0 = 10.0            # initial vel/acc sigma
sig_n0 = 30.0            # inital bias sigma

#  ---------------- Ambiguity resolution options ----------------
armode = 1               # 0:off, 1:continuos,2:instantaneous,3:fix-and-hold
thresar = 3              # AR threshold
thresar1 = 0.1           # max pos variation for AR
elmaskar = 15            # elevation mask for AR
var_holdamb = 0.1

# ----------- Single precision options ----------------------------
sing_p0 = 100
sing_v0 = 10
sing_dt = 1
sing_elmin = 10
sing_sq = 1e-2
sing_q5 = 1e-2
sing_err = [0.0, 0.3, 0.3]

# -------------  Base and Rover positions ------------------
# base position, set to zeros to use rinex header pos
rb = [0, 0, 0]

# initial rover position/velocity for alignment to RTKLIB solution for debug
# Set to zero to use standard precision computed starting position
rr_f = rr_b  = [0, 0, 0]


# ----------- Configure observation signals ----------------

skip_sig_tbl = {uGNSS.GPS: [],   # skip these obs
                uGNSS.GAL: [],
                uGNSS.QZS: []}
gnss_t = [uGNSS.GPS, uGNSS.GAL]
freq_ix0 = {uGNSS.GPS: 0, uGNSS.GAL: 0} # L1
freq_ix1 = {uGNSS.GPS: 1, uGNSS.GAL: 3} # L2/E5b


# ---------- Frequencies currently supported-------------
freq = [1.57542e9,   # L1/E1
        1.22760e9,   # L2
        1.17645e9,   # L5/E5a/B2a
        1.20714e9]   # E5b

