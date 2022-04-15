# configuration settings

from rtkcmn import uGNSS, rSIG

# ----------- PPK options -------------------------------
nf = 2                   # num frequencies ( 1 or 2)
pmode = 'kinematic'      # static, kinematic
filtertype = 'forward'  # forward, backward, combined, combined_noreset
use_sing_pos = False     # run initial single precision sol each epoch, not
                         # necessary unless receiever clock errors are large
elmin = 15
cnr_min = 35             # min signal strength
excsats = []             # excluded sats

maxinno = 1              # outlier threshold for phase
maxage = 0.6*2              # mag age of differential, set at half base sample rate for now
maxout = 20               # maximum outage [epoch]
thresdop = 6             # cycle slip detection by doppler method
thresslip = 0.05         # cycle slip detection by geom-free LC
interp_base = False       # interpret base observations

# ------------  Kalman Filter Statistics ------------------------
eratio = [300, 300]
efact = {uGNSS.GPS: 1.0, uGNSS.GLO: 1.5, uGNSS.GAL: 1.0} # relative weighting of each constellation
err = [0, 0.003, 0.003, 0.0, 0, 5e-12]  # err sigmas [-, base, el, rcvstd, bl, satclk]
#err = [0, 0.00, 0.00, 0.5, 0, 5e-12] 
accelh = 3
accelv = 1
prnbias = 1e-4
sig_p0 = 30.0            # initial pos sigma
sig_v0 = 10.0            # initial vel/acc sigma
sig_n0 = 30.0            # inital bias sigma

#  ---------------- Ambiguity resolution options ----------------
armode = 3              # 0:off, 1:continuos,3:fix-and-hold
thresar = 3              # AR threshold
thresar1 = 0.1           # max pos variation for AR
elmaskar = 15            # elevation mask for AR
var_holdamb = 0.1
minfix = 20           # min fix samples to set hold
minfixsats = 4        # min sat pairs to test for fix
minholdsats = 5       # min sat pairs to test for hold
mindropsats = 10

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
#rr_f =[-1276972.378274, -4717193.586414,  4087245.657488,
#       -0.010286,       -0.015413,        0.015250]
#rr_b = [-1276984.364211, -4717218.261086,  4087215.802648,
#        -0.005020,        0.000248,        0.008825]
# Set to zero to use standard precision computed starting position
rr_f = [0, 0, 0, 0, 0, 0]
rr_b  = [0, 0, 0, 0, 0, 0]


# ----------- Configure observation signals ----------------

gnss_t = [uGNSS.GPS, uGNSS.GLO, uGNSS.GAL]

sig_tbl = {'1C': rSIG.L1C, '1X': rSIG.L1X, '1W': rSIG.L1W,
           '2W': rSIG.L2W, '2C': rSIG.L2C, '2X': rSIG.L2X,
           '5Q': rSIG.L5Q, '5X': rSIG.L5X, '7Q': rSIG.L7Q,
           '7X': rSIG.L7X}

skip_sig_tbl = {uGNSS.GPS: [],   # skip these obs
                uGNSS.GLO: [],
                uGNSS.GAL: [],
                uGNSS.QZS: []}

# set these from table below
freq_ix0 = {uGNSS.GPS: 0, uGNSS.GLO: 4, uGNSS.GAL: 0}# L1
freq_ix1 = {uGNSS.GPS: 1, uGNSS.GLO: 5, uGNSS.GAL: 3} # L2/E5b

# ---------- Frequencies currently supported-------------
freq = [1.57542e9,   # L1/E1
        1.22760e9,   # L2
        1.17645e9,   # L5/E5a/B2a
        1.20714e9,   # E5b
        1.60200E9,   # G1
        1.24600E9]   # G2
dfreq_glo = [0.56250E6, 0.43750E6]


