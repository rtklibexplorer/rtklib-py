# configuration settings

from rtkcmn import uGNSS

# ----------- PPK options -------------------------------
nf = 2                   # num frequencies ( 1 or 2)
pmode = 'kinematic'      # static, kinematic
filtertype = 'forward'  # forward, backward, combined, combined_noreset
use_sing_pos = False     # run initial single precision sol each epoch, not
                         # necessary unless receiever clock errors are large
elmin = 15
cnr_min = 24             # min signal strength, not currently supported
excl_sat = []            # excluded sats

maxinno = 1              # outlier threshold for phase
maxage = 30              # mag age of differential
maxout = 4               # maximum outage [epoch]
thresdop = 5             # cycle slip detection by doppler method
thresslip = 0.10         # cycle slip detection by geom-free LC

# ------------  Kalman Filter Statistics ------------------------
eratio = [300, 100]
err = [0, 0.003, 0.003, 0.0, 0, 5e-12]  # err sigmas [-, base, el, rcvstd, bl, satclk]
#err = [0, 0.00, 0.00, 0.5, 0, 5e-12] 
accelh = 3
accelv = 1
prnbias = 1e-2
sig_p0 = 30.0            # initial pos sigma
sig_v0 = 10.0            # initial vel/acc sigma
sig_n0 = 30.0            # inital bias sigma

#  ---------------- Ambiguity resolution options ----------------
thresar = 3              # AR threshold
thresar1 = 0.05           # max pos variation for AR
armode = 0               # 0:off, 1:contunous,2:instantaneous,3:fix-and-hold
elmaskar = 15            # elevation mask for AR
var_holdamb = 0.01

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
rb = [-2703115.9211, -4291767.2078, 3854247.9065]

# initial rover position/velocity for alignment to RTKLIB
# rr_f = [-2694560.126,  -4296492.833,   3854817.084, # forwards
#         0.016,       0.043,         -0.027] 
# rr_b = [-2709836.020965, -4269097.509128, 3874376.664473, # backwards
#       0.018097,      0.051379,     -0.027851 ]
# Set to zero to use standard precision computed starting position
rr_f = rr_b  = [0, 0, 0]


# ----------- Configure observation signals ----------------

skip_sig_tbl = {uGNSS.GPS: [],   # skip these obs
                uGNSS.GAL: [],
                uGNSS.QZS: []}
gnss_t = [uGNSS.GPS, uGNSS.GAL]
freq_ix0 = {uGNSS.GPS: 0, uGNSS.GAL: 0}
freq_ix1 = {uGNSS.GPS: 2, uGNSS.GAL: 2}

# ---------- Frequencies currently supported-------------
freq = [1.57542e9,   # L1/E1
        1.22760e9,   # L2
        1.17645e9,   # L5/E5a/B2a
        1.20714e9]   # E5b

