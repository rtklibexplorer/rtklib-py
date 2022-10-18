# configuration settings

from rtkcmn import uGNSS, rSIG

# ----------- PPK options -------------------------------
nf = 2                   # num frequencies ( 1 or 2)
pmode = 'kinematic'      # static, kinematic
filtertype = 'forward'   # forward, backward, combined, combined_noreset
use_sing_pos = False     # run initial single precision sol each epoch, not
                         # necessary unless receiever clock errors are large
elmin = 15               # minimum elevation for float solution (degrees)
cnr_min = [35, 35]       # min signal strength [freq1, freq2] (dB-Hz)
excsats = []             # excluded sats

maxinno = 1              # outlier threshold for phase (m)
maxcode = 10             # outlier threshold for code (m)
maxage = 30              # mag age of differential
maxout = 20              # maximum outage [epoch]
thresdop = 6             # cycle slip detection by doppler method
thresslip = 0.05         # cycle slip detection by geom-free LC
interp_base = False      # interpolate base observations

# ------------  Kalman Filter Statistics ------------------------
eratio = [300, 300]      # ratio between psuedorange noise and carrier phase noise for L1, L2
efact = {uGNSS.GPS: 1.0, uGNSS.GLO: 1.5, uGNSS.GAL: 1.0} # relative weighting of each constellation
err = [0, 0.003, 0.003, 0.0, 0, 0, 5e-12]  # error sigmas [-, base, el, bl, snr, rcvstd, satclk]
snrmax = 52              # max signal strength for variance calc (dB-Hz)
accelh = 3               # horiz accel noise sigma (m/sec2)
accelv = 1               # vert accel noise sigma (m/sec2)
prnbias = 1e-4           # Carrier phase bias sigma ( cycles)
sig_p0 = 30.0            # initial pos sigma (m)
sig_v0 = 10.0            # initial vel/acc sigma (m/sec)
sig_n0 = 30.0            # inital bias sigma (m)

#  -------------Ambiguity resolution options ----------------
armode = 3               # 0:off, 1:continuos,3:fix-and-hold
thresar = 3              # AR threshold
minlock = 0              # min consecutive fix samples to include sat in AR 
glo_hwbias = 0.0         # GLONASS HW bias
thresar1 = 0.1           # max pos variation for AR (and accel update in kalman filter)
elmaskar = 15            # elevation mask for AR
var_holdamb = 0.1        # Hold ambiguity variance (m)
minfix = 20              # min fix samples to set hold
minfixsats = 4           # min sat pairs to test for fix
minholdsats = 5          # min sat pairs to test for hold
mindropsats = 10         # min sat pairs to drop sats from AR

# -----------  Single precision parameters ----------------------------
sing_p0 = 100            # initial pos sigma
sing_v0 = 10             # initial vel/acc sigma
sing_elmin = 10          # minimum elevation (degrees)

# -------------Base and Rover positions ------------------
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

# Valid signals
sig_tbl = {'1C': rSIG.L1C, '1X': rSIG.L1X, '1W': rSIG.L1W,
           '2W': rSIG.L2W, '2C': rSIG.L2C, '2X': rSIG.L2X,
           '5Q': rSIG.L5Q, '5X': rSIG.L5X, '7Q': rSIG.L7Q,
           '7X': rSIG.L7X}

skip_sig_tbl = {uGNSS.GPS: [],   # skip these obs
                uGNSS.GLO: [],
                uGNSS.GAL: [],
                uGNSS.QZS: []}

# set these from table below
freq_ix0 = {uGNSS.GPS: 0, uGNSS.GLO: 4, uGNSS.GAL: 0} # L1
freq_ix1 = {uGNSS.GPS: 1, uGNSS.GLO: 5, uGNSS.GAL: 3} # L2/E5b

# ---------- Frequencies currently supported-------------
freq = [1.57542e9,   # L1/E1
        1.22760e9,   # L2
        1.17645e9,   # L5/E5a/B2a
        1.20714e9,   # E5b
        1.60200E9,   # G1
        1.24600E9]   # G2
dfreq_glo = [0.56250E6, 0.43750E6]  # L1, L2


