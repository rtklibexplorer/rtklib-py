# configuration settings

from rtkcmn import uGNSS, rSIG

# ----------- PPK options -------------------------------
nf = 2                   # num frequencies ( 1 or 2)
pmode = 'kinematic'      # static, kinematic
filtertype = 'forward'   # forward, backward, combined, combined_noreset
use_sing_pos = False     # run initial single precision sol each epoch, not
                         # necessary unless receiever clock errors are large
elmin = 15               # minimum elevation for float solution (degrees)
cnr_min = [28, 20]   # min signal strength [freq1, freq2] (dB-Hz)
excsats = []             # excluded sats

maxinno = 1              # outlier threshold for phase (m)
maxcode = 10             # outlier threshold for code (m)
maxage = 30              # mag age of differential
maxout = 4               # maximum outage [epoch]
thresdop = 5             # cycle slip detection by doppler method
thresslip = 0.10         # cycle slip detection by geom-free LC
interp_base = False      # interpolate base observations

# ------------  Kalman Filter Statistics ------------------------
eratio = [300, 100]    # ratio between psuedorange noise and carrier phase noise for L1, L5
efact = {uGNSS.GPS: 1.0, uGNSS.GLO: 1.5, uGNSS.GAL: 1.0} # relative weighting of each constellation
err = [0, 0.003, 0.003, 0.0, 0, 0, 5e-12]  # err sigmas [-, base, el, bl, snr, rcvstd, satclk]
snrmax = 45              # max signal strength for variance calc (dB-Hz)
accelh = 3               # horiz accel noise sigma (m/sec2)
accelv = 1               # vert accel noise sigma (m/sec2)
prnbias = 1e-2           # Carrier phase bias sigma ( cycles)
sig_p0 = 30.0            # initial pos sigma (m)
sig_v0 = 10.0            # initial vel/acc sigma (m/sec)
sig_n0 = 30.0            # inital bias sigma (m)

#  -------------Ambiguity resolution options ----------------
armode = 0               # 0:off, 1:continuos,3:fix-and-hold
thresar1 = 0.05           # max pos variation for AR (and accel update in kalman filter)
# the rest of these aren't used since AR is disabled
minlock = 0              # min consecutive fix samples to include sat in AR 
thresar = 3              # AR threshold
glo_hwbias = 0.0         # GLONASS HW bias
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

# -------------  Base and Rover positions ------------------
# base position, set to zeros to use rinex header pos
#rb = [-2703115.9211, -4291767.2078, 3854247.9065] # SLAC
rb = [0, 0, 0]  # will set based on base file name

# initial rover position/velocity for alignment to RTKLIB
#rr_f = [-2694556.682853,  -4296492.691977,   3854819.048563, # forwards
#        -0.009033 ,       0.042808,         -0.027949] 
#rr_b = [-2709836.020965, -4269097.509128, 3874376.664473, # backwards
#      0.018097,      0.051379,     -0.027851 ]
# Set to zero to use standard precision computed starting position
rr_f = [0, 0, 0, 0, 0, 0]
rr_b  = [0, 0, 0, 0, 0, 0]


# ----------- Configure observation signals ----------------

gnss_t = [uGNSS.GPS, uGNSS.GLO, uGNSS.GAL]

# Valid signals
sig_tbl = {'1C': rSIG.L1C, '1X': rSIG.L1X, '1W': rSIG.L1W,
           '2W': rSIG.L2W, '2L': rSIG.L2L, '2X': rSIG.L2X,
           '5Q': rSIG.L5Q, '5X': rSIG.L5X, '7Q': rSIG.L7Q,
           '7X': rSIG.L7X}

skip_sig_tbl = {uGNSS.GPS: [],   # skip these obs
                uGNSS.GLO: [],
                uGNSS.GAL: [],
                uGNSS.QZS: []}

# set these from table below
freq_ix0 = {uGNSS.GPS: 0, uGNSS.GLO: 4, uGNSS.GAL: 0} # L1
freq_ix1 = {uGNSS.GPS: 2, uGNSS.GLO: 5, uGNSS.GAL: 2} # L5

# ---------- Frequencies currently supported-------------
freq = [1.57542e9,   # L1/E1
        1.22760e9,   # L2
        1.17645e9,   # L5/E5a/B2a
        1.20714e9,   # E5b
        1.60200E9,   # G1
        1.24600E9]   # G2
dfreq_glo = [0.56250E6, 0.43750E6]  # L1, L2


