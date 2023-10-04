# Example files fit RGB stars are in this directory
# Explanations on most of the parameters (there is more, but they are not critical):
# model_fullname : Model to be computed. Be aware that the model class io_asymptotic must be also set in the Config/default/config_default.cfg in order to work
# model_type:  MUST BE FIX. If set to 0, it will generate p modes according to the first order asymptotic. If set to 1, it will generate l=1 p modes by shifting the l=0 p mode frequencies (de facto using the l=0 curvatures)
# bias_type:  MUST BE FIX. If set to 0: No spline fitting on the top of the asymptotic (by pass the Extra parameters table). 
#                         If set to 1: Cubic spline fitting to attempt to catch departure from the asymptotic result
#                         If set to 2: Hermite spline fitting
# freq_smoothness: Smoothness condition on l=0,2,3 p modes frequencies: It is equivalent to use a cubic spline
# eta0_switch: MUST BE FIXED. If 0, it will not compute the centrifugal effects. This means that a2_measured = a2_CF + a2_ELSE
#                             If 1, it will compute the centrifugal effects and add them to l=0,2,3 p modes (l=1 are discarded because we would need to account for mixing). Then, a2_measured = a2_ELSE
# Asymetry : Lorentzian mode asymetry
# Inclination : Stellar inclination
# Visibility_lx : Bolometric Mode visibility for mode of degree x. For the sun Vl1 = 1.50, Vl2 = 0.53, Vl3 = 0.08
# Height : Prior applied to ALL heights
# Width : Prior applied to ALL width. Fix_Auto will set it according to Dnu.
# delta01 : small separation. Fix_Auto will set it according to Dnu
# DP1 : Period Spacing for l=1 mixed modes
# alpha_g : Phase offset for g modes. Often called epsilon_g in publications
# q : coupling strength
# rot_env : Average Envelope rotation
# rot_core : Average Core rotation
# a2_env : a2 coefficient for l=2,3 modes. This is an averaged quantity on l and n
# a2_core: [NOT FUNCTIONAL] Could be used in the future if we come up with a mixing model for l=1 a2 coefficients. LET IT FIX TO 0 FOR THE MOMENT
# a3_env: a3 coefficient for l=2,3 modes. This is an averaged quantity on l and n
# a4_env: a4 coefficient for l=2,3 modes. This is an averaged quantity on l and n
# a5_env: a5 coefficient for l=2,3 modes. This is an averaged quantity on l and n
# a6_env: a6 coefficient for l=2,3 modes. This is an averaged quantity on l and n
# Hfactor : Height correction factor from the Asymptotic theory and for the l=1 mixed modes
# Wfactor : Width correction factor from the Asymptotic theory and for the l=1 mixed modes
# trunc_c : Coefficient of truncation of the Lorentzian. The smaller trunc_c, the closer to the central frequencies we cut. This is used to allow faster computation BUT it will slightly bias the Height and noise background.
#           trunc_c ~ 30 lead to a "power leak" between noise and height of a few percent.
