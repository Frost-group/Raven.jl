# Raven.jl

Raven.jl: Small polaron mobility by hopping. (A work in progress.)

# --- Glossary ----------------------------------------------------------------
# a,b,c      = lattice dims
# N          = number of mobile ions (Nions)
# M          = number of defects (Ndef); occ = -1 means defect
# SITES      = a*b*c
# steps      = total attempted moves (N * sweeps)
# sweeps     = ceil(steps / N)
# lagtime    = τ used for MSD windows
# pos        = N×3 Int, ion coordinates (1-based, PBC)
# occ        = a×b×c Int, 0 empty, >0 ion id, -1 defect
# disp       = 3×N Int, cumulative displacement per ion
# dr_log     = (sweeps+1)×3×N Int, displacement snapshot by sweep
# acc_log    = (sweeps+1)×N Int, accepted hops per ion
# Dtr/Dbulk  = tracer / collective diffusion
# Haven      = Dtr / Dbulk
# ----------------------------------------------------------------------------- 
