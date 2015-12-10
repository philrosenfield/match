# run zcombine, hybridMC, zcombine, zcmerge, and diagnostics scripts (after calcsfh -mcdata finished).
# NB:
# 1) Need to make sure match_base is correct!
# 2) Check the extensions that are hard coded.
# Usage e.g.: bash theworks.sh kdg73_f475w_f814w.out
match_base="$HOME/match2.5/bin"
$match_base/zcombine $1 -bestonly > ${1/.out}.sfh
$match_base/hybridMC $1.dat ${1/.out}.mcmc -tint=2.0 -nmc=10000 -dt=0.015 > ${1/.out}.mcmc.scrn
$match_base/zcombine ${1/.out}.mcmc -unweighted -medbest -jeffreys -best=${1/.out}.sfh > ${1/.out}.mcmc.zc
$match_base/zcmerge -absolute ${1/.out}.sfh ${1/.out}.mcmc.zc > ${1/.out}.mcmc.zc.dat
python -m dweisz.match.scripts.diagnostics -n ${1/.out}*
