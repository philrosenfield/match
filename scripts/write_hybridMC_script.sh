###
# A script to print the commands for running hybridMC and MC on all *sfh files in a directory.
#
# run in a directory containing .sfh files (the output from calcsfh):
# $ bash write_hybridMC_script.sh > script.sh
# then edit script.sh as needed and
# $ bash script.sh
#
# NB: File extensions are hard coded.
###

# locations:
hybridMC="$HOME/research/match2.5/bin/hybridMC"
zcombine="$HOME/research/match2.5/bin/zcombine"

# if unix and want to run parallel this sets a command for each processor. 
# I then will edit the output to increase from 0 to however many processors I want.
preCMD="nice -n +19 taskset -c 0"

# flags for HMC and ZC
hmcflags="-tint=2.0 -nmc=10000 -dt=0.015"
zcflags="-unweighted -medbest -jeffreys"

# assumes extension is .sfh (could make it input command with *$1)
ext=".sfh"
for l in $(ls *$ext)
do 
p=${l/$ext}
echo "$preCMD $hybridMC $p.out.dat $p.mcmc $hmcflags > $p.mcmc.scrn &"
done

for l in $(ls *$ext)
do 
p=${l/$ext}
echo "$preCMD $zcombine $p.mcmc $zcflags -best=$p.sfh > $p.mcmc.zc &"
done

# I will cut and paste this between commands from taskset -c MAX and taskset -c 0 so it will wait until the jobs are done to start the next set.
echo "for job in \`jobs -p\`; do echo $job; wait $job; done"
