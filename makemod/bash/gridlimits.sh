#!/bin/bash
# find L and T extrema in "unprocessed" files.
# each sub directory
for dir in ov0.30 ov0.40 ov0.50 ov0.60
do
    cd $dir
    pwd
    for l in $(ls */*dat)
    do
        # take extrema of column 4 (Mbol) and 3 (logT)
        # head -2 to account for header
        sort -grk 4 $l | head -2 >> ../l.list
        sort -gk 4 $l | head -2 >> ../l.list
        sort -grk 3 $l | head -2 >> ../t.list
        sort -gk 3 $l | head -2 >> ../t.list
    done
    cd ../
done

# remove all empty lines and lines that start with #
sed -i ' ' -e 's/#.*$//' -e '/^$/d' l.list
sed -i ' ' -e 's/#.*$//' -e '/^$/d' t.list

# print extrema (delete | cut if you want the full line)
echo MOD_L0:
sort -grk 4 l.list | tail -1 | cut -d ' ' -f 4
echo MOD_LF:
sort -gk 4 l.list | tail -1 | cut -d ' ' -f 4

echo MOD_T0:
sort -grk 3 t.list | tail -1 | cut -d ' ' -f 3
echo MOD_TF:
sort -gk 3 t.list | tail -1 | cut -d ' ' -f 3

#rm l.list
#rm t.list
