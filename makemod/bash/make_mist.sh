cd ../; make MIST; cd MIST; ./makemod > makemod.out
cd ../; ./bin/calcbinmod -MIST
./bin/ascii2binmod -MIST
