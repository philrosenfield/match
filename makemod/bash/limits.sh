echo "age lmits:" > limits.txt 
ls mod1_* | sort -n -t _ -k 3 | cut -d _ -f 3 | tail -1  >> limits.txt
ls mod1_* | sort -rn -t _ -k 2 | cut -d _ -f 2 | tail -1 >> limits.txt
echo "logZ lmits:" >> limits.txt 
ls mod1_* | sort -n -t _ -k 4 | cut -d _ -f 4 | tail -1 >> limits.txt 
ls mod1_* | sort -rn -t _ -k 4 | cut -d _ -f 4 | tail -1 >> limits.txt 

