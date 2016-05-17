ls raw/*/*dat | cut -d '_' -f 4 | sed 's/.dat//' | sed 's/M//' | sort -ug | tr '\n' ','
