awk -F'\t' '$2=="species" {count++}END{print"Número de registros hasta el nivel de especie:" count}' ../exam1/data/metagen/infants_metagenome.txt | wc -c
