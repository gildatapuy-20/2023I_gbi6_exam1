N_archivo="denint.csv"
echo "Archivo,Filas,Columnas" > "$N_archivo"
Direc="./data/denvint"
for CS in "$Direc"/*.csv
do
Arch=$(basename "$CS" .csv)
filas=$(awk 'END{print NR}' "$CS")
columnas=$(awk -F ',' 'NR==1{print NF}' "$CS")
echo "$Arch,$filas,$columnas" >> "$N_archivo"
done
