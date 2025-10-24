continent="europe"
stn_ref_path="station_metadata"
filename="sample_station_metadata.csv"
full_path="$stn_ref_path/$filename"

while IFS="," read -r dc net sta
do
	#echo "$continent" "$dc" "$net" "$sta"
	sbatch dwnld_sngl_stn_slurm.sh "$continent" "$dc" "$net" "$sta"
	sleep 2

done < <(tail -n +2 "$full_path")
