# Impute2 Pipeline
# Matt Ravenhall
#
# Call as:
# bash imputeFlipped.sh <chr>
#
# Required files:
# Reference Files:
# 	genetic map: ./${refDir}/genetic_map_chrX_combined_b37.txt
# 	haplotype: ./${refDir}.hap
# 	legend: ./${refDir}.legend
# Input Files:
# 	${toImpFile}.gen

impute2dir='./'
refDir='./'
refName='PLACEHOLDER_CHR'
mapNamePre='genetic_map_chr'
mapNamePost='_combined_b37.txt'
toImpFile='PLACEHOLDER_CHR'
outprefix='./PLACEHOLDER'
psLimit=5

# loop each chromosome
if [ -z "$1" ]; then
	exit 'No chromosome number provided'
else
	i=$1
fi

	hapLegFile=${refDir}${refName}

	# set initial bp locations ≈5000000bp, also grab the last snp bp location
	bpStart=1
	bpEnd=$(($bpStart + 5000000 - 1))

	chrLen=$(tail -1 ${toImpFile}${i}.map | cut -f4)

	# Create working directory if non-existent
	if [ ! -d "$WorkDir" ]; then mkdir -p WorkDir; fi

	# Prep for multi-threading
	psHolder=()

	# Perform within-chromosome iterations
	while [ $bpStart -lt $chrLen ]; do

		if [ $bpEnd -gt $chrLen ]; then
			bpEnd=$chrLen
		fi

		# Hold max processes
		((i=i%$psLimit)); ((i++==0)) && wait

		${impute2dir}impute2 -m ${refDir}${mapNamePre}${i}${mapNamePost} -h ${hapLegFile}${i}.hap -l ${hapLegFile}${i}.legend -g ${toImpFile}${i}.gen -int ${bpStart} ${bpEnd} -o ${outprefix}_chr${i}_bp${bpStart}to${bpEnd}.segment &
		psHolder+=($!)

		bpStart=$(($bpEnd + 1))
		bpEnd=$(($bpStart + 5000000 - 1))

	done

	# Hold until all sub-processes are complete
	wait "${psHolder[@]}"

	# Create output directories if they don't exist
	if [ ! -d "$Imputed/info_by_sample" ]; then mkdir -p Imputed/info_by_sample; fi
	if [ ! -d "$Imputed/info" ]; then mkdir -p Imputed/info; fi
	if [ ! -d "$Imputed/warnings" ]; then mkdir -p Imputed/warnings; fi
	if [ ! -d "$Imputed/summary" ]; then mkdir -p Imputed/summary; fi

	# Move log/info files to sub-directories, keep imputed files in ./Imputed/
	mv ${outprefix}_chr${i}_bp*.segment_info_by_sample ./Imputed/info_by_sample/
	mv ${outprefix}_chr${i}_bp*.segment_info ./Imputed/info/
	mv ${outprefix}_chr${i}_bp*.segment_warnings ./Imputed/warnings/
	mv ${outprefix}_chr${i}_bp*.segment_summary ./Imputed/summary/

	# cat all output segments into one .impute2 file, tar.gz it and remove the segments
	ls ${outprefix}_chr${i}_bp*.segment | xargs cat > ${outprefix}_chr${i}.impute2
	tar -cvzf ./Imputed/chr${i}.tar.gz ${outprefix}_chr${i}.impute2
	rm ${outprefix}_chr${i}_bp*.segment

#done
