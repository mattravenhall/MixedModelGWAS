#########################################################################################
# EMMAX Pipeline - Matt Ravenhall		    					#
# - Start with a 12-coded tped, masking is automatic. 					#
# - Note that both the kinship matrix and index will be based on the default tped.	#
# - EMMAX can only apply the additive model, masking allows for Het, Dom & Rec models. 	#
#########################################################################################

## tped > model masking > kinship matrix > emmax > index > merge & plot

CHROMOSOME=$1
SUBTYPE=$2

EMMAXdir=./
TPED=./PLACEHOLDER_chr${CHROMOSOME}		# NB. TPED must be 12 format
PHENO=./PLACEHOLDER_${SUBTYPE}.pheno	# Pre-made phenotype file
COVARS=./PLACEHOLDER.covars				# Pre-made covariates file
OUTDIR=./								# Output directory
OUTPREFIX=PLACEHOLDER_${SUBTYPE}		# Output prefix
KINF=./PLACEHOLDER.kinf 				# Kinship matrix
SCRIPTDIR='./' 							# Supporting scripts directory

# Mask default input files to return Heterozygous, Recessive and Dominant model Ps when under the Additive model (as with EMMAX)
# Take SNP call count for masking loops
limit=$(head -1 ${TPED}.tped | wc -w)

# Mask for heterozygous (MM to mm, Mm stays Mm, mm stays mm) # Note that missing/Major calls have become missing/minor calls here.
echo 'Masking for heterozygous model'
cat ${TPED}.tped | awk -v limit=${limit} '{printf $1" "$2" "$3" "$4; for(i=5;i<limit;i=i+2) { j=i+1; if($i==2 && $j==2) {printf " 1 1";} else if($i==1 && $j==1) {printf " 1 1";} else if($i==1 && $j==2) {printf " 1 2";} else if($i==2 && $j==1) {printf " 2 1";} else if($i==0 && $j==0) {printf " 0 0";} else if($i==0 && $j==1) {printf " 0 1";} else if($i==0 && $j==2) {printf " 0 1";} else if($i==1 && $j==0) {printf " 1 0";} else if($i==2 && $j==0) {printf " 1 0";}} print""}' > ${OUTDIR}${OUTPREFIX}_Het.tped
cp ${TPED}.tfam ${OUTDIR}${OUTPREFIX}_Het.tfam

# Mask for recessive (MM to Mm, Mm to mm, mm stays mm)
echo 'Masking for recessive model'
cat ${TPED}.tped | awk -v limit=${limit} '{printf $1" "$2" "$3" "$4; for(i=5;i<limit;i=i+2) { j=i+1; if($i==2 && $j==2) {printf " 1 2";} else if($i==1 && $j==1) {printf " 1 1";} else if($i==2 && $j==1) {printf " 1 1";} else if($i==1 && $j==2) {printf " 1 1";} else if($i==0 && $j==0) {printf " 0 0";} else if($i==0 && $j==1) {printf " 0 1";} else if($i==0 && $j==2) {printf " 0 2";} else if($i==1 && $j==0) {printf " 1 0";} else if($i==2 && $j==0) {printf " 2 0";}} print""}' > ${OUTDIR}${OUTPREFIX}_Rec.tped
cp ${TPED}.tfam ${OUTDIR}${OUTPREFIX}_Rec.tfam

# Mask for dominant (MM to Mm, Mm stays Mm, mm stays mm)
echo 'Masking for dominant model'
cat ${TPED}.tped | awk -v limit=${limit} '{printf $1" "$2" "$3" "$4; for(i=5;i<limit;i=i+2) { j=i+1; if($i==2 && $j==2) {printf " 1 2";} else if($i==1 && $j==1) {printf " 1 1";} else if($i==1 && $j==2) {printf " 1 2";} else if($i==2 && $j==1) {printf " 2 1";} else if($i==0 && $j==0) {printf " 0 0";} else if($i==0 && $j==1) {printf " 0 1";} else if($i==0 && $j==2) {printf " 0 2";} else if($i==1 && $j==0) {printf " 1 0";} else if($i==2 && $j==0) {printf " 2 0";}} print""}' > ${OUTDIR}${OUTPREFIX}_Dom.tped
cp ${TPED}.tfam ${OUTDIR}${OUTPREFIX}_Dom.tfam

# Create Kinship Matrix, NB. kinship matrix is made from the default tped NOT the masked ones.
#echo 'Creating kinship matrix'
#${EMMAXdir}/emmax-kin -v -h -d 10 $TPED

# Run EMMAX
echo 'Running EMMAX under Additive model'
${EMMAXdir}/emmax -v -d 10 -t ${TPED} -p ${PHENO} -k ${KINF} -c ${COVARS} -o ${OUTDIR}${OUTPREFIX}_Add
echo 'Running EMMAX under Heterozygous model'
${EMMAXdir}/emmax -v -d 10 -t ${OUTDIR}${OUTPREFIX}_Het -p ${PHENO} -k ${KINF} -c ${COVARS} -o ${OUTDIR}${OUTPREFIX}_Het
echo 'Running EMMAX under Dominant model'
${EMMAXdir}/emmax -v -d 10 -t ${OUTDIR}${OUTPREFIX}_Dom -p ${PHENO} -k ${KINF} -c ${COVARS} -o ${OUTDIR}${OUTPREFIX}_Dom
echo 'Running EMMAX under Recessive model'
${EMMAXdir}/emmax -v -d 10 -t ${OUTDIR}${OUTPREFIX}_Rec -p ${PHENO} -k ${KINF} -c ${COVARS} -o ${OUTDIR}${OUTPREFIX}_Rec

## Convert EMMAX input to a Manhattan Plot
# Create Index
echo 'Creating SNP Index'
cut -f1,2,4 -d' ' ${TPED}.tped > ${OUTPREFIX}_MergeIndex.tmp

# Merge and Plot .ps
echo 'Identifying and plotting most significant Ps'
Rscript ${SCRIPTDIR}mergePlotEMMAX.r ${OUTDIR}${OUTPREFIX} # NB. If per-chr is passed, a duplication plot will be created

# Sort plots into output folder
echo 'Moving files to output directory'
if [ ! -d "$EMMAX_Output" ]; then
	mkdir ${OUTDIR}EMMAX_Output
fi

# Remove temporary files
rm ${OUTPREFIX}_MergeIndex.tmp

# Move output to final directory
mv ${OUTDIR}${OUTPREFIX}*_EMMAX_MPlot.png ${OUTDIR}EMMAX_Output/
mv ${OUTDIR}${OUTPREFIX}*_EMMAX_QQplot.png ${OUTDIR}EMMAX_Output/
mv ${OUTDIR}${OUTPREFIX}.ps ${OUTDIR}EMMAX_Output/
mv ${OUTDIR}${OUTPREFIX}*.log ${OUTDIR}EMMAX_Output/
mv ${OUTDIR}${OUTPREFIX}*.reml ${OUTDIR}EMMAX_Output/
mv ${OUTDIR}${OUTPREFIX}*.tped ${OUTDIR}EMMAX_Output/
mv ${OUTDIR}${OUTPREFIX}*.tfam ${OUTDIR}EMMAX_Output/
