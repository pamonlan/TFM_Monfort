BAM=$1
BEDPATHS=$2

BASE=$(basename $BAM .bam)
mkdir /storage/tbc/melanoma/Simulacion_artificial/${BASE}

for dir in $(ls $BEDPATHS)
do
	echo $dir
	mkdir /storage/tbc/melanoma/Simulacion_artificial/${BASE}/${dir}
	mkdir /storage/tbc/melanoma/Simulacion_artificial/${BASE}/${dir}/bam
	mkdir /storage/tbc/melanoma/Simulacion_artificial/${BASE}/${dir}/TruthVcf/
        mkdir /storage/tbc/melanoma/Simulacion_artificial/${BASE}/${dir}/TumorVcf/

	for bed in $(ls ${BEDPATHS}/${dir})
	do
		BED=$(realpath ${BEDPATHS}/${dir}/${bed})
		echo $bed
		sh /storage/tbc/melanoma/script_sbatch/AddSnv.1.1.batch  ${BAM} ${BED} /storage/tbc/melanoma/Simulacion_artificial/${BASE}/${dir}
done
done
