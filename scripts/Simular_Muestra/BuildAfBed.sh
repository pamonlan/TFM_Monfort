#!bin/bash

BED=$1

for i in 0.6
do         
	SNV_PATH="af."${i}
        echo $SNV_PATH
        mkdir $SNV_PATH
        
	for j in {1..10}         
	do           
		FILE="af."${i}"."${j}".bed"                 
		echo $FILE
		
		~/scripts/CreateMutationSnv 80 ${BED}
		
		cat "snv_80.bed" | awk -v var=$i '{print $1 "\t" $2 "\t" $3 "\t" var}' > ${SNV_PATH}/${FILE} 
		

done 
done

