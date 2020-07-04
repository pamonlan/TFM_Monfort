#!bin/bash
BED=$1

for i in 10 20 40 60 80 
do         
	SNV_PATH="snv."${i}
        echo $SNV_PATH
        mkdir $SNV_PATH
        
	for j in {1..10}         
	do           
		FILE="snv."${i}"."${j}".bed"                 
		echo $FILE
		
		~/scripts/CreateMutationSnv ${i} ${BED}
		mv "snv_"${i}".bed" ${SNV_PATH}/${FILE}          

done 
done

