load_samtools_1.9

cd /safs-data02/dennislab/Public/Collaborators/Phipps/all_genomes/milk_mags/variants/

for k in `ls *.sam`; do
  #loop through bins
  cd ../
  for i in `ls *.fasta.gz`; do
    #move back to sams
    cd variants/
    #Pull @SQ lines
    grep '@SQ' ${k} | grep "${i%.fasta.gz}" | sed "s/${i%.fasta.gz}_contig_//g" >> ${i%.fasta.gz}_${k#all_milk_mags_}

    #Pull @PG line and edit to genome
    grep '@PG' ${k} | sed "s/all_milk_mags/${i%.fasta.gz}/g" >> ${i%.fasta.gz}_${k#all_milk_mags_}

    #Pull remaining genome lines
    grep -v '@SQ' ${k} | grep -v '@PG' | grep "$(printf '\t')${i%.fasta.gz}" | sed "s/\t${i%.fasta.gz}_contig_/\t/g" | sed '/\tHF0/d' | sed '/\tisolate/d' >> ${i%.fasta.gz}_${k#all_milk_mags_}

    #move back to bins
    cd ../
  done
  #remove old sam
  cd variants/
  rm ${k}
done
