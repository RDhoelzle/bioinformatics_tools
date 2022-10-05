cd /safs-data02/dennislab/Public/Collaborators/Phipps/all_genomes/milk_mags/

#loop through bins
for i in `ls *.fasta.gz`; do
  #move back to sams
  cd variants/plots_all_milk_mags/
  #Pull header lines
  zcat all_milk_mags.vcf.gz | head -5 >> ${i%.fasta.gz}.vcf

  #Pull contig ids
  zcat all_milk_mags.vcf.gz | grep "##contig" | grep ${i%.fasta.gz} >> ${i%.fasta.gz}.vcf

  #Pull remaining header lines
  zcat all_milk_mags.vcf.gz | grep "#" | grep -v "contig" | tail -n +6 >> ${i%.fasta.gz}.vcf

  #Pull contig lines
  zcat all_milk_mags.vcf.gz | grep -v "#" | grep ${i%.fasta.gz} >> ${i%.fasta.gz}.vcf

  #move back to bins
  cd ../../
done
#move back to variants directory
cd variants/plots_all_milk_mags/
