threads=23  # set as appropriate for your system (no more than 128)
genome_size=4224919

mkdir SRL662_assemblies

for i in 01 02 03 04; do
    /home/nik_arapitsas/Documents/Bacillus_project/scripts/flye.sh SRL662_subsampled_reads/sample_"$i".fastq SRL662_assemblies/SRL662_flye_"$i" "$threads" "$genome_size"
done
