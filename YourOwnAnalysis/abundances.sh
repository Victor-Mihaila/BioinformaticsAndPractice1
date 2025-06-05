stringtie --merge -p 8 -G ./genes/gencode.gtf -o ./stringtie/stringtie_merged.gtf mergelist.txt
stringtie -e -B -p 8 -G ./stringtie/stringtie_merged.gtf -o ballgown/RNA-control/RNA-control.gtf ./binfo1-datapack1/RNA-control.bam
stringtie -e -B -p 8 -G ./stringtie/stringtie_merged.gtf -o ballgown/RNA-siLin28a/RNA-siLin28a.gtf ./binfo1-datapack1/RNA-siLin28a.bam
stringtie -e -B -p 8 -G ./stringtie/stringtie_merged.gtf -o ballgown/RNA-siLuc/RNA-siLuc.gtf ./binfo1-datapack1/RNA-siLuc.bam
stringtie -e -B -p 8 -G ./stringtie/stringtie_merged.gtf -o ballgown/RPF-siLin28a/RPF-siLin28a.gtf ./binfo1-datapack1/RPF-siLin28a.bam
stringtie -e -B -p 8 -G ./stringtie/stringtie_merged.gtf -o ballgown/RPF-siLuc/RPF-siLuc.gtf ./binfo1-datapack1/RPF-siLuc.bam
stringtie -e -B -p 8 -G ./stringtie/stringtie_merged.gtf -o ballgown/CLIP-35L33G/CLIP-35L33G.gtf ./binfo1-datapack1/CLIP-35L33G.bam

python prepDE.py3 -i prep_deseq.txt