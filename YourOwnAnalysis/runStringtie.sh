stringtie -p 8 -G genes/gencode.gtf -o stringtie/RNA-control.gtf -l RNA-control ./binfo1-datapack1/RNA-control.bam
stringtie -p 8 -G genes/gencode.gtf -o stringtie/RNA-siLin28a.gtf -l RNA-siLin28a ./binfo1-datapack1/RNA-siLin28a.bam
stringtie -p 8 -G genes/gencode.gtf -o stringtie/RNA-siLuc.gtf -l RNA-siLuc ./binfo1-datapack1/RNA-siLuc.bam
stringtie -p 8 -G genes/gencode.gtf -o stringtie/RPF-siLin28a.gtf -l RPF-siLin28a ./binfo1-datapack1/RPF-siLin28a.bam
stringtie -p 8 -G genes/gencode.gtf -o stringtie/RPF-siLuc.gtf -l RPF-siLuc ./binfo1-datapack1/RPF-siLuc.bam
stringtie -p 8 -G genes/gencode.gtf -o stringtie/CLIP-35L33G.gtf -l CLIP-35L33G ./binfo1-datapack1/CLIP-35L33G.bam