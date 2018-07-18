seq.check.pl -i Athaliana_167_TAIR9.fa -o 11_genome.fas
gff2gtb.pl -i Athaliana_167_TAIR10.gene.gff3 -o 41.gtb
gtb.phytozome.pl -i 41.gtb -o 42.gtb
gtb.pickalt.pl -i 42.gtb -o 51.gtb

gtb2gff.pl -i 51.gtb -o 51.gff
gtb2fas.pl -i 51.gtb -d 11_genome.fas -o 51.fas
