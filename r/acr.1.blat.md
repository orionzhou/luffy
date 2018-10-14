## Maize ACR analysis

1. Create a working directory and switch to it
```bash
mkdir -p /home/springer/zhoux379/data/misc1/maize.acr
cd /home/springer/zhoux379/data/misc1/maize.acr
```

2. Remove header and take only first 3 columns to create a BED file
```bash
cut -f1-3 /home/springer/nosha003/wgbs_schmitz/ACR/B73_peaks/B73L_final_ACR.bed > 01.bed
sed -i '1d' 01.bed
```

3. Extract sequences using the coordinates in the BED file from reference genome
```bash
seqret.pl -d /home/springer/zhoux379/data/genome/Zmays_v4/11_genome.fas -b 01.bed -o 02.fas
```
  * [link to seqret.pl](https://github.com/orionzhou/luffy/blob/master/perl/seqret.pl)

4. Load the BLAT toolkit (on MSI)
```bash
module load blat/34
```

5. Build genome index for BLAT, this step is optional since BLAT can also take a fasta file directly
```bash
faToTwoBit Zmays_v4.fasta Zmays_v4.2bit
blat Zmays_v4.2bit tmp.fas tmp.out -makeOoc=Zmays_v4.2bit.tile11.ooc
```
  * [link to faToTwoBit](https://genome.ucsc.edu/goldenpath/help/blatSpec.html#faToTwoBitUsage)
  * The \*.ooc file contains over-occurring 11-mers and will increase the speed of BLAT by a factor of 40 in many cases, but is not required.

6. Blat against the (already indexed) PH207 genome
```bash
blat /home/springer/zhoux379/data/genome/PH207/21.blat/db.2bit -ooc=/home/springer/zhoux379/data/genome/PH207/21.blat/db.2bit.tile11.ooc 02.fas 04.psl
```
  * [link to blat usage](https://genome.ucsc.edu/goldenpath/help/blatSpec.html#blatUsage)
  * blat output is a 21-column text file in [PSL format](https://useast.ensembl.org/info/website/upload/psl.html)

7. Convert Blat PSL output to more human-readable tab-separated file
```bash
psl2tsv.pl -i 04.psl -o 05.tsv
```
  * [link to psl2tsv.pl](https://github.com/orionzhou/luffy/blob/master/perl/psl2tsv.pl)

8. The Blat step is computationally intensive so we can put it in a PBS script and submit to job queue:
```bash
#PBS -l nodes=1:ppn=1,walltime=10:00:00
#PBS -m ae
#PBS -M zhoux379@umn.edu
#PBS -q small

cd /home/springer/zhoux379/data/misc1/maize.acr
blat /home/springer/zhoux379/data/genome/PH207/21.blat/db.2bit -ooc=/home/springer/zhoux379/data/genome/PH207/21.blat/db.2bit.tile11.ooc 02.fas 04.psl
psl2tsv.pl -i 04.psl -o 05.tsv
```

