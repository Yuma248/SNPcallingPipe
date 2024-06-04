package BWAB2;

my $LEVEL = 1;
sub align{
my @arg = @_;
foreach $ar (@arg){
	if ($ar =~ /^-rg/){our $refgenom = (split(/ /,$ar))[1];}
	elsif ($ar =~ /^-i/){our $inputfolder = (split(/ /,$ar))[1];}
        elsif($ar=~ /^-o/){our $outputfolder = (split(/ /, $ar))[1];}
        elsif ($ar=~ /^-ncp/){$ncores=(split(/ /, $ar))[1];}
        elsif ($ar=~ /^-ncr/){$nc=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-al/){$B=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-t/){$type=(split(/ /, $ar))[1];}
}
use File::Path qw( make_path );
if (not defined ($refgenom && $inputfolder && $outputfolder)){print "\nThis script will map reads of several samples to a reference genome using bwa-mem, bowtie2 or snap-aligner  in parallel. It requires files (after demultiplexing and trimming) of all samples storaged in the same folder\n\nUsage:\nBWAparallel.pl\n\t-rg <path to the reference genome fasta file>\n\t-i <path to the folder with the input fasta files>\n\t-o <path to the output folder>\n Optional:\n\t-ncp <number of runs in parallel, default 1>\n\t-nc <number cores for each run, default 4>\n\t-al <aligner to be used, BWA, B2 for bowtie2, or SNAP,  default BWA>\n\t-t <sequencing type single-end \"S\" or paired-end \"P\", default P> \n\nFor example:\nBWAparallel.pl -rg /home/Yumafan/genome/reference_genome.fasta -i /home/Yumafan/demultiplex/trimmed/ -o /home/Yumafan/bwaout/ -ncp 10 -ncr 4 -al BWA -t P\n\n"; exit;}
if ( !-d $outputfolder ) {
    make_path $outputfolder or die "Failed to create path: $outputfolder";
}
my $mergefolder="$outputfolder\/merged";
if ( !-d $mergefolder ) {
    make_path $mergefolder or die "Failed to create path: $mergefolder";
}

if (not defined ($ncp)){$ncp=1;}
if (not defined ($nc)){$nc=4;}
if (not defined ($B)){$B="BWA";}
if (not defined ($type)){$type="P";}
opendir(DIR, "$inputfolder") or die "Can not open folder, $!\n" ;
my @files = readdir(DIR);
closedir (DIR);
our @samplesnames=();
our @names=();
if ($type eq "P"){
#foreach my $file (@files){ next unless ($file =~ /\.pair1\.truncated\.gz$/); my @fileinf = split (/\./, $file); push @samplesnames, $fileinf[0]; my @nameinf = split (/\_/, $fileinf[0]); push @names, $nameinf[0];}
foreach my $file (@files){ next unless ($file =~ /\.pair1\.truncated\.gz$/); my @fileinf = split (/\./, $file); push @samplesnames, $fileinf[0]; push @names, $fileinf[0];}
#foreach my $file (@files){ next unless ($file =~ /\.trim\.1\.fq\.gz$/); my @fileinf = split (/\./, $file); push @samplesnames, $fileinf[0]; @names = @samplesnames;}
#foreach my $name (@samplesnames){my $forward = "$inputfolder\/$name\.trim\.1\.fq\.gz"; my $reverse = "$inputfolder\/$name\.trim\.2\.fq\.gz"; my $output = " $outputfolder\/$name\.sam"; print "$name\t$forward\t$reverse\t$output\n";}
if ($B eq "BWA"){
my $cmd ="parallel -j $ncores --link --results $outputfolder\/logbwa1 --noswap \"bwa mem -t $nc $refgenom $inputfolder\/{1}\.pair1\.truncated\.gz $inputfolder\/{1}\.pair2\.truncated\.gz >> $outputfolder\/{2}\.pairs.sam\" ::: @samplesnames ::: @names";
my $cmd2 ="parallel -j $ncores --link --results $outputfolder\/logbwa2 --noswap \"bwa mem -t $nc $refgenom $inputfolder\/{1}\.collapsed\.truncated\.gz  >> $outputfolder\/{2}\.collapsed\.truncated\.sam\" ::: @samplesnames ::: @names";
my $cmd3 ="parallel -j $ncores --link --results $outputfolder\/logbwa3 --noswap \"bwa mem -t $nc $refgenom $inputfolder\/{1}\.collapsed\.gz >> $outputfolder\/{2}\.collapsed\.sam\" ::: @samplesnames ::: @names";
my $cmd4 ="parallel -j $ncores --link --results $outputfolder\/logbwa4 --noswap \"bwa mem -t $nc $refgenom $inputfolder\/{1}\.singleton\.truncated\.gz >> $outputfolder\/{2}\.singleton\.truncated\.sam\" ::: @samplesnames ::: @names";
my $cmd5 ="parallel -j $ncores --link --results $outputfolder\/logmerge --noswap samtools merge $mergefolder/{1}.bam $outputfolder/{1}.pairs.sam $outputfolder/{1}.collapsed.truncated.sam $outputfolder/{1}.collapsed.sam $outputfolder/{1}.singleton.truncated.sam -O BAM ::: @names";
system ($cmd);
system ($cmd2);
system ($cmd3);
system ($cmd4);
system ($cmd5);
}
elsif ($B eq "B2"){
my $cmd ="parallel -j $ncores --link --results $outputfolder\/logbowtie1 --noswap bowtie2 -p $nc -x $refgenom -1 $inputfolder\/{1}\.pair1\.truncated\.gz -2 $inputfolder\/{1}\.pair2\.truncated\.gz -S $outputfolder\/{2}\.pairs.sam --no-contain -X 1000 --rg-id {1} --rg SM:{1} --rg LB:library --rg PU:lane ::: @samplesnames ::: @names";
my $cmd2 ="parallel -j $ncores --link --results $outputfolder\/logbowtie2 --noswap bowtie2 -p $nc -x $refgenom -U $inputfolder\/{1}\.collapsed\.truncated\.gz -S $outputfolder\/{2}\.collapsed\.truncated\.sam --rg-id {1} --rg SM:{1} --rg LB:library --rg PL:ILLUMINA --rg PU:lane ::: @samplesnames ::: @names";
my $cmd3 ="parallel -j $ncores --link --results $outputfolder\/logbowtie3 --noswap bowtie2 -p $nc -x $refgenom -U $inputfolder\/{1}\.collapsed\.gz -S $outputfolder\/{2}\.collapsed\.sam --rg-id {1} --rg SM:{1} --rg LB:library --rg PL:ILLUMINA --rg PU:lane ::: @samplesnames ::: @names";
my $cmd4 ="parallel -j $ncores --link --results $outputfolder\/logbowtie4 --noswap bowtie2 -p $nc -x $refgenom -U $inputfolder\/{1}\.singleton\.truncated\.gz -S $outputfolder\/{2}\.singleton\.truncated\.sam --rg-id {1} --rg SM:{1} --rg LB:library --rg PL:ILLUMINA --rg PU:lane ::: @samplesnames ::: @names";
my $cmd5 ="parallel -j $ncores --link --results $outputfolder\/logmerge --noswap samtools merge $mergefolder/{1}.bam $outputfolder/{1}.pairs.sam $outputfolder/{1}.collapsed.truncated.sam $outputfolder/{1}.collapsed.sam $outputfolder/{1}.singleton.truncated.sam -O BAM ::: @names";
system ($cmd);
system ($cmd2);
system ($cmd3);
system ($cmd4);
system ($cmd5);
}
elsif ($B eq "SNAP"){
my $cmd ="parallel -j $ncores --link --results $outputfolder\/logsnap1 --noswap snap-aligner paired $refgenom -compressedFastq $inputfolder\/{1}\.pair1\.truncated\.gz $inputfolder\/{1}\.pair2\.truncated\.gz -so -o $outputfolder\/{2}\.pairs\.bam -t $nc -rg {1} -R ID\:{1}\\\\\\\\tPL\:Illumina\\\\\\\\tPU\:pu\\\\\\\\tLB\:lb\\\\\\\\tSM\:{1} ::: @samplesnames ::: @names";
my $cmd2 ="parallel -j $ncores --link --results $outputfolder\/logsnap2 --noswap snap-aligner single $refgenom -compressedFastq $inputfolder\/{1}\.collapsed\.truncated\.gz -so -o $outputfolder\/{2}\.collapsed\.truncated\.bam -t $nc -rg {1} -R ID\:{1}\\\\\\\\tPL\:Illumina\\\\\\\\tPU\:pu\\\\\\\\tLB\:lb\\\\\\\\tSM\:{1} ::: @samplesnames ::: @names";
my $cmd3 ="parallel -j $ncores --link --results $outputfolder\/logsnap3 --noswap snap-aligner single $refgenom -compressedFastq $inputfolder\/{1}\.collapsed\.gz -so -o  $outputfolder\/{2}\.collapsed\.bam -t $nc -rg {1} -R ID\:{1}\\\\\\\\tPL\:Illumina\\\\\\\\tPU\:pu\\\\\\\\tLB\:lb\\\\\\\\tSM\:{1} ::: @samplesnames ::: @names";
my $cmd4 ="parallel -j $ncores --link --results $outputfolder\/logsnap4 --noswap snap-aligner single $refgenom -compressedFastq $inputfolder\/{1}\.singleton\.truncated\.gz -so -o $outputfolder\/{2}\.singleton\.truncated\.bam -t $nc -rg {1} -R ID\:{1}\\\\\\tPL\:Illumina\\\\\\tPU\:pu\\\\\\tLB\:lb\\\\\\tSM\:{1} ::: @samplesnames ::: @names";
my $cmd5 ="parallel -j $ncores --link --results $outputfolder\/logmerge --noswap samtools merge $mergefolder/{1}.bam $outputfolder/{1}.pairs.bam $outputfolder/{1}.collapsed.truncated.bam $outputfolder/{1}.collapsed.bam $outputfolder/{1}.singleton.truncated.bam -O BAM ::: @names";
system ($cmd);
system ($cmd2);
system ($cmd3);
system ($cmd4);
system ($cmd5);
}
else {print "Please use a propert aligner description BWA, B2 (for bowtie2), or SNAP\n\n";}
}
elsif ($type eq "S"){
foreach my $file (@files){ next unless ($file =~ /\.fastq\.gz$/ || $file =~ /\.fq\.gz$/); my @fileinf = split (/\./, $file); push @samplesnames, $fileinf[0]; my @nameinf = split (/\_/, $fileinf[0]); push @names, $nameinf[0]; if ($file =~ /\.fastq\.gz$/){our $format="fastq";}else{our $format="fq";}}
if ($B eq "BWA"){
if ($format eq "fq"){
our  $cmd ="parallel -j $ncores --link --results $outputfolder\/logbwaSE --noswap \"bwa mem -t $nc $refgenom $inputfolder\/{1}\.fq\.gz >> $outputfolder\/{1}\.sam\" ::: @names";
}
elsif ($format eq "fastq"){
our $cmd ="parallel -j $ncores --link --results $outputfolder\/logbwaSE --noswap \"bwa mem -t $nc $refgenom $inputfolder\/{1}\.fastq\.gz >> $outputfolder\/{1}\.sam\" ::: @names";
}
system ($cmd);
}
elsif ($B eq "B2"){
if ($format eq "fq"){
our $cmd = "parallel -j $ncores --link --results $outputfolder\/logbowtieSE --noswap bowtie2 -p $nc -x $refgenom -U $inputfolder\/{1}\.fq\.gz  -S $outputfolder\/{1}.sam --rg-id {1} --rg SM:{1} --rg LB:library --rg PL:ILLUMINA --rg PU:lane ::: @samplesnames ::: @names";
}
elsif ($format eq "fastq"){
our $cmd = "parallel -j $ncores --link --results $outputfolder\/logbowtieSE --noswap bowtie2 -p $nc -x $refgenom -U $inputfolder\/{1}\.fastq\.gz  -S $outputfolder\/{1}.sam--rg-id {1} --rg SM:{1} --rg LB:library --rg PL:ILLUMINA --rg PU:lane ::: @samplesnames ::: @names";
}
system ($cmd);
}
elsif ($B eq "SNAP"){
if ($format eq "fq"){
our $cmd = "parallel -j $ncores --link --results $outputfolder\/logbowtieSE --noswap snap-aligner single $refgenom -compressedFastq $inputfolder\/{1}\.fq\.gz  -so -s  $outputfolder\/{1}.bam -t $nc ::: @samplesnames ::: @names";
}
elsif ($format eq "fastq"){
our $cmd = "parallel -j $ncores --link --results $outputfolder\/logbowtieSE --noswap snap-aligner single $refgenom -compressedFastq $inputfolder\/{1}\.fastq\.gz  -S $outputfolder\/{1}.sam -t $nc ::: @samplesnames ::: @names";
}
system ($cmd);
}

}

}
1;
