package bcftoolSNP;
my $LEVEL = 1;
sub callSNPs {

use File::Path qw( make_path );
my @arg = @_;
foreach $ar (@arg){
        if ($ar=~ /^-i/){our $inputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-o/){our $outputfolder= (split(/ /,$ar))[1];}
        elsif ($ar=~ /^-nc/){our $nc=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-rg/){$refgenome=(split(/ /,$ar))[1];}
}
#while (@ARGV){
#        $_=shift @ARGV;
#        if ($_=~ /^-i$/){$inputfolder=shift @ARGV;}
#        elsif ($_=~ /^-o$/){$outputfolder=shift @ARGV;}
#        elsif ($_=~ /^-rg$/){$refgenome=shift @ARGV;}
#        elsif ($_=~ /^-nc$/){$nc=shift @ARGV;}
#my $inputfolder = $ARGV[0];
#my $outputfolder = $ARGV[1];
#my $refgenome = $ARGV[2];
#my $ncores =$ARGV[3];
#}
if (not defined ($inputfolder && $outputfolder)){print "\nThis script will create the bcftools command to call SNPs for several samples in parallel. Its need the mapped bam files with realigned indels in a folder.\n\nUsage:\nMPbcftoolscallSNPs\n\t-i <input folder with mapped bam files with realigned indels>\n\t-o <output folder to save vcf files>\n\t-rg <reference genome>\n\t-nc <number the cores to be used in parallel, recomend to use the # of Chromosomes, default 20>\n\nFor example:\nbcftoolscallSNPs -i  ./Yuma/indelrealigned/ -o ./Yuma/rawvcf/ -rg ./Yuma/genome/reference_genome.fasta -nc 23\n\n"; exit;}
if (not defined ($nc)){$nc = 20;}
opendir(DIR, "$inputfolder") or die "Can not open folder, $!\n" ;
our @files = readdir(DIR);
closedir (DIR);
if ( !-d $outputfolder ) {
    make_path $outputfolder or die "Failed to create path: $outputfolder";
}
our @samplesnames=();
our @pops=();
our @scaff=`grep \">\" $refgenome | perl -p -e \'s/>//\' | perl -p -e \'s/\\s\.*\\n/\\n/\'`;
chomp @scaff;
our @cores= (1..$nc);
our @scaffr=();
#print "$scaff[0]\t$cores[0]\n";
@lines=`parallel --link -j 5 echo \"{1}\"\$\'\\t\'\"{2}\" ::: @scaff ::: @cores`;
our %comb;
#@lines=`cat combinations`;
foreach $line (@lines){
        chomp $line;
        @inf=split /\s/,$line;
        $comb{$inf[1]}{$inf[0]}=1;
#       print "$inf[1] es $inf[0]\n";
}
for (my $count=1;$count<=$nc;$count++){
        my $mix="";
        foreach my $key (keys %{$comb{$count}}){
                $mix .= $key.",";
}
        $mix =~ s/,$//;
        push @scaffr, $mix;
}

foreach my $file (@files){ next unless ($file =~ /\.bam$/ && $file =~ /realigned/); my @fileinf = split (/\./, $file);  push @samplesnames, $fileinf[0];}
#foreach my $name (@samplesnames){my $forward = "$inputfolder\/$name\.trim\.1\.fq\.gz"; my $reverse = "$inputfolder\/$name\.trim\.2\.fq\.gz"; my $output = " $outputfolder\/$name\.sam"; print "$name\t$forward\t$reverse\t$output\n";}

my $cmd1 ="parallel --link -j $nc bcftools mpileup -Ou -r {1} -f $refgenome "; 
foreach $samplename (@samplesnames){
$cmd1=$cmd1."$inputfolder\/$samplename\.bam ";
}
$cmd1=$cmd1."-d 1000 -a \"FORMAT/DP,FORMAT/AD\" \'|\' bcftools call -vmO z -o $outputfolder/call.{2}.vcf.gz ::: @scaffr ::: @cores";
#print "$cmd1\n";
#`echo $cmd1 \>\> ./TestYuma`; 
system ($cmd1);

#print "$count1\t$count2\n";
`parallel -j $nc bcftools index $outputfolder/call.{1}.vcf.gz ::: @cores`;
`bcftools concat -a -O z --rm-dups all  $outputfolder/call.*.vcf.gz -o $outputfolder/all.vcf.gz`;
`tabix -p vcf $outputfolder/all.vcf.gz`;
#print "$cmd\n";
#system ($cmd);

}

1;
