package indexgenome;
my $LEVEL = 1;
sub indrg {

my @arg = @_;
foreach $ar (@arg){
        if ($ar=~ /^-i/){$inputfolder=(split(/ /,$ar))[1];}
        elsif ($ar=~ /^-rg/){$ref=(split(/ /,$ar))[1];}
	elsif ($ar=~ /^-pf/){$pf=(split(/ /,$ar))[1];}
}
if (not defined ($inputfolder && $ref)){print "\nThis script will indexs a reference genome using samtools, picard, bowtie2 and snap. The genome should have extention fna.\n\nUsage:\n\t-i <input folder where the reference is, and where all the indeces and dictionaries will be saved>\n\t-rg <the name of the reference genome incluiding the extention>\n\t-pf <path to the picar jar executable>\n\nExample:SNPcallpipe.pl -stp indref -i yume/genomes/ -rg Taust.Chr.fna -pf /local/SNOcallPipe/\n"; exit;}

$ref =~ s/$inputfolder//;
our $orgn = $ref;
unless ($ref =~ m/\.fna$/){
	$orgn =~ s/\.fasta$//;
	$orgn =~ s/\.fa$//;
	`mv $inputfolder\/$ref $inputfolder\/$orgn\.fna`; 
}
$orgn =~ s/\.fna$//;
$OUTSNAP ="$inputfolder\/$orgn\_SNAP";
`mkdir $OUTSNAP`;
our @cmds = ();
$cmds[1] = "samtools faidx $inputfolder\/$orgn\.fna";
$cmds[2] = "java -Xmx4096m -jar $pf\/picard.jar CreateSequenceDictionary R=$inputfolder\/$orgn\.fna O=$inputfolder\/$orgn.dict";
$cmds[3] = "bowtie2-build  $inputfolder\/$orgn\.fna $inputfolder\/$orgn\.fna";
$cmds[4] = "snap-aligner index $inputfolder\/$orgn\.fna $OUTSNAP";

use Parallel::Loops;
my $plraw = Parallel::Loops->new(4);
$plraw->foreach (\@cmds, sub{
my $cmd = $_ ;
print ("$cmd\n");
`$cmd`;
#`samtools faidx $inputfolder\/$orgn\.fna`;
#`java -Xmx4096m -jar $pf\/picard.jar CreateSequenceDictionary R=$inputfolder\/$orgn\.fna O=$inputfolder\/$orgn.dict`;
#`bowtie2-build  $inputfolder\/$orgn\.fna $inputfolder\/$orgn\.fna`;
#`bwa index $inputfolder\/$orgn\.fna $inputfolder\/$orgn\.fna`;
#`dragen-os --build-hash-table true --ht-reference $inputfolder\/$orgn\.fna --output-directory $OUTSNAP`;
});
}
1;
