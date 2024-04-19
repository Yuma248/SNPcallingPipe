package samdedup;

my $LEVEL = 1;
sub dedup{
my @arg = @_;
foreach $ar (@arg){
	if ($ar=~ /^-i/){our $inf= (split(/ /,$ar))[1];}
	elsif ($ar=~ /^-o/){our $outf= (split(/ /,$ar))[1];}
	elsif ($ar=~ /^-lnc/){our $nc=(split(/ /,$ar))[1];}
	elsif ($ar=~ /^-snc/){our $ncp=(split(/ /,$ar))[1];}
}

if (not defined ($inf)){
	print "\nThis script will convert from sam to bam, sort by name, fixmates, sort by coordinates and markduplicates using samtools in parallel. Requires your sam files after mapping (sample01.sam), all save in one folder.\n\nUsages: samtoolsdup.pl\n\t-i <path to inputfolder>\n\nOptional:\n\t-o <path to outputfolder, default samout>\n\t-nc <number the cores or samples to use in parallel, default 4>\n\t-ncp <number of cores per sample, default 1>\n\n";
 	exit;
}
$outf //= "./samout";
our $tmpdir=$outf."/tmp";
our $names=$outf."/name";
our $fix=$outf."/fix";
our $coor=$outf."/cdnt";
our $ddp=$outf."/dedup";
if ( !-d $outf ) {`mkdir $outf`;}
for my $ofn ($tmpdir, $names, $fix, $coor, $ddp){
 mkdir $ofn unless -d $ofn;
}

$nc //= 4;
$ncp //=1;
our $nct=$nc + $ncp;
my $ext = `ls $inf\/ | tail  -n 1 | awk -F '.' '{print \$NF}'`;
my @names=`ls $inf/*.$ext | sed 's/.$ext//g'`;
foreach $name (@names){chomp $name; $name=~ s/$inf\///g; push (@nms, $name);}
if ($ext eq "sam") {
	our $cmd="parallel -j $nc samtools view -bS $inf/{1}.sam '|' samtools sort -n -@ $ncp -o $names/{1}_namesort.bam - ::: @nms";
}
elsif ($ext eq "bam") {
	our $cmd="parallel -j $nc samtools sort -n -@ $ncp -o $names/{1}_namesort.bam  $inf/{1}.bam ::: @nms";
}

my $cmd2="parallel -j $nc samtools fixmate -m $names/{1}_namesort.bam $fix/{1}_fixmate.bam ::: @nms";
my $cmd3="parallel -j $nc samtools sort -\@ $ncp -o $coor/{1}_positionsort.bam $fix/{1}_fixmate.bam ::: @nms";
my $cmd4="parallel -j $nc samtools markdup -r -s $coor/{1}_positionsort.bam $ddp/{1}_markdup.bam \'\>\' $tmpdir/{1}_log ::: @nms";
print "$cmd\n\n$cmd2\n\n$cmd3\n\n$cmd4\n\n";

if ($ext eq "sam") {
`parallel -j $nc samtools view -bS $inf/{1}.sam \'\|\' samtools sort -n -\@ $ncp -o $names/{1}_namesort.bam - ::: @nms`;
}elsif ($ext eq "sam") {
`parallel -j $nc samtools sort -n -\@ $ncp -o $names/{1}_namesort.bam  $inf/{1}.bam ::: @nms`;
}
`parallel -j $nc samtools fixmate -m $names/{1}_namesort.bam $fix/{1}_fixmate.bam ::: @nms`;
`parallel -j $nc samtools sort -\@ $ncp -o $coor/{1}_positionsort.bam $fix/{1}_fixmate.bam ::: @nms`;
`parallel -j $nc --results $tmpdir samtools markdup -r -s $coor/{1}_positionsort.bam $ddp/{1}_markdup.bam ::: @nms`; 
}

1;
