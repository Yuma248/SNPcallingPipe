package dDocent;

my $LEVEL = 1; 
sub demul{
my @arg = @_;
foreach $ar (@arg){
        if ($ar =~ /^-i/){our $input = (split(/ /,$ar))[1];}
        elsif($ar=~ /^-o/){our $output = (split(/ /, $ar))[1];}
        elsif ($ar=~ /^-bf/){$barcode_file=(split(/ /, $ar))[1];}
        elsif ($ar=~ /^-rad/){$radtag=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-lnc/){$lnc=(split(/ /, $ar))[1];}
	elsif ($ar=~ /^-snc/){$snc=(split(/ /, $ar))[1];}
##        elsif ($ar=~ /^-stp/){(split(/ /, $ar))[1];}
}
if (not defined ($input && $output && $barcode_file)){print "\nThis script uses stacks's process_rad to demultiplex fasta files.\n\nUsage:\ndDocent_process_PE.pl\n\t-i <directory with raw sequencing files>\n\t-o <output folder, if it does no exist it will be created>\n\t-bf <barcode file, tab delimited (LaneName SampleName Barcode Single Popnumber)>\nOptional:\n\t-lnc <number of lanes in parallel, or number of R1 files in you folder. default 1>\n\t-snc <number of samples perl lane in parallel, optimum 58/number of R1 files in you folder. default 10>\n\t-rad <RAD_tag, default TGCAGG TAA>\n\nFor example:\nSNPcallPipe.pl -stp demul -i /yuma/rawread/ -o /yuma/demultiplex -bf /yuma/barcodefile -lnc 1 -snc 10 \n\n"; exit;}
my $format ='';
unless (-d $input){ print "No $input folder exists\n";exit;}
unless (-d $output){`mkdir $output`;}
$outdDocent="$output\/dDocent";
unless (-d $outdDocent){`mkdir $outdDocent`;}
unless (-f $barcode_file){print "$barcode_file file exists\n";exit;}
if (not defined $radtag){$radtag = "TGCAGG\tTAA\tsample";} ## default RAD-tag
if (not defined $lnc){$nc=1;}
if (not defined $snc){$nc=10;}
if (not defined $step){$step=1;}
#my @inputfiles = glob("$input/*");

open(FILE, $barcode_file);
my @DATA = <FILE>;
close FILE;

my $trace=0;
open(POPMAP, ">popmapseq");
foreach my $line (@DATA){
chomp ($line);
my @try = split(/\t/,$line);
#my $test=scalar @try;
#print "$test\n";
if ($#try !=4){print "check format of $barcode_file file\nData should be in this order:\nRaw_sequencing_file\tsample_name\tBarcode\ttype(Parent or Single or Progeney)\tpop\n"; exit;}
my ($raw_file,$sample_name, $barcode, $type,$popmap) = split(/\t/,$line);
$hash_raw{$raw_file}=1;
$hash{$sample_name} = $type;
$hash_pop{$sample_name}=$popmap;
$hash{$raw_file}{traceno}{$trace} = $barcode;
$hash{$raw_file}{code}{$barcode}=$sample_name;
$trace++;
print POPMAP "$sample_name\_1.RAD\t$popmap\n";
print POPMAP "$sample_name\_2.RAD\t$popmap\n";
}
close BARCODE;
open (RADFILE, '>RAD.txt');
print RADFILE $radtag;
close RADFILE;

`mkdir $output/log`;

	my $options="";



my @RAW=(sort keys%hash_raw);
##foreach my $raw (sort keys%hash_raw){
use Parallel::Loops;
my $plraw = Parallel::Loops->new($lnc);
##$plAC->share( \%GFI, \%GFIS, \%SMMPR, \%SMFPR,\%SMSPR);
$plraw->foreach (\@RAW, sub{
my $raw = $_ ;
chomp($raw);
 my @inputfile = glob("$input/$raw*.f*q*");
  if ($inputfile[0] =~ /R1|_1|F\./ && $inputfile[1] =~ /R2|_2|R\./){print "pair-end files found\n";}else{print "pair-end files not found or check file names (file name should contain R1 and R2 for pairs)\n";}
my $tmpbarcodef="$raw\_barcode_tmp.txt";
open (BARCODE , '>', $tmpbarcodef);
	foreach my $file_number (keys %{$hash{$raw}{'traceno'}}){
				print BARCODE $hash{$raw}{traceno}{$file_number},"\n";
				}
			close BARCODE;
	if (($inputfile[0] =~ /\.fastq$|\.fq$/i) && ($inputfile[1] =~ /\.fastq$|\.fq$/i)){$format = 'fastq';}
	if (($inputfile[0] =~ /\.gz$/i) && ($inputfile[1] =~ /\.gz$/i)){$format = 'gzfastq';}
	our $outtemp="$output\/$raw\_outtemp";
	our $radtemp="$output\/$raw\_radtemp";
        if ($step<2){
	`mkdir $outtemp $radtemp`;
	`process_radtags -P -1 $inputfile[0] -2 $inputfile[1] -o $outtemp -b $tmpbarcodef -e sbfI -E phred33 -r --disable_rad_check --barcode_dist_1 1 -i $format`;
	`mv $outtemp/*\.log $output/log/$raw\_barcode.log`;
	`mkdir $output/$raw`;
	`mkdir $output/$raw/remain_reads`;
	`mv $outtemp/*.rem.*.fq.gz $output/$raw/remain_reads`;
	print "first step  is working\n";
	}
my @namefiles = `ls $outtemp/*.1.fq.gz`;
s/$outtemp\///g for @namefiles;
chomp @namefiles;
use Parallel::Loops;
my $plsmpl = Parallel::Loops->new($snc);
$plsmpl->foreach (\@namefiles, sub{
my $file = $_ ;
chomp($file);
my $file2 = $file;
$file2 =~ s/.1.fq.gz/.2.fq.gz/;
my $f_name = $file;
$f_name =~ s/.1.fq.gz//;
my $bcode = $f_name;
$bcode =~ s/sample_//;
#print "This is the barcode $bcode\n\n"; 
my $s_name = $hash{$raw}{code}{$bcode};
`mkdir $radtemp/$s_name/`;
	    my @outfiles =`ls $outtemp/$f_name.*.fq.gz`;
            chomp @outfiles;
#            print "Second step is working too\n";
#	    s/$raw\_outtemp\///g for @outfiles;
#		print "This is the $outfiles[0] and $outfiles[1]\n\n\n";
                                                                                                `process_radtags -P -1 $outfiles[0] -2 $outfiles[1] -o $radtemp/$s_name -b RAD.txt -e sbfI -E phred33 -r --barcode_dist_1 2 --barcode_dist_2 2 --inline_inline --disable_rad_check -i gzfastq`;
                                                                                               `mv $radtemp/$s_name/*.log $output/log/$raw\_$s_name\_radtag.log`;
                                                                                                `mv $radtemp/$s_name/sample.1.fq.gz $outdDocent/pop$hash_pop{$s_name}\_$s_name.F.fq.gz`;
	                       	                                                                `mv $radtemp/$s_name/sample.2.fq.gz $outdDocent/pop$hash_pop{$s_name}\_$s_name.R.fq.gz`;
												#$mycmd = "zcat $outdDocent/pop$hash_pop{$s_name}\_$s_name\.R.fq.gz \| perl -pe  \'s/\\/2\\/2/\\/2/g\' \| gzip  -c > $outdDocent/pop$hash_pop{$s_name}\.Ra.fq.gz";
												#print "$mycmd\n";	
 	 											`zcat $outdDocent/pop$hash_pop{$s_name}\_$s_name\.R.fq.gz \| perl -pe  \'s/\\/2\\/2/\\/2/g\' \| gzip  -c > $outdDocent/pop$hash_pop{$s_name}\_$s_name.Ra.fq.gz`;												
												`zcat $outdDocent/pop$hash_pop{$s_name}\_$s_name\.F.fq.gz \| perl -pe  \'s/\\/1\\/1/\\/1/g\' \| gzip -c > $outdDocent/pop$hash_pop{$s_name}\_$s_name\.Fa.fq.gz`;
												#`rm $outdDocent/pop$hash_pop{$s_name}\_$s_name\.[FR].fq.gz`;
												`mv $outdDocent/pop$hash_pop{$s_name}\_$s_name\.Fa.fq.gz $outdDocent/pop$hash_pop{$s_name}\_$s_name\.Fb.fq.gz`;
												`mv $outdDocent/pop$hash_pop{$s_name}\_$s_name\.Ra.fq.gz $outdDocent/pop$hash_pop{$s_name}\_$s_name\.Rb.fq.gz`;
												

						});
`rm -r $tmpbarcodef`;


	});
}





1;

