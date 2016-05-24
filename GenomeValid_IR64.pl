#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"gff:s","bam:s","fa:s","genome:s","step:s","project:s","help");


my $help=<<USAGE;
perl $0 --step 1
--fa : genome sequence
--genome: genome sequence used to valid SV.
--gff: gff of SV to valid
  Chr1	SVpipe	Insertion	755001	755021	.	.	.	Size=-1;Method=pindel;
--bam: list of bam files for each library
  /rhome/cjinfeng/BigData/01.Rice_genomes/HEG4/00.Bam/HEG4_MSU7_BWA/FC52_7.MSU7_BWA.bam	500
--step: 1,2,3,4
--project: project dir, which will be created and store all files of this run.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{fa}  ||="/rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa";
$opt{genome} ||="/rhome/cjinfeng/BigData/01.Rice_genomes/HEG4/00.Genome/Chromosome/HEG4_ALLPATHLG_v1.chr.masked.fasta";
$opt{gff} ||="insertion.gff";
$opt{bam} ||="HEG4.bam.list";
$opt{project} ||= "Genome";
`mkdir $opt{project}` unless (-e "$opt{project}");

 
my $sv=readtable($opt{gff});
#my $bam=readlist($opt{bam});
my $seq=getfastaseq($opt{fa});
my $genome=getfastaseq($opt{genome});
my $len=getfastalen($opt{fa});
my $flank=2000;


my %type=(
   "0" => "No SV exists",
   "1" => "Local Assembly find SV and sequence: Insertion",
   "2" => "Local Assembly find SV and sequence: Deletion",
   "3" => "Need Manual verify if SV exists"
);

my $prefix=basename($opt{gff},".gff");
$prefix.=".0";

###generate target inf from genome assembly
`cut -f1,4,5 $opt{gff} > $prefix.table` unless (-e "$prefix.table");
`perl table2inf_laV2_IR64.pl --table $prefix.table --project $prefix` unless (-e "$prefix.inf");
#`perl table2inf_laV2_HEG4.pl --table $prefix.table --project $prefix` unless (-e "$prefix.inf");
`grep "HEG4" $prefix.inf | awk '{len=\$8-\$7;if(len < 20000){print}}' > $prefix.draw.inf`;
my $targetinf=readinf("$prefix.draw.inf");

###creat output gff file
writefile("","$prefix.Check.gff");
writefile("","$prefix.NoSV.gff");
writefile("","$prefix.NoMatch.gff");
writefile("","$prefix.LocalAssembly.gff");
writefile("","$prefix.Manual.gff");
###
my %summary;
foreach my $p (sort keys %$sv){
      `mkdir $opt{project}/$p` unless (-e "$opt{project}/$p");
      `mkdir $opt{project}/$p/genome` unless (-e "$opt{project}/$p/genome");
      my $start =$sv->{$p}->[3]-$flank >= 0 ? $sv->{$p}->[3]-$flank : 0;
      my $end   =$sv->{$p}->[4]+$flank <= $len->{$sv->{$p}->[0]} ? $sv->{$p}->[4]+$flank : $len->{$sv->{$p}->[0]};
      my $region="$sv->{$p}->[0]:$start-$end";
      #my @insert=values %$bam;
      print "$p\t$region\t$start\t$end\n";
      if ($opt{step}=~/1/){ 
      ####Region sequence
      my $subseq=substr($seq->{$sv->{$p}->[0]},$start,$flank*2); 
      $subseq   =formatseq($subseq,100);
      #print "$subseq\n";
      writefile(">$region\n$subseq\n","$opt{project}/$p/Region.fa");
      
      ####target sequence
      `mkdir $opt{project}/$p/genome` unless (-e "$opt{project}/$p/genome");
      my $targetseq=substr($genome->{$targetinf->{$p}->[0]},$targetinf->{$p}->[1],$targetinf->{$p}->[2]-$targetinf->{$p}->[1]+1); 
      $targetseq =formatseq($targetseq,100);
      writefile(">$targetinf->{$p}->[0]_$targetinf->{$p}->[1]_$targetinf->{$p}->[2]\n$targetseq\n","$opt{project}/$p/genome/scaffold.fa");
      }###step1
=pod
      ####Region reads
      foreach (@insert){
         writefile("","$opt{project}/$p/$_.raw.sam"); ### create new sam file
      }
      foreach (keys %$bam){
         `/usr/local/bin/samtools view $_ $region >> $opt{project}/$p/$bam->{$_}.raw.sam`;
      }
      foreach (@insert){
         `sort $opt{project}/$p/$_.raw.sam > $opt{project}/$p/$_.sam`;
      }
      }###step1
      
      if ($opt{step}=~/2/){
      ####Assembly
      `/rhome/cjinfeng/software/tools/Velvet/velvet/velveth $opt{project}/$p/assembly 31 -shortPaired -sam $opt{project}/$p/$insert[0].sam -shortPaired2 -sam $opt{project}/$p/$insert[1].sam`;
      `/rhome/cjinfeng/software/tools/Velvet/velvet/velvetg $opt{project}/$p/assembly -exp_cov 200 -ins_length $insert[0] -ins_length2 $insert[1] -min_contig_lgth 200 -scaffolding yes`;
      }###step2
=cut      
      if ($opt{step}=~/3/){
      ####Alignment
      my $query="$opt{project}/$p/genome/scaffold.fa";
      my $target="$opt{project}/$p/Region.fa";
      #`/opt/linux/centos/7.x/x86_64/pkgs/exonerate/2.2.0/bin/exonerate --querytype dna --targettype dna --gapextend -3 --query $query --bestn 50 --model affine:local --joinrangeext 300 --score 15 --target $target --gappedextension false --hspfilter 200 --dnahspdropoff 10 --showvulgar TRUE --showcigar TRUE --ryo "%S %pi %ql %C\n" > $opt{project}/$p/alignment.3`;
      `/rhome/cjinfeng/BigData/software/SVcaller/age_v0.4/src/age_align -indel -both $target $query > $opt{project}/$p/alignment.4`;
      }###step3

      if ($opt{step}=~/4/){
      ####Parse alignment and output refined SV
      print ">$p\n";
      ###indel=0,1,2: 0 mean non indel/sv, 1 mean insertion, 2 mean deletion, 3 mean not sure if exists a SV.
      my ($indel,$sequence)=parse_age("$opt{project}/$p/alignment.4","$opt{project}/$p/Region.fa","$opt{project}/$p/genome/scaffold.fa"); ###indel=0,1,2: 0 mean non indel/sv, 1 mean insertion, 2 mean deletion, 3 mean not sure if exists a SV.
      #my $prefix=basename($opt{gff},".gff");
      $summary{$indel}++;
      print "PARSE:\nindel: $indel\nSeq: $sequence\n";
      open CK, ">>$prefix.Check.gff" or die "$!";
         $temp=join("\t",@{$sv->{$p}});
         print CK "$temp\t$indel\t$type{$indel}\n";
      close CK;
      open NO, ">>$prefix.NoSV.gff" or die "$!";
      open NM, ">>$prefix.NoMatch.gff" or die "$!";
      open LOCAL, ">>$prefix.LocalAssembly.gff" or die "$!";
      open MANUAL, ">>$prefix.Manual.gff" or die "$!";
      if ($indel > 0){ ### if indel exists
         if ($indel == 1){
            $sv->{$p}->[8].="INDEL=Insertion;Seq=$sequence;";
            my $line=join("\t",@{$sv->{$p}});
            print LOCAL "$line\n";
         }elsif($indel == 2){
            $sv->{$p}->[8].="INDEL=Deletion;Seq=$sequence;";
            my $line=join("\t",@{$sv->{$p}});
            print LOCAL "$line\n";
         }elsif($indel == 3){ ### indel =3, need manual check
            my $line=join("\t",@{$sv->{$p}});
            print MANUAL "$line\n";
         }elsif($indel == 4){
            my $line=join("\t",@{$sv->{$p}});
            print NM "$line\n";
         }
      }else{
         my $line=join("\t",@{$sv->{$p}});
         print NO "$line\n";
      }
      close NO;
      close LOCAL;
      close MANUAL;
      }###step4 
}
foreach(sort {$a <=> $b} keys %summary){
      print STDERR "$_\t$summary{$_}\t$type{$_}\n";
}
sortfile("$prefix.NoSV.gff");
sortfile("$prefix.NoMatch.gff");
sortfile("$prefix.LocalAssembly.gff");
sortfile("$prefix.Manual.gff");
sortfile("$prefix.Check.gff");

##############################
sub sortfile
{
my ($file)=@_;
`sort -k1,1 -k4,4n $file > temp.sort`;
`mv temp.sort $file`;
}

sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    $hash{"$unit[0]_$unit[3]_$unit[4]"}=[@unit];
}
close IN;
return \%hash;
}

sub readinf
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    $hash{"$unit[1]_$unit[2]_$unit[3]"}=[$unit[5],$unit[6],$unit[7]];
}
close IN;
return \%hash;
}




sub readlist
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}



sub writefile
{
my ($lines,$file)=@_;
open WF, ">$file" or die "$!";
     print WF "$lines";
close WF;
}


sub formatseq
{
### format a single line sequence into lines with user specific length
my ($seq,$step)=@_;
my $length=length $seq;
my $run=int ($length/$step);
my $newseq;
for(my $i=0;$i<=$run;$i++){
   my $start=$i*$step;
   my $line=substr($seq,$start,$step);
   $newseq.="$line\n";
}
return $newseq;
}



 
sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}

sub getfastalen
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}= length $seq;
}
close IN;
$/="\n";
return \%hash;
}

###
#Alignment:
# first  seq =>  [  2,  9] EXCISED REGION [ 10,173]
# second seq =>  [174,168] EXCISED REGION [165,  2]
sub parse_age0
{
my ($align,$target,$query)=@_;

my $indel=3;
my $sequence="NA";
my $flank= 2000;
my $refseq=getfastaseq($query);
$/="\n\n\n";
open IN, "$align" or die "$!";
while(<IN>){
   my $block=$_;
   while($block=~/(MATCH.*Alignment time is .*? s)/sg){  ### match every alignment block for each contig
      #print "MATCH\n$1\nEND\n";
      my $match=$1;
      my $contig;
      if ($match=~/Second seq .* nucs \'(.*)\'\n/){  ### match contig name
         $contig=$1;
         print "$contig\n";
      }
      ### no excise
      if ($match=~/Alignment:\n first\s+seq =>\s+\[\s*(\d+)\,\s*(\d+)\]\n\s+second seq \=\> \s+\[\s*(\d+)\,\s*(\d+)\]\n/){
         print "$1\t$2\t$3\t$4\n";
         if ($1 < $flank-100 and $2 > $flank+100){ ### perfect alignment cover 200 bp of breakpoint
            $indel=0 unless ($indel == 1); ### indicate no SV exists
         }
      ### excise
      }elsif($match=~/Alignment:\n first\s+seq =>\s+\[\s*(\d+)\,\s*(\d+)\] EXCISED REGION \[\s*(\d+)\,\s*(\d+)\]\n\s+second seq \=\> \s+\[\s*(\d+)\,\s*(\d+)\] EXCISED REGION \[\s*(\d+)\,\s*(\d+)\]\n/){
         print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\n";
         if ($1 < $flank-100 and $4 > $flank+100){ ### excised alignment cover 200 bp of breakpoint
            if ($2 <= $flank+100 and $2 >= $flank-100 and $3 <= $flank+100 and $3 >= $flank-100){ ### excise region cover 100 bp of breakpoint
               print "SV Contig $contig: $6,$7\n";
               $indel=1; ### indicate SV exists
               my $seqlen=length $refseq->{$contig};
               my $start=$6 < $7 ? $6 : $7;
               my $len  =abs($7-$6+1);
               $sequence=substr($refseq->{$contig},$start,$len);
               my $gap;
               while($sequence=~/(N+)/g){
                  $gap+=length $1;
               }
               if ($gap > $len*0.3){
                  $indel=3 unless ($indel == 1 or $indel == 0); ### too much unknown sequence in SV sequence, set to not sure
               }
               if ($len < 10){
                  $indel=0 unless ($indel == 1); ### insertion sequence too small, probably not a SV
               }
            }elsif(($1 <= $flank-100 and $2 >= $flank+100) or ($3 <= $flank-100 and $4 >= $flank+100)){
               $indel=0 unless ($indel == 1); ### indicate no SV exists
            }else{
               $indel=3 unless ($indel == 1 or $indel == 0); ### indicate not sure
            }
         }
      }
   }
}
close IN;
$/="\n";
return ($indel,$sequence);
}


sub parse_age
{
my ($align,$target,$query)=@_;

my $indel=3;
my $sequence="NA";
my $flank = 2000;
my $jun   = 500;
my $refseq=getfastaseq($query);
if (exists $refseq->{'__'}){
   $indel = 4;
}
$/="\n\n\n";
open IN, "$align" or die "$!";
while(<IN>){
   my $block=$_;
   while($block=~/(MATCH.*Alignment time is .*? s)/sg){  ### match every alignment block for each contig
      #print "MATCH\n$1\nEND\n";
      my $match=$1;
      my $contig;
      if ($match=~/Second seq .* nucs \'(.*)\'\n/){  ### match contig name
         $contig=$1;
         print "$contig\n";
      }
      ### no excise
      if ($match=~/Alignment:\n first\s+seq =>\s+\[\s*(\d+)\,\s*(\d+)\]\n\s+second seq \=\> \s+\[\s*(\d+)\,\s*(\d+)\]\n/){
         print "$1\t$2\t$3\t$4\n";
         if ($1 < $flank-$jun and $2 > $flank+$jun){ ### perfect alignment cover 1000 bp of breakpoint
            $indel=0 unless ($indel == 1); ### indicate no SV exists
         }
      ### excise
      }elsif($match=~/Alignment:\n first\s+seq =>\s+\[\s*(\d+)\,\s*(\d+)\] EXCISED REGION \[\s*(\d+)\,\s*(\d+)\]\n\s+second seq \=\> \s+\[\s*(\d+)\,\s*(\d+)\] EXCISED REGION \[\s*(\d+)\,\s*(\d+)\]\n/){
         print "$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\n";
         if ($1 < $flank-$jun and $4 > $flank+$jun){ ### excised alignment cover 200 bp of breakpoint
            if ($2 <= $flank+$jun and $2 >= $flank-$jun and $3 <= $flank+$jun and $3 >= $flank-$jun){ ### excise region cover 100 bp of breakpoint
               print "SV Contig $contig: $6,$7\n";
               if (abs($7-$6) > 100){
                   $indel=1; ### indicate SV exists
               }
               my $seqlen=length $refseq->{$contig};
               my $start=$6 < $7 ? $6 : $7;
               my $len  =abs($7-$6+1);
               $sequence=substr($refseq->{$contig},$start,$len);
               my $gap;
               while($sequence=~/(N+)/g){
                  $gap+=length $1;
               }
               if ($gap > $len*0.7){
                  $indel=3 unless ($indel == 1 or $indel == 0); ### too much unknown sequence in SV sequence, set to not sure
               }elsif ($len < 20){
                  $indel=0 unless ($indel == 1); ### insertion sequence too small, probably not a SV
               }else{
                  $indel=1;
               }
            }elsif(($1 <= $flank-$jun and $2 >= $flank+$jun) or ($3 <= $flank-$jun and $4 >= $flank+$jun)){
               $indel=0 unless ($indel == 1); ### indicate no SV exists
            }else{
               $indel=3 unless ($indel == 1 or $indel == 0); ### indicate not sure
            }
         }
      }
   }
}
close IN;
$/="\n";
return ($indel,$sequence);
}

