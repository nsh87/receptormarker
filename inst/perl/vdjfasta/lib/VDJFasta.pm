package VDJFasta;
use strict;
use vars qw($rootdir);
no warnings 'recursion';
use DB_File;
use IO::Uncompress::AnyUncompress qw($AnyUncompressError);

# Author:  Jacob Glanville 
# Contact: jake@distributedbio.com
#
# Glanville J, Zhai W, Berka J et al. Precise determination of the diversity 
# of a combinatorial antibody library gives insight into the human immunoglobulin 
# repertoire. Proc Natl Acad Sci USA. 2009;106:20216–20221

BEGIN {
  use File::Basename;
  use Cwd qw(abs_path getcwd);
  $rootdir = dirname(dirname(abs_path($0)));
};

# Constructor ----------------------------------------------

# took these out of new - they were a major cpu hit - just do it once at package init time

my $check_blast;
my $check_hmmsearch;
my $check_hmmalign;

BEGIN {
  $check_blast     = `which blastn`;
  $check_hmmsearch = `which hmmsearch`;
  $check_hmmalign  = `which hmmalign`;
}

sub new {
  my ($class) = @_;
  my $self = {};
  $self->{filename}     = "";
  $self->{headerdb}     = "";
  $self->{sequencedb}   = "";
  $self->{headers}      = [];
  $self->{sequence}     = [];
  $self->{germline}     = [];
  $self->{nseqs}        =  0;
  $self->{mids}         = {};

  $self->{accVsegQstart}  = {}; # example: 124
  $self->{accVsegQend}    = {}; # example: 417
  $self->{accJsegQstart}  = {};
  $self->{accJsegQend}    = {};
  $self->{accDsegQstart}  = {};
  $self->{accDsegQend}    = {};
  $self->{accCsegQstart}  = {};
  $self->{accCsegQend}    = {};

  $self->{ranVseg}        =  0; # boolean
  $self->{ranDseg}        =  0;
  $self->{ranJseg}        =  0;
  $self->{ranCseg}        =  0;

  $self->{Vseg}           = {}; # example: IGHV1-46 293 3
  $self->{Dseg}           = {};
  $self->{Jseg}           = {};
  $self->{Cseg}           = {};

  $self->{aaVhseg}        = []; # example: IGHV1-46 98.9% 41  V1M 
  $self->{aaJhseg}        = []; # example: IGHJ6 91.6% 134 K104Q
  $self->{aaVhFWID}       = []; # example: 95% V1M K104Q

  $self->{aaVlseg}        = [];
  $self->{aaJlseg}        = [];
  $self->{aaVlFWID}       = []; # example: 88% Q130V I148V S155T S157K S158K

  $self->{V}              = []; # example: IGHV1-46 293 3
  $self->{D}              = [];
  $self->{J}              = [];
  $self->{C}              = [];
  $self->{H3}             = [];
  $self->{L3}             = [];
  $self->{coords}         = [];

  $self->{dnaVsegdb}    = "$rootdir/db/imgt.V.dna.fa"; 
  $self->{dnaJsegdb}    = "$rootdir/db/imgt.J.dna.fa";
  $self->{dnaDsegdb}    = "$rootdir/db/imgt.HD.dna.nr.fa";
  $self->{dnaCsegdb}    = "$rootdir/db/imgt.CH1.dna.nr.fa";

  $self->{VhVkHMM}      = "$rootdir/db/Vh-linker-Vk.hmm";

  $self->{blast}        = "blastn";
  $self->{tblastn}      = "tblastn";
  $self->{hmmsearch}    = "hmmsearch";
  $self->{hmmalign}     = "hmmalign";

  $self->{vbasefile}    = "$rootdir/db/vbase.imgt.map.txt";
  $self->{vhpssmfile}   = "$rootdir/db/human.equally.weighted.average.VH.pssm";
  $self->{vhpssm}       = {};  

  $self->{a2mref}       = "$rootdir/db/biologic.reference.c2m";

  $self->{blosum_file}  = "$rootdir/db/blosum62.txt";

  $self->{canon_file}   = "$rootdir/db/ighv.canonical.txt";
 
  $self->{verbose}      = 0;

  bless $self,'VDJFasta';

  # make sure blastn, hmmsearch and hmmalign are in the current path

  if($check_blast eq ""){
    print "VDJFasta cannot find blastn in \$PATH...\n";
    print "\tSolution: 	1. install NCBI blast 2.2.25+\n";
    print "\t		2. add the path to your blast bin dir in .bash_profile\n";
    print "\t		   (it will look something like /tools/apps/ncbi-blast-2.2.25+/bin/\n";
    print "\t		3. source ~/.bash_profile\n";
    exit;
  }
  if($check_hmmsearch eq ""){
    print "VDJFasta cannot find hmmsearch in \$PATH...\n";
    print "\tSolution:  1. install HMMER 3.0\n";
    print "\t           2. add the path to your hmmer bin dir in .bash_profile\n";
    print "\t		   (it will look something like /tools/apps/HMMER/binaries\n";
    print "\t           3. source ~/.bash_profile\n";
    exit;
  }
  if($check_hmmalign eq ""){
    print "VDJFasta cannot find hmmalign in \$PATH...\n";
    print "\tSolution:  1. install HMMER 3.0\n";
    print "\t           2. add the path to your hmmer bin dir in .bash_profile\n";
    print "\t              (it will look something like /tools/apps/HMMER/binaries\n";
    print "\t           3. source ~/.bash_profile\n";
    exit;
  }

  return $self;
}

# Methods --------------------------------------------------

sub convertVgene2Canonical {
  my($self)=@_;

  # first populate vgene2canon hash
  my %vgene2canon=();
  open(CANON,$self->{canon_file});
  my @lines=<CANON>;
  close(CANON);
  chomp(@lines);
  for(my $x=0;$x<scalar(@lines);$x++){
    my($vgene,$canonical_group)=split(/\t/,$lines[$x]);
    $vgene2canon{$vgene}=$canonical_group;
    #print "Setting $vgene to $canonical_group\n";
  } 

  # now go through each sequence and alter their V-gene assignment to a canon assignment
  # and alter their J-gene assignment to be generically IGHJ
  
  for(my $seqid=0;$seqid<$self->{nseqs};$seqid++){
    my $this_header = $self->getHeader($seqid);
    my @fields      = split(/;/,$this_header);

    # field 1 is the V-gene, field 3 is the J-gene
    my $vgene_field = $fields[1];
    my $jgene_field = $fields[3];
       $vgene_field =~s/ .*//;
    #print "got $vgene_field\n";
    if(defined($vgene2canon{$vgene_field})){
      $vgene_field=$vgene2canon{$vgene_field} . " 100 0";
    }else{
      #print "No match for $vgene_field\n";
      $vgene_field="IGHV 100 0";
    }
       $jgene_field = "IGHJ 30 0";
    $fields[1]=$vgene_field;
    $fields[3]=$jgene_field;

    my $new_header=join(";",@fields);

    $self->setHeader($seqid,$new_header); 
  }
}

sub parseHeaderAnnotations {
  my($self)=@_;
  print "stub\n";

  my $annotated=0;

  for(my $x=0;$x<scalar(@{$self->{headers}});$x++){
    my @fields = split(/;/,${$self->{headers}}[$x]);
    if(scalar(@fields)<6){
      print "Not annotated!\n";
    }else{
      ${$self->{V}}[$x]     = $fields[1];  
      ${$self->{D}}[$x]     = $fields[2];
      ${$self->{J}}[$x]     = $fields[3];
      ${$self->{H3}}[$x]    = $fields[4];
      ${$self->{L3}}[$x]    = $fields[5];
      ${$self->{C}}[$x]     = $fields[6];
      ${$self->{coords}}[$x]= $fields[7];
      ${$self->{H1}}[$x]    = $fields[8];
      ${$self->{H2}}[$x]    = $fields[9];
      ${$self->{L1}}[$x]    = $fields[10];
      ${$self->{L2}}[$x]    = $fields[11];
      $annotated=1;
    }
  }
  return $annotated;
}

sub printDNAPlasmidCoordinates {
  my($self)=@_;

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $accession=$self->getAccession($s);

    my($header,$seq)=$self->getSeq($s);
    # header fields
    my @header_fields=split(/;/,$header);

    # get the DNA sequence
    my $dnaseq = $header_fields[81];
    $seq = $self->translateSeq($dnaseq);

    # report on heavy chain features
    my $vh      = $header_fields[2];
    my $jh      = $header_fields[14];
    unless(defined($vh)){
      $vh="";
    }
    if(defined($jh)){
      $jh =~ s/ .*//;
    }else{
      $jh =  "";
    }

    my $vh_fw1  = $header_fields[17];
    my $vh_cdr1 = $header_fields[18];
    my $vh_fw2  = $header_fields[19];
    my $vh_cdr2 = $header_fields[20];
    my $vh_fw3  = $header_fields[21];
    my $vh_cdr3 = $header_fields[22];
    my $vh_fw4  = $header_fields[23];

    unless(defined($vh_fw1)){
      $vh_fw1="";
    }
    unless(defined($vh_cdr1)){
      $vh_cdr1="";
    }   
    unless(defined($vh_fw2)){
      $vh_fw2="";
    }    
    unless(defined($vh_cdr2)){
      $vh_cdr2="";
    } 
    unless(defined($vh_fw3)){
      $vh_fw3="";
    }    
    unless(defined($vh_cdr3)){
      $vh_cdr3="";
    } 
    unless(defined($vh_fw4)){
      $vh_fw4="";
    }    

    my($fw1_start,$fw1_stop)=$self->match_positions($seq,$vh_fw1);
    my($fw3_start,$fw3_stop)=$self->match_positions($seq,$vh_fw3);

    if($fw3_stop>$fw1_start){
      print $accession . "\t" . $vh . "\t" . ($fw1_start*3) . "\t" . ($fw3_stop*3) . "\n";
    }

    my($start,$stop)=$self->match_positions($seq,$vh_cdr1);
    if($stop>$start){
      print $accession . "\t" . "CDR-H1\t" . ($start*3) . "\t" . ($stop*3) . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vh_cdr2);
    if($stop>$start){
      print $accession . "\t" . "CDR-H2\t" . ($start*3) . "\t" . ($stop*3) . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vh_cdr3);
    if($stop>$start){
      print $accession . "\t" . "CDR-H3\t" . ($start*3) . "\t" . ($stop*3) . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vh_fw4);
    if($stop>$start){
      print $accession . "\t" . $jh . "\t" . ($start*3) . "\t" . ($stop*3) . "\n"; 
    }

    my $vl      = $header_fields[8];
    my $jl      = $header_fields[15];
    if(defined($jl)){
      $jl      =~ s/ .*//;
    }else{
      $jl="";
    }

    # now report on light chain features
    my $vl_fw1  = $header_fields[24];
    my $vl_cdr1 = $header_fields[25];
    my $vl_fw2  = $header_fields[26];
    my $vl_cdr2 = $header_fields[27];
    my $vl_fw3  = $header_fields[28];
    my $vl_cdr3 = $header_fields[29];
    my $vl_fw4  = $header_fields[30];


    ($fw1_start,$fw1_stop)=$self->match_positions($seq,$vl_fw1);
    ($fw3_start,$fw3_stop)=$self->match_positions($seq,$vl_fw3);

    if($fw3_stop>$fw1_start){
      print $accession . "\t" . $vl . "\t" . ($fw1_start*3) . "\t" . ($fw3_stop*3) . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vl_cdr1);
    if($stop>$start){
      print $accession . "\t" . "CDR-L1\t" . ($start*3) . "\t" . ($stop*3) . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vl_cdr2);
    if($stop>$start){
      print $accession . "\t" . "CDR-L2\t" . ($start*3) . "\t" . ($stop*3) . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vl_cdr3);
    if($stop>$start){
      print $accession . "\t" . "CDR-L3\t" . ($start*3) . "\t" . ($stop*3) . "\n";
    } 

    ($start,$stop)=$self->match_positions($seq,$vl_fw4);
    if($stop>$start){
      print $accession . "\t" . $jl . "\t" . ($start*3) . "\t" . ($stop*3) . "\n";
    }

    print $accession . "\tSequence\t" . $dnaseq . "\n";
  }
}


sub printAAPlasmidCoordinates {
  my($self)=@_;

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $accession=$self->getAccession($s);

    my($header,$seq)=$self->getSeq($s);
    #raw sequence
    $seq=~s/\.*//g;
    $seq=~s/\-*//g;
    $seq=uc($seq);
    # header fields
    my @header_fields=split(/;/,$header);

    # report on heavy chain features
    my $vh      = $header_fields[2];
    my $jh      = $header_fields[14];
    unless(defined($vh)){
      $vh="";
    }
    if(defined($jh)){
      $jh =~ s/ .*//;
    }else{
      $jh =  "";
    }

    my $vh_fw1  = $header_fields[17];
    my $vh_cdr1 = $header_fields[18];
    my $vh_fw2  = $header_fields[19];
    my $vh_cdr2 = $header_fields[20];
    my $vh_fw3  = $header_fields[21];
    my $vh_cdr3 = $header_fields[22];
    my $vh_fw4  = $header_fields[23];

    unless(defined($vh_fw1)){
      $vh_fw1="";
    }
    unless(defined($vh_cdr1)){
      $vh_cdr1="";
    }   
    unless(defined($vh_fw2)){
      $vh_fw2="";
    }    
    unless(defined($vh_cdr2)){
      $vh_cdr2="";
    } 
    unless(defined($vh_fw3)){
      $vh_fw3="";
    }    
    unless(defined($vh_cdr3)){
      $vh_cdr3="";
    } 
    unless(defined($vh_fw4)){
      $vh_fw4="";
    }    

    my($fw1_start,$fw1_stop)=$self->match_positions($seq,$vh_fw1);
    my($fw3_start,$fw3_stop)=$self->match_positions($seq,$vh_fw3);

    if($fw3_stop>$fw1_start){
      print $accession . "\t" . $vh . "\t" . $fw1_start . "\t" . $fw3_stop . "\n";
    }

    my($start,$stop)=$self->match_positions($seq,$vh_cdr1);
    if($stop>$start){
      print $accession . "\t" . "CDR-H1\t" . $start . "\t" . $stop . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vh_cdr2);
    if($stop>$start){
      print $accession . "\t" . "CDR-H2\t" . $start . "\t" . $stop . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vh_cdr3);
    if($stop>$start){
      print $accession . "\t" . "CDR-H3\t" . $start . "\t" . $stop . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vh_fw4);
    if($stop>$start){
      print $accession . "\t" . $jh . "\t" . $start . "\t" . $stop . "\n"; 
    }

    my $vl      = $header_fields[8];
    my $jl      = $header_fields[15];
    if(defined($jl)){
      $jl      =~ s/ .*//;
    }else{
      $jl="";
    }

    # now report on light chain features
    my $vl_fw1  = $header_fields[24];
    my $vl_cdr1 = $header_fields[25];
    my $vl_fw2  = $header_fields[26];
    my $vl_cdr2 = $header_fields[27];
    my $vl_fw3  = $header_fields[28];
    my $vl_cdr3 = $header_fields[29];
    my $vl_fw4  = $header_fields[30];


    ($fw1_start,$fw1_stop)=$self->match_positions($seq,$vl_fw1);
    ($fw3_start,$fw3_stop)=$self->match_positions($seq,$vl_fw3);

    if($fw3_stop>$fw1_start){
      print $accession . "\t" . $vl . "\t" . $fw1_start . "\t" . $fw3_stop . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vl_cdr1);
    if($stop>$start){
      print $accession . "\t" . "CDR-L1\t" . $start . "\t" . $stop . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vl_cdr2);
    if($stop>$start){
      print $accession . "\t" . "CDR-L2\t" . $start . "\t" . $stop . "\n";
    }

    ($start,$stop)=$self->match_positions($seq,$vl_cdr3);
    if($stop>$start){
      print $accession . "\t" . "CDR-L3\t" . $start . "\t" . $stop . "\n";
    } 

    ($start,$stop)=$self->match_positions($seq,$vl_fw4);
    if($stop>$start){
      print $accession . "\t" . $jl . "\t" . $start . "\t" . $stop . "\n";
    }

    print $accession . "\tSequence\t" . $seq . "\n";
  }
}

sub match_positions {
  my ($self,$string,$substring) = @_;
  if($string =~ /$substring/){
    return ($-[0], $+[0]);
  }else{
    return (0,0);
  }
}

sub sangerTrimBack {
  my($self,$window_size,$min_percent_coverage)=@_;
  
  for(my $s=0;$s<$self->{nseqs};$s++){
    my($header,$seq)=$self->getSeq($s);
   
    # trim forward
    my $start_trim_site=0;
    for(my $start=0;$start<length($seq);$start++){
      my $window=substr($seq,$start,$window_size);
      if($self->sangerPercentLegible($window)>$min_percent_coverage){
        $start_trim_site=$start;
        $start=length($seq);
      }
    }     

    # trim backward
    my $stop_trim_site=0;
    for(my $stop=(length($seq)-1);($stop-$window_size)>=0;$stop--){
      my $window=substr($seq,($stop-$window_size),$window_size);
      if($self->sangerPercentLegible($window)>$min_percent_coverage){
        $stop_trim_site=$stop;
        $stop=0;
      }
    } 
 
    my $trimmed_seq = substr($seq,$start_trim_site,($stop_trim_site-$start_trim_site));
    ${$self->{sequence}}[$s]=$trimmed_seq;
  }
}

sub sangerPercentLegible {
  my($self,$sequence)=@_;
  my $bases=0;
  my $legible=0;
  my @chars=split(/ */,$sequence);
  for(my $c=0;$c<scalar(@chars);$c++){
    $bases++;
    unless($chars[$c] =~ m/[NnXxZZ]/){
      $legible++;
    }
  }
  return ($legible/$bases);
}

sub lambdaNTermSmoothing {
  my($self,$print)=@_;
  # lambda chain FW1 regions are sometimes not aligned particularly well in the kappa/lambda scFv HMM. 
  # lamdba FW1 should look like
  #  QSVLTQPPSASGTPGQRVTISCSGSSSNIGSNTVNWYQQLPGTAPKLLIYSNNQRPSGVP
  #  QSALTQPRSVSGSPGQSVTISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGV
  #  SYELMQPPSVSVSPGQTARITCSGDALPKQYAYWYQQKPGQAPVLVIYKDSERPSGIPER

  # Error modes are               should be
  # 28 GsYELTQPP-SVSVSPGQTARITC   SYELTQPPSVSVSPGQTARITC
  # 38 GQsVLTQ-PPSVSGAPGQRVTISC   QSVLTQPPSVSGAPGQRVTISC
  # 40 GsYELTQPP-SVSVSPGQTASITC   SYELTQPPSVSVSPGQTASITC
  # 55 GQsVLTQ-PPSVSAAPGQKVTISC   QSVLTQPPSVSAAPGQKVTISC

  # or
  # 13 GQSALTQPASVSGSPGQSITISC    QSALTQPASVSGSPGQSITISC
  # 14 GQSVVTQPPSVSAAPGQKVTISC    QSVVTQPPSVSAAPGQKVTISC
  # 18 EIVLTQSPGTLSLSPGERATLSC    EIVLTQSPGTLSLSPGERATLSC
  # 59 GQSVVTQPPSVSGAPGQRVTISC    QSVVTQPPSVSGAPGQRVTISC

  # discover all correctable offenses
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $fw1 =  $self->getHMMCOLRange($s,128,140);
    my $xlinker =  $self->getHMMCOLRange($s,120,127);
    #print "$xlinker $fw1\n";

    if($fw1 =~ m/-/){
      my $linker =  $self->getHMMCOLRange($s,120,127);
      if( ($fw1 =~ m/[a-z]/) || ($linker =~ m/[a-z]/) ){
        my $only_match_states = $fw1;
           $only_match_states =~ s/[a-z-\.]//g; 
        if(length($only_match_states)>8){
          if($fw1=~m/[A-Z]$/){
            #print "$linker $fw1 $only_match_states\n";
            #  fix it: .GsYELTQPP-SVS
            # the fix: .G.SYELTQPPSVS
            $self->hmmBubbleRightShift($s,128,140);
          }
        }
      }
    }
  }
}

sub hmmBubbleRightShift {
  my($self,$s,$colstart,$colstop)=@_;
  # right-shift characters to fill in a missing gap state
  #  fix it: .GsYELTQPP-SVS
  # the fix: .G.SYELTQPPSVS

  my $nterm=$self->getHMMCOLRange($s,0,($colstart-1));
  my $sequence=$self->getHMMCOLRange($s,$colstart,$colstop);
  my $cterm=$self->getHMMCOLRange($s,($colstop+1),10000000000);
  my $chunk=$self->getHMMCOLRange($s,($colstart-20),($colstop+20));

  my $newseq="";
  my @chars=split(/ */,$sequence);
  my $activated_shift_cycle=0;
  for(my $c=(scalar(@chars)-1);$c>=0;$c--){
    if($c == 0){
      $newseq = $chars[$c] . $newseq;
    }elsif($chars[$c] =~ m/[A-Z]/){
      $newseq = $chars[$c] . $newseq;
    }elsif($chars[$c] =~ m/[a-z]/){
      $newseq = $chars[$c] . $newseq;
    }elsif($chars[$c] eq "."){
      $newseq = $chars[$c] . $newseq;
    }elsif($chars[$c] eq "-"){
      # grab the next available N-terminal character and move it over
      for(my $d=($c-1);$d>=0;$d--){
        if($chars[$d] =~ m/[A-Z]/){
          $activated_shift_cycle=1;
          $newseq = $chars[$d] . $newseq;           
          $chars[$d] = "-";
          $d=-1;
        }elsif($chars[$d] =~ m/[a-z]/){
          $activated_shift_cycle=1;
          $newseq = uc($chars[$d]) . $newseq;
          $chars[$d] = ".";
          $d=-1;
        }   
      }
    }
  }
  ${$self->{sequence}}[$s] = $nterm . $newseq . $cterm; 
  my $chunk2=$self->getHMMCOLRange($s,($colstart-20),($colstop+20));

  #print "\t$sequence\tinput\n";
  #print "\t$newseq\toutput\n";
  #print $nterm . "\t" . $newseq . "\t" .  $cterm . "\n";
  #print "\t" . $chunk . "\n";
  #print "\t" . $chunk2 . "\n";

}

sub cTermSmoothing {
  my($self,$print)=@_;
  # the last match state sometimes is left unoccupied, with the residue placed in an insert state prior
  # this code slides through looking for the pattern 
  # ("[a-z]-" or "[a-z][a-z]--") and corrects them (".[A-Z]" or "..[A-Z][A-Z]")
  # only relevant for a2m sequences

  for(my $s=0;$s<$self->{nseqs};$s++){
    my($header,$seq)=$self->getSeq($s);
    if( ($seq=~m/[A-Z]-[a-z][a-z\.]*$/) || ($seq=~m/[A-Z]--[a-z][a-z][a-z\.]*$/) || ($seq=~m/[A-Z]---[a-z][a-z][a-z][a-z\.]*$/) || ($seq=~m/[A-Z]----[a-z][a-z][a-z][a-z][a-z\.]*$/)){
      # .....nnnnnnSSSSSSSSSS--ccxxxx......

      my $cterm_insert=$seq;
         $cterm_insert=~s/.*[A-Z-]//; # stores ccxxxx.....

      my $nterm_section=$seq;
         $nterm_section=~s/[a-z\.]*$//; # ....nnnnnnSSSSSSSSSS--

      my $residues=0;
      if($seq=~m/[A-Z]-[a-z][a-z\.]*$/){
        $residues=1;
      }elsif($seq=~m/[A-Z]--[a-z][a-z][a-z\.]*$/){
        $residues=2;
      }elsif($seq=~m/[A-Z]---[a-z][a-z][a-z][a-z\.]*$/){
        $residues=3;
      }elsif($seq=~m/[A-Z]----[a-z][a-z][a-z][a-z][a-z\.]*$/){
        $residues=4;
      }
      # grab the last $residues characters from the nterm_insert and replace them with $residues "."s at the beginning
      # remove the first $residues dashes in the beginning of $cterm section and replace them with uc() residues from above
      my $stripped_cterm   = substr($cterm_insert,$residues) . $self->charRepeat(".",$residues);
      my $final_characters = uc(substr($cterm_insert, 0, $residues)); 
      my $stripped_nterm   = substr($nterm_section,0,0-$residues);

      # reset the value 
      $seq = $stripped_nterm . $final_characters . $stripped_cterm;
      ${$self->{sequence}}[$s] = $seq;
    }
    if($print){
      print $header . "\n";
      print $seq . "\n";
    }
  }
}

sub nTermSmoothing {
  my($self,$print)=@_;
  # the first match state often is left unoccupied, with the residue placed in an insert state prior
  # this code slides through looking for the pattern 
  # ("[a-z]-" or "[a-z][a-z]--") and corrects them (".[A-Z]" or "..[A-Z][A-Z]")
  # only relevant for a2m sequences

  for(my $s=0;$s<$self->{nseqs};$s++){
    my($header,$seq)=$self->getSeq($s);
    if( ($seq=~m/^[a-z\.]*[a-z]-[A-Z]/) || ($seq=~m/^[a-z\.]*[a-z][a-z]--[A-Z]/) || ($seq=~m/^[a-z\.]*[a-z][a-z][a-z]---[A-Z]/)){
      my $nterm_insert=$seq;
         $nterm_insert=~s/[A-Z-].*//;
         #$nterm_insert=~s/\.*//;
      my $cterm_section=$seq;
         $cterm_section=~s/[a-z\.*]*//;
      my $residues=0;
      if($seq=~m/^[a-z\.]*[a-z]-[A-Z]/){
        $residues=1;
      }elsif($seq=~m/^[a-z\.]*[a-z][a-z]--[A-Z]/){
        $residues=2;
      }elsif($seq=~m/^[a-z\.]*[a-z][a-z][a-z]---[A-Z]/){
        $residues=3;
      }
      # grab the last $residues characters from the nterm_insert and replace them with $residues "."s at the beginning
      # remove the first $residues dashes in the beginning of $cterm section and replace them with uc() residues from above
      my $stripped_nterm   = $self->charRepeat(".",$residues) . substr($nterm_insert,0,-$residues);
      my $final_characters = uc(substr($nterm_insert, -$residues)); 
      my $stripped_cterm   = substr($cterm_section,$residues);

      # reset the value 
      $seq = $stripped_nterm . $final_characters . $stripped_cterm;
      ${$self->{sequence}}[$s] = $seq;
      #print $stripped_nterm . $final_characters . $stripped_cterm . "\n";
    }
    if($print){
      print $header . "\n";
      print $seq . "\n";
    }
  }
}

sub charRepeat {
  my($self,$char,$repeat)=@_;
  my $string="";
  for(my $x=0;$x<$repeat;$x++){
    $string.=$char;
  }
  return $string;
}


sub getH3CoordsFromA2M {
  my($self,$seq)=@_;
  # length to ends, length of H3 
  # want do do the H3 curate correction here
  my $h3Begins = $self->getHMMCOLRange($seq,0,92);
  my $h3Seq    = $self->getHMMCOLRange($seq,92,104);
  my $h3Ends   = $self->getHMMCOLRange($seq,104,400); 

  $h3Begins =~ s/[\.\-]//g;
  $h3Seq =~ s/[\.\-]//g;
  $h3Ends =~ s/[\.\-]//g;
  $h3Begins =~ s/.$//;
  $h3Ends =~ s/^.//;

  my $orientation=$self->getHeader($seq);
     $orientation=~s/;.*//;
     $orientation=~s/..*frame_//;
  my $accession=$self->getAccession($seq);

  return($accession,$orientation,(length($h3Begins)+2),(length($h3Seq)-4),(length($h3Ends)+2));
}

sub getVariableCoordsFromA2M {
  my($self,$seq)=@_;
  # length to ends, length of H3 
  # want do do the H3 curate correction here
  my $h3Begins = $self->getHMMCOLRange($seq,0,92);
  my $h3Seq    = $self->getHMMCOLRange($seq,92,104);
  my $h3Ends   = $self->getHMMCOLRange($seq,104,400);

  my $jhSeq    = $self->getHMMCOLRange($seq,103,112);
  my $gsSeq    = $self->getHMMCOLRange($seq,113,128);
  my $vlSeq    = $self->getHMMCOLRange($seq,128,236);

  $h3Begins =~ s/[\.\-]//g;
  $h3Seq =~ s/[\.\-]//g;
  $h3Ends =~ s/[\.\-]//g;
  $h3Begins =~ s/.$//;
  $h3Ends =~ s/^.//;

  $jhSeq =~ s/[\.\-]//g;
  $gsSeq =~ s/[\.\-]//g;
  $gsSeq =~ s/.$//;
  $vlSeq =~ s/[a-z\.]*$//g;
  $vlSeq =~ s/[\.\-]//g;
  $vlSeq =~ s/^[a-z]*//;

  my $vhSeq = $h3Begins;
     $vhSeq =~ s/^[a-z\.]*//;

  my $orientation=$self->getHeader($seq);
     $orientation=~s/;.*//;
     $orientation=~s/..*frame_//;
  my $accession=$self->getAccession($seq);

  return($accession,$orientation,(length($h3Begins)+2),(length($h3Seq)-4),(length($h3Ends)+2),
         length($jhSeq),length($gsSeq),length($vlSeq),length($vhSeq));
}

sub getVariableSegmentsFromDNA {
  my ($self,$dna_seqid,$orientation,$BeginLen,$h3len,$EndLen,$jSeq,$gsSeq,$vlLength,$vhLength)=@_;
  #print "Values are $jSeq,$gsSeq,$vlLength\n";
  my $start=0;
  my $stop=0;
  my $aasequence_length = $BeginLen + $h3len + $EndLen;
  my $cdrh3_dna_length= 3 * $h3len; #($EndLen - $BeginLen);
  my $lvseg  = ""; # leader sequence and vseg: anything before cdr3
  my $h3_dna = "";
  my $jcseg  = ""; # anything after cdr3
  my $gsseg  = ""; # linker
  my $vlseg  = ""; # variable domain
  my $vhseg  = "";
  if($orientation>=0){                          #forward orientation
    $start = $orientation + ( 3 * $BeginLen);
    $stop  = $start + $cdrh3_dna_length;
    $h3_dna=$self->getSubSequence($dna_seqid,$start,$stop);
    $lvseg =$self->getSubSequence($dna_seqid,0,$start);
    #$vhseg =$self->getSubSequence($dna_seqid,($start-($vhLength*3)),($vhLength*3));
    $jcseg =$self->getSubSequence($dna_seqid,$stop,($stop+($jSeq*3)));
    $gsseg =$self->getSubSequence($dna_seqid,($stop+($jSeq*3)),($stop+(($jSeq+$gsSeq)*3)));
    $vlseg =$self->getSubSequence($dna_seqid,($stop+(($jSeq+$gsSeq)*3)),($stop+(($jSeq+$gsSeq)*3)+($vlLength*3)));

    #print "forward - " . length($jcseg) . "\t" . length($dna_seqid) . "\n";
  }else{                                        # reverse orientation
    my $reverse_seq=$self->getReverseStrand($dna_seqid);
    my $negative_offset = abs($orientation) - 1;
    $start = ( 3 * $BeginLen) + $negative_offset;
    $stop  = $start + $cdrh3_dna_length;

    # now we go from start to stop
    my @chars=split(/ */,$reverse_seq);
    for(my $c=$start;$c<$stop;$c++){
      $h3_dna.=$chars[$c];
    }
    #select the variable domain
    $lvseg =substr($reverse_seq,0,$start);
    #$vhseg =substr($reverse_seq,($start-($vhLength*3)),($vhLength*3)); 
    #$self->getSubSequence($dna_seqid,($start-($vhLength*3)),($vhLength*3));
    # select the jsegment. should always be a fixed length beyond the stop!
    $jcseg =substr($reverse_seq,$stop,($jSeq*3));  #  
    #print "reverse - " . length($jcseg) . "\t" . length($reverse_seq) . "\n";
    $gsseg =substr($reverse_seq,($stop+($jSeq*3)),($gsSeq*3)); #$self->getSubSequence($dna_seqid,($stop+30),($stop+75));
    $vlseg =substr($reverse_seq,($stop+(($jSeq+$gsSeq)*3)),($vlLength*3)); #$self->getSubSequence($dna_seqid,($stop+75),($stop+75+($vlLength*3)));
  }
  # make sure this stretch equals what is observed in -d";" -f5
  my $h3=VDJFasta->new();
     $h3->addSeq("header",$h3_dna);
  my ($h3_header,$h3_translation)=$h3->dna2aa(0,0);
  my $qc_h3 = $self->getHeaderField($dna_seqid,4);

  # now get heavy chain
  $vhseg = substr($lvseg,(length($lvseg)-($vhLength*3)-6));

  if ($qc_h3 ne $h3_translation){
    return ("","","");
  }else{
    return ($vhseg,$h3_dna,$jcseg,$gsseg,$vlseg);
  }
}


sub getH3fromDNA { 
  my ($self,$dna_seqid,$orientation,$BeginLen,$h3len,$EndLen)=@_;

  unless(defined($dna_seqid)){
    print "Error: getH3fromDNA passed a null id\n";
  #  exit;
  }

  my $vlLength = 90;
  my $start=0;
  my $stop=0;
  my $aasequence_length = $BeginLen + $h3len + $EndLen;
  my $cdrh3_dna_length= 3 * $h3len; #($EndLen - $BeginLen);
  my $lvseg  = ""; # leader sequence and vseg: anything before cdr3
  my $h3_dna = "";
  my $jcseg  = ""; # anything after cdr3
  my $gsseg  = ""; # linker
  my $vlseg  = ""; # variable domain
  if($orientation>=0){   			#forward orientation
    $start = $orientation + ( 3 * $BeginLen);
    $stop  = $start + $cdrh3_dna_length; 
    $h3_dna=$self->getSubSequence($dna_seqid,$start,$stop);
    $lvseg =$self->getSubSequence($dna_seqid,0,$start); #added +9
    $jcseg =$self->getSubSequence($dna_seqid,$stop,($stop+30)); # stop-3 stop+30
    $gsseg =$self->getSubSequence($dna_seqid,($stop+30),($stop+75));   
    $vlseg =$self->getSubSequence($dna_seqid,($stop+75),($stop+75+($vlLength*3))); 
 
  }else{     					# reverse orientation
    my $reverse_seq=$self->getReverseStrand($dna_seqid);
    my $negative_offset = abs($orientation) - 1;
    $start = ( 3 * $BeginLen) + $negative_offset;
    $stop  = $start + $cdrh3_dna_length;

    # now we go from start to stop
    my @chars=split(/ */,$reverse_seq);
    for(my $c=$start;$c<$stop;$c++){
      if(defined($chars[$c])){
      $h3_dna.=$chars[$c];
      }
    }
    my $reverse_seq_length=length($reverse_seq);

    #select the variable domain
    $lvseg =substr($reverse_seq,0,$start);
    # select the jsegment. should always be a fixed length beyond the stop!
    # Aiming for 0 272.83726 30 of (no sequence)
    #print "Aiming for $reverse_seq_length $stop 30 of $reverse_seq\n";
    if(($stop + 30) < $reverse_seq_length){
      $jcseg =substr($reverse_seq,$stop,30);  # 
    } 
    if(($stop + 75) < $reverse_seq_length){
      $gsseg =substr($reverse_seq,($stop+30),45); #$self->getSubSequence($dna_seqid,($stop+30),($stop+75));
    }
    if(($stop + 75 + ($vlLength*3)) < $reverse_seq_length){
      $vlseg =substr($reverse_seq,($stop+75),($vlLength*3)); #$self->getSubSequence($dna_seqid,($stop+75),($stop+75+($vlLength*3)));
    }
  }
  # make sure this stretch equals what is observed in -d";" -f5
  my $h3=VDJFasta->new();
     $h3->addSeq("header",$h3_dna);
  my ($h3_header,$h3_translation)=$h3->dna2aa(0,0);
  my $qc_h3 = $self->getHeaderField($dna_seqid,4);

  #$jcseg="~~~" . $jcseg . "|||";
  #print "MONKEY Got $qc_h3 vs $h3_translation\n";
  if ($qc_h3 ne $h3_translation){
    return ("","","");
  }else{
    return ($lvseg,$h3_dna,$jcseg,$gsseg,$vlseg);
  }
}

sub getSubSequence {
  my($self,$seq,$start,$stop)=@_;

  if($start eq ""){
    return "";
  }
  if($stop eq ""){
    return "";
  }
  unless(defined($seq)){
    return "";
  }
  if($stop<$start){
    my $temp=$start;
    $start=$stop;
    $stop=$temp;
  }
  if($start<0){
    $start=0;
  }
  if($stop<0){
    $stop=0;
  }
  if($stop > length(${$self->{sequence}}[$seq])){
    $stop=length(${$self->{sequence}}[$seq]);
  }
  if( ($stop - $start) < 0){
    return "";
  }
  my $substr=substr(${$self->{sequence}}[$seq],$start,($stop - $start));
  return $substr;
}


sub reportGermFreqs {
  my($self,$domain,$minlength,$min_shm,$max_shm,$isRedundant)=@_;
  print "stub\n";
}

sub getGermlines {
  my($self,$domain,$sequence_type,$custom_reference_db,$get_alleles)=@_;
  my $blast_type="blastn";
  if($sequence_type eq "aa"){
    $blast_type="tblastn";
  }
  my %segment_classifier_hash=();

  my $reference_db="";
  my $evalue="";
  if($domain eq "V"){
    $reference_db=$self->{dnaVsegdb};
    $evalue=10e-10;
  }elsif($domain eq "J"){
    $reference_db=$self->{dnaJsegdb};
    $evalue=0.001;
  }elsif($domain eq "D"){
    $reference_db=$self->{dnaDsegdb};
    $evalue=1;
  }elsif($domain eq "C"){
    $reference_db=$self->{dnaCsegdb};
    $evalue=0.001;
  }else{
    print " (V J D C)\n";
    exit;
  }
  if(defined($custom_reference_db) and ($custom_reference_db ne "")){
    if(-f $custom_reference_db){
      $reference_db=$custom_reference_db;
    }else{
      print "Attempted to load domain $domain custom reference db $custom_reference_db that does not exist.\n";
      exit(1);
    }
  }

  if($domain eq "D"){
    %segment_classifier_hash = $self->batchBlastParseDsegClassifier($reference_db,$evalue,$blast_type,$domain);
  }else{
    %segment_classifier_hash = $self->batchBlastParseClassifier($reference_db,$evalue,$blast_type,$domain,$get_alleles);
  }
  # handle ambiguous cases
  my $ambiguous_reduce=1;
  if($ambiguous_reduce){
    my @keys = keys %segment_classifier_hash;
    for(my $i=0;$i<scalar(@keys);$i++){
      my @fields = split(/ /,$segment_classifier_hash{$keys[$i]});
      if(scalar(@fields)>4){
        # stub for the advanced case resolver - to be added here later
        if($fields[0]=~m/IG[HKL]V/){
          $fields[0]=~s/-.*//;
        }elsif($fields[0]=~m/IG[HJL]J/){
          $fields[0]=~s/[0-9].*//;
        }else{
          $fields[0]=~s/[0-9]$//;
        }
        $segment_classifier_hash{$keys[$i]} = $fields[0] . " " . $fields[1] . " " . $fields[2];
      }
    }
  }
  return %segment_classifier_hash;
}

sub cleanHeaders {
  my($self)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    ${$self->{headers}}[$s]=~s/;/./g;
  }
}

sub reportH3charge {
  my($self,$max_shm)=@_;
  my %charge_counts=();
  my $total_counts=0;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my @vgene=split(/ */,$self->getHeaderField($s,1));
    if($vgene[2] <= $max_shm){
      my @h3=split(/ */,$self->getHeaderField($s,4));
      my $this_charge=0;
      for(my $p=0;$p<scalar(@h3);$p++){
        if($h3[$p]=~m/[RK]/){
          $this_charge++;
        }elsif($h3[$p]=~m/[DE]/){
          $this_charge--;
        }
      }
      $total_counts++;
      if(defined($charge_counts{$this_charge})){
        $charge_counts{$this_charge}++;
      }else{
        $charge_counts{$this_charge}=1;
      }
    }
  }
  for(my $charge=-10;$charge<10;$charge++){
    my $percentage=0;
    if(defined($charge_counts{$charge})){
      $percentage=$charge_counts{$charge}/$total_counts;
    }
    print $percentage . "\n";
  }
}

sub reportSHMDist {
  my($self,$min_vgene_span)=@_;

  my %shm_counts=();
  my $total_counts=0;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my @vgene_fields=split(/ /,$self->getHeaderField($s,1));
    if(scalar(@vgene_fields)==3){
      if($vgene_fields[1] >= $min_vgene_span){
        if(defined($shm_counts{$vgene_fields[2]})){
	  $shm_counts{$vgene_fields[2]}++;
	}else{
	  $shm_counts{$vgene_fields[2]}=1;
	}
	$total_counts++;
      }
    }
  }
  for(my $shm=0;$shm<50;$shm++){
    my $percentage=0;
    if(defined($shm_counts{$shm})){
      $percentage=$shm_counts{$shm}/$total_counts;
    }
    print $percentage . "\n";
  }
}

sub batchBlastTopHit {
  my($self,$db,$e,$outfile,$blast_type,$domain)=@_;
  # grabbing top hit without probabilistic classification
  my @lines=$self->getBlastLines($self->{filename},$db,$e,$blast_type);
  
  my %accessions=();
  my %acc_tophit=();
  for(my $x=0;$x<scalar(@lines);$x++){
    my($acc,$match,$pid,$alnlength,$mismatch,$gaps,$qstart,$qend,$mstart,$mend,$eval,$bit)=split(/\t/,$lines[$x]);
    my $hitname=$match; #$self->parse_hitname($match);
    if(!(defined($acc_tophit{$acc}))){
      $acc_tophit{$acc} = $hitname;
      $accessions{$acc} = $hitname . " " . $alnlength . " " . $mismatch . " ";
    }#else{
      # if this hit is equal to the top hit, then ambiguous
     # my($current_best_hitname,$current_best_alnlength,$current_best_mismatch)=split(/ /,$acc_tophit{$acc});
     # if($current_best_alnlength <= $alnlength){
     #   if($current_best_mismatch > $mismatch){
     #     print "Updating to a better match\n";
     #     $accessions{$acc} = $hitname . " " . $alnlength . " " . $mismatch . " ";
     #   }elsif($current_best_mismatch == $mismatch){
     #     print "Ambiguous! Will not call (length $current_best_alnlength <= $alnlength; mismatch $current_best_mismatch >= $mismatch\n";
     #     $accessions{$acc} = $current_best_hitname . "/" . $hitname . " " . $alnlength . " " . $mismatch . " ";
     #   }
     # }
    #}
  }
  my @keys=keys %accessions;
  open(OUT,">$outfile");
  for(my $k=0;$k<scalar(@keys);$k++){
    print OUT $keys[$k] . "\t" . $accessions{$keys[$k]} . "\n";
  }
  print OUT "WRITEDONEFLAG" . "\t" . "\n";
  close(OUT);
}

sub getAlleleCallHash {
  my($self,$db)=@_;
  my %allele_hash=();
  my $counts=$self->getSeqCount();
  
  for(my $seqid=0;$seqid<$counts;$seqid++){
     my $acc = $self->getAccession($seqid);
     if(defined(${$self->{"tophit" . $db}}{$acc})){
       #print ${$self->{"tophit" . $db}}{$acc} . "\n";
       $allele_hash{$acc}=${$self->{"tophit" . $db}}{$acc};
     }
  }
  return %allele_hash;
}

sub batchBlastParseClassifier {
  my($self,$db,$e,$blast_type,$domain,$get_alleles)=@_;
  my %segment_classifier_hash=();
  my @lines=$self->getBlastLines($self->{filename},$db,$e,$blast_type);

  my %accessions=();
  my %acc_tophit=();
  my $best_allele="";
  my $longest_alignment=0;
  my $lowest_mismatches=10000000;
  for(my $x=0;$x<scalar(@lines);$x++){
    my($acc,$match,$pid,$alnlength,$mismatch,$gaps,$qstart,$qend,$mstart,$mend,$eval,$bit)=split(/\t/,$lines[$x]);
    # stub for allowing allele level retention of hitname
    my $hitname=$self->parse_hitname($match);

    if(!(defined($acc_tophit{$acc}))){    
      # store_tophit   
      ${$self->{"tophit" . $db}}{$acc}=$match;  # altered compared with old form of code
      $acc_tophit{$acc} = $hitname;
      # set the best match allele
      $best_allele=$match;
      $longest_alignment=$alnlength;
      $lowest_mismatches=$mismatch;
      my $prob = 1;
      if($get_alleles){
        $accessions{$acc} = $best_allele . " " . $alnlength . " " . $mismatch . " ";
      }else{
        $accessions{$acc} = $hitname . " " . $alnlength . " " . $mismatch . " ";
      }
      if($domain eq "V"){
        ${$self->{accVsegQstart}}{$acc}  = $qstart;
        ${$self->{accVsegQend}}{$acc}    = $qend;
      }elsif($domain eq "J"){
        ${$self->{accJsegQstart}}{$acc}  = $qstart;
        ${$self->{accJsegQend}}{$acc}    = $qend;
      }elsif($domain eq "C"){
        ${$self->{accCsegQstart}}{$acc}  = $qstart;
        ${$self->{accCsegQend}}{$acc}    = $qend;
      }
    }elsif($acc_tophit{$acc} ne $hitname){
      my $prob = 1;
 
      if(defined($accessions{$acc})){
        my($best_match,$best_alnlength,$best_mismatch)=split(/ /,$accessions{$acc});
        $prob = $self->miscall_probability($best_alnlength,$alnlength,$best_mismatch,$mismatch);
        $prob = (int(1000 * $prob))/1000;
        if($mismatch>60){ 
          $prob=0;
        }
        # Add the allele classification
        if($get_alleles){
          if($alnlength>=$best_alnlength){
            if($mismatch<=$best_mismatch){
              # allele call is no longer valid
              #print $accessions{$acc} . "\n";
              $accessions{$acc}=~s/\*[^ ]*//;
              #print $accessions{$acc} . "\n";
            }
          }
        }
      }
      if(($blast_type eq "blastn") && ($prob > 0.01)){
        $accessions{$acc} .= $hitname . " " . $alnlength . " " . $mismatch . " ";
      }elsif(($blast_type eq "tblastn") && ($prob > 0.5)){
        $accessions{$acc} .= $hitname . " " . $alnlength . " " . $mismatch . " ";
      }
    }
  }

  my @keys=keys %accessions;
  #open(OUT,">$outfile");
  for(my $k=0;$k<scalar(@keys);$k++){
    #print OUT $keys[$k] . "\t" . $accessions{$keys[$k]} . "\n";
    $segment_classifier_hash{$keys[$k]}=$accessions{$keys[$k]};
    $segment_classifier_hash{$keys[$k]}=~s/ $//;
    ${$self->{$domain . "seg"}}{$keys[$k]}=$accessions{$keys[$k]};
  }
  #print OUT "WRITEDONEFLAG" . "\t" . "\n";
  #close(OUT);
  $self->{"ran" . $domain . "seg"}=1;
  return %segment_classifier_hash;
}

sub getDsegCoords {
  my($self)=@_;
  my %dseg_coords=();

  my @accs= keys %{$self->{accDsegQstart}};

  #open(OUTFILE,">$outfile");
  for(my $a=0;$a<scalar(@accs);$a++){
    $dseg_coords{$accs[$a]}=${$self->{accVsegQstart}}{$accs[$a]} . " " 
          . ${$self->{accVsegQend}}{$accs[$a]} . " "
          . ${$self->{accDsegQstart}}{$accs[$a]} . " "
          . ${$self->{accDsegQend}}{$accs[$a]} . " "
          . ${$self->{accJsegQstart}}{$accs[$a]} . " "
          . ${$self->{accJsegQend}}{$accs[$a]};

    #print OUTFILE  $accs[$a] . "\t"
    #      . ${$self->{accVsegQstart}}{$accs[$a]} . " "
    #      . ${$self->{accVsegQend}}{$accs[$a]} . " "
    #      . ${$self->{accDsegQstart}}{$accs[$a]} . " "
    #      . ${$self->{accDsegQend}}{$accs[$a]} . " "
    #      . ${$self->{accJsegQstart}}{$accs[$a]} . " "
    #      . ${$self->{accJsegQend}}{$accs[$a]} . "\n";
  }
  #close(OUTFILE);
  return %dseg_coords;
}

sub isDsegSandwiched {
  my($self,$acc,$dseg_qstart,$dseg_qend)=@_;

  my $answer=1;

  # if not defined, then exit
  unless(defined(${$self->{accVsegQstart}}{$acc}) && defined(${$self->{accVsegQend}}{$acc})
         && ${$self->{accJsegQstart}}{$acc} && ${$self->{accJsegQend}}{$acc} ){
    return 0;
  }

  # 1. get middle of each segment
  my $v_average= (${$self->{accVsegQstart}}{$acc} + ${$self->{accVsegQend}}{$acc}) /2;
  my $j_average= (${$self->{accJsegQstart}}{$acc} + ${$self->{accJsegQend}}{$acc}) /2;
  my $d_average= ($dseg_qstart + $dseg_qend) / 2;
  
  # 2. determine orientation
  my $orientation="forward";
  my $vseg_3prime=${$self->{accVsegQend}}{$acc};
  my $jseg_5prime=${$self->{accJsegQstart}}{$acc};
  if( $v_average > $j_average ){
    $orientation="reverse";
    $vseg_3prime=${$self->{accVsegQstart}}{$acc};
    $jseg_5prime=${$self->{accJsegQend}}{$acc};
  }

  # get lower
  my $lower_d=$dseg_qstart;
  my $higher_d=$dseg_qend;
  if($dseg_qend<$lower_d){
    $lower_d=$dseg_qend;
    $higher_d=$dseg_qend;
  }

  # make sure D-seg is sandwiched
  if( ($v_average < $d_average) && ($d_average < $j_average) ) {
    # state is forward
    my $v_overlap= ${$self->{accVsegQend}}{$acc} - $lower_d;
    my $j_overlap= $higher_d - ${$self->{accJsegQstart}}{$acc};
    if($v_overlap>6){
      $answer=0;
    }
    if($j_overlap>6){
      $answer=0;
    }
  }elsif( ($v_average > $d_average) && ($d_average > $j_average)){
    # state is reversed
    my $v_overlap= $higher_d - ${$self->{accVsegQstart}}{$acc};
    my $j_overlap= ${$self->{accJsegQend}}{$acc} - $lower_d;
    if($v_overlap>6){
      $answer=0;
    }
    if($j_overlap>6){
      $answer=0;
    }
  }else{
    $answer=0;
  }

  return $answer;
}

sub batchBlastParseDsegClassifier {
  my($self,$db,$e,$blast_type,$domain)=@_;
  my %segment_classifier_hash=();
  unless($self->{ranVseg} && $self->{ranJseg}){
    print "V-seg and J-seg need to run first. Exiting...\n";
    exit;
  }
  my @lines=$self->getBlastLines($self->{filename},$db,$e,$blast_type);
 
  my %accessions=();
  my %acc_tophit=();

  for(my $x=0;$x<scalar(@lines);$x++){
    my($acc,$match,$pid,$alnlength,$mismatch,$gaps,$qstart,$qend,$mstart,$mend,$eval,$bit)=split(/\t/,$lines[$x]);
    my $hitname=$self->parse_hitname($match);

    if($self->isDsegSandwiched($acc,$qstart,$qend)){
      if(!(defined($acc_tophit{$acc}))){
        $acc_tophit{$acc} = $hitname;
        my $prob = 1;
        $accessions{$acc} = $hitname . " " . $alnlength . " " . $mismatch . " ";
        ${$self->{accDsegQstart}}{$acc}  = $qstart;
        ${$self->{accDsegQend}}{$acc}    = $qend;
      }elsif($acc_tophit{$acc} ne $hitname){
        my $prob = 1;
        if(defined($accessions{$acc})){
          my($best_match,$best_alnlength,$best_mismatch)=split(/ /,$accessions{$acc});
          $prob = $self->miscall_probability($best_alnlength,$alnlength,$best_mismatch,$mismatch);
          $prob = (int(1000 * $prob))/1000;
        }
        if($prob > 0.01){
          $accessions{$acc} .= $hitname . " " . $alnlength . " " . $mismatch . " ";
        }
      }
    }
  }

  my @keys=keys %accessions;
  for(my $k=0;$k<scalar(@keys);$k++){
    $segment_classifier_hash{$keys[$k]}=$accessions{$keys[$k]};
    ${$self->{$domain . "seg"}}{$keys[$k]}=$accessions{$keys[$k]};
    $segment_classifier_hash{$keys[$k]}=~s/ $//;
  }
  $self->{"ran" . $domain . "seg"}=1;
  return %segment_classifier_hash;
}

sub parse_hitname {
  my($self,$hitname)=@_;
  $hitname=~s/\*.*//;
  if($hitname=~m/IGH[DGMEA]/){
  #  $hitname=~s/[0-9]*$//;
    $hitname=~s/P$//;
  }
  return $hitname;
}

sub miscall_probability{
  my($self,$sequence_length,$sequence_length2,$mutations1,$mutations2)=@_;
  my $length_diff=0;
  if($sequence_length2<$sequence_length){
    $length_diff=$sequence_length - $sequence_length2;
  }
  my $specific_mutations=$mutations2 - $mutations1 + $length_diff;
  my $total_mutations=$mutations2 + $specific_mutations;

  # PNAS 2009 original code
  #my $result=$self->nCk( ($sequence_length - $specific_mutations),($total_mutations - $specific_mutations) );
  #my $result2=$self->nCk($sequence_length,$total_mutations);
  #my $probability=($result/$result2);
  #return $probability;

  # Bob's replacement code
  # doesn't add a lot of value - just separation for Bob's benefit
  my $n=($sequence_length - $specific_mutations);
  my $k=($total_mutations - $specific_mutations);

  # return probability 0 if n<=k or sequence length < total_mutations, as this is uninterpretable
  my $probability=0;

  if( ($n<=$k) or ($sequence_length <= $total_mutations)){
    return $probability;
  }
  my $result=$self->nCk( ($sequence_length - $specific_mutations),($total_mutations - $specific_mutations) );
  my $result2=$self->nCk($sequence_length,$total_mutations);
  if($result2>0){
    $probability=($result/$result2);
  }
  return $probability;
}

# validated this nCk is performing equivalently to old vdjfasta version
# to digits limit of computers used
sub nCk {
  my($self,$n,$k)=@_;
  my $result; 

  #   n!
  # -----
  # k!(n-k)!

  my $numerator=1;
  my $denominator=1;
  my $nminusk=$n-$k;

  for (my $x=$n; $x > ($nminusk); $x--){
    $numerator *= $x;
  }

  for(my $x=$k;$x>1;$x--){
    $denominator *= $x;
  }

  $result=($numerator/$denominator);
  return $result;
}

sub aascFvc2m {
  my($self)=@_;
  my $name = $self->{filename};
     $name =~ s/.[Ff][Nn]*[Aa][Ss]*[Tt]*[Aa]*$//;
  my $unique_file = $name . ".unique.fa";
  my $stockfile   = $name . ".stock";
  my $a2m_file    = $name . ".Vh-gs-Vk.a2m";
  my $c2m_file    = $name . ".Vh-gs-Vk.c2m";
  $self->makeAccessionsUnique();
  if($self->{nseqs}>0){
    $self->writeSeqs($unique_file);
    $self->scFvAlign($unique_file,$stockfile);
    `rm $unique_file`;
    $self->stockholmToFasta($stockfile,$a2m_file);
    `rm $stockfile`;

    my $a2m=VDJFasta->new();
       $a2m->loadSeqs($a2m_file);
       $a2m->nTermSmoothing(0);
       $a2m->cTermSmoothing(0);
       $a2m->lambdaNTermSmoothing(0);
       $a2m->writeSeqs($a2m_file);

    $self->cropToModel($a2m_file,$c2m_file);
  }else{
    `touch $c2m_file $a2m_file`;
  }
  return $c2m_file;
}

sub getVHDNA {
  my($self,$dnafile,$a2mfile,$segment)=@_;

  my %region_hash=();

  my $dnafasta=VDJFasta->new();
    $dnafasta->loadSeqs($dnafile);
  my $a2mfasta=VDJFasta->new();
    $a2mfasta->loadSeqs($a2mfile);
  my %dnaAccSeqid=();
    $dnafasta->getAccession2seqidHash(\%dnaAccSeqid);

    # go through each sequence in the a2m alignment (as frameshifts can cause multiple appearances
    # of same accessions.
  my $a2m_counts=$a2mfasta->getSeqCount();
  for(my $a2m_seqid=0;$a2m_seqid<$a2m_counts;$a2m_seqid++){
     my $acc = $a2mfasta->getAccession($a2m_seqid);
     my $dna_seqid = $dnaAccSeqid{$acc};
     # get appropriate coordinates
     my($accession,$orientation,$BeginLen,$h3len,$EndLen,$jSeq,$gsSeq,$vlSeq,$vhSeq)=$a2mfasta->getVariableCoordsFromA2M($a2m_seqid);
     my ($vhseg,$dnaH3,$jcseg,$gsseg,$vlseg) = $dnafasta->getVariableSegmentsFromDNA($dna_seqid,$orientation,$BeginLen,$h3len,$EndLen,$jSeq,$gsSeq,$vlSeq,$vhSeq);
     
     unless($dnaH3 eq ""){
       if($segment eq "VH"){
         $region_hash{$acc}=$vhseg . $dnaH3 . $jcseg;
       }elsif($segment eq "VL"){
         $vlseg=~s/^GGA//;
         $region_hash{$acc}=$vlseg;
       }elsif($segment eq "H3"){
         $region_hash{$acc}=$dnaH3;
       }elsif($segment eq "GS"){
         if($vlseg=~m/^GGA/){
           $gsseg.="GGA";
         }
         $region_hash{$acc}=$gsseg;
       }elsif($segment eq "scFv"){
         $region_hash{$acc}=$vhseg . $dnaH3 . $jcseg . $gsseg . $vlseg;
       }
       #print ">scFv-$a2m_seqid\n";
       #print $vhseg . $dnaH3 . $jcseg . $gsseg . $vlseg . "\n";
       #$region_hash{$acc}=$lvseg . $dnaH3 . $jcseg;
       #print ">VH\n" . $lvseg . "\n>h3\n" . $dnaH3 . "\n>jseg\n" . $jcseg . "\n>gs\n" . $gsseg . "\n>light\n" . $vlseg . "\n";
     }
   }
  return %region_hash;
}


sub getH3dna {
  my($self,$dnafile,$a2mfile)=@_;

  my %h3_hash=();

  my $dnafasta=VDJFasta->new();
    $dnafasta->loadSeqs($dnafile);
  my $a2mfasta=VDJFasta->new();
    $a2mfasta->loadSeqs($a2mfile);
  my %dnaAccSeqid=();
    $dnafasta->getAccession2seqidHash(\%dnaAccSeqid);

    # go through each sequence in the a2m alignment (as frameshifts can cause multiple appearances
    # of same accessions.
  my $a2m_counts=$a2mfasta->getSeqCount();
  for(my $a2m_seqid=0;$a2m_seqid<$a2m_counts;$a2m_seqid++){
     my $acc = $a2mfasta->getAccession($a2m_seqid);
     my $dna_seqid = $dnaAccSeqid{$acc};
     my($accession,$orientation,$BeginLen,$h3len,$EndLen)=$a2mfasta->getH3CoordsFromA2M($a2m_seqid);
     my ($lvseg,$dnaH3,$jcseg) = $dnafasta->getH3fromDNA($dna_seqid,$orientation,$BeginLen,$h3len,$EndLen);
     unless($dnaH3 eq ""){
       $h3_hash{$acc}=$dnaH3;
     }
   }
  return %h3_hash;
}

sub cropToModel {
  my($self,$a2m,$c2m)=@_;

  open(OUT,">$c2m");
  my $seqname="";
  my $seqseq="";

  my $x=-1;
  open MSA, $a2m;
  while (<MSA>) {
   chomp;
   if (/^>/) {
     my $tempname = $_;

     if($seqname ne ""){
       print OUT $seqname . "\n";
       $seqseq=~s/^[\.a-z]*//;
       $seqseq=~s/[\.a-z]*$//;
       print OUT $seqseq  . "\n";
       $seqseq="";
     }

     $seqname=$tempname;
     $x++;
   } else{
     my $tempseq = $_;
     $seqseq .= $tempseq;
   }
  }
  close(MSA);
  $seqseq=~s/^[\.a-z]*//;
  $seqseq=~s/[\.a-z]*$//;
  print OUT $seqname . "\n";
  print OUT $seqseq  . "\n";
  close(OUT);
}

sub stockholmToFasta {
  my($self,$stockfile,$outfile)=@_;
  my %seqs = ();
  my $start=0;

  open(MSA,$stockfile);
  while(<MSA>){
    my $line=$_;
    chop($line);
    if($line =~ m/#/){
      $start=1;
    }
    if($start){
      if( !(($line =~ m/^#/) || ($line =~ m/^ *$/) || ($line =~ m/^\/\//)) ){
        my @fields = split(/  */,$line);
        if(exists($seqs{$fields[0]})){
          $seqs{$fields[0]} .= $fields[1];
        } else {
          $seqs{$fields[0]} = $fields[1];
        }
      }
    }
  }
  close(MSA);

  open(OUT,">$outfile");
  while ( my($key,$value) = each(%seqs)){
    print OUT ">$key\n";
    print OUT "$value\n";
  }
  close(OUT);
}

sub makeAccessionsUnique {
  my($self)=@_;
  my %accessions_already_seen=();
  
  for(my $s=0;$s<$self->{nseqs};$s++){
    ${$self->{headers}}[$s]=~s/ /_/g;
    if(defined($accessions_already_seen{${$self->{headers}}[$s]})){
      while(defined($accessions_already_seen{${$self->{headers}}[$s]})){
        ${$self->{headers}}[$s].= "." . int(rand(100000));
      }
    }
    $accessions_already_seen{${$self->{headers}}[$s]}=1;
  }
}

sub scFvAlign {
  my($self,$seqfile,$outfile)=@_;  
  my $cmd = $self->{hmmalign} 
          . " --allcol --amino " 
          . $self->{VhVkHMM} 
          . " $seqfile > $outfile";
  `$cmd`;
}

sub igScore {
  my($self,$evalue)=@_;

  # score all frames with hmm. add "-A 0" to store alignment of all hits
  my $score_file = $self->{filename} . ".$evalue.score.txt";
  my $cmd = $self->{hmmsearch} 
          . " -E $evalue " 
          . $self->{VhVkHMM} 
          . " " . $self->{filename} 
          . " > $score_file"; 
  `$cmd`;
  # generate hmmscore hits
  my @scored_hits=();
 
  open(FILE,$score_file);
  my @lines=<FILE>;
  close(FILE);
  `rm $score_file`;
 
  chomp(@lines);
  my $recording=0;
  for(my $x=0;$x<scalar(@lines);$x++){
    if($lines[$x]=~m/^>>/){
      $lines[$x]=~s/^>>* *//;
      $lines[$x]=~s/  *$//;
      push @scored_hits,$lines[$x];
    }
  }
  return @scored_hits;
}

sub addHeaderAnnotation {
  my($self,$input_feature)=@_;

  # input_feature could be an annotation_file, or a hash of features
  my $acc2cdr;
  if(-f $input_feature){
    my $annotation_file=$input_feature;
    $acc2cdr = {};
    # load file (accession	annotation)
    open(FILE,$annotation_file);
    my @lines=<FILE>;
    close(FILE);
    chomp(@lines);
    for(my $i=0;$i<scalar(@lines);$i++){
      my @fields=split(/\t/,$lines[$i]);
      if(defined($fields[1])){
        $fields[1]=~s/ *$//;
      }else{
        $fields[1]="";
      }
      $fields[0]=~s/^> *//;
      $fields[0]=~s/ *//;
      $acc2cdr->{$fields[0]}=$fields[1];
    }
  }else{
    $acc2cdr = $input_feature;
  }

  # apply annotations
  for(my $n=0;$n<$self->{nseqs};$n++){
    ${$self->{headers}}[$n] .= ";";
    if(defined($acc2cdr->{$self->getAccession($n)})){
      ${$self->{headers}}[$n] .= $acc2cdr->{$self->getAccession($n)};
    }
  }
}

sub getBlastLines {
  my($self,$seqfile,$blastdb,$eval,$blast_type)=@_;
  
  my $command="";
  if($blast_type eq "blastn"){
    $command = $self->{blast} . " -task blastn-short -query $seqfile -db $blastdb -outfmt 6 -evalue $eval";
  }elsif($blast_type eq "tblastn"){
    $command =  $self->{tblastn} . " -word_size 3 -query $seqfile -db $blastdb -outfmt 6 -evalue $eval";
  }
  my @lines=`$command`;
  return @lines;
}

sub getAccession2seqidHash {
  my($self,$hash)=@_;
  for(my $n=0;$n<$self->{nseqs};$n++){
    $$hash{$self->getAccession($n)}=$n;
  }
}

sub setHeader {
  my($self,$c,$new_header)=@_;
  unless(defined($c)){
    return "";
  }
  if(defined(${$self->{headers}}[$c])){
    ${$self->{headers}}[$c]=$new_header;
  }else{
    print "You are trying to setHeader for sequence $c that doesn't exist\n";
    exit;
  }
}

sub getHeader {
  my($self,$c)=@_;
  unless(defined($c)){
    return "";
  }
  if(defined(${$self->{headers}}[$c])){
    return ${$self->{headers}}[$c];
  }else{
    print "You are trying to getHeader for sequence $c that doesn't exist\n";
    exit;
  } 
}

sub getHeaderField {
  my($self,$seqid,$field)=@_;
  my $this_header = $self->getHeader($seqid);
  my @fields=split(/;/,$this_header);
  chomp(@fields);
  if(defined($fields[$field])){
    return $fields[$field];
  }else{
    return "";
  }
}


sub getAccession {
  my($self,$c)=@_;
  my $header=$self->getHeader($c);
  $header=~s/ .*//;
  $header=~s/;.*//;
  $header=~s/^>//;
  $header=~s/_frame_[0123-]* *$//;
  return $header;
}

sub loadSeqs {
  my ($self,$file)=@_;
  if(-e $file){
    $self->{filename}=$file;
  } else {
    print "File does not exist.\n";
    return 0;
  }

  # VDJFASTA_USE_TEMPFILES indicates that we should try to use file-backed
  # (i.e. tied) arrays, usually because the caller knows that the amount
  # of memory needed to hold the files in memory would be too much for
  # the system we're running on.
  #
  # If the tied arrays can't be created, we go ahead and try in memory anyway :-)
  #
  if (exists($ENV{VDJFASTA_USE_TEMPFILES})) {
    my $basefilename = basename($file);
    my $headerdb = $basefilename.'_'.$$.'_headers.db';
    my $sequencedb = $basefilename.'_'.$$.'_sequences.db';
    my @header_array;
    my @sequence_array;

    # put DB_File temporary files in the current working directory
    # otherwise they go into /var/tmp or /tmp, and can fill root filesystems
    my $cwd = getcwd();
    if (open(C, "> DB_CONFIG")) {
      print C "set_tmp_dir $cwd\n";
      close(C);
    }

    if (tie @header_array, 'DB_File', $headerdb, O_RDWR|O_CREAT, 0644, $DB_RECNO) {
      $self->{headerdb} = $headerdb; # so it gets removed by the destructor
      if (tie @sequence_array, 'DB_File', $sequencedb, O_RDWR|O_CREAT, 0644, $DB_RECNO) {
        $self->{sequencedb} = $sequencedb;

        # copy any existing sequences into the new array
        push(@header_array, @{$self->{headers}})
          if $self->{headers};
        push(@sequence_array, @{$self->{sequence}})
          if $self->{sequence};

        # replace the array references in the VDJFasta object
        $self->{headers} = \@header_array;
        $self->{sequence} = \@sequence_array;
      }
      else {
        print "could not create file-based array $sequencedb: $!\n";
      }
    }
    else {
      print "could not create file-based array $headerdb: $!\n";
    }
  }

  # use of AnyUncompress will allow us to transparently handle
  # files whether compressed or not
  my $filehandle = new IO::Uncompress::AnyUncompress $self->{filename};
  if (not $filehandle) {
    print "could not open input file $self->{filename}: $AnyUncompressError\n";
    return 0;
  }

  # load seqs into filename, headers and sequence
  # will check for fastq file input to transform to fasta headers
  # and throw away the quality scores
  my $x=-1;
  while (<$filehandle>) {
    chomp;
    if (/^\+/) {
      # preceeds the fastq quality line
      $filehandle->getline();
      next;
    }
    s/^@/>/; # transform fastq file header
    if (/^>/) {
      push @{$self->{headers}},$_;
      $x++;
    } else {
      ${$self->{sequence}}[$x] .= $_;
    }
  }
  $filehandle->close();
  $x++;
  $self->{nseqs}=$x;
}

sub addSeq {
  my($self,$header,$sequence)=@_;
  push @{$self->{headers}},$header;
  push @{$self->{sequence}},$sequence;
  $self->{nseqs}++;
}

sub getSeqLen {
  my($self,$c)=@_;
  my $seq=${$self->{sequence}}[$c];
  my $count=length($seq);
  return $count;
}

sub printSeqs {
  my($self)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    $self->printSeq($s);
  }
}

sub subsampleViableCloneHash {
  my($self,$depth)=@_;
  my @id_array=();
  my %clone_hash=();
  unless(defined($depth)){
    $depth=$self->{nseqs};
  }
  if($depth>$self->{nseqs}){
    $depth=$self->{nseqs};
  }
  for(my $s=0;$s<$self->{nseqs};$s++){
    push @id_array,$s;
  }
  $self->fisher_yates_shuffle(\@id_array);
  my $counts=0;
  for(my $s=0;($counts<$depth) and ($s<$self->{nseqs});$s++){
    # monkey
    my $vseg = $self->getHeaderField($id_array[$s],1);
    my $jseg = $self->getHeaderField($id_array[$s],3);
    my @vseg_fields = split(/ /,$vseg);
    my @jseg_fields = split(/ /,$jseg);
    if(scalar(@vseg_fields)>1){
      if(scalar(@jseg_fields)>1){
        $vseg=$vseg_fields[0];
        $jseg=$jseg_fields[0];
        my $cdr3_field=4;
        my $cdr3 = $self->getHeaderField($id_array[$s],$cdr3_field);
        unless(defined($clone_hash{$vseg . "_" . $jseg . "_" . $cdr3})){
          $counts++;
        }
        $clone_hash{$vseg . "_" . $jseg . "_" . $cdr3}=1;
      }
    }
  }
  return %clone_hash;
}

sub shuffle {
  my($self,$depth)=@_;
  my @id_array=();
  unless(defined($depth)){
    $depth=$self->{nseqs};
  }
  if($depth>$self->{nseqs}){
    $depth=$self->{nseqs};
  }
  for(my $s=0;$s<$self->{nseqs};$s++){
    push @id_array,$s;
  }
  $self->fisher_yates_shuffle(\@id_array);  
  for(my $s=0;$s<$depth;$s++){
    $self->printSeq($id_array[$s]);
  }
}

# randomly permutate @array in place
sub fisher_yates_shuffle {
  my ($self,$array) = @_;
  my $i = @$array;
  while ( --$i ) {
    my $j = int rand( $i+1 );
    @$array[$i,$j] = @$array[$j,$i];
  }
}

sub printClean {
  my($self)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    print ${$self->{headers}}[$s] . "\n";
    my $seq = ${$self->{sequence}}[$s];
       $seq = uc($seq);
       $seq =~ s/-//g;
       $seq =~ s/\.//g;
       $seq =~ s/_//g;
    print $seq . "\n";
  }
}

sub printReverse {
  my($self)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $reverse=$self->getReverseStrand($s); 
    print ${$self->{headers}}[$s] . "\n";
    print $reverse . "\n";
  }
}

sub filterPseudoGenes {
  # remove pseudogenes
  my($self)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my($trans_head,$aa)=$self->dna2aa($s,0);
    unless($aa=~m/[\-_X]/){
      print ${$self->{headers}}[$s] . "\n";
      print ${$self->{sequence}}[$s] . "\n";
    }
  }
}

sub printSeqRange {
  my($self,$start,$stop)=@_;
  for(my $s=$start;$s<$stop;$s++){
    $self->printSeq($s);
  }
}

sub printSeq {
  my($self,$s)=@_;
  print ${$self->{headers}}[$s] . "\n";
  print ${$self->{sequence}}[$s] . "\n";
}

sub printReverseSeq {
  my($self,$s)=@_;
  my $reverse=$self->getReverseStrand($s);
  print ${$self->{headers}}[$s] . "\n";
  print $reverse . "\n";
}

sub printUniqueSeqs {
  my($self)=@_;
  my $id=1;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $header=">$id-";
    $id++;
    my $header_part2 = ${$self->{headers}}[$s];
       $header_part2 =~ s/^>//;
       $header .= $header_part2;
    print $header . "\n";
    print ${$self->{sequence}}[$s] . "\n";
  }
}

sub getCDRHashFromHeaders {
  my($self,$cdr3)=@_;
  # assumes the headers have been encoded with H3 in field 5
  my $cdr3_field=4;
  if($cdr3 eq "L3"){
    $cdr3_field=5;
  }
  print "Fields is $cdr3_field\n";
  my %cdr_hash=();

  for(my $s=0;$s<$self->{nseqs};$s++){
    my @fields=split(/;/,${$self->{headers}}[$s]);
    if(defined($fields[$cdr3_field])){
      if(length($fields[$cdr3_field]) > 0){
        if(defined($cdr_hash{$fields[$cdr3_field]})){
          $cdr_hash{$fields[$cdr3_field]}.="," . $s;
        }else{
          $cdr_hash{$fields[$cdr3_field]} = $s;
        }
      }
    }
  }
  return %cdr_hash;
}

sub printSeqSubsetbyList {
  my($self,$accession_list,$outfile,$no_frames)=@_;

  my %selected_accessions=();
  my @lines=@$accession_list;

  open(OUTFILE,">$outfile");

  for(my $x=0;$x<scalar(@lines);$x++){
    
    my @data=split(/\t/,$lines[$x]);
    $data[0]=~s/ .*//;
    #print "monkey Setting up |" . $data[0] . "|\n";
    $selected_accessions{$data[0]}=1;
  } 
 
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $acc=${$self->{headers}}[$s];
    $acc=~s/ .*//;
    $acc=~s/^>//;
    if(defined($no_frames)){
      if($no_frames == 1){
        $acc=~s/_frame.*//;
        $acc=~s/_PID=\d+_EID=\d+.*//;
        $acc=~s/;.*//;
      }
    }
    #print "monkey - Checking on |$acc|\n";
    if($selected_accessions{$acc}){
      print OUTFILE ${$self->{headers}}[$s] . "\n";
      print OUTFILE ${$self->{sequence}}[$s] . "\n";
    }
  }
  close(OUTFILE);
}

sub getSeq{
  my($self,$seq)=@_;
  return (${$self->{headers}}[$seq],${$self->{sequence}}[$seq]);
}

sub getSeqCount{
  my($self)=@_;
  return $self->{nseqs};
}

sub getFreqvSHM {
  my($self)=@_;

  my %clone_counts=();
  my %clone_pids=();
  my $counts=0;

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $this_shmpid = $self->getPercentSHM($s,200);
    if($this_shmpid =~ m/[0-9]/){
      my $clone_string = $self->getSeqVJH3($s);
      unless($clone_string eq ""){
        if(defined($clone_counts{$clone_string})){
          $clone_counts{$clone_string}++;
        }
        if(defined($clone_pids{$clone_string})){
          $clone_pids{$clone_string}+=$this_shmpid;
        }
        $counts++;
      }
    }
  }

  my @keys = keys %clone_counts;
  my $output="";
  for(my $k=0;$k<scalar(@keys);$k++){
     my $this_pid        = int($clone_pids{$keys[$k]} / $clone_counts{$keys[$k]});
     my $this_frequency = ($clone_counts{$keys[$k]} / $counts);
    $output .= $this_frequency . "\t" . $this_pid . "\n";
  }
  return $output;
}

sub getSeqVJH3 {
  my($self,$s)=@_;
  
  my @header_fields=$self->splitSeqHeader($s);
  my $clone="";
  if(defined($header_fields[4])){
    my @vseg_data=split(/ /,$header_fields[1]);
    my @jseg_data=split(/ /,$header_fields[3]);
    my $h3=$header_fields[4];
    if(scalar(@vseg_data)==3){
      if(scalar(@jseg_data)==3){
        $vseg_data[0]=~s/\*.*//;
        $jseg_data[0]=~s/\*.*//;
        $clone=$vseg_data[0] . "_" . $jseg_data[0] . "_" . $h3;
      }
    }
  }
  return $clone;
}

sub getFrequencyHash {
  my($self)=@_;
  
  my @features = ( "5e-1", "1e-1", "5e-2", "1e-2", "5e-3", "1e-3", "5e-4", "1e-4", "5e-5", "1e-5", "5e-6", "1e-6");

  my %cloneCounts = $self->getCloneCountsHash();

  my @clones = keys %cloneCounts;

  my $sum=0;
  for(my $c=0;$c<scalar(@clones);$c++){
    $sum+= $cloneCounts{$clones[$c]};
  }
  
  my %outData=("5e-1",0, "1e-1",0, "5e-2",0, "1e-2",0, "5e-3",0, "1e-3",0, "5e-4",0, "1e-4",0, "5e-5",0, "1e-5",0, "5e-6",0, "1e-6",0);
  for(my $c=0;$c<scalar(@clones);$c++){
    my $freq=$cloneCounts{$clones[$c]} / $sum;
    my $log = $self->getDiversityCategory($freq); #int($self->log10($freq));
    if(defined($outData{$log})){
      $outData{$log}+= (1/scalar(@clones));  
    }else{
      $outData{$log}=1/scalar(@clones);
    }
  }
  return %outData;
}

sub getDiversityCategory {
  my($self,$value)=@_;
  if($value > 5e-1){
    return "5e-1";
  }elsif($value > 1e-1){
    return "1e-1";
  }elsif($value > 5e-2){
    return "5e-2";
  }elsif($value > 1e-2){
    return "1e-2";
  }elsif($value > 5e-3){
    return "5e-3";
  }elsif($value > 1e-3){
    return "1e-3";
  }elsif($value > 5e-4){
    return "5e-4";
  }elsif($value > 1e-4){
    return "1e-4";
  }elsif($value > 5e-5){
    return "5e-5";
  }elsif($value > 1e-5){
    return "1e-5";
  }elsif($value > 5e-6){
    return "5e-6";
  }else{
    return "1e-6";
  }
}

sub getSegmentHash {
 my($self,$segment,$mutation_status)=@_;

  my $min_mut=0;
  my $max_mut=60;
  if(defined($mutation_status)){
    if($mutation_status eq "naive"){
      $max_mut=4;
    }elsif($mutation_status eq "memory"){
      $min_mut=5;
    }elsif($mutation_status eq "redundant"){
      $min_mut=0;
      $max_mut=60;
    }elsif($mutation_status eq "expanded"){
      $min_mut=0;
      $max_mut=60;
    }
  }

  my %unique_clones=();
  my %count_clones=();

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $clone = $self->getSeqVJH3($s);
    if(defined($mutation_status)){
      my($this_vgene,$this_length,$this_mutations)=split(/ /,$self->getHeaderField($s,1));
      if($this_length>30){  # 100 is best
        if($this_mutations>=$min_mut){
          if($this_mutations<=$max_mut){
            $unique_clones{$clone}=1;
            # get clone counts for redundant measurements
            if(defined($count_clones{$clone})){
              $count_clones{$clone}++;
            }else{
              $count_clones{$clone}=1;
            }
          }
        }
      }
    }else{
      $unique_clones{$clone}=1;
      # get clone counts for redundant measurements
      if(defined($count_clones{$clone})){
        $count_clones{$clone}++;
      }else{
        $count_clones{$clone}=1;
      }
    }
  }

  # now get the total count for expanded clones
  my $expanded_total=0;
  my @clone_names=keys %unique_clones;
  for(my $a=0;$a<scalar(@clone_names);$a++){
    if($count_clones{$clone_names[$a]}>1){
      $expanded_total+=$count_clones{$clone_names[$a]};
    }
  }  



  my %outDir=();


  for(my $c=0;$c<scalar(@clone_names);$c++){
    my($vgene,$jgene,$h3)=split(/_/,$clone_names[$c]);

    my $gene_of_interest="";
    if($segment eq "IGHV"){
      $gene_of_interest=$vgene;
    }elsif($segment eq "IGHJ"){
      $gene_of_interest=$jgene;
    }
    if(defined($outDir{$gene_of_interest})){
      if($mutation_status eq "redundant"){
        $outDir{$gene_of_interest}+=($count_clones{$clone_names[$c]}/$self->{nseqs});
      }elsif($mutation_status eq "expanded"){
        if($count_clones{$clone_names[$c]}>1){
          $outDir{$gene_of_interest}+=($count_clones{$clone_names[$c]}/$expanded_total);
        } 
      }else{
        $outDir{$gene_of_interest}+=(1/scalar(@clone_names));
      }
    }else{
      if($mutation_status eq "redundant"){
        $outDir{$gene_of_interest}=($count_clones{$clone_names[$c]}/$self->{nseqs});
      }elsif($mutation_status eq "expanded"){
        if($count_clones{$clone_names[$c]}>1){                  
          $outDir{$gene_of_interest}=($count_clones{$clone_names[$c]}/$expanded_total);
        }
      }else{
        $outDir{$gene_of_interest}=(1/scalar(@clone_names));
      }
    }    
  }

  return %outDir;
}


sub log10 {
  my($self,$data)=@_;
  if($data eq 0){
    return -10;
  }else{
    return log($data)/log(10);
  }
}
sub getSHMFrequency {
  my($self)=@_;
  # bin from 100 down to 70
  my %hist=( 100,0,     99,0,   98,0,   97,0,   96,0,
              95,0,     94,0,   93,0,   92,0,   91,0,
              90,0,     89,0,   88,0,   87,0,   86,0,
              85,0,     84,0,   83,0,   82,0,   81,0,
              80,0,     79,0,   78,0,   77,0,   76,0,
              75,0,     74,0,   73,0,   72,0,   71,0);

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $this_shmpid = $self->getPercentSHM($s,200);
    if($this_shmpid =~ m/[0-9]/){
      my $shmpid = int(100 * $this_shmpid);
      if(defined($hist{$shmpid})){
        $hist{$shmpid}++;
      }
    }
  }
  return %hist;
}
sub getSHMHistogram {
  my($self)=@_;
  # bin from 100 down to 70
  my %hist=( 100,0,	99,0,	98,0,	97,0,	96,0,
	      95,0,	94,0,	93,0,	92,0,	91,0,
	      90,0,     89,0,   88,0,   87,0,   86,0,
              85,0,     84,0,   83,0,   82,0,   81,0,
              80,0,     79,0,   78,0,   77,0,   76,0,
              75,0,     74,0,   73,0,   72,0,   71,0);
  
  my $counts=0;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $this_shmpid = $self->getPercentSHM($s,200);
    if($this_shmpid =~ m/[0-9]/){
      my $shmpid = int(100 * $this_shmpid);
      if(defined($hist{$shmpid})){ 
        $hist{$shmpid}++;
        $counts++;
      }
    }
  }
  my @pids=keys %hist;
  for(my $p=0;$p<scalar(@pids);$p++){
    if($counts>0){
      $hist{$pids[$p]} /= $counts;
    }
  }
  return %hist;
}

sub getSeqLengthHistogram {
  my($self)=@_;
  # bin by 25, up to 800

  my @hist=( 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	     0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	     0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0);

  my $counts=0;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $count = $self->getSeqLen($s);
    my $position = int($count/25);
    #print "Got count $count and position $position\n";
    if($position>1800){
      $position=1800;
    }
    $hist[$position]++;
    $counts++;
  }
  my $values= 0;
  if($counts>0){
    $values=(int(($hist[0]/$counts)*1000))/10;
  }

  # now populate the hash
  my %dataHash=();

  for(my $p=0;$p<scalar(@hist);$p++){
    my $percent = 0;
    if($counts>0){
      if(defined($hist[$p])){
        $percent=(int(($hist[$p]/$counts)*1000))/10;
      }
    }
    $dataHash{(($p+1)*25)}=$percent;
    $values.="\t" . $percent;
  }

  return %dataHash;
}

sub selectOnlyH3 {
  my($self,$outfile)=@_;
  open(FILE,">$outfile");
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $h3 = $self->getHeaderField($s,4);
    if(length($h3)>2){
      print FILE ${$self->{headers}}[$s]  . "\n";
      print FILE ${$self->{sequence}}[$s] . "\n";
    }
  }
  close(FILE);
}

sub writeSeq {
  my($self,$seq,$file)=@_;
  open(FILE,">$file");
  print FILE ${$self->{headers}}[$seq]  . "\n";
  print FILE ${$self->{sequence}}[$seq] . "\n";
  close(FILE);
  return $file;
}

sub appendSeq {
  my($self,$seq,$file)=@_;
  open(FILE,">>$file");
  print FILE ${$self->{headers}}[$seq]  . "\n";
  print FILE ${$self->{sequence}}[$seq] . "\n";
  close(FILE);
  return $file;
}

sub writeSeqs {
  my($self,$file)=@_;
  open(FILE,">$file");
  for(my $s=0;$s<$self->{nseqs};$s++){
    print FILE ${$self->{headers}}[$s]  . "\n";
    print FILE ${$self->{sequence}}[$s] . "\n";
  }
  close(FILE);
}

sub writeSeqRange {
  my($self,$file,$start,$stop)=@_;
  open(FILE,">$file");
  if($stop>$self->{nseqs}){
    $stop=$self->{nseqs};
  }
  for(my $s=$start;$s<$stop;$s++){
    print FILE ${$self->{headers}}[$s]  . "\n";
    print FILE ${$self->{sequence}}[$s] . "\n";
  }
  close(FILE);
}

sub getMutatedSeq {
  my($self,$seq,$nmutations,$trim5prime,$trim3prime)=@_;
  my $mutseq=${$self->{sequence}}[$seq];

  # trime 5 prime and 3 prime
  if($trim5prime){
    $mutseq=substr($mutseq,$trim5prime);
  }
  if($trim3prime){
    $mutseq=substr($mutseq,0,(0 - $trim3prime));
  }

  # perform mutations
  my @chars=('C','A','T','G');
  my @res=split(/ */,$mutseq);
  my $mut_down=$nmutations;
  while($mut_down>0){
    my $position=int(rand(scalar(@res)));
    my $new=$chars[int(rand(4))];
    while($new eq uc($res[$position])){
      $new=$chars[int(rand(4))];
    }
    $res[$position]=$new;
    $mut_down--;
  }
  $mutseq=join('',@res);

  return $mutseq;
}

sub writeMutatedSeq {
  my($self,$seq,$file,$nmutations,$trim5prime,$trim3prime)=@_;

  my $mutseq=$self->getMutatedSeq($seq,$nmutations,$trim5prime,$trim3prime);
  
  open(FILE,">$file");
  print FILE ${$self->{headers}}[$seq]  . " with $nmutations mutations\n";
  print FILE $mutseq . "\n";
  close(FILE);
  return $file;
}



sub hmmSearch {
  my($self,$seqfile,$hmm,$eval)=@_;
  my $hmmsearch=$self->{hmmsearch};
  my @results=`$hmmsearch -E $eval $hmm $seqfile`;
  my @frames=();
  my $print=1;
  for(my $x=0;$x<scalar(@results);$x++){
    if($results[$x]=~m/^frame/){
      if($results[$x]=~m/domain/){
        my $line=$results[$x];
        $line=~s/^frame//;
        $line=~s/: .*//;
        chomp($line);
        push @frames,$line;
      }
    }
  }
  return @frames;
}

sub getPercentGc {
  my($self,$c)=@_;
  my ($header,$oligo)=$self->getSeq($c);
  my @chars=split(/ */,$oligo);
  my $gc=0;
  for(my $g=0;$g<scalar(@chars);$g++){
    if($chars[$g]=~m/[gcGC]/){
      $gc++;
    }
  }
  my $percent_gc=0;
  if(scalar(@chars)>0){
    $percent_gc=((  int(1000 * ($gc/scalar(@chars))) )/10);
  }
  return $percent_gc;
}

sub rounded {
  my ($self,$value)=@_;
  my $rounded=((  int(10 * $value) )/10);
  return $rounded;
}

sub getReverseStrand {
  my($self,$c)=@_;

  my $seq=uc($self->getSequence($c)); 
  my @nucs=split(/ */,$seq);
  my $len=scalar(@nucs);
  my $reverse="";
  for(my $x=($len-1);$x>=0;$x--){
    if($nucs[$x] eq "T"){
      $reverse.="A";
    }elsif($nucs[$x] eq "A"){
      $reverse.="T";
    }elsif($nucs[$x] eq "G"){
      $reverse.="C";
    }elsif($nucs[$x] eq "C"){
      $reverse.="G";
    }elsif($nucs[$x] eq "N"){
      $reverse.="N";
    }else{
      $reverse.="X";
    }
  }
  return $reverse;
}

sub getSequence {
  my($self,$seq)=@_;
  if(defined($seq)){
    if(defined(${$self->{sequence}}[$seq])){
      return ${$self->{sequence}}[$seq];
    }else{
      return "";
    }
  }else{
    print "ERROR: getSequence passed an undefined seq from somewhere\n";
    return "";
  }
}

sub setHMM {
  my($self,$newHMM)=@_;
  $self->{VhVkHMM}=$newHMM;
}

sub setSequence {
  my($self,$c,$sequence)=@_;
  return ${$self->{sequence}}[$c]=$sequence;
}

sub codon2aa {
  my($codon,$readthrough_amber_stop) = @_;

  if(length($codon)<3){
    return "";
  }

  if(!defined($readthrough_amber_stop)){
    $readthrough_amber_stop=1;
  }

  if    ( $codon =~ /TC[ACGTN]/i ) { return 'S' } # Serine
  elsif ( $codon =~ /TT[CT]/i ) { return 'F' } # Phenylalanine
  elsif ( $codon =~ /TT[AG]/i ) { return 'L' } # Leucine
  elsif ( $codon =~ /TA[CT]/i ) { return 'Y' } # Tyrosine
  elsif ( $codon =~ /TAA/i ) { return 'X' } # Stop
  elsif ( ($codon =~ /TAG/i) && ($readthrough_amber_stop == 0) ) { return 'X' } # Stop
  elsif ( ($codon =~ /TAG/i) && ($readthrough_amber_stop == 1) ) { return 'Q' } # Stop
  elsif ( $codon =~ /TG[CT]/i ) { return 'C' } # Cysteine
  elsif ( $codon =~ /TGA/i ) { return 'X' } # Stop
  elsif ( $codon =~ /TGG/i ) { return 'W' } # Tryptophan
  elsif ( $codon =~ /CT[ACGTN]/i ) { return 'L' } # Leucine
  elsif ( $codon =~ /CC[ACGTN]/i ) { return 'P' } # Proline
  elsif ( $codon =~ /CA[CT]/i ) { return 'H' } # Histidine
  elsif ( $codon =~ /CA[AG]/i ) { return 'Q' } # Glutamine
  elsif ( $codon =~ /CG[ACGTN]/i ) { return 'R' } # Arginine
  elsif ( $codon =~ /AT[ACT]/i ) { return 'I' } # Isoleucine
  elsif ( $codon =~ /ATG/i ) { return 'M' } # Methionine
  elsif ( $codon =~ /AC[ACGTN]/i ) { return 'T' } # Threonine
  elsif ( $codon =~ /AA[CT]/i ) { return 'N' } # Asparagine
  elsif ( $codon =~ /AA[AG]/i ) { return 'K' } # Lysine
  elsif ( $codon =~ /AG[CT]/i ) { return 'S' } # Serine
  elsif ( $codon =~ /AG[AG]/i ) { return 'R' } # Arginine
  elsif ( $codon =~ /GT[ACGTN]/i ) { return 'V' } # Valine
  elsif ( $codon =~ /GC[ACGTN]/i ) { return 'A' } # Alanine
  elsif ( $codon =~ /GA[CT]/i ) { return 'D' } # Aspartic Acid
  elsif ( $codon =~ /GA[AG]/i ) { return 'E' } # Glutamic Acid
  elsif ( $codon =~ /GG[ACGTN]/i ) { return 'G' } # Glycine
  elsif ( $codon =~ m/N/ ) { return 'Z' }   # N - unknown character
  else { return "X"; } # bad codon
}

sub dna2aa {
  my($self,$seq,$frame,$bypass_amber_stop)=@_;

  unless(defined($bypass_amber_stop)){
    $bypass_amber_stop=0;
  }

  # ($header,$translation)=$h3->($seq,$frame);
  # header

  my $header=$self->getHeader($seq);
     $header=~s/ .*//;
     $header.="_frame_$frame";

  # frame options: 
  my @nucs=();
  my $startframe=$frame;
  if($frame>=0){
    @nucs=split(/ */,$self->getSequence($seq));
  }else{
    @nucs=split(/ */,$self->getReverseStrand($seq));
    $startframe=-1 - $frame;
  }

  # get translation
  my $translation="";
  for(my $n=$startframe;$n<scalar(@nucs);$n+=3){
    if(defined($nucs[($n+2)])){
      my $codon=$nucs[$n] . $nucs[($n+1)] . $nucs[($n+2)];
      $translation.=codon2aa($codon,$bypass_amber_stop);
    }
  }
  # deal with short sequences
  if(length($translation)<1){
    $translation="Z";
  }
  return($header,$translation);
}

sub translateSeq {
  my($self,$dnaseq)=@_;
  my @nucs=split(/ */,$dnaseq);
  my $translation="";
  for(my $n=0;$n<scalar(@nucs);$n+=3){
    if(defined($nucs[($n+2)])){
      my $codon=$nucs[$n] . $nucs[($n+1)] . $nucs[($n+2)];
      $translation.=codon2aa($codon);  
    }
  }
  return $translation;
}

sub translateAllFramesToFile {
  my($self,$outfile)=@_;
  open(FILE,">$outfile");
  for(my $s=0;$s<$self->{nseqs};$s++){ 
    for(my $f=-3;$f<3;$f++){
      my($header,$seq)=$self->dna2aa($s,$f);
      print FILE $header .  "\n";
      print FILE $seq . "\n";  
    }
  }
  close(FILE);  
}

sub restackRange {
  my($self,$range_start,$range_stop)=@_;
  # restack the residues within this range
  # n_right_stack_residues: the number of residues to stack to the end of the range
  # all others are stacked from left to right

  # first determine the longest segment in the range
  my $longest=0;
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $segment=$self->getHMMCOLRange($n,$range_start,$range_stop);
       $segment=~s/\.*//g;
       $segment=~s/^[A-Z-]//;
       $segment=~s/[A-Z-]$//;
    if(length($segment)>$longest){
      $longest=length($segment);
    }
  }
  
  # next proceed through all sequences, restacking the range and filling in the extra
  # space with enough dots to fit the $longest sequence
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $segment_before=$self->getHMMCOLRange($n,0,$range_start);
    my $segment_after =$self->getHMMCOLRange($n,$range_stop,length($self->getSeq($n)));
       $segment_before=~s/[a-z\.]*$//;
       $segment_after=~s/^[a-z\.]*//;

    # obtain segment
    my $segment=$self->getHMMCOLRange($n,$range_start,$range_stop);
       $segment=~s/\.*//g;
       $segment=~s/^[A-Z-]//; # remove boundary position, as it is represented in segment_before
       $segment=~s/[A-Z-]$//; # remove boundary position, as it is represented in segment_after
    # get hmmcol count
    my $hmmcols=$segment;
       $hmmcols=~s/[a-z]//g;
    my $hmmcol_count=length($hmmcols);
    # split in half
    my $halfway = int ( (length($segment) + 1)/2 ); #half way, round up
    my $segment_first_half=substr($segment,0,$halfway);
    my $segment_second_half=substr($segment,$halfway); # +1 I think
    # create middle gap
    my $middle_gap="";
    my $length_difference=$longest - length($segment);
    for(my $x=0;$x<$length_difference;$x++){
      $middle_gap.=".";
    }
    # create the new string
    my $structurally_correct_loop=lc($segment_first_half . $middle_gap . $segment_second_half);
    # reassign the HMM columns
       $structurally_correct_loop=~s/\-/./g;
    my @chars=split(/ */,$structurally_correct_loop);
    my $hmmcol_halfway=int ($hmmcol_count/2);
    for(my $p=(scalar(@chars)-1);(($hmmcol_halfway>0)&&($p>=0));$p--){
      $chars[$p]=uc($chars[$p]);
      $chars[$p]=~s/\./-/;
      $hmmcol_halfway--;
      $hmmcol_count--;
    }
    for(my $p=0;(($hmmcol_count>0)&&($p<scalar(@chars)));$p++){
      $chars[$p]=uc($chars[$p]);
      $chars[$p]=~s/\./-/;
      $hmmcol_count--;
    }
    $structurally_correct_loop="";
    for(my $x=0;$x<scalar(@chars);$x++){
      $structurally_correct_loop.=$chars[$x];
    }
    #print $segment_before . "\t" .  $structurally_correct_loop . "\t" . $segment_after . "\n";
    $self->setSequence($n,($segment_before . $structurally_correct_loop . $segment_after)); 
  }
}

sub getHMMCOLRange{
  my($self,$s,$startcol,$stopcol)=@_;
  my $sequence=$self->getSeq($s);
  my $subseq="";
  my @residues=split(/ */,$sequence);
  my $current_hmmcol=0;
  for(my $x=0;( ($x<scalar(@residues)) && ($current_hmmcol<=$stopcol));$x++){
    if( ($current_hmmcol >= $startcol) && ($current_hmmcol <= $stopcol) ){
      $subseq.=$residues[$x];
    }
    if($residues[$x]=~m/[A-Z-]/){
      $current_hmmcol++;
    }
  }
  return $subseq;
}

sub getLightChainInterfaceResidues {
  my($self,$seqid)=@_;

  my $seq  = $self->getHMMCOLRange($seqid,158,161);
     $seq .= $self->getHMMCOLRange($seqid,163,163);
     $seq .= $self->getHMMCOLRange($seqid,165,165);
     $seq .= $self->getHMMCOLRange($seqid,169,173);
     $seq .= $self->getHMMCOLRange($seqid,176,178);

     $seq .= $self->getHMMCOLRange($seqid,182,182);

     $seq .= $self->getHMMCOLRange($seqid,214,214);
     $seq .= $self->getHMMCOLRange($seqid,226,228);

  return $seq;
}

sub getHeavyChainInterfaceResidues {
  my($self,$seqid)=@_;
 
  my $seq  = $self->getHMMCOLRange($seqid,32,32);
     $seq .= $self->getHMMCOLRange($seqid,34,34);
     $seq .= $self->getHMMCOLRange($seqid,36,36);
     $seq .= $self->getHMMCOLRange($seqid,38,38);
     $seq .= $self->getHMMCOLRange($seqid,42,46);
     $seq .= $self->getHMMCOLRange($seqid,49,49);
     $seq .= $self->getHMMCOLRange($seqid,57,60);
     $seq .= $self->getHMMCOLRange($seqid,93,93);
     $seq .= $self->getHMMCOLRange($seqid,102,103);

  return $seq;
}

sub parseHeavyChainFrameworks {
  my($self,$seqid,$form)=@_;

#                                        ***** (kabat)                      ************************ (kabat)
#                                  * ********* (Aho)                        *************************** (Aho)
#                                  * ******    (chothia)                      ***** (chothia)
#(DVQLLQSGGG.LVQP.GG..AL..GLSCAA.S.G,FTFRNFAMS,W....VRQ.A..P..G.K..G.LE.W,VSAITGSGFNTYYA,.DFA.K....G.RF..TVSRDSSTN..TVHLHLNSLRAEDTALYYC,ARSPkeyy......dvyggPDY,WGQGTLVTVSS)
#(EVQLVESGGG.LVQP.GG..SL..RLSCAA.S.G,FTFSSYAMS,W....VRQ.A..P..G.K..G.LE.W,VSAISGSGGSTYYA,.DSV.K....G.RF..TISRDNSKN..TLYLQMNSLRAEDTAVYYC,AKRRtaa.........wyyFDY,WGQGTLVTVSS)
#(EVQLVDSGGG.LVQP.GG..SL..RLSCAA.S.G,FTFSTFAMS,W....VRQ.G..P..G.K..G.LE.W,VSTLSGSARGTYYA,.DSV.K....G.RF..TISRDNSRN..TLYLQMNSLRAEDTAVYYC,AKLSglrl.......geyfFDK,WGQGTLVTVSS)

# sidestack first, then run the command
# H1 move back by one extra residue on N-terminal side
# kabat - scoot 5 N-terminal residues forward
# chothia - scott 3 C-terminal residues back
                                               
# H2 move forward by two residues on N-terminal side
# H2 move forward 7 residues on C-terminal side
# kabat trim back 2 residues on C-terminal side
# chothia - move 2 residues forward on N-term and back 17 residues 

  if($form eq "aho"){
    my $fw1   =  $self->getHMMCOLRange($seqid,0,24);
     $fw1   =~ s/\.//g;
     $fw1   =~ s/^[a-z\.]*//g;
    my $cdr1  =  uc($self->getHMMCOLRange($seqid,25,34));
     $cdr1  =~ s/\.//g;
    my $fw2   =  $self->getHMMCOLRange($seqid,35,48);
     $fw2   =~ s/\.//g;
    my $cdr2  =  uc($self->getHMMCOLRange($seqid,49,66));
     $cdr2  =~ s/\.//g;
    my $fw3   =  $self->getHMMCOLRange($seqid,67,94);
     $fw3   =~ s/\.//g;
    my $cdr3  =  $self->getHMMCOLRange($seqid,94,102);
     $cdr3  =~ s/^.//;
     $cdr3  =~ s/.$//;
     $cdr3  =~ s/\.//g;
     
    my $fw4   =  $self->getHMMCOLRange($seqid,102,112);
       $fw4   =~ s/[a-z\.]*$//g;
       $fw4   =~ s/\.//g;
    return(uc($fw1),uc($cdr1),uc($fw2),uc($cdr2),uc($fw3),uc($cdr3),uc($fw4));
  }elsif($form eq "kabat"){
    my $fw1   =  $self->getHMMCOLRange($seqid,0,29);
     $fw1   =~ s/\.//g;
     $fw1   =~ s/^[a-z\.]*//g;
    my $cdr1  =  uc($self->getHMMCOLRange($seqid,30,34));
     $cdr1  =~ s/\.//g;
    my $fw2   =  $self->getHMMCOLRange($seqid,35,48);
     $fw2   =~ s/\.//g;
    my $cdr2  =  uc($self->getHMMCOLRange($seqid,49,64));
     $cdr2  =~ s/\.//g;
    my $fw3   =  $self->getHMMCOLRange($seqid,65,94);
     $fw3   =~ s/\.//g;
    my $cdr3  =  $self->getHMMCOLRange($seqid,94,102);
     $cdr3  =~ s/^.//;
     $cdr3  =~ s/.$//;
     $cdr3  =~ s/\.//g;
    my $fw4   =  $self->getHMMCOLRange($seqid,102,112);
       $fw4   =~ s/[a-z\.]*$//g;
       $fw4   =~ s/\.*//g;
    return(uc($fw1),uc($cdr1),uc($fw2),uc($cdr2),uc($fw3),uc($cdr3),uc($fw4));
  }elsif($form eq "chothia"){
    my $fw1   =  $self->getHMMCOLRange($seqid,0,24);
     $fw1   =~ s/\.//g;
     $fw1   =~ s/^[a-z\.]*//g;
    my $cdr1  =  uc($self->getHMMCOLRange($seqid,25,31));
     $cdr1  =~ s/\.//g;
    my $fw2   =  $self->getHMMCOLRange($seqid,32,50);
     $fw2   =~ s/\.//g;
    my $cdr2  =  uc($self->getHMMCOLRange($seqid,51,55));
     $cdr2  =~ s/\.//g;
    my $fw3   =  $self->getHMMCOLRange($seqid,56,94);
     $fw3   =~ s/\.//g;
    my $cdr3  =  $self->getHMMCOLRange($seqid,94,102);
     $cdr3  =~ s/^.//;
     $cdr3  =~ s/.$//;
     $cdr3  =~ s/\.//g;
    my $fw4   =  $self->getHMMCOLRange($seqid,102,112);
       $fw4   =~ s/[a-z\.]*$//g;
    return(uc($fw1),uc($cdr1),uc($fw2),uc($cdr2),uc($fw3),uc($cdr3),uc($fw4));
  }
}

sub parseLightChainFrameworks {
  my($self,$seqid)=@_; 

  my $fw1   =  $self->getHMMCOLRange($seqid,128,150);
     $fw1   =~ s/\.//g;
     $fw1   =~ s/^[a-z\.]*//g;
  my $cdr1  =  uc($self->getHMMCOLRange($seqid,151,161));
     $cdr1  =~ s/\.//g;
  my $fw2   =  $self->getHMMCOLRange($seqid,162,178);
     $fw2   =~ s/\.//g;
  my $cdr2  =  uc($self->getHMMCOLRange($seqid,179,183));
     $cdr2  =~ s/\.//g;
  my $fw3   =  $self->getHMMCOLRange($seqid,184,215);
     $fw3   =~ s/\.//g;
  my $cdr3  =  $self->getHMMCOLRange($seqid,215,226);
     $cdr3  =~ s/^.//;
     $cdr3  =~ s/.$//;
     $cdr3  =~ s/\.*//g;
  my $fw4   =  $self->getHMMCOLRange($seqid,226,236);
     $fw4   =~ s/\.//g;
     $fw4   =~ s/[a-z\.]*$//g;
  return($fw1,$cdr1,$fw2,$cdr2,$fw3,$cdr3,$fw4);
}

sub scoreMotif {
  my($self,$seq,$motif)=@_;
  # a sequence "TAVVYC" and a motif array "[ATY]","[CV]","[MTH]"
  $seq=~s/[a-z]//g;
  # next score states
  my @residues=split(/ */,$seq);

  my $score=0;
  for(my $m=0;$m<scalar(@$motif);$m++){
    if(defined($residues[$m])){
      my $position_motif=$$motif[$m];
      if($residues[$m]=~m/$position_motif/){
        $score++;
      }
    }
  }
  return $score;
}

sub registerAAGermlineSHM {
  my($self)=@_;
  my $a2m=VDJFasta->new();
     $a2m->loadSeqs($self->{a2mref});
  my $a2m_count = $a2m->getSeqCount(); #
  for(my $n=0;$n<$self->{nseqs};$n++){
    my($this_header,$this_sequence)=$self->getSeq($n);

    my %closest_domain_score  = ();
    my %closest_domain_id     = ();
    my %closest_mutation_list = ();

    # IGHV IG[KL]V IGHJ IG[KL]J    
    for(my $r=0;$r<$a2m_count;$r++){
      my($a2m_header,$a2m_sequence)=$a2m->getSeq($r);
      my $domain=$a2m_header;
         $domain=~s/[0-9].*//;
         $domain=~s/^>//;
      $this_sequence =~ s/[a-z\.]//g;
      $a2m_sequence  =~ s/[a-z\.]//g;
      my ($distance,$length,$mutation_list)=$self->getHMMSeqDistance($a2m_sequence,$this_sequence);
      if($length>0){
        if(defined($closest_domain_score{$domain})){
          if($closest_domain_score{$domain}<(($length-$distance)/$length)){
            $closest_domain_score{$domain}=(($length-$distance)/$length);
            $closest_domain_id{$domain}=$r;
            $closest_mutation_list{$domain}=$mutation_list;
          }
        }else{
          $closest_domain_score{$domain}=(($length-$distance)/$length);
          $closest_domain_id{$domain}=$r;
          $closest_mutation_list{$domain}=$mutation_list;
        }        
      }
    }
    # Light chain collapse
    if(defined($closest_domain_score{"IGLV"})){
      if(defined($closest_domain_score{"IGKV"})){
        if($closest_domain_score{"IGKV"} > $closest_domain_score{"IGLV"}){
          $closest_domain_score{"IGLV"}  = $closest_domain_score{"IGKV"};
          $closest_domain_id{"IGLV"}     = $closest_domain_id{"IGKV"};
          $closest_mutation_list{"IGLV"} = $closest_mutation_list{"IGKV"};
        }
      }
    }
    if(defined($closest_domain_score{"IGLJ"})){
      if(defined($closest_domain_score{"IGKJ"})){
        if($closest_domain_score{"IGKJ"} > $closest_domain_score{"IGLJ"}){
          $closest_domain_score{"IGLJ"}  = $closest_domain_score{"IGKJ"};
          $closest_domain_id{"IGLJ"}     = $closest_domain_id{"IGKJ"};
          $closest_mutation_list{"IGLJ"} = $closest_mutation_list{"IGKJ"};
        }
      }
    }
    # Heavy chain collapse
    if(defined($closest_domain_score{"IGHV"})){
      my($a2m_header,$a2m_sequence)=$a2m->getSeq($closest_domain_id{"IGHV"});
      $a2m_header=~s/>//;
      my $pid = ((int(1000 * $closest_domain_score{"IGHV"}))/10);
      ${$self->{aaVhseg}}[$n]   = "$a2m_header $pid " 
                                . $closest_domain_id{"IGHV"} 
                                . $closest_mutation_list{"IGHV"}; # example: IGHV1-46 98.9% 41  V1M
    }
    if(defined($closest_domain_score{"IGHJ"})){
      my($a2m_header,$a2m_sequence)=$a2m->getSeq($closest_domain_id{"IGHJ"});
      $a2m_header=~s/>//;
      my $pid = ((int(1000 * $closest_domain_score{"IGHJ"}))/10);
      ${$self->{aaJhseg}}[$n]   = "$a2m_header $pid " 
                                . $closest_domain_id{"IGHJ"} 
                                . $closest_mutation_list{"IGHJ"}; # example: IGHV1-46 98.9% 41  V1M
    }
    if(defined($closest_domain_score{"IGLV"})){
      my($a2m_header,$a2m_sequence)=$a2m->getSeq($closest_domain_id{"IGLV"});
      $a2m_header=~s/>//;
      my $pid = ((int(1000 * $closest_domain_score{"IGLV"}))/10);
      ${$self->{aaVlseg}}[$n]   = "$a2m_header $pid " 
                                . $closest_domain_id{"IGLV"} 
                                . $closest_mutation_list{"IGLV"}; # example: IGHV1-46 98.9% 41  V1M
    }
    if(defined($closest_domain_score{"IGLJ"})){
      my($a2m_header,$a2m_sequence)=$a2m->getSeq($closest_domain_id{"IGLJ"});
      $a2m_header=~s/>//;
      my $pid = ((int(1000 * $closest_domain_score{"IGLJ"}))/10);
      ${$self->{aaJlseg}}[$n]   = "$a2m_header $pid " 
                                . $closest_domain_id{"IGLJ"} 
                                . $closest_mutation_list{"IGLJ"}; # example: IGHV1-46 98.9% 41  V1M
    }

    # for each
    #my @domains=("IGHV","IGHJ","IGLV","IGLJ");
    #for(my $d=0;$d<scalar(@domains);$d++){
    #  if(defined($closest_domain_score{$domains[$d]})){
    #    my($a2m_header,$a2m_sequence)=$a2m->getSeq($closest_domain_id{$domains[$d]});
    #    print $domains[$d] . "\t" . ((int(1000 * $closest_domain_score{$domains[$d]}))/10) 
    #          . "\t" . $closest_domain_id{$domains[$d]} 
    #          . "\t" . $a2m_header 
    #          . "\t" . $closest_mutation_list{$domains[$d]} . "\n";
    #  }
    #}
  }
}

sub getAAVHGermHash {
  my($self)=@_;
  my %output=();

  for(my $n=0;$n<$self->{nseqs};$n++){
    if(defined(${$self->{aaVhseg}}[$n])){
      $output{$self->getAccession($n)}=${$self->{aaVhseg}}[$n];
    }
  }
  return %output;
}

sub getAAVLGermHash {
  my($self)=@_;
  my %output=();

  for(my $n=0;$n<$self->{nseqs};$n++){
    if(defined(${$self->{aaVlseg}}[$n])){
      $output{$self->getAccession($n)}=${$self->{aaVlseg}}[$n];
    }
  }
  return %output;
}

sub getAAJHGermHash {
  my($self)=@_;
  my %output=();

  for(my $n=0;$n<$self->{nseqs};$n++){
    if(defined(${$self->{aaJhseg}}[$n])){
      $output{$self->getAccession($n)}=${$self->{aaJhseg}}[$n];
    }
  }
  return %output;
}

sub getAAJLGermHash {
  my($self)=@_;
  my %output=();

  for(my $n=0;$n<$self->{nseqs};$n++){
    if(defined(${$self->{aaJlseg}}[$n])){
      $output{$self->getAccession($n)}=${$self->{aaJlseg}}[$n];
    }
  }
  return %output;
}

sub getHydrophobicityClusters {
  my($self)=@_;

  # first get all the sequences with numbering
  my %kabat_numbered_hash_by_acc=$self->assignKabatNumbering();

  # now check the cardinal patches in all CDRs
  my %cdr_patches = (	'H1 N-term', '1,2,3,24,25,26,27,28,29,97,106',
			'H1 central','27,28,29,30,31,32,33,52,52A,52B,52C,71,73,74,77',
			'H1 C-term', '34',

			'H2 N-term', '51',
                        'H2 central','54',
                        'H2 C-term', '58',

                        'H3 N-term', '95',
                        'H3 C-term', '104',

                        'L1 N-term', '27',
                        'L1 central','31',
                        'L1 C-term', '34',

                        'L2 N-term', '51',
                        'L2 C-term', '55',

                        'L3 N-term', '91',
                        'L3 C-term', '95'
		    );
  # this is a stub
  my %output=();
  return %output;
}

sub getChargeClusters {
  my($self)=@_;
  # this is a stub
  my %output=();
  return %output;
}

sub assignKabatNumbering {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    
    my($header,$sequence)=$self->getSeq($n);
    my @chars=split(/ */,$sequence);

    # FW1
    my $fw1=$self->getHMMCOLRange($n,0,25);
       $fw1=~s/[a-z\.]*([a-z])/$1/;
    if($fw1=~m/^[a-z]\-/){
      $fw1=~s/\-//;
      $fw1=uc($fw1);
    }else{
      $fw1=~s/^[a-z]//;
    }

    # CDR1. overlap to handle insert states at boundaries: CpsnYDFwsgycnW wsgycnWGQGTLVTVSS
    my $cdr1=$self->getHMMCOLRange($n,25,35);
       $cdr1=~s/\.*//g;
       $cdr1=~s/[A-Z\-]$//;
       $cdr1=~s/^[A-Z\-]//;

    # FW2
    my $fw2=$self->getHMMCOLRange($n,35,46);
       $fw2=~s/\.*//g;

    # CDR2
    my $cdr2=$self->getHMMCOLRange($n,46,65);
       $cdr2=~s/\.*//g;
       $cdr2=~s/[A-Z\-]$//;
       $cdr2=~s/^[A-Z\-]//;

    # FW3
    my $fw3=$self->getHMMCOLRange($n,65,94);
       $fw3=~s/\.*//g;

    # CDR3
    my $cdr3=$self->getHMMCOLRange($n,94,102);
       $cdr3=~s/\.*//g;
       # overlap to handle insert states at boundaries: CpsnYDFwsgycnW wsgycnWGQGTLVTVSS
       $cdr3=~s/[A-Z\-]$//;
       $cdr3=~s/^[A-Z\-]//;

    # FW4
    my $fw4=$self->getHMMCOLRange($n,102,112);
       $fw4=~s/\.*//g;
       $fw4=~s/^[a-z]*//;  

    #-------------------------- VL ---------------------------

    # FW1 VL
    my $vl_fw1=$self->getHMMCOLRange($n,128,150);
       $vl_fw1=~s/\.*//g;
       $vl_fw1=~s/^[a-z]*//;

    # CDR1 VL
    my $vl_cdr1=$self->getHMMCOLRange($n,150,162);
       $vl_cdr1=~s/\.*//g;
       $vl_cdr1=~s/^[A-Z]//;  
       $vl_cdr1=~s/[A-Z]$//; 

    # FW2 VL
    my $vl_fw2=$self->getHMMCOLRange($n,162,176);
       $vl_fw2=~s/\.*//g;

    # CDR2 VL
    my $vl_cdr2=$self->getHMMCOLRange($n,176,184);
       $vl_cdr2=~s/\.*//g;
       $vl_cdr2=~s/^[A-Z]//; 
       $vl_cdr2=~s/[A-Z]$//;

    # FW3 VL
    my $vl_fw3=$self->getHMMCOLRange($n,184,215);
       $vl_fw3=~s/\.*//g;

    # CDR3 VL
    my $vl_cdr3=$self->getHMMCOLRange($n,215,226);
       $vl_cdr3=~s/\.*//g;
       $vl_cdr3=~s/^[A-Z]//; 
       $vl_cdr3=~s/[A-Z]$//;
       $vl_cdr3=~s/\-*//g;
       $vl_cdr3=uc($vl_cdr3);

    # FW4 VL
    my $vl_fw4=$self->getHMMCOLRange($n,226,236);
       $vl_fw4=~s/\.*//g;

   #print $vl_fw1 . " " . $vl_cdr1 . " " . $vl_fw2 . " " . $vl_cdr2 . " " . $vl_fw3 . " " . $vl_cdr3 . " " . $vl_fw4 . "\n";

    #---------------------- OK, assign numbers ------------------------

    # FW1
    @chars=split(/ */,uc($fw1));
    my $kabat_num=1;
    for(my $c=0;$c<scalar(@chars);$c++){
      $output{$self->getAccession($n)} .=  "H" .  $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
    }

    # CDR1 27 - 35 35A 35B 35C
    $kabat_num=27;
    @chars=split(/ */,uc($cdr1));
    my $c=0;
    while(($c<scalar(@chars))&&($kabat_num<36)){
      $output{$self->getAccession($n)} .=  "H" .  $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }
    $kabat_num=0;
    my @h1_num=("35A","35B","35C","35D","35E","35F");
    while($c<scalar(@chars)){
      $output{$self->getAccession($n)} .=  "H" .  $h1_num[$kabat_num] . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }

    # FW2    
    @chars=split(/ */,uc($fw2));
    $kabat_num=36;
    for(my $c=0;$c<scalar(@chars);$c++){
      $output{$self->getAccession($n)} .=  "H" .  $kabat_num . "." .  $chars[$c] . " ";    
      $kabat_num++;
    }

    # CDR2 48 - 52 52A 52B 52C 
    $kabat_num=48;
    @chars=split(/ */,uc($cdr2));
    $c=0;
    while(($c<scalar(@chars))&&($kabat_num<53)){
      $output{$self->getAccession($n)} .=  "H" .  $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }
    my $h2_length=length($cdr2);
    $kabat_num=0;
    my @h2_num=("52A","52B","52C","52D","52E","52F");
    while(($c<scalar(@chars))&&($h2_length>18)){
      $output{$self->getAccession($n)} .= "H" .  $h2_num[$kabat_num] . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
      $h2_length--;
    }
    $kabat_num=53;
    while($c<scalar(@chars)){
      $output{$self->getAccession($n)} .=  "H" .  $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }

    # FW3 66 - 82 82A 82B 82C - 83 - 92
    @chars=split(/ */,uc($fw3));
    $kabat_num=66;
    $c=0;
    while(($c<scalar(@chars))&&($kabat_num<83)){
      $output{$self->getAccession($n)} .=  "H" .  $kabat_num . "." .  $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }
    my $fw3_length = scalar(@chars);
    my @fw3_num=("82A","82B","82C","82D","82E","82F","82G","82H","82I","82J");
    $kabat_num=0;
    while($fw3_length>27){
      if(defined($fw3_num[$kabat_num])){
        $output{$self->getAccession($n)} .=  "H" .  $fw3_num[$kabat_num] . "." . $chars[$c] . " ";
      }else{
        $output{$self->getAccession($n)} .=  "H" . "82?" . "." . $chars[$c] . " ";
      }
      $kabat_num++;
      $c++;
      $fw3_length--;
    } 
    $kabat_num=83;      
    while($c<scalar(@chars)){
      $output{$self->getAccession($n)} .=  "H" .  $kabat_num . "." .  $chars[$c] . " ";
      $c++;
      $kabat_num++;
    } 

    # CDR3
    @chars=split(/ */,uc($cdr3));
    $kabat_num=93;
    my %h3_label=(	93,93,	94,94,	95,95,	96,96,	97,97,
			98,98,	99,99,	100,100,
			101,	"100A",	102,	"100B",
			103,	"100C",	104,	"100D",
                        105,    "100E", 106,    "100F",
                        107,    "100G", 108,    "100H",
                        109,    "100I", 110,    "100J",
                        111,    "100K", 112,    "100L",
                        113,    "100M", 114,    "100N",
                        115,    "100O", 116,    "100P",

                        117,    "100Q", 118,    "100R",
                        119,    "100S", 120,    "100T",
                        121,    "100U", 122,    "100V",
                        123,    "100W", 124,    "100X",
                        125,    "100Y", 126,    "100Z");
    for(my $c=0;$c<(scalar(@chars)-2);$c++){
      $output{$self->getAccession($n)} .=  "H" . $h3_label{$kabat_num} . "." .  $chars[$c] . " ";
      $kabat_num++;
    }
   $output{$self->getAccession($n)} .= "H101." . $chars[(scalar(@chars)-2)] . " ";
   $output{$self->getAccession($n)} .= "H102." . $chars[(scalar(@chars)-1)] . " ";

    # FW4
    @chars=split(/ */,uc($fw4));
    $kabat_num=103;
    for(my $c=0;$c<scalar(@chars);$c++){
      $output{$self->getAccession($n)} .= "H" .  $kabat_num . "." .  $chars[$c] . " ";    
      $kabat_num++;
    }

    #----- light chains ------

    # VL-FW1
    @chars=split(/ */,uc($vl_fw1));
    $kabat_num=1;
    for(my $c=0;$c<scalar(@chars);$c++){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
    }

    # VL-CDR1
    $kabat_num=24;
    @chars=split(/ */,uc($vl_cdr1));
    $c=0;
    while(($c<scalar(@chars))&&($kabat_num<28)){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }
    my $vl_cdr1_length=length($vl_cdr1);
    my @l1_num=("27A","27B","27C","27D","27E","27F","27G","27H","27I","27J","27K","27L");
    $kabat_num=0;
    while($vl_cdr1_length>11){
      $vl_cdr1_length--;
      #print $l1_num[$kabat_num] . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $l1_num[$kabat_num] . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }
    $kabat_num=28;
    while($c<scalar(@chars)){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }

    # VL-FW2
    @chars=split(/ */,uc($vl_fw2));
    $kabat_num=35;
    for(my $c=0;$c<scalar(@chars);$c++){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
    }

    # VL-CDR2
    @chars=split(/ */,uc($vl_cdr2));
    $kabat_num=50;
    for(my $c=0;$c<scalar(@chars);$c++){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
    }

    # VL-FW3
    @chars=split(/ */,uc($vl_fw3));
    $kabat_num=57;
    for(my $c=0;$c<scalar(@chars);$c++){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
    }

    # VL-CDR3
    $kabat_num=89;
    @chars=split(/ */,uc($vl_cdr3));
    $c=0;
    while(($c<scalar(@chars))&&($kabat_num<95)){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }
    my $vl_cdr3_length=length($vl_cdr3);
    my @l3_num=("95","95A","95B","95C","95D","95E","95F","95G","95H","95I","95J","95K","95L");
    $kabat_num=0;
    while($vl_cdr3_length>8){
      $vl_cdr3_length--;
      #print $l3_num[$kabat_num] . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $l3_num[$kabat_num] . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }
    $kabat_num=96;
    while($c<scalar(@chars)){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
      $c++;
    }

    # VL-FW4
    @chars=split(/ */,uc($vl_fw4));
    $kabat_num=98;
    for(my $c=0;$c<scalar(@chars);$c++){
      #print $kabat_num . "L." . $chars[$c] . " ";
      $output{$self->getAccession($n)} .= "L" . $kabat_num . "." . $chars[$c] . " ";
      $kabat_num++;
    }

    #print "\n";


    #print $fw1 . "\t" . $cdr1 . "\t" . $fw2 . "\t" . $cdr2 . "\t" . $fw3 . "\t" . $cdr3 . "\t" . $fw4 . "\n";
  }
  return %output;
}

sub loadVHPSSM {
  my($self)=@_;
  open(PSSM,$self->{vhpssmfile});
  my @lines=<PSSM>;
  chomp(@lines);
  close(PSSM);

  my @title_fields=split(/\t/,$lines[0]);

  for(my $x=1;$x<scalar(@lines);$x++){
    my @fields=split(/\t/,$lines[$x]);
    for(my $f=2;$f<scalar(@fields);$f++){
      ${$self->{vhpssm}}{$title_fields[$f]}{$fields[0]}=$fields[$f];
      #print "Setting " . $title_fields[$f] . "-" . $fields[0] . " to " . $fields[$f] . "\n";
    }
  }
}


sub qcVH {
  my($self)=@_;
  my $window_size=10;
  my $lowest_score_allowed=4;
  my $range_start=0;
  my $range_stop=95;

  # generate motif
  my @VH_motif=();
  my @aas=('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
  for(my $p=0;$p<$range_stop;$p++){
    my $position_pattern="[";
    for(my $a=0;$a<scalar(@aas);$a++){
      if(${$self->{vhpssm}}{$aas[$a]}{$p} > 0){
        $position_pattern .= $aas[$a];
      }
    }
    $position_pattern.="]";
    #print $position_pattern;
    push @VH_motif,$position_pattern;
  }
  #print "\n";
     
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $sequence_passes=1;
    for(my $window_start=$range_start;$window_start<($range_stop-$window_size);$window_start++){
       my @window_motif= @VH_motif[$window_start..($window_start+$window_size)];
       my $range_subseq=$self->getHMMCOLRange($n,$window_start,($window_start+$window_size));
          $range_subseq=~s/[a-z\.]//g;
          $range_subseq=uc($range_subseq);
       my $range_score = $self->scoreMotif($range_subseq,\@window_motif);
       if($range_score<$lowest_score_allowed){
         $sequence_passes=0;
         #print "score: $range_score position $window_start sequence $n\n";
       }
    }
    if($sequence_passes){
      $self->printSeq($n);
      #print "Sequence failed: " . $self->getHeader($n) . "\n";
    }
  }
}

sub qcH3TCR {
  my($self)=@_;
  my $range1_start=81;
  my $range1_stop=91;
  my $range2_start=104;
  my $range2_stop=114;
  my $cdr_start=90; 
  my $cdr_stop=104;

  my @range1_motif=("[AEKLQSGYVR]","[KLPQSET]","[EGNRS]","[EDHQ]","[AMSTE]","[ASG]","[FLMVSTK]","Y","[FLY]","C");
  my @range2_motif=("[FV]","G","[ADENPQSKT]","[EG]","[MSTI]","[QWRYK]","[LV]","[IEFLST]","[IV]","[QELTV]");

  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub qcL3TCR {
  my($self)=@_;
  my $range1_start=209;  
  my $range1_stop=219;
  my $range2_start=232; 
  my $range2_stop=242;
  my $cdr_start=218; 
  my $cdr_stop=232; 
  
  my @range1_motif=("[LHQVRSEI]","[MEILHPTWK]","[EGSTQRN]","D","[ASTIE]","[ATG]","[ISELTVP]","[YF]","[FLIY]","C");
  my @range2_motif=("[FW]","[G]","[LADKQSTHP]","G","[TI]","[KLQIRST]","[LV]","[ILQSTVF]","[IV]","[EIKLNQRST]");

  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub qcH3 {
  my($self)=@_;
  my $range1_start=85;
  my $range1_stop=94;
  my $range2_start=102;
  my $range2_stop=112;
  my $cdr_start=94;
  my $cdr_stop=102;
  my @range1_motif=("[IGRKTDQ]","[ASTVP]","[EDAS]","D","[TSA]","[AG]","[VTILF]","Y","[YF]","C");
  my @range2_motif=("W","[GS]","[QPR]","G","[TASM]","[LTM]","[VLIA]","[TAINS]","[VAI]","[SAP]","[SAP]");

  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub csabaMatch {
  my($self,$outfile)=@_;
  # "TA[CT][CT][AT][CT]TG[TC][GA][AGCT].*[AGCT]TGGGG[GCT][CA]"
  # -9 +7
  # my($header,$seq)=$self->dna2aa($s,$f);
  my %outhash=();

  #open(OUT,">$outfile");

  for(my $n=0;$n<$self->{nseqs};$n++){
    my $acc = $self->getAccession($n);
    my $forward_seq=uc($self->getSequence($n));
    my $reverse_seq=uc($self->getReverseStrand($n));
    if($forward_seq=~m/TA[CT][CT][AT][CT]TG[TC][GA][AGCT].*[AGCT]TGGGG[GCT][CA]/){
      $forward_seq=~s/^.*TA[CT][CT][AT][CT](TG[TC][GA][AGCT].*)/$1/;
      $forward_seq=~s/(.*[AGCT]TGG)GG[GCT][CA].*$/$1/;
      if( (length($forward_seq) % 3) == 0 ){
        my $fasta=VDJFasta->new();
           $fasta->addSeq(">forward",$forward_seq);
        my($new_header,$new_aa)=$fasta->dna2aa(0,0);
        $outhash{$acc}=$new_aa;
      }
    }elsif($reverse_seq=~m/TA[CT][CT][AT][CT]TG[TC][GA][AGCT].*[AGCT]TGGGG[GCT][CA]/){
      $reverse_seq=~s/^.*TA[CT][CT][AT][CT](TG[TC][GA][AGCT].*)/$1/;
      $reverse_seq=~s/(.*[AGCT]TGG)GG[GCT][CA].*$/$1/;
      if( (length($reverse_seq) % 3) == 0 ){
        my $fasta=VDJFasta->new();
           $fasta->addSeq(">forward",$reverse_seq);
        my($new_header,$new_aa)=$fasta->dna2aa(0,0);
        $outhash{$acc}=$new_aa;
      }
    }
  }
  return %outhash;
}

sub qcH3ionTorrent {
  my($self,$outfile)=@_;
  my $range1_start=89;
  my $range1_stop=94;
  my $range2_start=102;
  my $range2_stop=108;
  my $cdr_start=94;
  my $cdr_stop=102;
  my @range1_motif=("[TSA]","[AG]","[VTILF]","Y","[YF]","C");
  my @range2_motif=("W","[GS]","[QPR]","G","[TASM]","[LTM]","[VLIA]");

  my %outhash = $self->qcCDR($outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub qcL3 {
  my($self)=@_;
  my $range1_start=207;
  my $range1_stop=215;
  my $range2_start=226;
  my $range2_stop=237;
  my $cdr_start=215;
  my $cdr_stop=226;
  my @range1_motif=("[ADGEPSTV]","[DEG]","[DNT]","[AEFILRSTV]","[AGVLT]","[DILMSTV]","[YL]","[YFHL]","[C]");
  my @range2_motif=("[FI]","[G]","[AGQPT]","[G]","[T]","[TKQR]","[LV]","[DET]","[IV]","[KL]");
 
  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub qcAhoL3 {
  my($self)=@_;
  my %outhash = $self->qcL3();
  my @accs = keys %outhash;
  for(my $a=0;$a<scalar(@accs);$a++){
    $outhash{$accs[$a]}=~s/^C//;
    $outhash{$accs[$a]}=~s/.$//;
  }
  return %outhash;
}

sub qcKabatH1 {
  my($self)=@_;
  my %outhash = $self->qcH1();
  my @accs = keys %outhash;
  for(my $a=0;$a<scalar(@accs);$a++){
    $outhash{$accs[$a]}=~s/^....//;
  }
  return %outhash;
}

sub qcChothiaFW2 {
  my($self)=@_;
  my %outhash_fw2  = $self->qcFW_H2();
  my %outhash_cdr1 = $self->qcH1();
  my %outhash_cdr2 = $self->qcH2();

  my @accs = keys %outhash_fw2;
  for(my $a=0;$a<scalar(@accs);$a++){
    if(defined($outhash_cdr1{$accs[$a]})){
    $outhash_cdr1{$accs[$a]}=~s/..*(...)$/$1/;
    }else{
      $outhash_cdr1{$accs[$a]}="";
    }
    if(defined($outhash_cdr2{$accs[$a]})){
      $outhash_cdr2{$accs[$a]}=~s/^(..).*/$1/;
    }else{
      $outhash_cdr2{$accs[$a]}="";
    }
    $outhash_fw2{$accs[$a]}= $outhash_cdr1{$accs[$a]} . $outhash_fw2{$accs[$a]} . $outhash_cdr2{$accs[$a]};
  }
  
  return %outhash_fw2;
}

sub qcChothiaH2 {
  my($self)=@_;
  my %outhash_cdr2 = $self->qcH2();
  my @accs = keys %outhash_cdr2;
  for(my $a=0;$a<scalar(@accs);$a++){
    $outhash_cdr2{$accs[$a]}=~s/^..//;
    $outhash_cdr2{$accs[$a]}=~s/..$//;
  }
  return %outhash_cdr2;
}

sub qcChothiaFW3 {
 my($self)=@_;
  my %outhash_fw3  = $self->qcFW_H3();
  my %outhash_cdr2 = $self->qcH2();
  my %outhash_cdr3 = $self->qcH3();

  my @accs = keys %outhash_fw3;
  for(my $a=0;$a<scalar(@accs);$a++){
    if(defined($outhash_cdr2{$accs[$a]})){
      $outhash_cdr2{$accs[$a]}=~s/.*(..)$/$1/;
    }else{
      $outhash_cdr2{$accs[$a]}="";
    }
    if(defined($outhash_cdr3{$accs[$a]})){
      $outhash_cdr3{$accs[$a]}=~s/^.(..).*$/$1/;
      $outhash_fw3{$accs[$a]}= $outhash_cdr2{$accs[$a]} . $outhash_fw3{$accs[$a]} . $outhash_cdr3{$accs[$a]};
    }else{
      $outhash_fw3{$accs[$a]}= $outhash_cdr2{$accs[$a]} . $outhash_fw3{$accs[$a]} . "XX";
    }
  }
  return %outhash_fw3;
}

sub qcAhoCDR3 {
 my($self)=@_;
  my %outhash_cdr3 = $self->qcH3();

  my @accs = keys %outhash_cdr3;
  for(my $a=0;$a<scalar(@accs);$a++){
    $outhash_cdr3{$accs[$a]}=~s/^..//;
    $outhash_cdr3{$accs[$a]}=~s/.$//;
  }
  return %outhash_cdr3;
}


sub qcChothiaH1 {
  my($self)=@_;
  my %outhash = $self->qcH1();
  my @accs = keys %outhash;
  for(my $a=0;$a<scalar(@accs);$a++){
    $outhash{$accs[$a]}=~s/...$//;
  }
  return %outhash;
}

sub qcKabatFW1 {
  my($self)=@_;
  my %outhash_fw1 = $self->qcFW_H1();
  my %outhash_h1  = $self->qcH1();
  my @accs = keys %outhash_fw1;
  
  for(my $a=0;$a<scalar(@accs);$a++){
    if(defined($outhash_h1{$accs[$a]})){
      $outhash_h1{$accs[$a]}=~s/^(....).*/$1/;
      $outhash_fw1{$accs[$a]}.=$outhash_h1{$accs[$a]};
    }else{
      $outhash_h1{$accs[$a]}="";
      $outhash_fw1{$accs[$a]}="";
    }
  }
  return %outhash_fw1;
}

sub qcKabatH2 {
  my($self)=@_;
  my %outhash_fw3 = $self->qcFW_H3();
  my %outhash_h2  = $self->qcH2();
  my @accs = keys %outhash_fw3;
  for(my $a=0;$a<scalar(@accs);$a++){
    $outhash_fw3{$accs[$a]}  =~ s/^(.......).*/$1/;
    $outhash_h2{$accs[$a]}  .= $outhash_fw3{$accs[$a]};
  }
  return %outhash_h2;
}

sub qcKabatFW3 {
  my($self)=@_;
  my %outhash_fw3 = $self->qcFW_H3();
  my %outhash_cdr3 = $self->qcH3();

  my @accs = keys %outhash_fw3;
  for(my $a=0;$a<scalar(@accs);$a++){
    if(defined($outhash_cdr3{$accs[$a]})){
      $outhash_cdr3{$accs[$a]}=~s/^.(..).*$/$1/;
      $outhash_fw3{$accs[$a]}  =~ s/^.......//;
      $outhash_fw3{$accs[$a]} .= $outhash_cdr3{$accs[$a]};
    }else{
      $outhash_fw3{$accs[$a]}  =~ s/^.......//;
      $outhash_fw3{$accs[$a]} .= "XX";
    }
  }
  return %outhash_fw3;
}



sub qcKabatH3 {
  my($self)=@_;
  my %outhash = $self->qcH3();
  my @accs = keys %outhash;
  for(my $a=0;$a<scalar(@accs);$a++){
    $outhash{$accs[$a]}=~s/^...//;
    $outhash{$accs[$a]}=~s/.$//;
  }
  return %outhash;
}

sub qcFW_H1 {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $VH_FW1=$self->getHMMCOLRange($n,0,21);
       $VH_FW1=~s/^[a-z\.]*//;
       $VH_FW1=~s/\.*//g;
    $output{$self->getAccession($n)}=$VH_FW1;
  }
  return %output;
}

sub qcFW_H2 {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $VH_FW2=$self->getHMMCOLRange($n,35,46);
       $VH_FW2=~s/^[a-z\.]*//;
       $VH_FW2=~s/\.*//g;
    $output{$self->getAccession($n)}=$VH_FW2;
  }
  return %output;
}

sub qcFW_H3 {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $VH_FW3=$self->getHMMCOLRange($n,60,94);
       $VH_FW3=~s/^[a-z\.]*//;
       $VH_FW3=~s/\.*//g;
    $output{$self->getAccession($n)}=$VH_FW3;
  }
  return %output;
}

sub qcFW_Aho_H3 {
  my($self)=@_;
  my %outhash_fw3  = $self->qcFW_H3();
  my %outhash_cdr3 = $self->qcH3();

  my @accs = keys %outhash_fw3;
  for(my $a=0;$a<scalar(@accs);$a++){
    if(defined($outhash_cdr3{$accs[$a]})){
      $outhash_cdr3{$accs[$a]}=~s/^.(..).*$/$1/;
      $outhash_fw3{$accs[$a]}.=$outhash_cdr3{$accs[$a]};
    }else{
      $outhash_fw3{$accs[$a]}.="XX";
    }
  }
  return %outhash_fw3;
}

sub qcFW_H4 {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $VH_FW4=$self->getHMMCOLRange($n,102,112);
       $VH_FW4=~s/[a-z\.]*$//;
       $VH_FW4=~s/\.*//g;
       $VH_FW4=~s/^[a-z\.]*//;
    $output{$self->getAccession($n)}=$VH_FW4;
  }
  return %output;
}

sub qcFW_L4 {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $VL_FW4=$self->getHMMCOLRange($n,226,236);
       $VL_FW4=~s/[a-z\.]*$//;
       $VL_FW4=~s/\.*//g;
    $output{$self->getAccession($n)}=$VL_FW4;
  }
  return %output;
}

sub qcFW_L3 {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $VL_FW3=$self->getHMMCOLRange($n,184,215);
       $VL_FW3=~s/\.*//g;
    $output{$self->getAccession($n)}=$VL_FW3;
  }
  return %output;
}

sub qcFW_L2 {
  my($self)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $VL_FW2=$self->getHMMCOLRange($n,162,176);
       $VL_FW2=~s/\.*//g;
    $output{$self->getAccession($n)}=$VL_FW2;
  }
  return %output;
}

sub qcFW_L1 {
  my($self,$VLFamilies)=@_;
  my %output=();
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $accession=$self->getAccession($n);
    my $VL_FW1=$self->getHMMCOLRange($n,128,150);
       $VL_FW1=~s/^[a-z\.]*//;
       $VL_FW1=~s/\.*//g;
    if(defined($$VLFamilies{$accession})){
      if($$VLFamilies{$accession} =~ m/IGLV/){
        $VL_FW1=~s/^.//;
      }
    }
    $output{$accession}=$VL_FW1;
  }
  return %output;
}

sub qcH1 {
  my($self)=@_;

  my $range1_start=13;
  my $range1_stop=23;
  my $range2_start=35;
  my $range2_stop=45;
  my $cdr_start=26;
  my $cdr_stop=34;
  my @range1_motif=("[P]","[GST]","[ADEGQRST]","[ST]", "[MLV]",    "[QFKRST]","[ILMV]","[ST]", "[C]",    "[AEKSTV]");
  my @range2_motif=("[W]","[AFIV]","[CHKQR]",   "[HKQ]","[AKPRSTVN]","[HPST]", "[EGS]",  "[KNQ]","[AGKRS]","[FLP]");

  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub qcH2 {
  my($self)=@_;

  my $range1_start=35;
  my $range1_stop=45;
  my $range2_start=65;
  my $range2_stop=75;
  my $cdr_start=48;
  my $cdr_stop=59;
  my @range1_motif=("[W]","[AFIV]","[CHKQR]",   "[HKQ]","[AKPRSTVN]","[HPST]", "[EGS]",  "[KNQ]","[AGKRS]","[FLP]");
  my @range2_motif=("[KQR]","[AFILTV]","[AISTKV]","[FILMV]", "[STFDN]", "[AKLRV]","[DE]","[DKNT]","[APS]","[IKST]");

  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub qcL1 {
  my($self)=@_;

  my $range1_start=141;
  my $range1_stop=151;
  my $range2_start=162;
  my $range2_stop=172;
  my $cdr_start=151;
  my $cdr_stop=161;
  my @range1_motif=("[LAGST]","[AILPSTV]","[GRADEQ]","[DEKPQRST]", "[PKTQSR]", "[AIV]",    "[KRST]","[ILM]","[ST]", "[C]");
  my @range2_motif=("[W]","[FHLYV]","[FLEQ]","[HEQS]","[HKQR]","[ANPQTS]", "[RDHEG]", "[PGHKQSTE]","[NALPSTV]","[FIPV]");

  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}

sub qcL2 {
  my($self)=@_;

  my $range1_start=162;
  my $range1_stop=172;
  my $range2_start=184;
  my $range2_stop=193;
  my $cdr_start=177;
  my $cdr_stop=183;
  my @range1_motif=("[W]","[FHLYV]","[FLEQ]","[HEQS]","[HKQR]","[ANPQTS]", "[RDHEG]", "[EPGHKQST]","[NALPSTV]","[FIPV]");
  my @range2_motif=("[GDE]","[ITV]","[LPS]","[ADESV]","[RPM]","[F]","[ST]","[AGSV]","[S]");

  my %outhash = $self->qcCDR($range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
  return %outhash;
}


sub qcCDR {
  my($self,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,$range1_motif,$range2_motif)=@_;
  # confirms boundaries of CDRs are passing motif requirements
  my $qcCDRcutoff=0.5;
  my %outhash=();
  my $range1_min_score= int(($range1_stop - $range1_start) * $qcCDRcutoff);
  my $range2_min_score= int(($range2_stop - $range2_start) * $qcCDRcutoff);
  #open(OUT,">$outfile");
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $range1_subseq=$self->getHMMCOLRange($n,$range1_start,$range1_stop);
       $range1_subseq=~s/[a-z\.]//g;
       $range1_subseq=uc($range1_subseq);
    my $range2_subseq=$self->getHMMCOLRange($n,$range2_start,$range2_stop);
       $range2_subseq=~s/[a-z\.]//g;
       $range2_subseq=uc($range2_subseq);
    my $cdr_subseq=$self->getHMMCOLRange($n,$cdr_start,$cdr_stop);
       $cdr_subseq=~s/[\.\-]//g;
       $cdr_subseq=uc($cdr_subseq);
    my $range1_score = $self->scoreMotif($range1_subseq,\@$range1_motif);
    my $range2_score = $self->scoreMotif($range2_subseq,\@$range2_motif);
    #print "testing $outfile\t" . $range1_subseq . "\t" . $range2_subseq . "\t" . $cdr_subseq . "\t" . $range1_score . "\t" . $range2_score . "\n";
    if($range1_score>$range1_min_score){
      if($range2_score>$range1_min_score){
         unless($cdr_subseq=~m/ZZZ/){
           #print OUT $self->getAccession($n) . "\t" . $cdr_subseq . "\n";
           $outhash{$self->getAccession($n)}=$cdr_subseq;
         }
       }
     }
   }
  #close(OUT);
  return %outhash;
}

sub printDNAseqs {
  my($self)=@_;
  for(my $n=0;$n<$self->{nseqs};$n++){
    if($self->isDNA($n)){
      $self->printSeq($n);
    }
  }
}

sub printAAseqs {
  my($self)=@_;
  for(my $n=0;$n<$self->{nseqs};$n++){
    unless($self->isDNA($n)){
      $self->printSeq($n);
    }
  }
}

sub isFastaDNA {
  my($self)=@_;

  my $isdna=1;

  for(my $n=0;$n<$self->{nseqs};$n++){
    unless($self->isDNA($n)){
      $isdna=0;
    }
  }
  return $isdna;
}

sub isDNA {
  my($self,$c)=@_;
  my $seq=uc($self->getSequence($c));
     $seq=~s/\.*//g;
     $seq=~s/\-*//g;
     $seq=~s/X//g;
     $seq=~s/Z//g;
  my $chars=length($seq);
     $seq=~s/[ACGTN]*//g;
  my $non_nuc=length($seq);
  my $non_nuc_fraction=$non_nuc/$chars;
  if($non_nuc_fraction<0.10){
    return 1;
  }else{
    return 0;
  }
}

sub deMultiplexReads {
  my($self)=@_;

  my @midseqs   = keys %{$self->{mids}};
  my @midnames  = values %{$self->{mids}};
  my $midlength = length($midseqs[0]);		# there is an expectation that all mids are of the same length
  my %outfiles  = ();
  my %outfile_handles = ();
  my %mid_counts  = ();

  my $input_file_prefix = $self->{filename};
     $input_file_prefix =~ s/.fasta$//;
     $input_file_prefix =~ s/.fa$//;
     $input_file_prefix =~ s/.fna$//;

  # open all files for writing
  for(my $n=0;$n<scalar(@midnames);$n++){
    $outfiles{$midnames[$n]}=$input_file_prefix . "-" . $midnames[$n] . ".fa";
    open($outfile_handles{$midnames[$n]},">" . $outfiles{$midnames[$n]} );
    $mid_counts{$midnames[$n]}=0;
  } 
  open($outfile_handles{"unknown"},">$input_file_prefix-unknown.MID.fa");
  $mid_counts{"unknown"}=0;
  push @midseqs,"unknown";
  push @midnames,"unknown";

  # write all files out
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $this_mid = $self->getMid($s,$midlength);
    if(defined($outfile_handles{$this_mid})){
      print {$outfile_handles{$this_mid}} ${$self->{headers}}[$s] . "\n";
      print {$outfile_handles{$this_mid}} ${$self->{sequence}}[$s] . "\n";
      $mid_counts{$this_mid}++;
    }else{
      print "No file for |$this_mid|\n";
    }
  }
   
  #close all files for writing, and return counts
  for(my $n=0;$n<scalar(@midnames);$n++){
    close($outfile_handles{$midnames[$n]});
    my $percent = (int(1000 * ($mid_counts{$midnames[$n]} / $self->{nseqs})))/10 ;
    if($percent>0.1){
      print $midnames[$n] . "\t" . $mid_counts{$midnames[$n]} . "\t$percent%\n";
    }else{
      my $command="rm " . $outfiles{$midnames[$n]};     
      `$command`;
    }
  }
}

sub getMid {
  my($self,$seq,$midlength)=@_;
  my($header,$sequence)=$self->getSeq($seq);
  my $midrange=substr($sequence,0,$midlength);

  my %mids = %{$self->{mids}}; #$self->getMIDhash($midlength);
  if(defined(${$self->{mids}}{$midrange})){
    return ${$self->{mids}}{$midrange};
  }else{
    return "unknown";
  }
}

sub loadMIDs {
  my($self,$midfile)=@_;

  unless(-f $midfile){
    print "Error! No midfile $midfile was detected!\n";
    exit;
  }else{
    open(FILE,$midfile);
    my @lines=<FILE>;
    chomp(@lines);
    close(FILE);
    
    for(my $x=0;$x<scalar(@lines);$x++){
      my @fields=split(/\t/,$lines[$x]);
      unless(scalar(@fields) == 2){
        print "Error! Wrong number of fields on line $x in midfile $midfile: " . $lines[$x] . "\n";
        print "Format should be: midname<tab>midseq\n";
        exit;
      }else{
        my $midname=$fields[0];
        my $midseq =$fields[1];
        ${$self->{mids}}{$midseq}=$midname;
      }
    }
  }
}

sub splitSeqHeader {
  my($self,$seq)=@_;
  my @header_fields=split(/;/,${$self->{headers}}[$seq]);
  return @header_fields;
}

sub printInSilicoSort {
  my($self,$minlength,$minshm,$maxshm)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my @header_fields=$self->splitSeqHeader($s);
    # >GOIECX001C9I8F...;IGHV4-30-4/31 44 0 1;IGHD3-10 18 1 1;IGHJ4 41 0 1;CARTLSYASGSYDYW;;IGHM 88 0 1 
    if(defined($header_fields[1])){
      my @vseg_data=split(/ /,$header_fields[1]);
      if(defined($vseg_data[2])){
        if($vseg_data[1]>$minlength){
          if($vseg_data[2]>=$minshm){
            if($vseg_data[2]<=$maxshm){
              if(scalar(@vseg_data)<5){
                print ${$self->{headers}}[$s] . "\n";
                print ${$self->{seqs}}[$s] . "\n";
              }
            }
          }
        }
      }
    }
  }
}

sub getGeneticClone {
  my($self,$x)=@_;
  my $vseg = $self->getHeaderField($x,1);
  my $jseg = $self->getHeaderField($x,3);
  my @vseg_fields = split(/ /,$vseg);
  my @jseg_fields = split(/ /,$jseg);

  my $clone="";

  if(scalar(@vseg_fields)<5){ 
    if(scalar(@vseg_fields)>1){
      if(scalar(@jseg_fields)<5){
        if(scalar(@jseg_fields)>1){
          $vseg=$vseg_fields[0];
          $jseg=$jseg_fields[0];
          my $cdr3_field=5;
          if($vseg=~m/IGH/){
            $cdr3_field=4;
          }
          my $cdr3 = $self->getHeaderField($x,$cdr3_field);
          if(length($cdr3)>1){
            unless($cdr3=~m/[XZ]/i){
              # add this clone

              $clone=$vseg . "_" . $jseg . "_" . $cdr3;
            }
          }
        }
      }
    }
  }
  return $clone;
}

sub getCloneFreqsHash {
  my($self,$domain,$max_sampling_depth)=@_;

  my %counts_hash=$self->getCloneCountsHash($domain,$max_sampling_depth);
  my $depth=0;
  my @clones=keys %counts_hash;
  for(my $x=0;$x<scalar(@clones);$x++){
    $depth+=$counts_hash{$clones[$x]};
  }
  for(my $x=0;$x<scalar(@clones);$x++){
    $counts_hash{$clones[$x]} /= $depth;
  }
  return %counts_hash;
}

sub getViableCloneDepth {
  my($self)=@_;
  my $max_sampling_depth=$self->{nseqs};
  my $viable_clones=0;

  for(my $x=0;$x<$max_sampling_depth;$x++){
    my $vseg = $self->getHeaderField($x,1);
    my $jseg = $self->getHeaderField($x,3);
    my @vseg_fields = split(/ /,$vseg);
    my @jseg_fields = split(/ /,$jseg);
    if(scalar(@vseg_fields)<5){
      if(scalar(@vseg_fields)>1){
        if(scalar(@jseg_fields)<5){
          if(scalar(@jseg_fields)>1){
            $vseg=$vseg_fields[0];
            $jseg=$jseg_fields[0];
            my $cdr3_field=4;
            my $cdr3 = $self->getHeaderField($x,$cdr3_field);
            if(length($cdr3)>1){
              unless($cdr3=~m/[XZ]/i){
                $viable_clones++;
              }
            }
          }
        }
      }
    }
  }
  return $viable_clones;
}

sub getCloneCountsHash {
  my($self,$domain,$max_sampling_depth)=@_;
  my %clone_counts_hash=();

  if(!defined($domain)){
    $domain="H";
  }
  if(!defined($max_sampling_depth)){
    $max_sampling_depth=$self->{nseqs};
  }
  if($max_sampling_depth>$self->{nseqs}){
    $max_sampling_depth=$self->{nseqs};
  }
  for(my $x=0;$x<$max_sampling_depth;$x++){
    my $vseg = "";
    my $jseg = "";
    if($domain eq "H"){
      $vseg = $self->getHeaderField($x,1);
      $jseg = $self->getHeaderField($x,3);
    }elsif($domain eq "K" or $domain eq "L"){
      $vseg = $self->getHeaderField($x,12);
      $jseg = $self->getHeaderField($x,13);
    }
    my @vseg_fields = split(/ /,$vseg);
    my @jseg_fields = split(/ /,$jseg);
    #print $vseg_fields[0] . " " . $jseg_fields[0] . "\n";
    if(scalar(@vseg_fields)<5){
      if(scalar(@vseg_fields)>1){ 
        if(scalar(@jseg_fields)<5){
          if(scalar(@jseg_fields)>1){
            $vseg=$vseg_fields[0];
            $jseg=$jseg_fields[0];
            my $cdr3_field=5;
            if($vseg=~m/IGH/){
              $cdr3_field=4;
            }
            my $cdr3 = $self->getHeaderField($x,$cdr3_field);
            if(length($cdr3)>1){
              unless($cdr3=~m/[XZ]/i){
                # add this clone
                my $clone=$vseg . "_" . $jseg . "_" . $cdr3;
                if(defined($clone_counts_hash{$clone})){
                  $clone_counts_hash{$clone}++;
                }else{
                  $clone_counts_hash{$clone}=1;
                }
              }
            }
          }
        }
      }
    }
  }
  return %clone_counts_hash;
}


sub getCloneIdsHash {
  my($self,$domain,$max_sampling_depth)=@_;
  my %clone_ids_hash=();

  if(!defined($max_sampling_depth)){
    $max_sampling_depth=$self->{nseqs};
  }
  if($max_sampling_depth>$self->{nseqs}){
    $max_sampling_depth=$self->{nseqs};
  }
  for(my $x=0;$x<$max_sampling_depth;$x++){
    my $vseg = "";
    my $jseg = "";
    if($domain eq "H"){
      $vseg = $self->getHeaderField($x,1);
      $jseg = $self->getHeaderField($x,3);
    }elsif($domain eq "K" or $domain eq "L"){
      $vseg = $self->getHeaderField($x,12);
      $jseg = $self->getHeaderField($x,13);
    }
    my @vseg_fields = split(/ /,$vseg);
    my @jseg_fields = split(/ /,$jseg);
    #print $vseg_fields[0] . " " . $jseg_fields[0] . "\n";
    if(scalar(@vseg_fields)<5){
      if(scalar(@vseg_fields)>1){
        if(scalar(@jseg_fields)<5){
          if(scalar(@jseg_fields)>1){
            $vseg=$vseg_fields[0];
            $jseg=$jseg_fields[0];
            my $cdr3_field=5;
            if($vseg=~m/IGH/){
              $cdr3_field=4;
            }
            my $cdr3 = $self->getHeaderField($x,$cdr3_field);
            if(length($cdr3)>1){
              unless($cdr3=~m/[XZ]/i){
                # add this clone
                my $clone=$vseg . "_" . $jseg . "_" . $cdr3;
                if(defined($clone_ids_hash{$clone})){
                  $clone_ids_hash{$clone} .=",$x";
                }else{
                  $clone_ids_hash{$clone}="$x";
                }
              }
            }
          }
        }
      }
    }
  }
  return %clone_ids_hash;
}

sub writeFullClones {
  my($self,$domain,$output)=@_;
  open(FILE,">$output");

  for(my $s=0;$s<$self->{nseqs};$s++){
    my @header_fields=$self->splitSeqHeader($s);
    my $v_header="";
    my $j_header="";
    my $cdr3     ="";
    if($domain eq "H" and (defined($header_fields[1]) or defined($header_fields[3]))){
      $v_header=$header_fields[1];
      $j_header=$header_fields[3];
      $cdr3    =$header_fields[4];
    }elsif($domain eq "K" or $domain eq "L" and (defined($header_fields[12]) or defined($header_fields[13]))){
      $v_header=$header_fields[12];
      $j_header=$header_fields[13];
      $cdr3    =$header_fields[5];
    }
    if($v_header ne "" or $j_header ne ""){
        my @vseg_data=split(/ /,$v_header);
        my @jseg_data=split(/ /,$j_header);
        if(scalar(@vseg_data)==3){
          if($vseg_data[0]=~m/IG$domain/){
            if(scalar(@jseg_data)==3){
              if($jseg_data[0]=~m/IG$domain/){
                if( length($cdr3)>0 ){
                  my $cdr_data = $cdr3;
                  if($cdr_data =~ m/^[A-WY][A-WY]*$/){
                    print FILE ${$self->{headers}}[$s] . "\n";
                    print FILE ${$self->{sequence}}[$s] . "\n";
                  }
                }
              }
            }
          }
        }
     }
  }
  close(FILE);
}

sub printUniqueClones {
  my($self)=@_;
  # take the subset of unique clones
  my %clones_already_seen=();
  for(my $s=0;$s<$self->{nseqs};$s++){
    my @header_fields=$self->splitSeqHeader($s);
    if(defined($header_fields[4])){
      my @vseg_data=split(/ /,$header_fields[1]);
      my @jseg_data=split(/ /,$header_fields[3]);
      my $h3=$header_fields[4];
      my $clone=$vseg_data[0] . $h3 . $jseg_data[0];
      unless(defined($clones_already_seen{$clone})){
        print ${$self->{headers}}[$s] . "\n";
        print ${$self->{seqs}}[$s] . "\n";
      }
      $self->addCloneVariants(\%clones_already_seen,$vseg_data[0],$h3,$jseg_data[0]);
    }
  }
}

sub addCloneVariants {
  my($self,$clones_hash,$vseg,$h3,$jseg)=@_;
  # registers all distance 1 variants in a hash of CDR3's that have been previously seen
  my($hash,$germline,$sequence)=@_;
  my @aa=("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
  my @chars=split(/ */,$h3);
  #print $sequence . "\n";
  for(my $position=1;$position<scalar(@chars);$position++){
    my $left_string=substr($h3,0,$position);
    my $right_string=substr($h3,($position+1),(scalar(@chars)-$position));
    for(my $a=0;$a<scalar(@aa);$a++){
      # altering position $position to residue $a
      $$clones_hash{$vseg . $left_string . $aa[$a]  . $right_string . $jseg}=1;
    }
  }
}

sub singleLinkageByCDRLength {
  my($self,$header_field,$max_distance,$bool_print_cluster,$bool_print_network,$outFile)=@_;
  
  my %cdr_lengths=();

  for(my $x=0;$x<$self->{nseqs};$x++){
    my $seq_length = length($self->getHeaderField($x,$header_field));
    if(defined($cdr_lengths{$seq_length})){
      $cdr_lengths{$seq_length}++;
    }else{
      $cdr_lengths{$seq_length}=1;
    }
  }

  my @cdr_lengths = keys %cdr_lengths;

  for(my $j=0;$j<scalar(@cdr_lengths);$j++){
      $self->singleLinkageCDRLength($header_field,$max_distance,$cdr_lengths[$j],$bool_print_cluster,$bool_print_network,$outFile);
  }
}

sub getGeneticGroup {
  my($self,$seqid)=@_;
      
  my $this_vgene  = $self->getHeaderField($seqid,1);
  my $this_jgene  = $self->getHeaderField($seqid,3);
  my $this_h3     = $self->getHeaderField($seqid,4);
  my $this_l3     = $self->getHeaderField($seqid,5);

  my $final_vgene = "";
  my $final_jgene = "";
  my $final_cdr3  = "";

  # eliminate allele calls
  $this_vgene=~s/\*.*//;
  $this_jgene=~s/\*.*//;
 
  # Vgene
  unless($this_vgene =~ m/ IG/){
    unless($this_vgene =~ m/ mIG/){  # mTRB
      $this_vgene=~s/ .*//;
      if(($this_vgene=~m/IG/)||($this_vgene=~m/TRB/)){
        $final_vgene=$this_vgene;       
      }
    }
  }
  # Jgene
  unless($this_jgene =~ m/ IG/){
    unless($this_jgene =~ m/ mIG/){  # mTRB
      $this_jgene=~s/ .*//;
      if(($this_jgene=~m/IG/)||($this_jgene=~m/TRB/)){
        $final_jgene=$this_jgene;
      }
    }
  }
  # CDR=H3
  if(length($this_h3)>0){
    $final_cdr3=length($this_h3);
  }
  if(length($this_l3)>0){
    $final_cdr3=length($this_l3);
  }
  unless( (length($final_vgene)>0) && (length($final_jgene)>0) && (length($final_cdr3)>0) ){
    $final_vgene="V";
    $final_jgene="J";
    $final_cdr3="0";
  }
  # genetic group
  my $genetic_group = $final_vgene . "-" . $final_jgene . "-" . $final_cdr3;
     $genetic_group =~ s/\//~/g;

  return $genetic_group;
}

sub singleLinkageByGermline {
  my($self,$header_field,$max_distance,$bool_print_cluster,$bool_print_network)=@_;

  # first get germline list
  my %vgenes=();
  my %jgenes=();
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_vgene  = $self->getHeaderField($x,1);
    my $this_jgene  = $self->getHeaderField($x,3);
    unless($this_vgene =~ m/ IG/){
      unless($this_vgene =~ m/ mIG/){  # mTRB
        $this_vgene=~s/ .*//;
        if(($this_vgene=~m/IG/)||($this_vgene=~m/TR/)){
          $vgenes{$this_vgene}=1;
        }
      }
    }
    unless($this_jgene =~ m/ IG/){
      unless($this_jgene =~ m/ mIG/){
        $this_jgene=~s/ .*//;
        if(($this_jgene=~m/IG/)||($this_jgene=~m/TR/)){
          $jgenes{$this_jgene}=1;
        }
      }
    }
  }
  my @vgene_list=sort {$vgenes{$b} cmp $vgenes{$a}} keys %vgenes;
  my @jgene_list=sort {$jgenes{$b} cmp $jgenes{$a}} keys %jgenes;

  # now go through all combinations
  for(my $v=0;$v<scalar(@vgene_list);$v++){
    for(my $j=0;$j<scalar(@jgene_list);$j++){
      $self->singleLinkageVJ($header_field,$max_distance,$vgene_list[$v],$jgene_list[$j],$bool_print_cluster,$bool_print_network);
    }
  }
}

sub reportClosestParatopeAcrossFiles {
  my($self,$file2,$scoring_mode)=@_;

  # scoring mode: paratope_pid or paratope_blosum62

  # prepare file 2, for the comparison
  my $fasta2=VDJFasta->new();
     $fasta2->loadSeqs($file2);
  my $fasta2_nseqs = $fasta2->getSeqCount();
  my @f2h1_array=();
  my @f2h2_array=();
  my @f2h3_array=();
  for(my $s2=0;$s2<$fasta2_nseqs;$s2++){
    my $f2h1 = $self->getHeaderField($s2,8);
    my $f2h2 = $self->getHeaderField($s2,9);
    my $f2h3 = $self->getHeaderField($s2,4);
    if($f2h1 ne ""){
      if($f2h2 ne ""){
        if($f2h3 ne ""){
          push @f2h1_array,$f2h1;
          push @f2h2_array,$f2h2;
          push @f2h3_array,$f2h3;
        }
      }
    }
  } 

  for(my $s1=0;$s1<$self->{nseqs};$s1++){
    # check each CDR
    my $h1 = $self->getHeaderField($s1,8);
    my $h2 = $self->getHeaderField($s1,9);
    my $h3 = $self->getHeaderField($s1,4);

    # check that all CDRs are defined
    if($h1 ne ""){
      if($h2 ne ""){
        if($h3 ne ""){
          # compare this to file2 seqs
          my $best_score=0;
          my $best_sequence=0;
          for(my $s2=0;$s2<scalar(@f2h3_array);$s2++){
            my $this_score=0;
            my $this_sequence_length=0;
            $this_sequence_length += length($h1);
            $this_sequence_length += length($h2);
            $this_sequence_length += length($h3);

            if(length($h1) == length($f2h1_array[$s2])){
              if(length($h2) == length($f2h2_array[$s2])){
                if(length($h3) == length($f2h3_array[$s2])){
                #if( (length($h3) < (length($f2h3_array[$s2])+2)) && (length($h3) > (length($f2h3_array[$s2])-2))){
                  $this_score           += $self->getSeqDistance($h1,$f2h1_array[$s2]);
                  $this_score           += $self->getSeqDistance($h2,$f2h2_array[$s2]);
                  $this_score           += $self->getSeqDistance($h3,$f2h3_array[$s2]);
                  #$this_sequence_length	+= $self->max(length($h3),$f2h3_array[$s3]);
                  # for H3, allow some length variation at edges of CAR|SDFSDFSD|DYW
                  $this_score = 1 - ($this_score / $this_sequence_length);
                  #print "\t$this_score\n";
                  if($this_score>$best_score){
                    $best_score=$this_score;
                    $best_sequence=$s2;
                  }
                }
              }
            }
          }
          print $s1 . "\t" . $best_score . "\n";
        }
      }
    }
  }
}

sub reportClosestVJClonesAcrossFiles {
  my($self,$file2,$header_field)=@_;
  #$self->reportClosestClones($header_field);
  # first get germline list
  my %vgenes=();
  my %jgenes=();
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_vgene  = $self->getHeaderField($x,1);
    my $this_jgene  = $self->getHeaderField($x,3);
    unless($this_vgene =~ m/ IG/){
      $this_vgene=~s/ .*//;
      $vgenes{$this_vgene}=1;
    }
    unless($this_jgene =~ m/ IG/){
      unless($this_jgene =~ m/ mIG/){
        $this_jgene=~s/ .*//;
        $jgenes{$this_jgene}=1;
      }
    }
  }
  my @vgene_list=sort {$vgenes{$b} cmp $vgenes{$a}} keys %vgenes;
  my @jgene_list=sort {$jgenes{$b} cmp $jgenes{$a}} keys %jgenes;

  my $fasta2=VDJFasta->new();
     $fasta2->loadSeqs($file2);
      
  # now go through all combinations, generating vj files for file2
  my %vjfiles=();
  for(my $v=0;$v<scalar(@vgene_list);$v++){
    for(my $j=0;$j<scalar(@jgene_list);$j++){
      my($subfile,$reads_found)=$fasta2->writeVJsubset($vgene_list[$v],$jgene_list[$j]);
      my $vj = $vgene_list[$v] . $jgene_list[$j];
      #print "$vj Got $subfile and $reads_found\n";
      $vjfiles{$vj}=$subfile;
    }
  }

  # now go through sequences in $self->{nseqs}, and compare to appropriate secondary file
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_vgene  = $self->getHeaderField($x,1);
    my $this_jgene  = $self->getHeaderField($x,3);
    my $this_cdr    = $self->getHeaderField($x,$header_field);
    unless($this_vgene =~ m/ IG/){
      unless($this_jgene =~ m/ IG/){
        unless($this_jgene =~ m/ mIG/){
          $this_vgene=~s/ .*//;
          $this_jgene=~s/ .*//;
          
          #now run on that file
          my $runfile=$this_vgene . $this_jgene;
          my $file2=VDJFasta->new();
             $file2->loadSeqs($vjfiles{$runfile});

             $file2->reportClosestClone($header_field,$this_cdr);
        }
      }
    }
  }

  #ok, now remove files
}


sub reportClosestVJClones {
  my($self,$header_field)=@_;
  #$self->reportClosestClones($header_field);
  # first get germline list
  my %vgenes=();
  my %jgenes=();
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_vgene  = $self->getHeaderField($x,1);
    my $this_jgene  = $self->getHeaderField($x,3);
    unless($this_vgene =~ m/ IG/){
      $this_vgene=~s/ .*//;
      $vgenes{$this_vgene}=1;
    }
    unless($this_jgene =~ m/ IG/){
      unless($this_jgene =~ m/ mIG/){
        $this_jgene=~s/ .*//;
        $jgenes{$this_jgene}=1;
      }
    }
  }
  my @vgene_list=sort {$vgenes{$b} cmp $vgenes{$a}} keys %vgenes;
  my @jgene_list=sort {$jgenes{$b} cmp $jgenes{$a}} keys %jgenes;

  # now go through all combinations
  for(my $v=0;$v<scalar(@vgene_list);$v++){
    for(my $j=0;$j<scalar(@jgene_list);$j++){
      my($subfile,$reads_found)=$self->writeVJsubset($vgene_list[$v],$jgene_list[$j]);
      if($reads_found>0){
        my $vjsub=VDJFasta->new();
           $vjsub->loadSeqs($subfile);
           $vjsub->reportClosestClones($header_field);
      }
      `rm $subfile`;
    }
  }
}

sub singleLinkageVJ {
  my($self,$header_field,$max_distance,$vgene,$jgene,$bool_print_cluster,$bool_print_network)=@_;
  my($subfile,$reads_found)=$self->writeVJsubset($vgene,$jgene);
  if($reads_found>0){
    my $vjsub=VDJFasta->new();
       $vjsub->loadSeqs($subfile);
       $vjsub->singleLinkageAllClones($header_field,$max_distance,$vgene . "_" . $jgene,$bool_print_cluster,$bool_print_network);
  }
  `rm $subfile`;
}

sub singleLinkageH3 {
  my($self,$h3,$max_distance,$vgene_pattern)=@_;
  my $h3_length = length($h3);
    
}

sub writeIDsubset {
  my($self,$id_list)=@_;
  my $outfile=int(rand(10000000000));
     $outfile.=".fa";
  open(OUT,">$outfile");
  for(my $i=0;$i<scalar(@$id_list);$i++){
    print OUT ${$self->{headers}}[$$id_list[$i]]  . "\n";
    print OUT ${$self->{sequence}}[$$id_list[$i]] . "\n";
  }
  close(OUT);
  return $outfile;
}

sub writeH3subset {
  my($self,$vgene,$h3)=@_;
  my $outfile=int(rand(10000000000));
     $outfile.=".fa";
  open(OUT,">$outfile");
  my $reads_matching=0;
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_vgene  = $self->getHeaderField($x,1);
    my $this_h3     = $self->getHeaderField($x,4);
    if($this_vgene=~m/$vgene/){
      if(length($h3) eq length($this_h3)){
        print OUT ${$self->{headers}}[$x]  . "\n";
        print OUT ${$self->{sequence}}[$x] . "\n";
        $reads_matching+=1;
      }
    }
  }
  close(OUT);
  return($outfile,$reads_matching);
}

sub writeVJsubset {
  my($self,$vgene,$jgene)=@_;
  my $outfile=int(rand(10000000000));
     $outfile.=".fa";
  open(OUT,">$outfile");
  my $reads_matching=0;
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_vgene  = $self->getHeaderField($x,1);
    my $this_jgene  = $self->getHeaderField($x,3);
    if($this_vgene=~m/^$vgene [0-9][0-9]* [0-9][0-9]*$/){
      if($this_jgene=~m/^$jgene [0-9][0-9]* [0-9][0-9]*$/){
        print OUT ${$self->{headers}}[$x]  . "\n";
        print OUT ${$self->{sequence}}[$x] . "\n";
        $reads_matching+=1;
      }
    }
  }
  close(OUT);
  return($outfile,$reads_matching);
}

sub singleLinkageCDRLength {
  my($self,$header_field,$max_distance,$length,$bool_print_cluster,$bool_print_network,$outFile)=@_;
  my($subfile,$reads_found)=$self->writeCDRlengthSubset($header_field,$length);
  if($reads_found>0){
    my $vjsub=VDJFasta->new();
       $vjsub->loadSeqs($subfile);
    my $outdata = $vjsub->singleLinkageAllClones($header_field,$max_distance,$length,$bool_print_cluster,$bool_print_network);
    open(FILE,">>$outFile");
    print FILE $outdata;
    close(FILE);
  }
  `rm $subfile`;
}



sub writeCDRlengthSubset {
  my($self,$header_field,$length)=@_;
  my $outfile=int(rand(10000000000));
     $outfile.=".len.$header_field.$length.fa";
  open(OUT,">$outfile");
  my $reads_matching=0;
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_seq  = $self->getHeaderField($x,$header_field);
    if(length($this_seq) == $length){
      print OUT ${$self->{headers}}[$x]  . "\n";
      print OUT ${$self->{sequence}}[$x] . "\n";
      $reads_matching+=1;
    }
  }
  close(OUT);
  return($outfile,$reads_matching);
}

sub reportClosestClones {
  my($self,$header_field)=@_;
  

  for(my $seed_id=0;$seed_id<$self->{nseqs};$seed_id++){
     my $closest_distance=100;
     for(my $search_id=0;$search_id<$self->{nseqs};$search_id++){
       unless($seed_id eq $search_id){  
         my $this_distance = $self->getDistance($seed_id,$search_id,$header_field);
        
         if( $this_distance < $closest_distance){
           #print "New distance: $this_distance vs $closest_distance\n";
           $closest_distance=$this_distance;
         }
       }
     }
     print $seed_id . "\t" . $closest_distance . "\n"; 
  }
}

sub writeCluster {
  my($self,$cluster_line)=@_;

  # parse cluster line
  my ($counts,$name,$accessions)=split(/\t/,$cluster_line);
  my @accs=split(/ /,$accessions);

  # create output file
  my $filename = $self->{filename};
     $filename =~ s/.fa$//;
     $filename =~ s/.fasta$//;
  my $output_file=$filename . "_" . $name . ".fa"; 
     $output_file=~s/\//_/g;
     $output_file=~s/\*/~/g;

  # print all sequences to this file
  my %acc2id = ();
  $self->getAccession2seqidHash(\%acc2id);

  open(OUT,">$output_file");
  for(my $a=0;$a<scalar(@accs);$a++){
    my $id = $acc2id{$accs[$a]};
    print OUT ${$self->{headers}}[$id] . "\n";
    print OUT ${$self->{sequence}}[$id] . "\n";
  }     
  close(OUT);
  return $output_file;
}

sub reportClosestClone {
  my($self,$header_field,$query)=@_;

  #print " I received $query and $header_field. i have " . $self->{nseqs} . "\n";
  my $closest_distance=100;
  for(my $search_id=0;$search_id<$self->{nseqs};$search_id++){
    my $this_distance = $self->getSeqDistance($query,$self->getHeaderField($search_id,$header_field));
    #print "Got distance $this_distance for $query and " . $self->getHeaderField($search_id,$header_field) . "\n";
    if( $this_distance < $closest_distance){
      $closest_distance=$this_distance;
    }
  }
  print $closest_distance . "\n";
}

sub singleLinkageAllClones {
  my($self,$header_field,$max_distance,$optional_clone_description,$bool_print_cluster,$bool_print_network)=@_;
  my %already_seen=();
  my $output="";

  # speed optimization: if more than 750 sequences, render unique and rank by frequency
  # store redundant_ids to add back into the final output so as to cause this operation to be
  # transparent
  my %cdr3_to_uniqid=();
  my %cdr3_counts=();
  my %cdr3_to_acc=();
  my %uniqacc_to_redundantaccs=();
  my @final_unique_id_list=();

  if($self->{nseqs}>0){
    for(my $seed_id=0;$seed_id<$self->{nseqs};$seed_id++){
      my $clustering_seq=$self->getHeaderField($seed_id,$header_field);  
      if($clustering_seq ne ""){
        my $this_acc = $self->getAccession($seed_id);
        if(defined($cdr3_to_acc{$clustering_seq})){
          $cdr3_counts{$clustering_seq}++;
          $uniqacc_to_redundantaccs{$cdr3_to_acc{$clustering_seq}}.= $this_acc . " ";
        }else{
          # new representative
          $cdr3_to_uniqid{$clustering_seq} = $seed_id;
          $cdr3_counts{$clustering_seq}=1;
          $cdr3_to_acc{$clustering_seq}=$this_acc;
          $uniqacc_to_redundantaccs{$this_acc}="";
        }
      }
    }
    # ok, now write a unique file containing ranked representatives by frequency
    my @ranked_list=sort {$cdr3_counts{$b} <=> $cdr3_counts{$a}} keys %cdr3_counts;
    for(my $i=0;$i<scalar(@ranked_list);$i++){
      push @final_unique_id_list,$cdr3_to_uniqid{$ranked_list[$i]};
      #print "Items in the ranked list:\t" . $ranked_list[$i] . "\t" . $cdr3_counts{$ranked_list[$i]} . "\tin a file with\t" . $self->{nseqs} . "\n";
    }  
  }
  my $nonredundant_file = $self->writeIDsubset(\@final_unique_id_list);
  my $uniqfasta=VDJFasta->new();
     $uniqfasta->loadSeqs($nonredundant_file);
  my $unique_count=$uniqfasta->getSeqCount();

  for(my $seed_id=0;$seed_id<$unique_count;$seed_id++){   # just these seqid 
    unless(defined($already_seen{$seed_id})){
      my @cluster_member_id_list=();
      my $seed_seq=$uniqfasta->getHeaderField($seed_id,$header_field);  # uniqfasta
      if(length($seed_seq)>1){
        $uniqfasta->getThisCloneCluster($seed_seq,$header_field,$max_distance,\@cluster_member_id_list,$bool_print_network); # uniqfasta
        my @unique_cluster_ids=sort { $a <=> $b } (keys %{{ map{$_=>1}@cluster_member_id_list}});

        # print the clone cluster: |253 IGHV_IGHJ_CARRYFDLW |
        if($bool_print_cluster){
          my $line = scalar(@unique_cluster_ids) . "\t$optional_clone_description" . "_$seed_seq\t";
          #print $line;
          $output .= $line;
        }
        my $expanded_accession_list="";
        for(my $n=0;$n<scalar(@unique_cluster_ids);$n++){
          $already_seen{$unique_cluster_ids[$n]}=1;
          my $seed_acc=$uniqfasta->getAccession($unique_cluster_ids[$n]); # uniqfasta
          if($bool_print_cluster){
            $output .= "$seed_acc ";
            $expanded_accession_list .= " $seed_acc " . $uniqacc_to_redundantaccs{$seed_acc};
            #print "$seed_acc ";
            #print $uniqacc_to_redundantaccs{$seed_acc} . " ";
          }
        }
        if($bool_print_cluster){
          $output .= "\n";
        }
        $expanded_accession_list=~s/  */ /g;
        $expanded_accession_list=~s/^ *//;
        $expanded_accession_list=~s/ *$//;
        #  my $line = scalar(@unique_cluster_ids) . "\t$optional_clone_description" . "_$seed_seq\t";
        #  print $line;
        my @expanded_accessions=(split(/ /,$expanded_accession_list));
        print scalar(@expanded_accessions) . "\t$optional_clone_description" . "_$seed_seq\t" . $expanded_accession_list . "\n";
      }
    }
  }
  `rm $nonredundant_file`;
  return $output;
}

sub blosumAverageSimilarity {
  my($self,$seq1,$seq2,$blosum)=@_;

  my @chars1=split(/ */,$seq1);
  my @chars2=split(/ */,$seq2);

  if( scalar(@chars1) ne scalar(@chars2) ){
    return -10;
  }
  if( scalar(@chars1) < 1 ){
    return -10;
  }
  if( scalar(@chars2) < 1 ){
    return -10;
  }

  my $total_blosum=0;
  my $total_positions=0;

  for(my $c=0;$c<scalar(@chars1);$c++){
    $total_blosum+=$$blosum{$chars1[$c]}{$chars2[$c]};
    $total_positions++;
  }
  my $average_blosum=$total_blosum/$total_positions;
  return $average_blosum;
}

sub blosum2hash {
  my($self,$file)=@_;
  my %blosum=();

  unless(defined($file)){
    $file = $self->{blosum_file};
  }

  my @lines=$self->getFileLines($file);
  my @headers=();

  for(my $x=0;$x<scalar(@lines);$x++){
    if($lines[$x]=~m/^#/){

    }elsif($lines[$x]=~m/^ /){
      @headers=split(/  */,$lines[$x]);
    }else{
      my @fields=split(/  */,$lines[$x]);
      for(my $f=1;$f<scalar(@fields);$f++){
        $blosum{$headers[$f]}{$fields[0]}=$fields[$f];
      }
    }
  }
  return %blosum;
}

sub getFileLines {
  my($self,$file)=@_;
  open(FILE,$file);
  my @lines=<FILE>;
  close(FILE);
  chomp(@lines);
  return @lines;
}

sub getThisStereotypyCluster {
  # a clever little devil of a subroutine
  my($self,$seed_h1,$seed_h2,$seed_h3,$max_distance,$cluster_member_id_list,$bool_print_network)=@_;

  # blosum62 distance between clones, splitting H3, length restricting H1 and H2
  # first get a list of already seen unique sequences (H1-H2-H3's)
  #my %already_seen_unique_seqs=();
  #for(my $id=0;$id<scalar(@$cluster_member_id_list);$id++){    
  #  my $clone = $self->getHeaderField($$cluster_member_id_list[$id],4) . "---" .
  #              $self->getHeaderField($$cluster_member_id_list[$id],4) . "---" .
  #              $self->getHeaderField($$cluster_member_id_list[$id],4);
  #  $already_seen_unique_seqs{$seed_seq}=1;
  #}
}

sub getThisCloneCluster {
  # inputs:  takes a seed sequence (H3, L3 etc), a header field (for H3, L3 etc), a max distance, and "cluster member" sequence id list
  # outputs: returns the list of all variants (passed by reference during recursion)
  my($self,$seed_seq,$header_field,$max_distance,$cluster_member_id_list,$bool_print_network)=@_;
 
  # first get a list of already seen unique sequences (H3s or whatever)
  my %already_seen_unique_seqs=();
  for(my $id=0;$id<scalar(@$cluster_member_id_list);$id++){
    $already_seen_unique_seqs{$self->getHeaderField($$cluster_member_id_list[$id],$header_field)}=1;
  }
  $already_seen_unique_seqs{$seed_seq}=1;

  # now get the observed counts of each clone
  my %clone_counts=();
  my %clone_vgene=();
  my %clone_jgene=();
  if($bool_print_network){
    for(my $x=0;$x<$self->{nseqs};$x++){
      my $this_seq  = $self->getHeaderField($x,$header_field);
      #my $this_iso  = $self->getHeaderField($x,$header_field);

      if(defined($clone_counts{$this_seq})){
        $clone_counts{$this_seq}++;
        # lowest common denom for Vgene and Jgene?   
      }else{
        $clone_counts{$this_seq} = 1;
        $clone_vgene{$this_seq}  = $self->getHeaderField($x,1);
        $clone_jgene{$this_seq}  = $self->getHeaderField($x,3);
        #$clone_isotype{$this_seq}= $self->getHeaderField($x,12);  

        $clone_vgene{$this_seq} =~s/ .*//;
        $clone_jgene{$this_seq} =~s/ .*//;
      }
    }
  }

  # now identify clustering partners
  my @new_seqs=();
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $this_seq  = $self->getHeaderField($x,$header_field);
    my $this_distance = $self->getSeqDistance($seed_seq,$this_seq);
    if($this_distance<=$max_distance){
      if($bool_print_network){
        #my $seqname_1 = $clone_vgene{$seed_seq} . " " . $clone_jgene{$seed_seq} . " " . $seed_seq . " " . $clone_counts{$seed_seq};
        #my $seqname_2 = $clone_vgene{$this_seq} . " " . $clone_jgene{$this_seq} . " " . $this_seq . " " . $clone_counts{$this_seq};
        print "$seed_seq $this_seq $this_distance\n";
      }
      push @$cluster_member_id_list,$x;
      push @new_seqs,$this_seq;
    }
  }

  # recurse for each new CDR3 found
  my @new_uniq_seqs=();
  for(my $id=0;$id<scalar(@new_seqs);$id++){
    unless(defined($already_seen_unique_seqs{$new_seqs[$id]})){
      push @new_uniq_seqs,$new_seqs[$id];
      $already_seen_unique_seqs{$new_seqs[$id]}=1;
    }
  }
  for(my $id=0;$id<scalar(@new_uniq_seqs);$id++){
    $self->getThisCloneCluster($new_uniq_seqs[$id],$header_field,$max_distance,\@$cluster_member_id_list,$bool_print_network);
  }
}

sub getSingleClusterRepresentative {
  my($self,$seq_ids,$acc2seqid)=@_;
  if(scalar(@$seq_ids)<1){
    print "";
  }else{
    my($header,$seq)=$self->getSeq($$acc2seqid{$$seq_ids[0]});
    print "$header\n";
    print "$seq\n";
  }
}

sub getPID {
  my($self,$id1,$id2)=@_;
  my $pid="";
  my $ids=0;
  my $align_length=0;

  my ($header1,$seq1)=$self->getSeq($id1);
  my ($header2,$seq2)=$self->getSeq($id2);
  my @chars1=split(/ */,$seq1);
  my @chars2=split(/ */,$seq2);
 
  for(my $c=0;$c<scalar(@chars1);$c++){
    if( ($chars1[$c] =~ m/[ACGTacgt]/) and ($chars2[$c] =~ m/[ACGTacgt]/)){
      $align_length++;
      if($chars1[$c] eq $chars2[$c]){
        $ids++;
      }
    }
  }
  if($align_length>0){
    $pid = $ids / $align_length;
  }
  return($pid,$align_length);
}

sub getDNAConsensus {
  my($self)=@_;

  my @matrix=();
  my $consensus="";

  for(my $x=0;$x<$self->{nseqs};$x++){
    my($header,$seq)=$self->getSeq($x);
    my @chars = split(/ */,$seq); 
    for(my $c=0;$c<scalar(@chars);$c++){
      if($chars[$c] =~ m/[ACGTacgt]/){
        if(defined($matrix[$c]{$chars[$c]})){
          $matrix[$c]{$chars[$c]}++;
        }else{
          $matrix[$c]{$chars[$c]}=1;
        }
      }
    }
  }

  # ok, now report the consensus sequence
  my($header,$seq)=$self->getSeq(0);
  for(my $c=0;$c<length($seq);$c++){
    my %options=("A",0,"C",0,"G",0,"T",0);
    if(defined($matrix[$c]{A})){
      $options{A}=$matrix[$c]{A};
    }
    if(defined($matrix[$c]{C})){
      $options{C}=$matrix[$c]{C};
    }
    if(defined($matrix[$c]{G})){
      $options{G}=$matrix[$c]{G};
    }
    if(defined($matrix[$c]{T})){
      $options{T}=$matrix[$c]{T};
    }
    my @sorted = sort { $options{$b} cmp $options{$a} } keys %options;
    if($options{$sorted[0]}>0){
      $consensus.=$sorted[0];
    }
  }
  #print $consensus . "\n";
  return $consensus;
}

sub getMuscleConsensus {
  my($self,$seq_ids,$acc2seqid)=@_;
  if(scalar(@$seq_ids)<1){
    print "";
  }elsif(scalar(@$seq_ids)==1){
    my($header,$seq)=$self->getSeq($$acc2seqid{$$seq_ids[0]});
    print "$header\n";
    print "$seq\n";
  }else{
    # return first member
    my $consensus_input_file=int(rand(10000000000));
    my $consensus_alignment=$consensus_input_file;

    $consensus_input_file .= ".consensus.fa";
    $consensus_alignment  .= ".consensus.msa";

    open(CONSENSUS,">$consensus_input_file");

    for(my $s=0;$s<scalar(@$seq_ids);$s++){
      my($header,$seq)=$self->getSeq($$acc2seqid{$$seq_ids[$s]});
      print CONSENSUS "$header\n";
      print CONSENSUS "$seq\n";
    }

    `muscle -in $consensus_input_file -out $consensus_alignment -maxiters 1 -diags`;

     `rm $consensus_input_file`;
     print "Created $consensus_alignment\n";

    # output header is the first sequence,and a list of other ids used, and the total number of sequences, and %PID

    # 1. reverse orientation, if in the reverse frame
    # 2. align, using muscle (or some other alignment algorithm)
    # 3. generate consensus, by coverage in each position
    # only sequences that have started by at least that column are counted in the calculation for that position
    # therefore, long consensus sequences are generated
    # a minimum pid to the consensus is required for individual sequences.
    # if not in the threshold, the sequence is ejected back into the "clustering pool"
  }
}

sub getThisCloneTopology {
  my($self,$cluster_member_id_list,$header_field)=@_;
  # take this list of clones, and generate a tree topology
  # output is size of clones, connectivity of clones, average shm of clones

  # count the number of times each clone is seen
  my %unique_seq_counts=();
  for(my $id=0;$id<scalar(@$cluster_member_id_list);$id++){
    my $this_seq=$self->getHeaderField($$cluster_member_id_list[$id],$header_field);
    if(defined($unique_seq_counts{$this_seq})){
      $unique_seq_counts{$this_seq}++;
    }else{
      $unique_seq_counts{$this_seq}=1;
    }
  }

  # print the clone sizes to one file

  # now rotate through each clone, and assign connectivity to each other clone
  # at first do distance 1, then distance 2
  my %topology_cloneA_cloneB_distance=();
  my @unique_seqs=keys %unique_seq_counts;

  for(my $x=0;$x<scalar(@unique_seqs);$x++){
            
    #my $this_seq  = $self->getHeaderField($x,$header_field);
    #my $this_distance = $self->getSeqDistance($seed_seq,$this_seq);
    #if($this_distance<=$max_distance){
    #  push @$cluster_member_id_list,$x;
    #  push @new_seqs,$this_seq;
    #}
  }

  # print back the clone connectivity, and the clone sizes.

}
 
sub getDistance {
  my($self,$seq_id1,$seq_id2,$header_field)=@_;
  my $seq1 = $self->getHeaderField($seq_id1,$header_field);
  my $seq2 = $self->getHeaderField($seq_id2,$header_field);
  return $self->getSeqDistance($seq1,$seq2);
}

sub getSeqDistance {
  my($self,$seq1,$seq2)=@_;
  if(length($seq1) ne length($seq2)){
    if(length($seq1)>0){
      return length($seq1); #$self->max(length($seq1),length($seq2));
    }else{
      return 10;
    }
  }else{
    my @chars1=split(/ */,$seq1);
    my @chars2=split(/ */,$seq2);
    my $distance=0;
    for(my $d=0;$d<scalar(@chars1);$d++){
      if($chars1[$d] ne $chars2[$d]){
        $distance++;
      }
    }
    return $distance;
  }
}

sub getHMMSeqDistance {
  my($self,$seq1,$seq2)=@_;

  my %hmmcol2kabat_map = ( '0', 'H1', '1', 'H2', '2', 'H3', '3', 'H4', '4', 'H5', '5', 'H6', '6', 'H7', '7', 'H8', 
'8', 'H9', '9', 'H10', '10', 'H11', '11', 'H12', '12', 'H13', '13', 'H14', '14', 'H15',
'15', 'H16', '16', 'H17', '17', 'H18', '18', 'H19', '19', 'H20', '20', 'H21', '21', 'H22', 
'22', 'H23', '23', 'H24', '24', 'H25', '25', 'H26', '26', 'H27', '27', 'H28', '28', 'H29', 
'29', 'H30', '30', 'H31', '31', 'H32', '32', 'H33', '33', 'H34', '34', 'H35', '35', 'H36', 
'36', 'H37', '37', 'H38', '38', 'H39', '39', 'H40', '40', 'H41', '41', 'H42', '42', 'H43',
'43', 'H44', '44', 'H45', '45', 'H46', '46', 'H47', '47', 'H48', '48', 'H49', '49', 'H50', 
'50', 'H51', '51', 'H52', '52', 'H53', '53', 'H54', '54', 'H55', '55', 'H56', '56', 'H57',
'57', 'H58', '58', 'H59', '59', 'H60', '60', 'H61', '61', 'H62', '62', 'H63', '63', 'H64',
'64', 'H65', '65', 'H66', '66', 'H67', '67', 'H68', '68', 'H69', '69', 'H70', '70', 'H71',
'71', 'H72', '72', 'H73', '73', 'H74', '74', 'H75', '75', 'H76', '76', 'H77', '77', 'H78',
'78', 'H79', '79', 'H80', '80', 'H81', '81', 'H82', '82', 'H82A', '83', 'H82B', '84', 'H82C', 
'85', 'H83', '86', 'H84', '87', 'H85', '88', 'H86', '89', 'H87', '90', 'H88', '91', 'H89', 
'92', 'H90', '93', 'H91', '94', 'H92', '95', 'H93', '96', 'H94', '97', 'H95', '98', 'H96',
'99', 'H97', '100', 'H101', '101', 'H102', '102', 'H103', '103', 'H104', '104', 'H105', 
'105', 'H106', '106', 'H107', '107', 'H108', '108', 'H109', '109', 'H110', '110', 'H111',
'111', 'H112', '112', 'H113', '113', 'G1', '114', 'G2', '115', 'G3', '116', 'G4', '117', 'G5', 
'118', 'G6', '119', 'G7', '120', 'G8', '121', 'G9', '122', 'G10', '123', 'G11', '124', 'G12',
'125', 'G13', '126', 'G14', '127', 'G15', '128', 'L1', '129', 'L2', '130', 'L3', '131', 'L4',
'132', 'L5', '133', 'L6', '134', 'L7', '135', 'L8', '136', 'L9', '137', 'L10', '138', 'L11',
'139', 'L12', '140', 'L13', '141', 'L14', '142', 'L15', '143', 'L16', '144', 'L17', '145', 'L18',
'146', 'L19', '147', 'L20', '148', 'L21', '149', 'L22', '150', 'L23', '151', 'L24', '152', 'L25',
'153', 'L26', '154', 'L27', '155', 'L28', '156', 'L29', '157', 'L30', '158', 'L31', '159', 'L32',
'160', 'L33', '161', 'L34', '162', 'L35', '163', 'L36', '164', 'L37', '165', 'L38', '166', 'L39',
'167', 'L40', '168', 'L41', '169', 'L42', '170', 'L43', '171', 'L44', '172', 'L45', '173', 'L46',
'174', 'L47', '175', 'L48', '176', 'L49', '177', 'L50', '178', 'L51', '179', 'L52', '180', 'L53',
'181', 'L54', '182', 'L55', '183', 'L56', '184', 'L57', '185', 'L58', '186', 'L59', '187', 'L60',
'188', 'L61', '189', 'L62', '190', 'L63', '191', 'L64', '192', 'L65', '193', 'L66', '194', 'L67',
'195', 'L68', '196', 'L69', '197', 'L70', '198', 'L71', '199', 'L72', '200', 'L73', '201', 'L74',
'202', 'L75', '203', 'L76', '204', 'L77', '205', 'L78', '206', 'L79', '207', 'L80', '208', 'L81',
'209', 'L82', '210', 'L83', '211', 'L84', '212', 'L85', '213', 'L86', '214', 'L87', '215', 'L88',
'216', 'L89', '217', 'L90', '218', 'L91', '219', 'L92', '220', 'L93', '221', 'L94', '222', 'L95',
'223', 'L95A', '224', 'L96', '225', 'L97', '226', 'L98', '227', 'L99', '228', 'L100', '229', 'L101',
'230', 'L102', '231', 'L103', '232', 'L104', '233', 'L105', '234', 'L106', '235', 'L107');

  if(length($seq1) ne length($seq2)){
    if(length($seq1)<length($seq2)){
      return (length($seq1),length($seq1),""); 
    }else{
      return (length($seq2),length($seq1),"");
    }
  }else{
    my @chars1=split(/ */,$seq1);
    my @chars2=split(/ */,$seq2);
    my $distance=0;
    my $non_empty_columns=0;
    my $mutation_list="";
    for(my $d=0;$d<scalar(@chars1);$d++){
      if( ($chars1[$d] ne "-") && ($chars2[$d] ne "-")){
        $non_empty_columns++;
        if($chars1[$d] ne $chars2[$d]){
          # d must be mapped to hmmcol positions
          $mutation_list .= " " . $hmmcol2kabat_map{$d} . "(" . $chars1[$d] . "-" . $chars2[$d] . ")";
          $distance++;
        }
      }
    }
    return ($distance,$non_empty_columns,$mutation_list);
  }
}

sub max {
  my($self,$x,$y)=@_;
  if($x<$y){
    return $y;
  }else{
    return $x;
  }
}

sub writeSHMUsage {
  my($self,$outfile)=@_;

  open(OUT,">$outfile");

  for(my $x=0;$x<$self->{nseqs};$x++){
    if(defined($self->getHeaderField($x,1))){
      my $vgene         = $self->getHeaderField($x,1);
      my @vgene_fields  = split(/ /,$vgene);
      if(scalar(@vgene_fields) == 3){
        if($vgene_fields[1] > 249){
          print OUT $vgene_fields[2] . "\n";
        }
      }
    }
  }
  close(OUT);
}

sub getSHMvCounts {
  my($self,$outfile)=@_;

  my %clones=();

  for(my $x=0;$x<$self->{nseqs};$x++){
    # a clone is a unique Vgene-Jgene-CDRH3 combination
    if(defined($self->getHeaderField($x,4))){
      my $vgene         = $self->getHeaderField($x,1);
      my $jgene         = $self->getHeaderField($x,3);
      my $h3            = $self->getHeaderField($x,4);
      if($h3=~m/^C..*W$/){
        my @vgene_fields  = split(/ /,$vgene);
        my @jgene_fields  = split(/ /,$jgene);
        if(scalar(@vgene_fields) == 3){
	  if(scalar(@jgene_fields) == 3){
            if($vgene_fields[1] > 249){
              my $clone = $vgene_fields[0] . "_" . $jgene_fields[0] . "_" . $h3;
              if(defined($clones{$h3})){
                $clones{$h3} .= " " . ($vgene_fields[2] + $jgene_fields[2]);
              }else{
                $clones{$h3} = ($vgene_fields[2] + $jgene_fields[2]);
              }
            }
          }
        }
      }
    }
  }
  # ok, now for each clone, report the counts, and the average SHM
  open(OUT,">$outfile"); 
  my @clone_names = keys %clones;
  for(my $c=0;$c<scalar(@clone_names);$c++){
    my @shm_scores = split(/ /,$clones{$clone_names[$c]});
    my $counts = scalar(@shm_scores);
    my $average=0;
    if($counts>0){
      for(my $x=0;$x<scalar(@shm_scores);$x++){
        $average+=$shm_scores[$x];
      }   
      $average/=$counts;
    }
    print OUT $counts . "\t" . $average . "\t" . $clone_names[$c] . "\n";
  } 
  close(OUT);
}

sub getGermlineUsage {
  my($self,$outfile)=@_;
 
  my %vGeneHash = ();

  open(OUT,">$outfile");
  my $total_counts=0;

  for(my $x=0;$x<$self->{nseqs};$x++){
    if(defined($self->getHeaderField($x,1))){
      my $vgene         = $self->getHeaderField($x,1);
      my @vgene_fields  = split(/ /,$vgene);
      if(scalar(@vgene_fields) == 3){
        if($vgene_fields[1] > 99){
          $total_counts++;
          if(defined $vGeneHash{$vgene_fields[0]}){
            $vGeneHash{$vgene_fields[0]}++;
          }else{
            $vGeneHash{$vgene_fields[0]}=1;
          } 
        }
      }
    }
  }
  foreach my $key (sort keys %vGeneHash) {
    if($key=~m/^IGHV/){
      print OUT $key . " " . (int(($vGeneHash{$key}/$total_counts)*10000)/100) . "\n";
    }
  }
  close(OUT);
}
 
sub getGermlineFreqs {
  my($self,$max_shm)=@_;
  my %vGeneHash = ();
  my $total_counts=0;

  for(my $x=0;$x<$self->{nseqs};$x++){
    if(defined($self->getHeaderField($x,1))){
      my $vgene         = $self->getHeaderField($x,1);
      my @vgene_fields  = split(/ /,$vgene);
      if(scalar(@vgene_fields) == 3){
        if($vgene_fields[1] > 99){
          if($vgene_fields[2] <= $max_shm){
            $total_counts++;
            if(defined $vGeneHash{$vgene_fields[0]}){
              $vGeneHash{$vgene_fields[0]}++;
            }else{
              $vGeneHash{$vgene_fields[0]}=1;
            }
          }
        }
      }
    }
  }
  # normalize
  foreach my $key (sort keys %vGeneHash) {
    $vGeneHash{$key}= $vGeneHash{$key} / $total_counts;
  }
  return %vGeneHash;
}

sub detectMissingScFvDomains {
  my($self,$nterm_vector_dna,$linker_vector_dna,$cterm_vector_dna)=@_;
  my %output=();
  return %output;
}


sub detectDNALiabilities {
  my($self)=@_;
  my %output_hash=();

  my %dna_liabilities = $self->getLiabilityHash("$rootdir/db/liabilities.nucleotide.txt");
  for(my $seqid=0;$seqid<$self->getSeqCount();$seqid++){
    my ($header,$dna) = $self->getSeq($seqid);
    my $reverse_seq   = $self->getReverseStrand($seqid);
    $dna .= "xxxxxxxxxxx" . $reverse_seq;
    my $output_liability_line="";
    while ( my ($liability_name, $liability_pattern) = each(%dna_liabilities) ) {
      if ($dna =~ m/$liability_pattern/) {
       $output_liability_line.= "$liability_name ($liability_pattern),";
      }
    }
    $output_liability_line=~s/,$//;
    if(length($output_liability_line)>0){
      my $acc=$self->getAccession($seqid);
      $output_hash{$acc}=$output_liability_line;
    }    
  } 
  return %output_hash;
}

sub getVHSequence {
  my($self)=@_;
  my $nseqs = $self->getSeqCount();
  my %output=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $seq=$self->getHMMCOLRange($seqid,0,112);
       $seq=~s/\.*//g;
       $seq=~s/^[a-z]*([a-z])/$1/g;
    unless($seq=~m/^[a-z]-/){
      $seq=~s/^[a-z]//;
    }
    $seq=uc($seq);
    $seq=~s/\-*//g;
    if(length($seq)<20){
      $seq="";
    }
    $output{$self->getAccession($seqid)}=$seq;
  } 
  return %output; 
}

sub getVLSequence {
  my($self,$VLFamilies)=@_;
  my $nseqs = $self->getSeqCount();
  my %output=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $seq=$self->getHMMCOLRange($seqid,128,235);
       $seq=~s/\.*//g;
       $seq=~s/^[a-z]*//;
       $seq=~s/\-*//g;
       $seq=uc($seq);
    my $accession=$self->getAccession($seqid);
    if(defined($$VLFamilies{$accession})){
      if($$VLFamilies{$accession} =~ m/IGLV/){
        $seq=~s/^.//;
      }
    }

    if(length($seq)<20){
      $seq="";
    }
    $output{$self->getAccession($seqid)}=$seq;
  }
  return %output;
}

sub getStructuralProblems {
  my($self)=@_;
  my $nseqs = $self->getSeqCount();
  my %output=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $char_c1=$self->getHMMCOLRange($seqid,21,21);
    my $char_c2=$self->getHMMCOLRange($seqid,94,94);
    unless($char_c1 =~ m/[C-]/){
      $output{$self->getAccession($seqid)}="VH c-c";
    }
    unless($char_c2 =~ m/[C-]/){
      $output{$self->getAccession($seqid)}="VH c-c";
    }
    my $char_c3=$self->getHMMCOLRange($seqid,150,150);
    my $char_c4=$self->getHMMCOLRange($seqid,215,215);
    unless($char_c3 =~ m/[C-]/){
      if(defined($output{$self->getAccession($seqid)})){
        $output{$self->getAccession($seqid)}.=", VL c-c";
      }else{
        $output{$self->getAccession($seqid)}="VL c-c";
      }
    } 
    unless($char_c3 =~ m/[C-]/){
      if(defined($output{$self->getAccession($seqid)})){
        $output{$self->getAccession($seqid)}.=", VL c-c";
      }else{
        $output{$self->getAccession($seqid)}="VL c-c";
      }
    }
    #my $char_w103 = $self->getHMMCOLRange($seqid,102,102); 
    #my $char_f100 = $self->getHMMCOLRange($seqid,226,226);
    #print $char_c1 . "\t" . $char_c2 . "\t" . $char_c3 . "\t" . $char_c4 . "\t" . $char_w103 . "\t" . $char_f100 . "\n";
  }
  return %output;
}

sub getStopCodons {
  my($self)=@_;
  my $nseqs = $self->getSeqCount();
  my %output=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $seq=$self->getHMMCOLRange($seqid,0,235);
       $seq=~s/\.*//g;
       $seq=~s/^[a-z]*//;
       $seq=uc($seq);
    
    if($seq=~m/X/){
      my $count = ($seq =~ tr/X//);
      $output{$self->getAccession($seqid)}="$count stop codons";
    }elsif($seq=~m/Z/){
      my $count = ($seq =~ tr/Z//);
      $output{$self->getAccession($seqid)}="$count unreliable codons";
    }else{
      $output{$self->getAccession($seqid)}="";
    }
  }
  return %output;
}

sub assembleCommonAccessions {
  my($self,$outfile)=@_;

  open(OUT,">$outfile");
  #my %dnaAccSeqid=();
  #   $self->getAccession2seqidHash(\%dnaAccSeqid);

  my %assembled_sequences=();
 
  my $nseqs = $self->getSeqCount();
  my %output=();

  # collect unique accession counts
  my %uniqueAccessions=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc = $self->getAccession($seqid);
    if(defined($uniqueAccessions{$acc})){
      $uniqueAccessions{$acc}.=",$seqid";
    }else{
      $uniqueAccessions{$acc}=$seqid;
    }
  }

  # now deal with any frameshift cases
  my @uniqAccs=keys %uniqueAccessions;
  for(my $acc=0;$acc<scalar(@uniqAccs);$acc++){
    my @ids=split(/,/,$uniqueAccessions{$uniqAccs[$acc]});
    my $acc = $self->getAccession($ids[0]);
    print OUT ">$acc\n";
    my $assembled_seq=$self->ids2assembly(@ids);
       $assembled_seq=~s/\.*//g;
       $assembled_seq=~s/-*//g;
       $assembled_seq=uc($assembled_seq);
    print OUT $assembled_seq . "\n";
    if(scalar(@ids)>1){
      #print "AN assembly!\n";
      $assembled_sequences{$acc}="frameshift";
    }
  }
  close OUT;
  return %assembled_sequences;
}

sub seq2Chars {
  my($self,$seqid)=@_;
  my($header,$seq)=$self->getSeq($seqid);
  $seq=~s/^[a-z\.]*//;
  $seq=~s/[a-z\.]*$//;
  my @chars=split(/ */,$seq);
  return @chars;
}

sub ids2assembly {
  my($self,@ids)=@_;
  # set the final character string to the first one observed
  my($header,$seq)=$self->getSeq($ids[0]);
  my @initial_characters=$self->seq2Chars($ids[0]);
  my @final_chars=@initial_characters;

  # now process all of the new characters
  for(my $i=1;$i<scalar(@ids);$i++){
    #my($this_header,$this_seq)=$self->getSeq($ids[$i]);
    my @char1=@final_chars;
    my @char2=$self->seq2Chars($ids[$i]);
    #ok, now progress through all characters
    for(my $c=0;$c<scalar(@char1);$c++){
      my $conf1=$self->getSlidingWindowConfidence($c,\@char1);
      my $conf2=$self->getSlidingWindowConfidence($c,\@char2);
      my $final_char="z";
      if($char1[$c] eq $char2[$c]){
        $final_char=$char1[$c];
      } elsif ($char1[$c] =~ m/[ZzXx]/){
        $final_char=$char2[$c]; 
      } elsif ($char2[$c] =~ m/[ZzXx]/){
        $final_char=$char1[$c];
      } elsif ($conf2 > ( $conf1 + 0.6) ){
        $final_char=$char2[$c];
      } elsif ($conf1 > ( $conf2 + 0.6) ){
        $final_char=$char1[$c];
      } elsif ($char1[$c] =~ m/[a-z\.]/) {
        $final_char="z";
      } elsif ($char1[$c] =~ m/[A-Z\-]/) {
        $final_char="Z";
      } else {
        $final_char="?";
      }
      #print $conf1 . "\t" . $conf2 . "\t" . $chars[$c] . "\t" . $char2[$c];
      #print "\t" . $final_char . "\n";
      $final_chars[$c]=$final_char;
    }
    #print join('',@char1) . "\n";
    #print join('',@char2) . "\n";
    #print "=\n";
  }
  my $final_seq=join('',@final_chars);
  #print $final_seq . "\n";
  #print "\n==============\n";
  $final_seq=~s/\.*//g;
  #print $final_seq . "\n";
  return($final_seq);
}

sub getSlidingWindowConfidence {
  my($self,$position,$chars)=@_;

  my $total_positions=0;
  my $total_good_positions=0;

  # working forward
  my $forward_positions_checked=0;
  my $i=$position;
  while(($forward_positions_checked<4)&&($i<scalar(@$chars))){
    unless($$chars[$i]=~m/\./){
      $forward_positions_checked++;
      $total_positions++;
      if($$chars[$i]=~m/[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]/){
        $total_good_positions++;
      }
    }
    $i++;
  }
  # working backward
  $i=$position-1;
  my $reverse_positions_checked=0;
  while(($reverse_positions_checked<4)&&($i>=0)){
    unless($$chars[$i]=~m/\./){
      $reverse_positions_checked++;
      $total_positions++;
      if($$chars[$i]=~m/[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]/){
        $total_good_positions++;
      }
    }
    $i--;
  }
  my $score=int(10 * ($total_good_positions / $total_positions))/10;
  return $score;
}

sub getLinkerSequences {
  my($self)=@_;
  my $nseqs = $self->getSeqCount();
  my %output=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $seq=$self->getHMMCOLRange($seqid,110,130);
    $seq=~s/\.*//g;
    unless($seq=~m/^[A-Z]/){
      $seq="";
    }
    unless($seq=~m/[A-Z]$/){
      $seq="";
    }
    $seq=~s/^...//;
    if($seq=~m/G..$/){
      $seq=~s/..$//;
    }else{
      $seq=~s/...$//;
    }
    $seq=uc($seq);
    $seq=~s/\-*//g;
    $output{$self->getAccession($seqid)}=$seq;
    #print $seq . "\n";
  }
  return %output;
}

sub getLeaderSequences {
  my($self)=@_;

  my $nseqs = $self->getSeqCount();
  my %output=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my($header,$seq)=$self->getSeq($seqid);
    $seq=~s/\.*//g;
    $seq=~s/[A-Z].*//;
    if($seq=~m/\-\-/){
      $seq="";
    }
    $seq=~s/[a-z]\-$//;
    $seq=~s/[^m]*//;
    if(length($seq)<15){
      $seq="";
    }
    unless($seq=~m/m/){
      $seq="";
    }
           
    if(length($seq)>40){
      $seq=~s/.*\(........................................\)$/$1/;
    } 
    #print $seq . "\n";
    $output{$self->getAccession($seqid)}=$seq;
  }
  return %output;
}

sub getVBaseNames {
  my($self,$domain)=@_;

  my $domain_field=1;
  if($domain eq "VH"){
    $domain_field=1;
  }elsif ($domain eq "VL"){
    $domain_field=12;
  }

  open(FILE,$self->{vbasefile});
  my @lines=<FILE>;
  chomp(@lines);
  close(FILE);
  
  # collect imgt->vbase conversion
  my %vbase_names=();
  for(my $i=0;$i<scalar(@lines);$i++){
    my($imgt,$vbase)=split(/\t/,$lines[$i]);
    $vbase_names{$imgt}=$vbase;	
  }  

  # perform vbase map assignment
  my $nseqs = $self->getSeqCount();
  my %output=();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $vh=$self->getHeaderField($seqid,$domain_field);
       $vh=~s/ .*//;
    #print $vh . "\n";
    if(defined($vbase_names{$vh})){
      $output{$self->getAccession($seqid)}=$vbase_names{$vh};
      #print $vbase_names{$vh} . "\n";
    }
  }  
  return %output;
}

sub getVHFamilies {
  my($self)=@_;
  my %output=(); 

  my $nseqs = $self->getSeqCount(); 
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc=$self->getAccession($seqid);
    my $vh=$self->getHeaderField($seqid,1);
       $vh=~s/\-.*//;
    $output{$acc}=$vh;
  }
  return %output;
}

sub getVLFamilies {
  my($self)=@_;
  my %output=();

  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc=$self->getAccession($seqid);
    my $vl=$self->getHeaderField($seqid,12);
       $vl=~s/\-.*//;
    $output{$acc}=$vl;
  }
  return %output;
}

sub getVHGermlines {
  my($self)=@_;
  my %output=();

  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc=$self->getAccession($seqid);
    my $vh=$self->getHeaderField($seqid,1);
       $vh=~s/ .*//;
    $output{$acc}=$vh;
  }
  return %output;
}

sub getVLGermlines {
  my($self)=@_;
  my %output=();
  
  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc=$self->getAccession($seqid);
    my $vl=$self->getHeaderField($seqid,12);
       $vl=~s/ .*//;
    $output{$acc}=$vl;
  }
  return %output;
}   

sub getVHGermlinePID {
  my($self)=@_;
  my %output=();
  
  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc=$self->getAccession($seqid);
    my @field_vh=split(/ /,$self->getHeaderField($seqid,1));
    my $pid="NA";
    if(scalar(@field_vh)>2){
      my $mismatch=$field_vh[2];
      my $alnlen=$field_vh[1];
      if($alnlen>0){
         $pid = (int(1000 * (($alnlen - $mismatch)/$alnlen)))/10;
      }
    }
    $output{$acc}=$pid;
  }
  return %output;
}

sub getVLGermlinePID {
  my($self)=@_;
  my %output=();
  
  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc=$self->getAccession($seqid);
    my @field_vl=split(/ /,$self->getHeaderField($seqid,12));
    my $pid="NA";
    if(scalar(@field_vl)>2){
      my $mismatch=$field_vl[2];
      my $alnlen=$field_vl[1];
      if($alnlen>0){
         $pid = (int(1000 * (($alnlen-$mismatch)/$alnlen)))/10;
      }
    }
    $output{$acc}=$pid;
  }
  return %output;
}


sub detectLiabilities {
  my($self)=@_;
  my %sequence_liability_assignments=();

  my %aa_liabilities = $self->getLiabilityHash("$rootdir/db/liabilities.aminoacid.txt");
  my $nseqs = $self->getSeqCount();
  my @cdr_name = ("H1","H2","H3","L1","L2","L3");

  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $paratope = $self->getParatope($seqid);
    my @paratope_split =split(/-/,$paratope);
    my $comma_count=0;

    if(defined($paratope_split[0])){
      $paratope_split[0] =~ s/M(.)$/A$1/;
    }
    if(defined($paratope_split[2])){
      $paratope_split[2] =~ s/^C//;
    }
    if(defined($paratope_split[5])){
      $paratope_split[5] =~ s/^C//;
    }

    my $liability_report="";
    #For loop for CDR Heavy Chain
    for(my $i=0;$i<6;$i++){
      if(defined($paratope_split[$i])){
        if($paratope_split[$i] ne ""){
          while ( my ($liability_name, $liability_pattern) = each(%aa_liabilities) ) {
            if ($paratope_split[$i] =~ m/$liability_pattern/) {
              if($comma_count!=0){
                $liability_report.=",";
              }
              $comma_count++;
              $liability_report.= $cdr_name[$i] . ":$liability_name ($liability_pattern)";
            }
          }
        }
      }
    }
    my $acc=$self->getAccession($seqid);
    $sequence_liability_assignments{$acc}=$liability_report;
  }
  return %sequence_liability_assignments;
}

sub getCanonicalHash {
  my($self,$cdr)=@_;
  my %output=();
  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc=$self->getAccession($seqid);

    if($cdr eq "H1"){
      my $h1 = $self->getHeaderField($seqid,8); 
      if(length($h1) == 0){
        $output{$acc}="";
      }elsif(length($h1) == 9){
        $output{$acc}="H1-MLT-canon1";
      }elsif(length($h1) == 10){
        $output{$acc}="H1-MLT-canon2";
      }elsif(length($h1) == 11){
        $output{$acc}="H1-MLT-canon3";
      }else{
        $output{$acc}="unknown";
      }
      #$output{$acc}.= " $h1\t" . length($h1);

    }elsif($cdr eq "H2"){
      my $h2 = $self->getHeaderField($seqid,9);
      if(length($h2) == 0){
        $output{$acc}="";
      }elsif(length($h2) == 12){
        $output{$acc}="H2-MLT-canon1";
      }elsif(length($h2) == 13){
        $output{$acc}="H2-MLT-canon2";
      }elsif(length($h2) == 14){
        $output{$acc}="H2-MLT-canon3";
      }elsif(length($h2) == 15){
        $output{$acc}="H2-MLT-canon4";
      }else{
        $output{$acc}="unknown";
      }
    }elsif($cdr eq "H3"){
      my $h3 = $self->getHeaderField($seqid,4);
      if(length($h3) == 0){
        $output{$acc}="";
      }elsif(($h3=~m/^..[RK]/)&&($h3=~m/D..$/)){
        $output{$acc}="H3-kinked";
      }else{
        $output{$acc}="H3-unkinked";
      }
    }elsif($cdr eq "L1"){
      my $l1 = $self->getHeaderField($seqid,10);
      if(length($l1) == 0){
        $output{$acc}="";
      }elsif(length($l1) == 18){
        $output{$acc}="L1-MLT-canon3";
      }elsif(length($l1) == 17){
        $output{$acc}="L1-MLT-canon4";
      }elsif(length($l1) == 16){
        $output{$acc}="L1-MLT-canon5";
      }elsif(length($l1) == 12){
        $output{$acc}="L1-MLT-canon6";
      }elsif(length($l1) == 11){
        $output{$acc}="L1-MLT-canon2";
      }elsif(length($l1) == 10){
        $output{$acc}="L1-MLT-canon1";
      }else{
        $output{$acc}="unknown";
      }
    }elsif($cdr eq "L2"){
      my $l2 = $self->getHeaderField($seqid,11);
      if(length($l2) == 0){
        $output{$acc}="";
      }elsif(length($l2) == 7){
        $output{$acc}="L2-MLT-canon1";
      }else{
        $output{$acc}="unknown";
      }
    }elsif($cdr eq "L3"){
      my $l3 = $self->getHeaderField($seqid,5);
      $output{$acc}="";
      if(length($l3) == 0){
        $output{$acc}="";
      }elsif(length($l3) == 10){
        $output{$acc}="L3-MLT-canon3";
      }elsif(length($l3) == 11){
        if($l3=~m/^C......P/){
	  $output{$acc}="L3-MLT-canon3";
        }else{
          $output{$acc}="L3-MLT-canon6";
        }
      }elsif(length($l3) == 12){
        if($l3=~m/^C.......P/){
          $output{$acc}="L3-MLT-canon1";
        }else{
          $output{$acc}="L3-MLT-canon2";
        }
      }elsif(length($l3) == 13){
        if($l3=~m/^C........P/){
          $output{$acc}="L3-MLT-canon5";
        }else{
          $output{$acc}="unknown";
        }
      }else{
        $output{$acc}="unknown";
      }    
    }
  }
  return %output;
}

sub getLiabilityHash {
  my($self,$file)=@_;
  my %liabilities=();
  if(-f $file){
    open(FILE,$file);
    my @lines=<FILE>;
    chomp(@lines);
    close(FILE);
    for(my $x=0;$x<scalar(@lines);$x++){
      my @fields=split(/\t/,$lines[$x]);
      $liabilities{$fields[1]}=$fields[0];
    }
  }
  return %liabilities;
}

sub getHydrophobicity {
  my($self)=@_;

  my %output_hash=();
  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc = $self->getAccession($seqid);
    
    # confirm that all 6 cdrs are present
    my $paratope="";
    my $h3 = $self->getHeaderField($seqid,4);
    my $h1 = $self->getHeaderField($seqid,8);
    my $h2 = $self->getHeaderField($seqid,9);
    my $l3 = $self->getHeaderField($seqid,5);
    my $l1 = $self->getHeaderField($seqid,10);
    my $l2 = $self->getHeaderField($seqid,11);
 
    if(length($h1)>0){
      if(length($h2)>0){
        if(length($h3)>0){
          if(length($l1)>0){
            if(length($h2)>0){
              if(length($l3)>0){
                my $sequence = $h1 . $h2 . $h3 . $l1 . $l2 . $l3;

                my @chars=split(/ */,$sequence);
                my $oil=0;
                my %isOily=('A',1,'F',1,'G',1,'I',1,'L',1,'V',1,'W',1,'Y',1);
                for(my $c=0;$c<scalar(@chars);$c++){
                  if(defined($isOily{$chars[$c]})){
                    $oil+=1;
                  }
                }
                $output_hash{$acc}=(int(1000 * ($oil/scalar(@chars))))/10;
              }
            }
          }
        }
      }
    }
  }
  return %output_hash;
}

sub getCharge {
  my($self)=@_;

  my %output_hash=();
  my $nseqs = $self->getSeqCount();
  for(my $seqid=0;$seqid<$nseqs;$seqid++){
    my $acc = $self->getAccession($seqid);

    # confirm that all 6 cdrs are present
    my $paratope="";
    my $h3 = $self->getHeaderField($seqid,4);
    my $h1 = $self->getHeaderField($seqid,8);
    my $h2 = $self->getHeaderField($seqid,9);
    my $l3 = $self->getHeaderField($seqid,5);
    my $l1 = $self->getHeaderField($seqid,10);
    my $l2 = $self->getHeaderField($seqid,11);

    if(length($h1)>0){
      if(length($h2)>0){
        if(length($h3)>0){
          if(length($l1)>0){
            if(length($h2)>0){
              if(length($l3)>0){
                my $sequence = $h1 . $h2 . $h3 . $l1 . $l2 . $l3;

                my @chars=split(/ */,$sequence);
                my $formalCharge=0;
                my %thisCharge=('D',-1,'E',-1,'K',1,'R',1);

                for(my $c=0;$c<scalar(@chars);$c++){
                  if(defined($thisCharge{$chars[$c]})){
                    $formalCharge+=$thisCharge{$chars[$c]};
                  }
                }
                $output_hash{$acc}=$formalCharge;
              }
            }
          }
        }
      }
    }
  }
  return %output_hash;
}

sub getParatope {
  my($self,$seqid)=@_;

  my $h3 = $self->getHeaderField($seqid,4);
  my $h1 = $self->getHeaderField($seqid,8);
  my $h2 = $self->getHeaderField($seqid,9);
  my $l3 = $self->getHeaderField($seqid,5);
  my $l1 = $self->getHeaderField($seqid,10);
  my $l2 = $self->getHeaderField($seqid,11);

  #print "$h1-$h2-$h3-$l1-$l2-$l3\n";
  return "$h1-$h2-$h3-$l1-$l2-$l3";
}

sub getParatopeGroup {
  my($self,$seqid)=@_;

  my $h3 = length($self->getHeaderField($seqid,4));
  my $h1 = length($self->getHeaderField($seqid,8));
  my $h2 = length($self->getHeaderField($seqid,9));
  my $l3 = length($self->getHeaderField($seqid,5));
  my $l1 = length($self->getHeaderField($seqid,10));
  my $l2 = length($self->getHeaderField($seqid,11));

  # blank CDRs unless the whole chain can be seen
  #unless( ($h1>0) && ($h2>0) && ($h3>0) ){
  #  $h1=0;
  #  $h2=0;
  #  $h3=0;
  #}
  #unless( ($l1>0) && ($l2>0) && ($l3>0) ){
  #  $l1=0;
  #  $l2=0;
  #  $l3=0;
  #}
  return "$h1-$h2-$h3-$l1-$l2-$l3";
}

sub getSHMPIDhash {
  my($self)=@_;
  my %shm_hash=();
 
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $pid = $self->getPercentSHM($x,180);
    my $acc = $self->getAccession($x);
    $shm_hash{$acc}=$pid;
  } 
  return %shm_hash;
}

sub getPercentSHM {
  my($self,$seqid,$minlength)=@_;
  my $pid="";
  my $length=0;
  my $mismatch=0;  
  my @vseg_fields = split(/ /,$self->getHeaderField($seqid,1));
  my @jseg_fields = split(/ /,$self->getHeaderField($seqid,3));
  my @vlseg_fields = split(/ /,$self->getHeaderField($seqid,12));
  my @jlseg_fields = split(/ /,$self->getHeaderField($seqid,13));
 
  if(scalar(@vseg_fields)==3){
    if($vseg_fields[1]>$minlength){
      $length=$vseg_fields[1];
      $mismatch=$vseg_fields[2];
      if(scalar(@jseg_fields)==3){
        $length   += $jseg_fields[1];
        $mismatch += $jseg_fields[2];
      }
    }
  }
  if(scalar(@vlseg_fields)==3){
    if($vlseg_fields[1]>$minlength){
      $length+=$vlseg_fields[1];
      $mismatch+=$vlseg_fields[2];
      if(scalar(@jlseg_fields)==3){
        $length   += $jlseg_fields[1];
        $mismatch += $jlseg_fields[2];
      }
    }
  }
  if($length != 0){
    $pid=( ($length - $mismatch) / $length);
    $pid = int(100 * $pid)/100;
  }
  return $pid;
}

sub getSuperQuickCloneMode {
  my($self,$domain,$id_array)=@_;
  # get list of unique H1,H2,H3
  # store PID SHM for each unique combo
  # collapse in sequences that don't have all CDRs but could be part of a clone (missing CDR1 or CDR2)

  my %clone_counts=();
  my %clone_to_aligned_bases=();
  my %clone_to_mismatched_bases=();
  my %clone_to_dsegments=();
  my %clone_to_jsegments=();
  my %clone_accessions=();

  my @vgene=();
  my @jgene=();
  my @dgene=();
  my $cdr1="";
  my $cdr2="";
  my $cdr3="";

  # first gather all the unique clones into a massive hash
  for(my $x=0;$x<scalar(@$id_array);$x++){
    # map x->s
    my $s=$$id_array[$x];

    my $acc = $self->getAccession($s);
    if(($domain eq "H") or ($domain eq "VH")){
      @vgene=split(/  */,$self->getHeaderField($s,1));
      unless(defined($vgene[0])){
        $vgene[0]="";
      }
      @jgene=split(/  */,$self->getHeaderField($s,3));
      unless(defined($jgene[0])){
        $jgene[0]="";
      }
      @dgene=split(/  */,$self->getHeaderField($s,2));
      unless(defined($dgene[0])){
        $dgene[0]="";
      }
      $cdr1=$self->getHeaderField($s,8);
      $cdr2=$self->getHeaderField($s,9);
      $cdr3=$self->getHeaderField($s,4);
    }elsif(($domain eq "K") or ($domain eq "VK") or ($domain eq "L") or ($domain eq "VL")){
      @vgene=split(/  */,$self->getHeaderField($s,12));
      unless(defined($vgene[0])){
        $vgene[0]="";
      }
      @jgene=split(/  */,$self->getHeaderField($s,13));
      unless(defined($jgene[0])){
        $jgene[0]="";
      }
      $cdr1=$self->getHeaderField($s,10);
      $cdr2=$self->getHeaderField($s,11);
      $cdr3=$self->getHeaderField($s,5);
    }
   
    my $clone_hash=$vgene[0] . ";" . $cdr1 . ";" . $cdr2 . ";" . $cdr3;

    # deal with counts
    if(defined($clone_counts{$clone_hash})){
      $clone_counts{$clone_hash}++;
    }else{
      $clone_counts{$clone_hash}=1;
    }

    # store SHM information
    if(defined($vgene[2])){
      if(defined($clone_to_aligned_bases{$clone_hash})){
        $clone_to_aligned_bases{$clone_hash}    += $vgene[1];
        $clone_to_mismatched_bases{$clone_hash} += $vgene[2];
      }else{
        $clone_to_aligned_bases{$clone_hash}     = $vgene[1];
        $clone_to_mismatched_bases{$clone_hash}  = $vgene[2];
      }
    }

    # store d-segments
    if(defined($dgene[2])){
      if(defined($clone_to_dsegments{$clone_hash})){
        $clone_to_dsegments{$clone_hash}+=" " . $dgene[0];
      }else{
        $clone_to_dsegments{$clone_hash}=$dgene[0];
      }
    }

    # store j-segments
    if(defined($jgene[2])){
      if(defined($clone_to_jsegments{$clone_hash})){
        $clone_to_jsegments{$clone_hash}+=" " . $jgene[0];
      }else{
        $clone_to_jsegments{$clone_hash}=$jgene[0];
      }
    }

    # store accessions
    if(defined($clone_accessions{$clone_hash})){
      $clone_accessions{$clone_hash}+=" " . $acc;
    }else{
      $clone_accessions{$clone_hash}=$acc;
    }
  }

  # get the sorted version by frequency
  my @ranked_clone_list=sort {$clone_counts{$b} <=> $clone_counts{$a}} keys %clone_counts;
  for(my $x=0;$x<scalar(@ranked_clone_list);$x++){
    print $x . "\t" . $ranked_clone_list[$x] . "\t" . $clone_counts{$ranked_clone_list[$x]} . "\n";
  } 
  # collapse clone groups that might belong together

  # all collected. now collapse SHM d-segments j-segments

}

sub getQuickCloneMode {
  my($self,$id_array)=@_;

  #my($self,$this_v,$this_j,$this_h3)=@_;

  my ($dominant_h1,$dominant_h2,$avg_shm,$dominant_dseg)=("","","","");

  my %cdr_dgene = ();
  my %cdr1      = ();
  my %cdr2      = ();
  my $matches   = 0;
  my $mutations = 0;
  my $all_accs  = "";

  for(my $x=0;$x<scalar(@$id_array);$x++){
    my $s = $$id_array[$x];
    my @vgene=split(/  */,$self->getHeaderField($s,1));
    my @vkgene=split(/  */,$self->getHeaderField($s,12));
    if(scalar(@vgene)==3 or scalar(@vkgene)==3){
      if(1 eq 1){
        my @jgene=split(/  */,$self->getHeaderField($s,3));
        my @jkgene=split(/  */,$self->getHeaderField($s,13));
        if(scalar(@jgene)==3 or scalar(@jkgene)==3){
          if(1 eq 1){
            my $h3=$self->getHeaderField($s,4);
            my $l3=$self->getHeaderField($s,5);
            if($h3 || $l3 ){
              $all_accs .= $self->getAccession($s) . ",";
              # d-segs
              my @dgene=split(/  */,$self->getHeaderField($s,2));
              if(!defined($dgene[0])){
                $dgene[0]="";
              }
              if(defined($cdr_dgene{$dgene[0]})){
                $cdr_dgene{$dgene[0]}++;
              }else{
                $cdr_dgene{$dgene[0]}=1;
              }

              # cdr1 8,10
                my $h1 = $self->getHeaderField($s,8);
                my $l1 = $self->getHeaderField($s,10);
                my $this_cdr1 = "";
                if(length($h1)>0){
                  $this_cdr1=$h1;
                }elsif(length($l1)>0){
                  $this_cdr1=$l1;
                }
                if(defined($cdr1{$this_cdr1})){
                  $cdr1{$this_cdr1}++;
                }else{
                  $cdr1{$this_cdr1}=1;
                }

              # cdr2 9,11
                my $cdr2 = "";
                my $h2 = $self->getHeaderField($s,9);
                my $l2 = $self->getHeaderField($s,11);
                my $this_cdr2="";
                if(length($h2)>0){
                  $this_cdr2=$h2;
                }elsif(length($l2)>0){
                  $this_cdr2=$l2;
                }
                if(defined($cdr1{$this_cdr1})){
                  $cdr2{$this_cdr2}++;
                }else{
                  $cdr2{$this_cdr2}=1;
                }

              # SHM
              
              if(defined($vgene[1])){
                $matches   += $vgene[1];
                $mutations += $vgene[2];
              }elsif(defined($vkgene[1])){
                $matches   += $vkgene[1];
                $mutations += $vkgene[2];
              }
              if(defined($jgene[1])){
                $matches   += $jgene[1];
                $mutations += $jgene[2];
              }elsif(defined($jkgene[1])){
                $matches   += $jkgene[1];
                $mutations += $jkgene[2];
              }
              #print $vgene[0] . "\t" . $dgene[0] . "\t" . $jgene[0] . "\t" . $this_h3 . "\t" . "\n";
            }
          }
        }
      }
    }
  }

  $all_accs=~s/,$//;

  # calculate average shm
  my $pid="";
  if($matches>0){
    $pid= (int(1000 * (($matches - $mutations) /$matches)))/10;
  }

  # get most common cdr1
  my @sorted_cdr1 = sort { $cdr1{$b} cmp $cdr1{$a} } keys %cdr1;
  if(defined($sorted_cdr1[0])){
    $dominant_h1 = $sorted_cdr1[0];
    if($dominant_h1 eq ""){
      if(defined($sorted_cdr1[1])){
        $dominant_h1 = $sorted_cdr1[1];
      }
    }
  }
  #print "Top: " . scalar(@sorted_cdr1) . " " . $sorted_cdr1[0] . " " . $sorted_cdr1[(scalar(@sorted_cdr1)-1)] . "\n";
  # get most common cdr2
  my @sorted_cdr2 = sort { $cdr2{$b} cmp $cdr2{$a} } keys %cdr2;
  if(defined($sorted_cdr2[0])){
    $dominant_h2 = $sorted_cdr2[0];
    if($dominant_h2 eq ""){
      if(defined($sorted_cdr2[1])){
        $dominant_h2 = $sorted_cdr2[1];
      }
    }
  }
  # get most common d-segment
  my @cdr_dgene = sort { $cdr_dgene{$b} cmp $cdr_dgene{$a} } keys %cdr_dgene;
  if(defined($cdr_dgene[0])){
    $dominant_dseg=$cdr_dgene[0];
    if($dominant_dseg eq ""){
      if(defined($cdr_dgene[1])){
        $dominant_dseg = $cdr_dgene[1];
      }
    }
  }

  #print "($dominant_h1,$dominant_h2,$pid,$dominant_dseg)\n";
  return($dominant_h1,$dominant_h2,$pid,$dominant_dseg,$all_accs);
}

sub getCloneMode {
  my($self,$this_v,$this_j,$this_h3)=@_;
  my ($dominant_h1,$dominant_h2,$avg_shm,$dominant_dseg)=("","","","");

  my %cdr_dgene = (); 
  my %cdr1      = ();
  my %cdr2      = ();
  my $matches   = 0;
  my $mutations = 0;
  my $all_accs  = "";

  for(my $s=0;$s<$self->{nseqs};$s++){
    my @vgene=split(/  */,$self->getHeaderField($s,1));
    if(scalar(@vgene)==3){
      if($vgene[0] eq $this_v){
        my @jgene=split(/  */,$self->getHeaderField($s,3));
        if(scalar(@jgene)==3){
          if($jgene[0] eq "$this_j"){
            my $h3=$self->getHeaderField($s,4);
            my $l3=$self->getHeaderField($s,5);
            if(($h3 eq $this_h3) || ($l3 eq $this_h3)){
              $all_accs .= $self->getAccession($s) . ",";
              # d-segs
	      my @dgene=split(/  */,$self->getHeaderField($s,2));
              if(!defined($dgene[0])){
                $dgene[0]="";
              }
              if(defined($cdr_dgene{$dgene[0]})){
                $cdr_dgene{$dgene[0]}++;
              }else{
                $cdr_dgene{$dgene[0]}=1;
              }

              # cdr1 8,10
                my $h1 = $self->getHeaderField($s,8);
                my $l1 = $self->getHeaderField($s,10);
                my $this_cdr1 = "";
                if(length($h1)>0){
                  $this_cdr1=$h1;
                }elsif(length($l1)>0){
                  $this_cdr1=$l1;
                }
                if(defined($cdr1{$this_cdr1})){
                  $cdr1{$this_cdr1}++;
                }else{
                  $cdr1{$this_cdr1}=1;
                }

              # cdr2 9,11
                my $cdr2 = "";
                my $h2 = $self->getHeaderField($s,9);
                my $l2 = $self->getHeaderField($s,11);
                my $this_cdr2="";
                if(length($h2)>0){
                  $this_cdr2=$h2;
                }elsif(length($l2)>0){
                  $this_cdr2=$l2;
                }
                if(defined($cdr1{$this_cdr1})){
                  $cdr2{$this_cdr2}++;
                }else{
                  $cdr2{$this_cdr2}=1;
                }

              # SHM
              if($vgene[1]>200){
                $matches   += $vgene[1];
                $mutations += $vgene[2];
              }
              if($jgene[1]>30){
                $matches   += $jgene[1];
                $mutations += $jgene[2];
              }
              #print $vgene[0] . "\t" . $dgene[0] . "\t" . $jgene[0] . "\t" . $this_h3 . "\t" . "\n";
            }
          }
        }
      }
    }
  }

  $all_accs=~s/,$//;

  # calculate average shm
  my $pid="";
  if($matches>0){
    $pid= (int(1000 * (($matches - $mutations) /$matches)))/10;
  }

  # get most common cdr1
  my @sorted_cdr1 = sort { $cdr1{$b} cmp $cdr1{$a} } keys %cdr1;
  $dominant_h1 = $sorted_cdr1[0];
  if($dominant_h1 eq ""){
    if(defined($sorted_cdr1[1])){
      $dominant_h1 = $sorted_cdr1[1];
    }
  }
  #print "Top: " . scalar(@sorted_cdr1) . " " . $sorted_cdr1[0] . " " . $sorted_cdr1[(scalar(@sorted_cdr1)-1)] . "\n";
  # get most common cdr2
  my @sorted_cdr2 = sort { $cdr2{$b} cmp $cdr2{$a} } keys %cdr2;
  $dominant_h2 = $sorted_cdr2[0];
  if($dominant_h2 eq ""){
    if(defined($sorted_cdr2[1])){
      $dominant_h2 = $sorted_cdr2[1];
    }
  }
  # get most common d-segment
  my @cdr_dgene = sort { $cdr_dgene{$b} cmp $cdr_dgene{$a} } keys %cdr_dgene;
  $dominant_dseg=$cdr_dgene[0];
  if($dominant_dseg eq ""){
    if(defined($cdr_dgene[1])){
      $dominant_dseg = $cdr_dgene[1];
    }
  }

  #print "($dominant_h1,$dominant_h2,$pid,$dominant_dseg)\n";
  return($dominant_h1,$dominant_h2,$pid,$dominant_dseg,$all_accs);
}

sub getNRCDRLengthDistribution {
  my($self,$header_field)=@_;
  
  my %cdrs=();

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $cdr=$self->getHeaderField($s,$header_field);
    if(length($cdr)>0){
      $cdrs{$cdr}=1;
    }
  }
 
  # ok, now count unique items by length
  my %cdr_lengths = ();

  my @unique_cdrs = keys %cdrs;
  for(my $u=0;$u<scalar(@unique_cdrs);$u++){
    my $length=length($unique_cdrs[$u]);
    if(defined($cdr_lengths{$length})){
      $cdr_lengths{$length}++;    
    }else{
      $cdr_lengths{$length}=1;
    }     
  }
  return %cdr_lengths;
}

sub getCDRLengthDistribution {
  my($self,$header_field)=@_;

  my %already_seen=();
  my %cdr_lengths = ();
  my $total=0;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $v=$self->getHeaderField($s,1);
       $v=~s/ .*//;
    my $j=$self->getHeaderField($s,3);
       $j=~s/ .*//;
    my $h3l3=$self->getHeaderField($s,4);
       $h3l3.=$self->getHeaderField($s,5);
    my $cdr=$self->getHeaderField($s,$header_field);
   
    if(length($cdr)>0){
      unless(defined($already_seen{$v . $j . $h3l3})){
        $already_seen{$v . $j . $h3l3}=1;
        $total++;
        my $length=length($cdr);
        if(defined($cdr_lengths{$length})){
          $cdr_lengths{$length}++;
        }else{
          $cdr_lengths{$length}=1;
        }
      }      
    }
  }
  # now normalize
  my @unique_cdrs = keys %cdr_lengths;
  for(my $u=0;$u<scalar(@unique_cdrs);$u++){
    $cdr_lengths{$unique_cdrs[$u]} /= $total;
  }
  return %cdr_lengths;
}

sub getReadingFrameDistribution {
  my($self)=@_;

  my %frames=('-3',0,'-2',0,'-1',0,'0',0,1,0,2,0);

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $header_field1=$self->getHeaderField($s,0);
    $header_field1=~s/..*frame_//;
    if(defined($frames{$header_field1})){
      $frames{$header_field1}++;
    }
  }
  my %dataHash = ("Reverse-3",$frames{'-3'},
                  "Reverse-2",$frames{'-2'},
                  "Reverse-1",$frames{'-1'},
                  "Forward1",$frames{'0'},
                  "Forward2",$frames{'1'},
                  "Forward3",$frames{'2'});
  return %dataHash;
}

sub getOnTargetCountsHash {
  my($self)=@_;
  my %data_hash=();

  $data_hash{'total'}    = $self->{nseqs};
  $data_hash{'VH'}       = 0;
  $data_hash{'VK'}       = 0;
  $data_hash{'VL'}       = 0;
  $data_hash{'Heterodimer'} = 0;
  $data_hash{'VHVJaaH3'} = 0;
  $data_hash{'VKJKaaL3'} = 0;
  $data_hash{'VLJLaaL3'} = 0;
  my $nohetero_jh3          = 0;
  my $nohetero_jl3          = 0;
  my $nohetero_jk3          = 0;

  for(my $s=0;$s<$self->{nseqs};$s++){
    my $vgene_hit = $self->getHeaderField($s,1);  
    my $jgene_hit = $self->getHeaderField($s,3);
    my $vlgene_hit = $self->getHeaderField($s,12);  
    my $jlgene_hit = $self->getHeaderField($s,13);
    my $h3        = $self->getHeaderField($s,4);
    my $l3        = $self->getHeaderField($s,5);
    
    #for heterodimers
    if(($vlgene_hit =~ m/IGLV/ or $vlgene_hit =~ m/IGKV/) and $vgene_hit =~ m/IGHV/){
      $data_hash{'Heterodimer'} += 1;
      #count heterodimers without J or cdr aa to be subtracted later for exclusivity
      if($jgene_hit !~ m/IGHJ/ or $h3 !~ m/^C[A-WY]*W$/){
          $nohetero_jh3 += 1;
      }
      if(($jlgene_hit !~ m/IGKJ/ or $l3 !~ m/^C[A-WY]*$/) and $vlgene_hit =~ m/IGKV/){
          $nohetero_jk3 += 1;
      }
      if(($jlgene_hit !~ m/IGLJ/ or $l3 !~ m/^C[A-WY]*$/) and $vlgene_hit =~ m/IGLV/){
          $nohetero_jl3 += 1;
      }

      #do not count V*J*aa*3 if heterdimer
      
    }else{
    if($vgene_hit =~ m/IGHV/){
      $data_hash{'VH'}       += 1;
      if($jgene_hit=~m/IGHJ/){
        if($h3=~m/^C[A-WY]*W$/){
          $data_hash{'VHVJaaH3'} += 1;
        }
      }
    }

    if($vlgene_hit =~ m/IGKV/){
      $data_hash{'VK'}       += 1;
      if($jlgene_hit=~m/IGKJ/){
        if($l3=~m/^C[A-WY]*$/){
          $data_hash{'VKJKaaL3'} += 1;
        }
      }
    }

    if($vlgene_hit =~ m/IGLV/){
      $data_hash{'VL'}       += 1;
      if($jlgene_hit=~m/IGLJ/){
        if($l3=~m/^C[A-WY]*$/){
          $data_hash{'VLJLaaL3'} += 1;
        }
      }
    }
    }

  }

  $data_hash{'total'}     = $self->{nseqs};
  $data_hash{'offtarget'} = $self->{nseqs} - ($data_hash{'VH'} + $data_hash{'VK'} + $data_hash{'VL'} + $data_hash{'Heterodimer'});

  $data_hash{'VH'}       -= $data_hash{'VHVJaaH3'};
  $data_hash{'VK'}       -= $data_hash{'VKJKaaL3'};
  $data_hash{'VL'}       -= $data_hash{'VLJLaaL3'};
#  $data_hash{'VH'}       -= $nohetero_jh3;
#  $data_hash{'VK'}       -= $nohetero_jk3;
#  $data_hash{'VL'}       -= $nohetero_jl3;
 
  return %data_hash;
}

sub getOnTargetHash {
  my($self)=@_;

  my %data_hash = $self->getOnTargetCountsHash();

  $data_hash{'VH'}       /= $self->{nseqs};
  $data_hash{'VK'}       /= $self->{nseqs};
  $data_hash{'VL'}       /= $self->{nseqs};
  $data_hash{'VHVJaaH3'} /= $self->{nseqs};
  $data_hash{'VKJKaaL3'} /= $self->{nseqs};
  $data_hash{'VLJLaaL3'} /= $self->{nseqs};
 
  $data_hash{'Heterodimer'} /= $self->{nseqs}; #1 - ( $data_hash{'VH'} + $data_hash{'VK'} + $data_hash{'VL'} );
  $data_hash{'offtarget'} /= $self->{nseqs}; #1 - ( $data_hash{'VH'} + $data_hash{'VK'} + $data_hash{'VL'} );

  # offtarget VH VK VL VHVJaaH3 VKJKaaL3 VLJLaaL3 
  return %data_hash; 
}

# we'll use the class destructor to clean up any files
# created by tied arrays
sub DESTROY {
  my $self = shift;
  if ($self->{headerdb}) {
    undef $self->{headers};
    unlink $self->{headerdb};
  }
  if ($self->{sequencedb}) {
    undef $self->{sequence};
    unlink $self->{sequencedb};
  }
  unlink "DB_CONFIG" if -e "DB_CONFIG";
}

1;
