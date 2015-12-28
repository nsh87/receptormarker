package Lineages;
use VDJFasta;
use strict;
use vars qw($rootdir);

# Author:  Jacob Glanville 

BEGIN {
  use File::Basename;
  use Cwd 'abs_path';
  $rootdir = dirname(dirname(abs_path($0)));
};

# Constructor ----------------------------------------------
sub new {
  my ($class) = @_;
  my $self = {};

  # raw features
  $self->{vdjfasta}     =    {};
  $self->{accession}    =    [];
  $self->{counts}       =    [];
  $self->{reads}        =    [];

  # processed features


  bless $self,'Lineages';
  return $self;
}

# Methods --------------------------------------------------

sub loadVDJFasta {
  my($self,$vdjfasta_file)=@_;

  $self->{vdjfasta}=VDJFasta->new();
  $self->{vdjfasta}->loadSeqs($vdjfasta_file);
}

sub loadLineages {
  my($self,$file)=@_;
  open(FILE,$file);
  my @lines=<FILE>;
  close(FILE);
  chomp(@lines);

  for(my $x=0;$x<scalar(@lines);$x++){
    if($lines[$x]=~m/^[0-9]/){
      my($counts,$accession,$reads)=split(/\t/,$lines[$x]);
      push @{$self->{counts}},$counts;
      push @{$self->{accession}},$accession;
      push @{$self->{reads}},$reads;   

      # Boyd format: process all sequences in that lineage
      # $self->getLineageFeatures($x); 
      # (${$self->{counts}}[$lineage],${$self->{accession}}[$lineage],${$self->{reads}}[$lineage])=split(/\t/,$lines[$x]);
      # my @acc_list = split(/ /,$read_accessions);
      #$lineage++;
    }
  }
}

sub getLineageCount {
  my($self)=@_;
  return scalar(@{$self->{accession}});
}

sub getCount {
  my($self,$x)=@_;
  return ${$self->{counts}}[$x];
}

sub getAccession {
  my($self,$x)=@_;
  return ${$self->{accession}}[$x];
}

sub getReads {
  my($self,$x)=@_;
  return ${$self->{reads}}[$x];
}

sub printLineages {
  my($self)=@_;
  
  my %acc2seqid=();
  $self->{vdjfasta}->getAccession2seqidHash(\%acc2seqid);
  my @keys=keys %acc2seqid;
  #for(my $k=0;$k<scalar(@keys);$k++){
  #  print $keys[$k] . "\t" . $acc2seqid{$keys[$k]} . "\n";
  #}
  for(my $x=0;$x<scalar(@{$self->{accession}});$x++){
    my ($vgene,$jgene,$h3)      = split(/_/,${$self->{accession}}[$x]);
    print $vgene . "\t" . $jgene . "\t" . $h3 . "\n";
    my @reads = split(/ /,${$self->{reads}}[$x]);
    for(my $r=0;$r<scalar(@reads);$r++){
      print $reads[$r] . "\n";
      #my($header,$seq)=$self->{vdjfasta}->getSeq($acc2seqid{$reads[$r]});
      #print $header . "\n";
      #print $seq . "\n";
    }  
  }
}

sub exportAlignments {
  my($self,$max_depth)=@_;
 
  my %acc2seqid=();
  $self->{vdjfasta}->getAccession2seqidHash(\%acc2seqid);
  my @keys=keys %acc2seqid;
  # lineages
  for(my $x=0;$x<scalar(@{$self->{accession}});$x++){
    my $lineage_file="lineage-" . ${$self->{accession}}[$x] . ".fa";
       $lineage_file=~s/\//./;
    open(FILE,">$lineage_file");
    if($self->countDonors($x)>0){
      my ($vgene,$jgene,$h3)      = split(/_/,${$self->{accession}}[$x]);
      my $total_donors=$self->countDonors($x);
      my $total_seqs=$self->countTotalSeqs($x);
      #print "$total_donors\t$total_seqs\t$vgene\t$jgene\t$h3\n";
      my @reads = split(/ /,${$self->{reads}}[$x]);
      for(my $r=0;$r<scalar(@reads);$r++){
        my($header,$seq)=$self->{vdjfasta}->getSeq($acc2seqid{$reads[$r]});
        print FILE $header . "\n";
        print FILE $seq . "\n";
        if(defined($max_depth)){
          if($r>=$max_depth){
            $r=scalar(@reads); # you are done
          }
        }
      }
    }
    close(FILE);
  }
}

sub countDonors {
  my($self,$lineage_id)=@_;
  my @reads = split(/ /,${$self->{reads}}[$lineage_id]);
  my %unique=();
  for(my $x=0;$x<scalar(@reads);$x++){
    $reads[$x]=~s/.*donor-//;
    $reads[$x]=~s/[ab]-//;
    $reads[$x]=~s/4/3/;
    $reads[$x]=~s/g[di].*//;
    $unique{$reads[$x]}=1;
  }
  my @keys = keys %unique;
  return scalar(@keys);
}

sub countTotalSeqs {
  my($self,$lineage_id)=@_;
  my @reads = split(/ /,${$self->{reads}}[$lineage_id]);
  my $total=0;
  for(my $x=0;$x<scalar(@reads);$x++){
    my @fields = split(/_/,$reads[$x]);
    $total+=$fields[1];
  }
  return $total;
}


sub printEllisonLineages {
  my($self)=@_;

  # prepare tabular output for ggplot2
  # Counts YearVisit Vfam Vgene Jgene H3 Isotype Replicates
  print "Counts\tYearVisit\tTimeTransect\tVfam\tVgene\tJgene\tIsotype\tReplicates\tH3\n";

  # prepare the vdjfasta acc2seqid hash
  my %acc2seqid=();
  $self->{vdjfasta}->getAccession2seqidHash(\%acc2seqid);
  my @keys=keys %acc2seqid;
  for(my $k=0;$k<scalar(@keys);$k++){
    my $short_acc=$keys[$k];
       $short_acc=~s/.*_//;
    $acc2seqid{$short_acc}=$acc2seqid{$keys[$k]};
  }
 
  # process each lineage
  for(my $x=0;$x<scalar(@{$self->{accession}});$x++){
    my $yearvisit		= "";
    my ($vgene,$jgene,$h3)	= split(/_/,${$self->{accession}}[$x]);
    my $vfam			= $vgene;
       $vfam			=~ s/\-.*//;
    my $isotype			= "";
    my $replicates		= "";
    my $chosen_time		= 0;

    # parse all data and get year/visit/isotype/replicate info for lineage
    my @data=split(/ /,${$self->{reads}}[$x]);
    my %boyd_time_transect = (  "Y2-V1",0,"Y2-V2",0,"Y2-V3",0,
				"Y3-V1",0,"Y3-V2",0,"Y3-V3",0); 
    my %boyd_time_shm      = (  "Y2-V1",0,"Y2-V2",0,"Y2-V3",0,
                                "Y3-V1",0,"Y3-V2",0,"Y3-V3",0);
    my %boyd_time_shm_count= (  "Y2-V1",0,"Y2-V2",0,"Y2-V3",0,
                                "Y3-V1",0,"Y3-V2",0,"Y3-V3",0);
    for(my $d=0;$d<scalar(@data);$d++){
      my($this_acc,$this_bfi,$this_date,$this_num,
         $this_Group,$this_pcr,$this_year,
         $this_visit,$this_replicate,$this_isotype)=$self->boyd_cloneparse($data[$d]);
      
      # dealing with this_date vs yearvisit
      $boyd_time_transect{$this_date}=1;
      if(($yearvisit eq "") || ($this_date eq $yearvisit)){
        $yearvisit=$this_date;
      }else{
        $yearvisit="mixed";
      }

      # dealing with isotype
      if(($isotype eq "") || ($isotype eq $this_isotype)){
        $isotype=$this_isotype;
      }else{
        $isotype="mixed";
      }
      
      # dealing with replicates
      if(($replicates eq "") || ($replicates eq $this_replicate)){
        $replicates=$this_replicate;
      }else{
        $replicates="mixed";
      }

      # dealing with SHM
      
      my $this_shm = "";
      if(defined($acc2seqid{$this_acc})){
        $this_shm=$self->{vdjfasta}->getPercentSHM($acc2seqid{$this_acc},100);
        if($this_shm ne ""){
          $boyd_time_shm{$this_date}+=$this_shm;
          $boyd_time_shm_count{$this_date}++;
        }
     }

      #print $this_acc . "\t" . $acc2seqid{$this_acc} . "\t" . $self->{vdjfasta}->getPercentSHM($acc2seqid{$this_acc},100) . "\n";

    }

    # prepare tabular output for ggplot2
    # Counts YearVisit Vfam Vgene Jgene H3 Isotype Replicates

    if($isotype eq ""){
      $isotype="?";
    }

    my $boyd_time_transect_line = 	$boyd_time_transect{"Y2-V1"} . $boyd_time_transect{"Y2-V2"} .
					$boyd_time_transect{"Y2-V3"} . $boyd_time_transect{"Y3-V1"} .
					$boyd_time_transect{"Y3-V2"} . $boyd_time_transect{"Y3-V3"};
    if($boyd_time_shm_count{"Y2-V1"}>0){
      $boyd_time_shm{"Y2-V1"}/=$boyd_time_shm_count{"Y2-V1"};
    }
    if($boyd_time_shm_count{"Y2-V2"}>0){
      $boyd_time_shm{"Y2-V2"}/=$boyd_time_shm_count{"Y2-V2"};
    }
    if($boyd_time_shm_count{"Y2-V3"}>0){
      $boyd_time_shm{"Y2-V3"}/=$boyd_time_shm_count{"Y2-V3"};
    }
    if($boyd_time_shm_count{"Y3-V1"}>0){
      $boyd_time_shm{"Y3-V1"}/=$boyd_time_shm_count{"Y3-V1"};
    }
    if($boyd_time_shm_count{"Y3-V2"}>0){
      $boyd_time_shm{"Y3-V2"}/=$boyd_time_shm_count{"Y3-V2"};
    }
    if($boyd_time_shm_count{"Y3-V3"}>0){
      $boyd_time_shm{"Y3-V3"}/=$boyd_time_shm_count{"Y3-V3"};
    }
    if($replicates=~m/mixed/){
      if($boyd_time_transect_line=~m/011[01]11/){
        # print all sequences to trees/${vgene}-${jgene}-${h3}-${counts}.fa
        my $outfile="trees/$vgene-$jgene-$h3-" . ${$self->{counts}}[$x] . ".fa";
        for(my $d=0;$d<scalar(@data);$d++){
          my($this_acc,$this_bfi,$this_date,$this_num,
           $this_Group,$this_pcr,$this_year,
           $this_visit,$this_replicate,$this_isotype)=$self->boyd_cloneparse($data[$d]);
          print "Writing $this_acc to $outfile...\n";
          $self->{vdjfasta}->appendSeq($acc2seqid{$this_acc},$outfile);
        }
        #print ${$self->{counts}}[$x] . "\t" . $yearvisit . "\t" . $boyd_time_transect_line . "\t" . $vfam . "\t" . $vgene . "\t" . $jgene . "\t" . $isotype . "\t" . $replicates . "\t" . $h3 . "\t";
        #if($boyd_time_shm{"Y2-V2"} > ($boyd_time_shm{"Y3-V2"} + 0.02)){
        #  print "matured\t";
        #}elsif(($boyd_time_shm{"Y2-V2"} + 0.02) < $boyd_time_shm{"Y3-V2"}){
        #  print "regressed\t";
        #}else{
        #  print "static\t";
        #}
        #print $boyd_time_shm{"Y2-V1"} . "_" . $boyd_time_shm{"Y2-V2"} . "_" . $boyd_time_shm{"Y2-V3"} . "_"
        #. $boyd_time_shm{"Y3-V1"} . "_" . $boyd_time_shm{"Y3-V2"} . "_" . $boyd_time_shm{"Y3-V3"} . "\n";
      }
    }
  }
}

sub summarizeLineage {
  my($self)=@_;

  my %acc2seqid=();
  $self->{vdjfasta}->getAccession2seqidHash(\%acc2seqid);
  my @keys=keys %acc2seqid;
  #for(my $k=0;$k<scalar(@keys);$k++){
  #  print $keys[$k] . "\t" . $acc2seqid{$keys[$k]} . "\n";
  #}

  # process each lineage
  for(my $x=0;$x<scalar(@{$self->{accession}});$x++){
    my %isotypes=();
    my $shm_pid_max=0;
    my $shm_pid_min=1000;
    my %unique_seqs=();
    my $total_seqs=0;
    #my %timeframes=("Y2-V1",0,"Y2-V2",0,"Y2-V3",0,  # optional timeframe tracking
    #                "Y3-V1",0,"Y3-V2",0,"Y3-V3",0); # should be totally abstracted
    my %timeframes=( "00:00:00",0,"56-days",0,"365-days",0);

    my ($vgene,$jgene,$h3)      = split(/_/,${$self->{accession}}[$x]);
    my @reads = split(/ /,${$self->{reads}}[$x]);
    for(my $r=0;$r<scalar(@reads);$r++){
      my($header,$seq)=$self->{vdjfasta}->getSeq($acc2seqid{$reads[$r]});
      my @header_fields=split(/;/,$header);
      # timeframe check
      if($header=~m/00:00:00/){
        $timeframes{"00:00:00"}++;
      }
      if($header=~m/56-days/){
        $timeframes{"56-days"}++;
      }
      if($header=~m/365-days/){
        $timeframes{"365-days"}++;
      }
      #if($header=~m/Y3-[0-9][0-9][0-9]-V1/){
      #  $timeframes{"Y3-V1"}++;
      #}
      #if($header=~m/Y3-[0-9][0-9][0-9]-V2/){
      #  $timeframes{"Y3-V2"}++;
      #}
      #if($header=~m/Y3-[0-9][0-9][0-9]-V3/){
      #  $timeframes{"Y3-V3"}++;
      #}
      # isotypes
      my $iso=$header_fields[6];
         $iso=~s/ .*//;
         $isotypes{$iso}=1; 
      # shm
      my @vgene_mut = split(/ /,$header_fields[1]);
      my @jgene_mut = split(/ /,$header_fields[3]);
      my $mismatches=0;
      my $alignment=0;
      if(defined($vgene_mut[1])){
        $alignment   = $vgene_mut[1];
        $mismatches  = $vgene_mut[2];
      } 
      if(defined($jgene_mut[1])){
        $alignment  += $jgene_mut[1];
        $mismatches += $jgene_mut[2];
      }
      if($alignment>0){
        my $shm= (int(1000 * (($alignment - $mismatches)/$alignment)))/1000;
        if($shm<$shm_pid_min){
          $shm_pid_min=$shm;
        }
        if($shm>$shm_pid_max){
          $shm_pid_max=$shm;
        }
      }
      # total and unique sequences
      $unique_seqs{$seq}=1;
      $total_seqs++;
    }
    # isotype
    my $isotype_list=join('-',keys %isotypes);
    # SHM range
    my $affinity_maturation_range = $shm_pid_max - $shm_pid_min;
    # now print out the lineage summary
    print $vgene . "\t" . $jgene . "\t" . $h3 
          . "\t" . $isotype_list . "\t" . $shm_pid_min . "-" . $shm_pid_max 
          . "\t" . $affinity_maturation_range . "\t" . $total_seqs 
          . "\t" . scalar(keys %unique_seqs) 
          . "\t" . $timeframes{"00:00:00"} . "-" . $timeframes{"56-days"} . "-" . $timeframes{"365-days"}
          #. "\t" . $timeframes{"Y2-V1"} . "-" . $timeframes{"Y2-V2"} . "-" . $timeframes{"Y2-V3"}
          #. "-"  . $timeframes{"Y3-V1"} . "-" . $timeframes{"Y3-V2"} . "-" . $timeframes{"Y3-V3"}
          . "\n";
  }
}

sub boyd_cloneparse {
  my($self,$data)=@_;

  # Potential formats
  # GYAI2TT03EW5HF_BFI-0000374_Y3-004-V3_11_S138_PCR_R02
  # G9L3WKS01A69TL_BFI-0000374_Y3-004-V3_18_S138_PCR_igM_R01
      
  my ($this_acc,$this_bfi,$this_date,$this_num,$this_Group,$this_pcr,@this_isotype)=split(/_/,$data);
     $this_date=~s/-[0-9][0-9][0-9]-/-/;
  my ($this_year,$this_visit)=split(/\-/,$this_date);
  my($this_replicate,$this_isotype)=$self->boyd_replicate_isotype(@this_isotype);
  return($this_acc,$this_bfi,$this_date,$this_num,$this_Group,$this_pcr,$this_year,$this_visit,$this_replicate,$this_isotype);
}

sub boyd_replicate_isotype {
  my($self,@this_isotype)=@_;

  my $this_replicate=$this_isotype[0];
  my $this_isotype="";

  if(scalar(@this_isotype)>1){
    $this_replicate="iso"; 
    $this_isotype=$this_isotype[0];
  }
  
  return($this_replicate,$this_isotype);
}

sub writeLineages {
  my($self,$file)=@_;
  open(FILE,">$file");
  for(my $x=0;$x<scalar(@{$self->{accession}});$x++){
    print FILE ${$self->{counts}}[$x] . "\t";
    print FILE ${$self->{accession}}[$x] . "\t";
    print FILE ${$self->{reads}}[$x] . "\n";
  }
  close(FILE);
}

 
sub getGermlineFrequencies {
  my($self)=@_;

  my %vcount_hash=();
  my %jcount_hash=();
  my %h3len_hash=();

  my $counts=0;
  for(my $x=0;$x<scalar(@{$self->{accession}});$x++){
    my($vseg,$jseg,$h3)=split(/_/,${$self->{accession}}[$x]);
    $vcount_hash{$vseg}++;
    $jcount_hash{$jseg}++;
    $h3len_hash{length($h3)}++;    
  }

}

sub findLongestRepresentative {
  my($self,$cluster)=@_;
    print "Finding longest sequence in lineage $cluster\n";
}

1;
