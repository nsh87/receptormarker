package MSA;
use strict;

# Author:  Jacob Glanville 
# Contact: jake@distributedbio.com
#
# Glanville J, Zhai W, Berka J et al. Precise determination of the diversity 
# of a combinatorial antibody library gives insight into the human immunoglobulin 
# repertoire. Proc Natl Acad Sci USA. 2009;106:20216–20221

# Constructor ----------------------------------------------
sub new {
  my ($class) = @_;
  my $self = {};
  $self->{filename}     = "";
  $self->{headers}      = [];
  $self->{sequence}     = [];
  $self->{seqnames}     = {};
  $self->{nseqs}        =  0;
  $self->{cols}         =  0;
  $self->{msa}          = [];

  bless $self,'MSA';
  return $self;
}

# Methods --------------------------------------------------

sub getMSALength {
  my($self)=@_;
  return $self->{cols};
}

sub getChar {
  my($self,$seq,$char)=@_;
  return ${$self->{msa}}[$seq][$char];
}

sub getCharFromSeqName {
  my($self,$seqname,$position)=@_;
  if(defined(${$self->{seqnames}}{$seqname})){
    if(defined(${$self->{msa}}[${$self->{seqnames}}{$seqname}][$position])){
      return ${$self->{msa}}[${$self->{seqnames}}{$seqname}][$position];
    }
  }
  print "Error! $seqname $position not found!\n";
  exit;
}

sub getSeqCount {
  my($self)=@_;
  return $self->{nseqs};
}

sub getHeader {
  my($self,$s)=@_;
  return ${$self->{headers}}[$s]; 
}

sub loadMSA {
  my ($self,$file)=@_;
  if(-e $file){
    $self->{filename}=$file;
  } else {
    print "File does not exist.\n";
    return 0;
  }

  # load msa into filename, headers and sequence
  my $x=-1;
  open MSA, $self->{filename};
  while (<MSA>) {
    chomp;
    if (/^>/) {
      push @{$self->{headers}},$_;
      $x++;
    } elsif ( /^$/ ) {
    } else {
      ${$self->{sequence}}[$x] .= $_;
    }
  }
  $x++;
  $self->{nseqs}=$x;
  $self->{cols}=length(${$self->{sequence}}[0]);
  # extract residues into msa
  for(my $i=0;$i<$self->{nseqs};$i++){
    my @residues = split(//,${$self->{sequence}}[$i]);
    for(my $j=0;$j<$self->{cols};$j++){
      ${$self->{msa}}[$i][$j] = $residues[$j];
    }
  }
  # set header hash
  for(my $id=0;$id<scalar(@{$self->{headers}});$id++){
    my $this_header = ${$self->{headers}}[$id];
       $this_header =~ s/^>//;
    ${$self->{seqnames}}{$this_header}=$id;
  }
}

sub printMSA{
  my($self)=@_;
  $self->printMSArange(0,$self->{cols});
}

sub printMSArange{
  my($self,$start,$stop)=@_;
  if($stop>$self->{cols}){
    $stop=$self->{cols};
  }

  for(my $x=0;$x<$self->{nseqs};$x++){
    print ${$self->{headers}}[$x] . "\n";
    for(my $j=$start;$j<=$stop;$j++){
      print ${$self->{msa}}[$x][$j];
    }
    print "\n";
  }
}

sub printConservedPeptides {
  my ($self,$window_size)=@_;
  my $seq=0;
  for(my $framestart=0;($framestart < ($self->{cols} - $window_size));$framestart++){
    my $pid=$self->getPIDRange($framestart,($framestart + $window_size));
    
    print "$pid\t";

    my $window_viewed=0;
    for(my $c=$framestart;(($c<$self->{cols})&&($window_viewed<$window_size));$c++){
      print ${$self->{msa}}[0][$c];
      $window_viewed++;
    }
    print " " . $framestart . "-" . ($framestart + $window_size) . "\n";
    print "\n";
  }
}


sub printHMMqc {
  my($self,$window_size,$min_match_states,$max_gap_states)=@_;
  # for each sequence
  for(my $x=0;$x<$self->{nseqs};$x++){
    my $seq_unverified=1;
    # for each position in the sequence, evalue the thirty frame window of residues (some are skipped: see below)
    for(my $framestart=0;(($framestart<$self->{cols}) && $seq_unverified);$framestart++){ 
      # if this position is a character, evaluate it. Ignore gappy positions.
      if(${$self->{msa}}[$x][$framestart]=~m/[A-Za-z]/){
        my $window_viewed=0;
        my $match_states=0;
        my $gap_states=0;
        #scan forward from this point until you see thirty residues or reach the end of the sequence
        for(my $c=$framestart;(($c<$self->{cols})&&($window_viewed<$window_size));$c++){
          if(${$self->{msa}}[$x][$c]=~m/[A-Z]/){
            $window_viewed++;
            $match_states++;
          }elsif(${$self->{msa}}[$x][$c]=~m/[a-z]/){
            $window_viewed++;
          }elsif(${$self->{msa}}[$x][$c]=~m/\-/){
            $gap_states++;
          }
        }
        if($window_viewed==$window_size){
          if($match_states>=$min_match_states){
            if($gap_states<=$max_gap_states){
              $seq_unverified=0;
            }
          }else{
            # move forward a number of residues if window is mostly inserts
            my $shift=$min_match_states - $match_states;
            $shift--;
            if($shift>0){
              $framestart+=$shift;
            }
          }
        }else{
          # not enough window left to evaluation: finish.
          $framestart=$self->{cols};
        }
      }
    }
    if($seq_unverified==0){
      print ${$self->{headers}}[$x] . "\n";
      print ${$self->{sequence}}[$x] . "\n";
    }
  }
}

sub isColIdentical {
  my($self,$col)=@_;
  my $id=1;
  for(my $x=0;$x<$self->{nseqs};$x++){
    for(my $y=($x+1);$y<$self->{nseqs};$y++){
      if(${$self->{msa}}[$x][$col] ne ${$self->{msa}}[$y][$col]){
        $id=0;
        return $id;
      }
    }
  }
  return $id;
}

sub getPIDRange {
  my($self,$start,$stop)=@_;
  my $ids=0;
  my $observations=0;
  for(my $j=$start;$j<$stop;$j++){
    $ids+=$self->isColIdentical($j);
    $observations++;
  }
  my $pid=0;
  if($observations>0){
    $pid=((int(($ids * 1000)/$observations))/10);
  }
  return $pid;
}

sub getPID {
  my($self)=@_;
  return $self->getPIDRange($self,0,$self->{cols});
}

sub printIDProfile{
  my($self)=@_;
  for(my $j=0;$j<$self->{cols};$j++){
    print $self->isColIdentical($j);
  }
  print "\n";
}

sub crop2model {
  my($self)=@_;

  # find boundaries
  my $hmmstart=0;
  my $hmmstop=$self->{cols};

  for(my $j=0;$j<$self->{cols};$j++){
    if((${$self->{msa}}[0][$j] eq "-")||(${$self->{msa}}[0][$j]=~m/[A-Z]/)){
      $hmmstart=$j;
      $j=$self->{cols}; # break
    }
  }
  for(my $j=($self->{cols} - 1);$j>=0;$j--){
    if((${$self->{msa}}[0][$j] eq "-")||(${$self->{msa}}[0][$j]=~m/[A-Z]/)){
      $hmmstop=$j;
      $j=-1; # break
    }
  }

  # print result
  $self->printMSArange($hmmstart,$hmmstop);
}


sub printConsensus {
  my ($self)=@_;
  my $consensus=$self->getConsensus();
  print ">consensus\n";
  print "$consensus\n";
}

sub getConsensus {
  my($self)=@_;
  # Output
  my $quality="";
  my $nextqual="";
  my $consensus_seq="";

  for(my $col=0;$col<$self->{cols};$col++){
    my %counts=();
    for(my $r=0;$r<$self->{nseqs};$r++){
      if(!exists($counts{${$self->{msa}}[$r][$col]})){
        $counts{${$self->{msa}}[$r][$col]}=1;
      } else {
        $counts{${$self->{msa}}[$r][$col]}++;
      }
    }
    # return highest
    my $highest_number=0;
    my $highest_name=0;
    while (my ($key, $value) = each(%counts)){
      #print $key . " " . $value . " \n";
      if($value>$highest_number){

        $highest_number=$value;
        $highest_name=$key;
      }
    }
    #print "highest $highest_number " . $self->{nseqs} . "\n";
    my $score = int(($highest_number*10)/($self->{nseqs}+1));
    if($highest_name ne "-"){
      #print " old score $score ";
      if(defined($counts{"-"})){
        #print "int(($highest_number*10)/( " . ($self->{nseqs}+1) . " - " . $counts{"-"} . " )) . \n";
        $score = int(($highest_number*10)/( ($self->{nseqs}+1) - $counts{"-"} ));
      }
    }
    # only add non-gap sequences
    unless($highest_name eq "-"){
      if($score>7){
        $quality .= ".";
      }else{
        $quality .= $score;
      }
      $consensus_seq.=$highest_name;
    }
  }
  #print "\n";
  print ">Quality:\n" . $quality . "\n";
  #print "Consensus:\t" . $consensus_seq . "\n";
  return $consensus_seq;
}


sub getAverageLength {
  my($self)=@_;
  my $total_length=0;

  for(my $l=0;$l<$self->{nseqs};$l++){
    my $residues=${$self->{sequence}}[$l];
    $residues=~s/\-*//g;
    $residues=~s/\.*//g;
    $total_length+=length($residues);
  }
  my $length=int($total_length/$self->{nseqs});
  return $length;
}


1;
