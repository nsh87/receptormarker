package AminoAcidRelativeAdaptiveness;
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
  $self->{aaName} 	    = undef;	# a single amino acid
  $self->{aaThreeLetter}    = undef;	# the three letter to represent this amino acid
  $self->{aaChar} 	    = undef;	# the one letter representation
  $self->{codons} 	    = [];	# all codons that represent this amino acid
  $self->{codonScores} 	    = [];	# all codon scores for this amino acid
  $self->{pickRecord} 	    = [];	# number of times codon was selected	
  $self->{cutoff} 	    = 0; 	# minimum frequency for selection
  $self->{scoresAreCounts}  = 1;	# yes
  $self->{bestScore} 	    = 0;	# the highest frequency codon for this aa. Used to establish relative adaptiveness
  bless $self,'AminoAcidRelativeAdaptiveness';
  return $self;
}

# Methods --------------------------------------------------

sub setName {
  my ($self, $name) = @_;
  $self->{aaName} = $name if defined($name);
  return 1;
}

sub setThreeLetter {
  my ($self, $aaThreeLetter) = @_;
  $self->{aaThreeLetter} = $aaThreeLetter if defined($aaThreeLetter);
  return 1;
}

sub setChar {
  my ($self, $aaChar) = @_;
  $self->{aaChar} = $aaChar if defined($aaChar);
  return 1;
}

sub setCutoff {
  my ($self, $aCutoff) = @_;
  if(defined($aCutoff)){
    if($aCutoff<100){
      $self->{cutoff} = $aCutoff if defined($aCutoff);
    }
    return 1;
  } else {
    return 0;
  }
}

sub getCodons {
  my $self = shift;
  my @list=();
  for(my $x=0;$x<scalar(@{$self->{codons}});$x++){
    push @list,${$self->{codons}}[$x];
  }
  return @list;
}

sub getScore {
  my ($self,$codon) = @_;
  for(my $x=0;$x<scalar(@{$self->{codons}});$x++){
    if(${$self->{codons}}[$x] eq $codon){
      return ${$self->{codonScores}}[$x];
    }
  }
  print "Error: no score available: codon not found\n";
  return 0;
}

sub getFractionOfBestScore {
  my ($self,$codon) = @_;
  for(my $x=0;$x<scalar(@{$self->{codons}});$x++){
    if(${$self->{codons}}[$x] eq $codon){
      return int((100 * ${$self->{codonScores}}[$x]) / $self->{bestScore});
    }
  }
}

sub addCodon {
  my ($self,$aCodon) = @_;
  if(!$self->hasCodon($aCodon)){
    push(@{ $self->{codons} }, $aCodon);
    push(@{ $self->{codonScores} }, 0);
    push(@{ $self->{pickRecord} }, 0);
  }
  return @{ $self->{codons} };
}

sub addCodonScore {
  my ($self,$aCodon,$aScore) = @_;
  if($self->hasCodon($aCodon)){
    $self->resetCodonScore($aCodon,$aScore);
  } else {
    push(@{ $self->{codons} }, $aCodon);
    push(@{ $self->{codonScores} }, $aScore);
    push(@{ $self->{pickRecord} }, 0);
    if($aScore > $self->{bestScore}){
      $self->{bestScore}=$aScore;
    }
  }
  return @{ $self->{codonScores} };
}

sub resetCodonScore {
  my ($self,$aCodon,$aScore) = @_;
  if($self->hasCodon($aCodon)){
    for(my $x=0;$x<scalar(@{$self->{codons}});$x++){
      if(${$self->{codons}}[$x] eq $aCodon){
        ${$self->{codonScores}}[$x]=$aScore;
        if($aScore > $self->{bestScore}){
          $self->{bestScore}=$aScore;
        }
      }
    }
  }
  return 0; 
}

sub hasCodon {
  my ($self,$aCodon) = @_;
  for(my $x=0;$x<scalar(@{$self->{codons}});$x++){
    if(${$self->{codons}}[$x] eq $aCodon){
      return 1;
    }
  }
  return 0;
}

sub printCodons {
  my $self = shift;
  for(my $x=0;$x<scalar(@{$self->{codons}});$x++){
    print ${$self->{codons}}[$x] . " " . ${$self->{codonScores}}[$x] .  "\n";
  }
}

sub codonHash {
  my $self = shift;
  my %encode = ();
  for(my $x=0;$x<scalar(@{$self->{codons}});$x++){
    $encode{ ${$self->{codons}}[$x] } = $self->{aaChar};
  }
  return %encode;
}

sub print {
  my $self = shift;
  print $self->{aaName} . ": " . $self->{aaThreeLetter} . " " . $self->{aaChar} . "\n";
  $self->printCodons();
}

sub normalize {
  my $self = shift;
  my $total=0;
  for(my $x=0;$x<scalar(@{$self->{codonScores}});$x++){
    $total+=${$self->{codonScores}}[$x];
  }
  if($total==0){
    print "For some reason " . $self->{aaChar} . "has no scores\n";
    exit;
  }
  $self->{bestScore}=0;
  for(my $x=0;$x<scalar(@{$self->{codonScores}});$x++){
    ${$self->{codonScores}}[$x]=int((100 * ${$self->{codonScores}}[$x])/$total);
    if( $self->{bestScore} < ${$self->{codonScores}}[$x] ) {
      $self->{bestScore}=${$self->{codonScores}}[$x];
    }
  }
  $self->{scoresAreCounts}=0;
}

sub stochasticCodonSample {
  my ($self,$attempts) = @_;
  my $total_reached=0;
  my $roll = int(rand(100));

  $attempts++;
  if($attempts>100){
    return "Error: codon could never be selected...\n";
    exit;
  }
  if($self->{scoresAreCounts}){
    $self->normalize();
  }
  for(my $x=0;$x<scalar(@{$self->{codonScores}});$x++){
    $total_reached+=${$self->{codonScores}}[$x];
    if($total_reached>=$roll){
      if(${$self->{codonScores}}[$x] < $self->{cutoff}){
        $self->stochasticCodonSample($attempts);
      } else {
        ${$self->{pickRecord}}[$x]++;
        return ${$self->{codons}}[$x];
      }
    }
  }
  $self->stochasticCodonSample($attempts);
}

1;
