package VirtualExpressionSytem;
use AminoAcidRelativeAdaptiveness;
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
  $self->{aminoacids}          	= {};
  $self->{cutoff}              	=  0;
  $self->{aaCodes}              = {};

  %{$self->{aaCodes}} 		= (    	'Ala','A','Arg','R','Asn','N','Asp','D','Asx','B','Cys','C','Glu','E',
                        		'Gln','Q','Glx','Z','Gly','G','His','H','Ile','I','Leu','L','Lys','K','Met','M',
                        		'Phe','F','Pro','P','Ser','S','Thr','T','Try','W','Trp','W','Tyr','Y','Val','V','End','.');

  bless $self,'VirtualExpressionSytem';
  return $self;
}

# Methods --------------------------------------------------

sub getExpressionSystemFile {
  my($self,$host)=@_;
  my $file="";
  if($host eq "human"){
    $file="$rootdir/db/codon-usage-homo.sapiens.txt";
  }elsif($host eq "ecoli"){
    $file="$rootdir/db/codon-usage-e.coli.txt";
  }elsif($host eq "yeast"){
    $file="$rootdir/codon-usage-s.cerevisiae.txt";
  }elsif(-f $host){
    $file=$host;
  }else{
    print "Unknown host!\n";
  }
  return $file;
}

sub loadCodonUsageBiases{
  my ($self,$expression_system_file) = @_;
  unless(-f $expression_system_file){
    print "Error! Expression system file $expression_system_file not found\n";
    exit;
  }
  # open expression system file
  open FILE,$expression_system_file;
  my @lines=<FILE>;
  close FILE;
  chomp(@lines);

  for(my $x=1;$x<scalar(@lines);$x++){
    if(!($lines[$x]=~m/^ *$/)){
      my ($aaTripleCode,$dnaCodon,$counts)=split(/  */,$lines[$x]); # monkey
      my $aaChar=${$self->{aaCodes}}{$aaTripleCode};
      if(!defined(${$self->{aminoacids}}{$aaChar})){
        ${$self->{aminoacids}}{$aaChar} = AminoAcidRelativeAdaptiveness->new();
        ${$self->{aminoacids}}{$aaChar}->setName($aaTripleCode);
        ${$self->{aminoacids}}{$aaChar}->setChar($aaChar);
        ${$self->{aminoacids}}{$aaChar}->setThreeLetter($aaTripleCode);
        ${$self->{aminoacids}}{$aaChar}->setCutoff($self->{cutoff});
      }
      ${$self->{aminoacids}}{$aaChar}->addCodonScore($dnaCodon,$counts);
    }
  }
}

sub stochasticCodonSample {
  my($self,$aaChar)=@_;
  if($aaChar eq "X"){
    return "NNN";
  }elsif($aaChar eq "Z"){
    return "NNN";
  }else{
    return ${$self->{aminoacids}}{$aaChar}->stochasticCodonSample();
  }
}

sub translateCodon {
  my($self,$codon)=@_;
if ( $codon =~ /^TC/i ) { return 'S' } # Serine
elsif ( $codon =~ /TT[CT]/i ) { return 'F' } # Phenylalanine
elsif ( $codon =~ /TT[AG]/i ) { return 'L' } # Leucine
elsif ( $codon =~ /TA[CT]/i ) { return 'Y' } # Tyrosine
elsif ( $codon =~ /TA[AG]/i ) { return 'X' } # Stop
elsif ( $codon =~ /TG[CT]/i ) { return 'C' } # Cysteine
elsif ( $codon =~ /TGA/i ) { return 'X' } # Stop
elsif ( $codon =~ /TGG/i ) { return 'W' } # Tryptophan
elsif ( $codon =~ /^CT/i ) { return 'L' } # Leucine
elsif ( $codon =~ /^CC/i ) { return 'P' } # Proline
elsif ( $codon =~ /CA[CT]/i ) { return 'H' } # Histidine
elsif ( $codon =~ /CA[AG]/i ) { return 'Q' } # Glutamine
elsif ( $codon =~ /^CG/i ) { return 'R' } # Arginine
elsif ( $codon =~ /AT[ACT]/i ) { return 'I' } # Isoleucine
elsif ( $codon =~ /ATG/i ) { return 'M' } # Methionine
elsif ( $codon =~ /^AC/i ) { return 'T' } # Threonine
elsif ( $codon =~ /AA[CT]/i ) { return 'N' } # Asparagine
elsif ( $codon =~ /AA[AG]/i ) { return 'K' } # Lysine
elsif ( $codon =~ /AG[CT]/i ) { return 'S' } # Serine
elsif ( $codon =~ /AG[AG]/i ) { return 'R' } # Arginine
elsif ( $codon =~ /^GT/i ) { return 'V' } # Valine
elsif ( $codon =~ /^GC/i ) { return 'A' } # Alanine
elsif ( $codon =~ /GA[CT]/i ) { return 'D' } # Aspartic Acid
elsif ( $codon =~ /GA[AG]/i ) { return 'E' } # Glutamic Acid
elsif ( $codon =~ /^GG/i ) { return 'G' } # Glycine
elsif ( $codon =~ m/N/ ) { return 'X' }   # N - unknown character
else { return "X"; } # bad codon
}

sub getRelativeAdaptiveness {
  my($self,$codon)=@_;
  my $aa = $self->translateCodon($codon);
  if($aa eq "X"){
    return 0;
  }else{
    return ${$self->{aminoacids}}{$aa}->getFractionOfBestScore($codon);
  }
}

1;
