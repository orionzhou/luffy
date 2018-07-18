package Snp;
use strict; 
use Path::Class;
use Data::Dumper;
use Common; 
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw//;
@EXPORT_OK = qw//;

sub ssp_read {
  my ($fi) = @_;
  open(FH, "<$fi");
  my ($nInd, $nPos) = (0, 0);
  my (@inds, @poss, @refs, @snps);
  my $i = 1;
  while( <FH> ) {
    chomp;
    next unless $_;
    my @ps = split(/\W+/);
    if($i == 1) {
      ($nInd, $nPos) = @ps;
    } elsif($i == 2) {
      die "line 2 not $nPos poss: ".@ps."\n" unless @ps == $nPos;
      @poss = @ps;
    } elsif($i == 3) {
      @ps = @ps[1..$#ps] if $ps[0] eq "anc";
      die "line 3 not $nPos poss: ".@ps."\n" unless @ps == $nPos;
      @refs = @ps;
    } else {
      die "line $i not $nPos poss: ".@ps."\n" unless @ps == $nPos + 1;
      push @inds, $ps[0];
      for my $j (1..$#ps) {
        my $allele = $ps[$j];
        die "\t[$i, $j]: $allele [invalid]\n" if $allele =~ /[^ATCGURYMKWSBDHVN]/i;
      }
      push @snps, [@ps[1..$#ps]];
    }
    $i ++;
  }
  close FH;
  die "not $nInd inds: ".@snps."\n" unless $nInd == @snps;
  printf "  %d snps in %d inds\n", $nPos, $nInd;
  my $h = {nInd=>$nInd, nPos=>$nPos, inds=>\@inds, poss=>\@poss, refs=>\@refs, snps=>\@snps};
  return $h;
}
sub ssp_write {
  my ($h, $fo) = @_;
  open(FH, ">$fo");
  printf FH "%d  %d\n", $h->{nInd}, $h->{nPos};
  print FH join("\t", @{$h->{poss}})."\t\n";
  print FH join("\t", @{$h->{refs}})."\t\n";
  for my $i (0..$h->{nInd}-1) {
    my $snps = $h->{snps}->[$i];
    print FH join("\t", $h->{inds}->[$i], @$snps)."\t\n";
  }
  close FH;
}
sub ssp_sample {
  my ($dirI, $dirO, $chrs) = @_;
  system("mkdir -p $dirO") unless -d $dirO;
  my $ssp;
  
  my ($n, $m) = (10000, 1250);
  my @idxs = sample_serial($n, $m);
  for my $chr (@$chrs) {
    my $f01 = "$dirI/$chr.ssp";
    my $h = ssp_read($f01);
    my $chrN = substr($chr, length($chr)-1, 1);
    $ssp->{nInd} ||= $h->{nInd};
    $ssp->{inds} ||= $h->{inds};
    $ssp->{nPos} ||= 0;
    $ssp->{poss} ||= [];
    $ssp->{refs} ||= [];
    my @poss = map {$h->{poss}->[$_]} @idxs;
    my @refs = map {$h->{refs}->[$_]} @idxs;
    $ssp->{nPos} += @poss;
    push @{$ssp->{poss}}, map {$chrN*100_000_000 + $_} @poss;
    push @{$ssp->{refs}}, @refs;
    for my $i (0..$ssp->{nInd}-1) {
      $ssp->{snps}->[$i] ||= [];
      my @snps = map {$h->{snps}->[$i]->[$_]} @idxs;
      push @{$ssp->{snps}->[$i]}, @snps;
    }
  }
  my $f_ssp = "$dirO/01.ssp";
  ssp_write($ssp, $f_ssp);
  print "sspConvert -i $dirO/01.ssp -o $dirO/02.structure -f structure\n";
  print "sspConvert -i $dirO/01.ssp -o $dirO/03.phylip -f phylip\n";
}

1;
__END__
