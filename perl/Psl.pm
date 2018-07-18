package Psl;
use strict;
use Data::Dumper;
use Common;
use Location;
use Seq;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/psl2gal net2gal/;
@EXPORT_OK = qw//;

sub net2gal {
  my ($fhi, $fho) = @_;
  print $fho join("\t", qw/level id type qId qBeg qEnd qSrd qLen 
    hId hBeg hEnd hSrd hLen
    match misMatch baseN dq dh score qSeq hSeq/)."\n";
  
  my ($hId, $hSize, @lines);
  while(<$fhi>) {
    chomp;
    if(/^\#/ || /^\s*$/) {
      next;
    } elsif(/^net (\w+) (\d+)/) {
      @lines == 0 || write_net( parse_net(\@lines, $hId, $hSize), $fho );
      ($hId, $hSize) = ($1, $2);
      @lines = ();
    } else {
      push @lines, $_;
    }
  }
  @lines == 0 || write_net( parse_net(\@lines, $hId, $hSize), $fho );
  close $fhi;
  close $fho;
}
sub write_one_net {
  my ($af, $idx, $fho) = @_;
  my ($lev, $pa, $stat, $statE, $gaps) = @{$af->[$idx]};
  my ($hId, $hb, $hl, $qId, $qSrd, $qb, $ql) = @$stat;
  my @keys = grep {exists($statE->{$_})} qw/id score ali qOver qFar qDup type/;
  my @strs = map {$_." ".$statE->{$_}} @keys;
  my ($id, $score, $type) = map {$statE->{$_}} qw/id score type/;
#    print $fho join('', (' ')x(2*$lev-1)).join(" ", 'fill', $hb, $hl, $qId, $qSrd, $qb, $ql, @strs)."\n";
  for (@$gaps) {
    my ($ghb, $ghl, $gqb, $gql, $chi) = @$_;
    my ($hBeg, $hEnd, $hLen) = ($hb+1, $ghb, $ghb-$hb);
    my ($qBeg, $qEnd, $qLen) = $qSrd eq "-" ?
      ($gqb+$gql+1, $qb+$ql, $qb+$ql-$gqb-$gql) : ($qb+1, $gqb, $gqb-$qb);
    $qLen >= 0 || print "$qb:$ql $gqb:$gql [$qLen] $qSrd $hb:$hl $ghb:$ghl [$hLen]\n"; 
    $hLen >= 0 || print "$qb:$ql $gqb:$gql [$qLen] $qSrd $hb:$hl $ghb:$ghl [$hLen]\n"; 
#      print $fho join('', (' ')x($lev*2)).join(" ", 'gap', $ghb, $ghl, $qId, $qSrd, $gqb, $gql)."\n";
    print $fho join("\t", $lev, $id, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
      $hId, $hBeg, $hEnd, "+", $hLen, ('')x3, $gql, $ghl, $score, '', '')."\n";
    for my $idC (@$chi) { 
      write_one_net($af, $idC, $fho);
    }
    ($hb, $hl) = ($ghb+$ghl, $hl-($ghb+$ghl-$hb));
    ($qb, $ql) = $qSrd eq "-" ? 
      ($qb, $gqb-$qb) : ($gqb+$gql, $ql-($gqb+$gql-$qb));
  }
  my ($hBeg, $hEnd, $hLen) = ($hb+1, $hb+$hl, $hl);
  my ($qBeg, $qEnd, $qLen) = ($qb+1, $qb+$ql, $ql); 
  $qLen >= 0 || print "$qb:$ql [$qLen] $qSrd $hb:$hl [$hLen]\n"; 
  $hLen >= 0 || print "$qb:$ql [$qLen] $qSrd $hb:$hl [$hLen]\n"; 
  print $fho join("\t", $lev, $id, $type, $qId, $qBeg, $qEnd, $qSrd, $qLen,
    $hId, $hBeg, $hEnd, "+", $hLen, ('')x3, 0, 0, $score, '', '')."\n";
}
sub write_net {
  my ($af, $fho) = @_;
  my @ids = grep {$af->[$_]->[0] == 1} (0..@$af-1);
  for my $id (@ids) {
    write_one_net($af, $id, $fho);
  }
}
sub parse_net {
  my ($lines, $hId, $hSize) = @_;
  my (@af, @al);
  my $idx = 0;
  for (@$lines) {
    if(/^(\s+)fill (.+)$/) {
      my $level = (length($1)+1) / 2;
      my @ps = split(" ", $2);
      my ($hb, $hl, $qId, $qSrd, $qb, $ql) = @ps[0..5];
      my $stat = [$hId, $hb, $hl, $qId, $qSrd, $qb, $ql];
      my @pps = @ps[6..$#ps];
      my $statE = {map {$pps[$_*2] => $pps[$_*2+1]} (0..@pps/2-1)};
      
      if($level > @al) {
        push @al, $idx;
      } elsif($level == @al) {
        $al[-1] = $idx;
      } else {
        for ($level+1..@al) { pop @al; };
        $al[-1] = $idx;
      }
      my $pa = @al>=2 ? $al[-2] : "-1";
      push @af, [$level, $pa, $stat, $statE, []];
      if($pa >= 0) {
        push @{$af[$pa]->[4]->[-1]->[4]}, $idx;
      }
      $idx ++;
    } elsif(/^(\s+)gap (.+)$/) {
      my $level = length($1) / 2;
      my @ps = split(" ", $2);
      $level <= @al || die "gap level[$level] wrong: $al[-1]\n";
      if($level < @al) {
        for ($level+1..@al) { pop @al; };
      }
      my $id = $al[-1];
      my ($hb, $hl, $qId, $qSrd, $qb, $ql) = @ps;
      my ($qId0, $qSrd0) = @{$af[$id]->[2]}[3..4];
      $qId eq $qId0 || die "fill[$id] $qId0 <> gap[$hb-$hl] $qId\n";
      $qSrd eq $qSrd0 || die "fill[$id] $qSrd0 <> gap[$hb-$hl] $qSrd\n";
      push @{$af[$id]->[4]}, [$hb, $hl, $qb, $ql, []]; 
    } else {
      die "unknonw line: $_\n";
    }
  }
  return \@af;
}

1;
__END__
