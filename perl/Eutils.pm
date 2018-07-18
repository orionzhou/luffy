package Eutils;
use strict; 
use Bio::Seq;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/gi2Taxid annotate_taxid/;
@EXPORT_OK = qw//;

my $dir = "$ENV{'data'}/db/ncbi_taxon";
our $fnodes = "$dir/nodes.dmp";
our $fnames = "$dir/names.dmp";

sub get_gi_note {
  my @ids = @_;
  my $fac = Bio::DB::EUtilities->new(-eutil=>'esummary', -email=>'zhoupenggeni@gmail.com', -db=>'nuccore', -id=>\@ids);
  my @notes;
  while(my $ds = $fac->next_DocSum) {
    my @items = $ds->get_Items_by_name("Title");
    my $note = $items[0]->get_content;
    push @notes, $note;
  }
  return @notes;
}
sub gi2Taxid {
  my @gis = @_;
  my $h = { map {$_=>''} @gis };
  @gis = keys(%$h);

  my $batch = ceil(@gis/10000);
  print "split task in $batch batches\n" if $batch > 1;

  for my $i (1..$batch) {
    my $idxb = ($i - 1) * 10000;
    my $idxe = min($i * 10000 - 1, $#gis);
    my @sgis = @gis[$idxb..$idxe];
    my $fac = my $factory = Bio::DB::EUtilities->new(-eutil=>'esummary', 
      -email=>'zhoupenggeni@gmail.com',
      -db=>'nucleotide', -id=>\@sgis);
    while( my $ds = $fac->next_DocSum) {
      my $gi;
      while(my $it = $ds->next_Item("flattened")) {
        $gi = $it->get_content if $it->get_name eq "Gi";
        if($it->get_name eq "TaxId") {
          die "no GI error\n" unless $gi;
          $h->{$gi} = $it->get_content;
        }
      }
    }
  }
  
  my $nt = scalar(keys(%$h));
  my $ns = scalar(grep {$_ ne ''} values(%$h));
  printf "%6d / %6d (unique) GIs searched\n", $ns, $nt;
  return $h;
}
sub gi2Taxid_entrez {
  my @gis = @_;
  my $db = Bio::DB::Taxonomy->new(-source=>'entrez');
  my $h1 = { map {$_=>''} @gis };
  my $h2 = { map {$_=>''} @gis };
  my @ids = keys(%$h1);
  for my $i (0..$#ids) {
    my $gi = $ids[$i];
    my $node = $db->get_Taxonomy_Node(-gi=>$gi, -db=>'nucleotide');
    $h2->{$gi} = $node->scientific_name;
    my %l;
    my @ary;
    for my $j (1..30) {
      my $pa = $node->ancestor();
      last if(!$pa);
      if($pa->rank eq "no rank") {
        push @ary, $pa->scientific_name;
      } else {
        $l{$pa->rank} = $pa->scientific_name;
      }
      $node = $pa;
    }
    my $cat1;
    if(!exists($l{'kingdom'}) && !exists($l{'superkingdom'})) {
      $cat1 = $ary[-1];
    } elsif(exists($l{'kingdom'})) {
      $cat1 = $l{'kingdom'};
      if($l{'kingdom'} =~ /(Viridiplantae)|(Metazoa)/i) {
        $cat1 = $l{"family"} if exists $l{"family"};
      }
    } elsif(exists($l{'superkingdom'})) {
      $cat1 = $l{'superkingdom'};
    }
    $h1->{$gi} = $cat1;
    printf "%4d | %4d\n", $i+1, scalar(@ids) if ($i+1) % 10 == 0;
  }
#    my @cats1 = map {$h1->{$_}} @gis;
#    my @cats2 = map {$h2->{$_}} @gis;
#    return (\@cats1, \@cats2);
  return ($h1, $h2);
}
sub annotate_taxid {
  my @ids = @_;
  my $h = { map {$_=>''} @ids };
  @ids = keys(%$h);

  my $db = Bio::DB::Taxonomy->new(
    -source => 'flatfile', 
    -nodesfile => $fnodes, 
    -namesfile => $fnames, 
    -directory => $dir);
  for my $i (0..$#ids) {
    my $id = $ids[$i];
    my @cats = ('') x 4;
    my $node = $db->get_taxon(-taxonid=>$id);
    my @ary;
    while($node) {
      if($node->rank eq "no rank") {
        push @ary, $node->scientific_name;
      }
      $cats[0] = $node->scientific_name if $node->rank eq "superkingdom";
      $cats[1] = $node->scientific_name if $node->rank eq "kingdom";
      $cats[2] = $node->scientific_name if $node->rank eq "family";
      $cats[3] = $node->scientific_name if $node->rank eq "species";
      
      $node = $node->ancestor();
    }
    if($cats[1] eq "") {
        $cats[1] = "Bacteria" if $cats[0] eq "Bacteria";
    }
    $h->{$id} = \@cats;
  }
  return $h;
}




1;
__END__
        my $cat1;
        if(!exists($l{'kingdom'}) && !exists($l{'superkingdom'})) {
            $cat1 = $ary[-1];
        } elsif(exists($l{'kingdom'})) {
            $cat1 = $l{'kingdom'};
            if($l{'kingdom'} =~ /(Viridiplantae)|(Metazoa)/i) {
                $cat1 = $l{"family"} if exists $l{"family"};
            }
        } elsif(exists($l{'superkingdom'})) {
            $cat1 = $l{'superkingdom'};
        }

