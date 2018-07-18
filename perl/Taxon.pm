package Taxon;
use strict; 
use Bio::DB::Taxonomy;
use Data::Dumper;
use List::Util qw/min max sum/; 
use POSIX qw/ceil floor/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/gi2Taxid annotate_taxid/;
@EXPORT_OK = qw//;

my $dir = "$ENV{'data'}/ncbi_taxon";
my $fnodes = "$dir/nodes.dmp";
my $fnames = "$dir/names.dmp";
my $nbin = "$dir/gi_taxid_nucl.bin";
my $pbin = "$dir/gi_taxid_prot.bin";

sub gi2Taxid {
  my @gis = @_;
  my $dict = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict => $nbin); 
  
  my $h = { map {$_=>''} @gis };
  @gis = keys(%$h);
  for my $gi (@gis) {
    my $taxid = $dict->get_taxid($gi);
  }
  printf "%6d (unique) GIs searched\n", scalar(@gis);
  return $h;
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

