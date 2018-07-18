package Qry;
use strict;
use InitPath;
use Common;
use Localdb; 
use Path::Class;
use Bio::AlignIO;
use Bio::Seq;
use Data::Dumper;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/sync storeExp storeLocIdx/;
@EXPORT_OK = qw//;
sub new {
    my $self = shift;
    my ($dbh, $table, $cols1, $cols2) = 
        rearrange(["dbh", "table", "cols1", "cols2"], @_);
    $dbh ||= getBioDb("mt30");
    $cols1 = ref($cols1) eq "ARRAY" ? $cols1 : [$cols1];
    $cols2 = ref($cols2) eq "ARRAY" ? $cols2 : [$cols2];
    my $condStr = join(" AND ", map {$_."=?"} @$cols1);
    my $qry1 = qq/SELECT * FROM $table WHERE $condStr/;
    my $setStr = join(", ", map {$_."=?"} @$cols2);
    my $qry2 = qq/UPDATE $table SET $setStr WHERE $condStr/;
    my $colStr = join(", ", @$cols1, @$cols2);
    my $colValueStr = join(", ", ("?") x (@$cols1+@$cols2));
    my $qry3 = qq/INSERT INTO $table ($colStr) VALUES ($colValueStr)/;
    my $sth1 = $dbh->prepare($qry1);
    my $sth2 = $dbh->prepare($qry2);
    my $sth3 = $dbh->prepare($qry3);
    return bless {table=>$table, h1=>$sth1, h2=>$sth2, h3=>$sth3, cols1=>$cols1, cols2=>$cols2}, 
        ref($self) || $self;
}
sub sync {
    my $self = shift;
    my ($keys, $values) = rearrange(['keys', 'values'], @_);
    $keys = ref($keys) eq "ARRAY" ? $keys : [$keys];
    $values = ref($values) eq "ARRAY" ? $values : [$values];
    my $cntCols1 = scalar(@{$self->{cols1}});
    my $cntCols2 = scalar(@{$self->{cols2}});
    die "key cols [".@$keys."] not $cntCols1" unless @$keys == $cntCols1;
    die "value cols [".@$values."] not $cntCols2" unless @$values == $cntCols2;
    my $tag;
    my ($h1, $h2, $h3) = map {$self->{$_}} qw/h1 h2 h3/;
    $h1->execute(@$keys);
    if($h1->rows == 0) {
        $h3->execute(@$keys, @$values) or die "Error query: " . $h3->errstr;
        $tag = 1;
    } elsif($h1->rows == 1) {
        $h2->execute(@$values, @$keys) or die "Error query: " . $h2->errstr;
        $tag = 2;
    }
    return $tag;
}
=store expression data in db
my $table = "exp_rnaseq";
my $DIR_Work = dir($DIR_Misc2, "rnaseq");
$table = "exp_mtgea";
$DIR_Work = dir($DIR_Misc2, "mtgea");
my $f01 = file($DIR_Work, "01_exp.txt");
my $f02 = file($DIR_Work, "02_exp.txt");
#processExp(-in=>$f01, -out=>$f02);
my $qry = Qry->new(-table=>$table, -cols1=>"id", -cols2=>"value");
#$qry->storeExp($f02);
my ($exp, $condS, $condD) = getExp('Mtr.35108.1.S1_at', "mtgea");
print "\t".join(" ", @$condS)."\n";
for my $id (sort(keys(%$exp))) {
    print join("\t", $id, join(" ", @{$exp->{$id}}))."\n";
}
=cut
sub storeExp {
    my $self = shift;
    my ($fIn) = @_;
    my $fInH = new IO::File $fIn, "r";
    my ($cntI, $cntU) = (0, 0);
    while( <$fInH> ) {
        chomp;
        next if /^\#/;
        my @ps = split("\t");
        my ($id, $value) = ($ps[0], join(" ", @ps[1..$#ps]));
        my $tag = $self->sync($id, $value);
        $cntI ++ if $tag == 1;
        $cntU ++ if $tag == 2;
    }
    print $self->{table}." syn completed:\n";
    print "\tinserted $cntI records\n";
    print "\tupdateed $cntU records\n";
}
sub storeLocIdx {
    my $self = shift;
    my ($chr, $locPos) = @_;
    my ($cntI, $cntU) = (0, 0);
    for my $loc (sort {$a<=>$b} keys(%$locPos)) {
        my $filePos = $locPos->{$loc};
        #print join(" - ", $gene, $genePos, $filePos)."\n";
        my $tag = $self->sync([$chr, $loc], $filePos);
        $cntI ++ if $tag == 1;
        $cntU ++ if $tag == 2;
    }
    print "[Insert]$cntI [Update]$cntU\n";
}


1;
__END__
