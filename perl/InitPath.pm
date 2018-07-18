package InitPath;
use strict;
use File::Path qw/make_path remove_tree/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/$DIR_work $DIR_data $DIR_src $DIR_code $DIR_tmp
  $DIR_in $DIR_out $DIR_misc1 $DIR_misc2 $DIR_misc3 $DIR_misc4
  $DIR_genome $DIR_db/;
@EXPORT_OK = qw//;
our $DIR_src = $ENV{'src'};
our $DIR_work = $ENV{'work'};
our $DIR_data = $ENV{'data'};
our $DIR_code = $ENV{'code'};
our $DIR_genome = "$DIR_data/genome";
our $DIR_db = "$DIR_data/db";
our $DIR_in = "$DIR_data/in";
our $DIR_out = "$DIR_data/out";
our $DIR_misc1 = "$DIR_data/misc1";
our $DIR_misc2 = "$DIR_data/misc2";
our $DIR_misc3 = "$DIR_data/misc3";
our $DIR_misc4 = "$DIR_data/misc4";
our $DIR_tmp = "/scratch/zhoup/tmp";
-d $DIR_tmp || make_path($DIR_tmp);

1;
__END__
