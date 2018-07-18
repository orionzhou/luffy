package Wiggle;
=head NAME
    my $fWib = file($DIR_Out, "tmp.wib");
    my $wig = Wiggle->new(-file=>$fWib, -write=>1);
    $wig->load(file($DIR_Out, "tmpCov.txt"));
    print join("=", @{$wig->read(-start=>10, -end=>20)});
=cut
use strict;
use warnings; 
use IO::File;
use Carp qw/croak confess/; 
use Common;
use List::Util qw/min max sum/; use POSIX qw/ceil floor/;
use constant BODY => 'n';
use constant MAXVALUE => 65535;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw//;
@EXPORT_OK = qw//;
sub new {
    my $self = shift;
    my ($fpath, $write) = rearrange(["file", "write"], @_);
    $fpath ||= '';
    my $mode = $write ? -e $fpath   # if file already exists...
                                                  ? '+<'    # ...open for read/write
                                                  : '+>'    # ...else clobber and open a new one
                                        : '<';       # read only
    my $fh = $fpath ? IO::File->new($fpath,$mode) : IO::File->new_tmpfile;
    $fh or die (($fpath||'temporary file').": $!");
    return bless { fh => $fh, write => $write }, ref $self || $self;
}
sub load {
    my $self = shift;
    my ($fIn) = @_;
    my $fInH = new IO::File $fIn, "r";
    my $fOutH = $self->{fh};
    my $packFormat = BODY."*";
    my (@buffer, @tmp);
    my $cnt = 0;
    my $sum = 0;
    while( <$fInH> ) {
        next if /^#/;
        next unless $_;
        $cnt ++;
        chomp;
        my $num = $_;
        if($num > MAXVALUE) {
            push @tmp, $cnt;
            $num = MAXVALUE;
        }
        push @buffer, $num;
        $sum += $num;
        if(@buffer >= 500_000) {
            print $fOutH pack($packFormat, @buffer);
            @buffer = ();
        }
    }
    print $fOutH pack($packFormat, @buffer) unless @buffer == 0;
    return $sum;
}
sub read {
    my $self = shift;
    my ($start, $end) = rearrange(['start', 'end'], @_);
    my $fh = $self->{fh};
    seek($fh, 0, 0);
    my $offSet = ($start-1) * 2;
    seek($fh, $offSet, 0);
    my $span = $end - $start + 1;
    my @values;
    for my $i (1..$span) {
        my $buffer;
        read($fh, $buffer, 2);
        my $value = unpack( BODY, $buffer );
        push (@values, $value);
    }
    return \@values;
}


1;
__END__
