#!/usr/bin/perl
# ===========================================================================
# return codes:
#     0 - success
#     1 - invalid arguments
#     2 - no hits found
#     3 - rid expired
#     4 - search failed
#     5 - unknown error
# ===========================================================================

use Getopt::Long;
use Pod::Usage;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
use strict;

my $fi = '';
my $fo = '';
my $fhi;
my $fho;
my $help_flag;

#----------------------------------- MAIN -----------------------------------#
GetOptions(
    "help|h"  => \$help_flag,
    "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my $ua = LWP::UserAgent->new;
my @ids;
($fi)= @ARGV;
if(!$fi) {
    pod2usage(2);
} elsif ($fi eq '-' || $fi eq "stdin") {
    $fhi = \*STDIN;
} else {
    open ($fhi, $fi) || die "Can't open file $fi: $!\n";
}

my $qry = "";
while (<$fhi>) {
    $qry .= uri_escape($_);
}

my ($program, $db) = ("blastn", "nr");
my $args = "CMD=Put&PROGRAM=$program&DATABASE=$db&QUERY=".$qry;

my $req = new HTTP::Request POST => 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
$req->content_type('application/x-www-form-urlencoded');
$req->content($args);

# get the response
my $response = $ua->request($req);

# parse out the request id
$response->content =~ /^    RID = (.*$)/m;
my $rid=$1;

# parse out the estimated time to completion
$response->content =~ /^    RTOE = (.*$)/m;
my $rtoe=$1;
print "  RID:$rid  waiting $rtoe seconds\n";
sleep $rtoe;

while(1) {
    sleep 3;

    $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
    $response = $ua->request($req);

    if ($response->content =~ /\s+Status=WAITING/m) {
        # print STDERR "Searching...\n";
        next;
    }
    if ($response->content =~ /\s+Status=FAILED/m) {
        print STDERR "Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
        exit 4;
    }
    if ($response->content =~ /\s+Status=UNKNOWN/m) {
        print STDERR "Search $rid expired.\n";
        exit 3;
    }
    if ($response->content =~ /\s+Status=READY/m) {
        if ($response->content =~ /\s+ThereAreHits=yes/m) {
            #  print STDERR "Search complete, retrieving results...\n";
            last;
        } else {
            print STDERR "No hits found.\n";
            exit 2;
        }
    }

    # if we get here, something unexpected happened.
    exit 5;
}

$req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid";
$response = $ua->request($req);

if(!$fo || $fo eq "stdout") {
    $fho = \*STDOUT;
} else {
    open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}
print $fho $response->content;
close $fho;

exit 0;
