#!/usr/bin/perl

my $help = "Usage: zcat test.vcf.gz | vcf-liftover chainfile [unmappedfile]\n";
my $chain = shift or die $help;
my $unmap = shift || '/dev/null';

my $LIFT;
my @F;

my $bodystart = 0;
while(<>){
  if(substr($_, 0, 1) eq "#"){
    print $_;
  } else {
    if($bodystart == 0){
      my $cmd = "| liftOver -minMatch=1 -bedPlus=3 /dev/stdin " . $chain . " /dev/stdout " . $unmap . " | cut -f1,3- ";
      open $LIFT, $cmd or die;
    }
    @F = split '\t';
    print $LIFT join "\t", $F[0],$F[1]-1, @F[1..$#F];
    $bodystart = 1;
  }
}

close $LIFT;
