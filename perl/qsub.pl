#!/usr/bin/perl -w
use strict;
use Common;
use Text::Template;
use Getopt::Long;

my ($opt, $prog) = ('') x 2;
GetOptions('opt|t=s'=>\$opt, 'prog|p=s'=>\$program);

my $dir = "$ENV{'code'}/pbs";
chdir($dir) or die "cannot chdir to $dir\n";
my $f_tmpl = "$dir/perl.tmpl";
my $tmpl = Text::Template->new( SOURCE=>"perl.tmpl" );

my ($queue, $nodes, $ppn, $mem, $walltime, $chdir_str, $module_str, $cmd) = ("") x 8;
my $ht = {
    'lab' => 72,
    'mirror' => 96,
    'lab-long' => 150
};

if($opt eq "run") {
    if($program eq "sv") {
        run_sv();
    } elsif($program eq "seqRec") {
        run_seqRec();
    } elsif($program eq "mpileup") {
        run_mpileup();
    } elsif($program eq "pipe_bwa") {
        run_pipe_bwa();
    } elsif($program eq "pipe_bam") {
        run_pipe_bam();
    } else {
        die("unknown program: $program\n");
    }
} elsif($opt eq "qclear") {
    clear_queue();
} elsif($opt eq "qprint") {
    get_job_stats();
} else {
    die("unknonw opt: $opt\n");
}
sub run_sv {
    ($queue, $nodes, $ppn, $mem, $walltime) = qw/lab 1 2 4GB 24:00:00/;
    my $cmd_fmt = "sv.pl --program hydra --beg %d --end %d";
    my @jobNums = ([0,9], [10,19], [20,29]);
    for my $i (0..$#jobNums) {
        my ($jobB, $jobE) = @{$jobNums[$i]};
        $cmd = sprintf $cmd_fmt, $jobB, $jobE;
        my $text = $tmpl->fill_in( HASH=> {
            nodes=>$nodes, ppn=>$ppn, mem=>$mem, walltime=>$walltime, 
            chdir_str=>$chdir_str, module_str=>$module_str, cmd=>$cmd
        } );
        my $f_job = sprintf "job%03d", $i+1;
        system("echo \"$text\" > $f_job");
        system("qsub $f_job");
        system("rm -rf $f_job");
    }
    get_job_stats();
}
sub run_seqRec {
    ($queue, $nodes, $ppn, $mem, $walltime) = qw/lab 1 2 4GB 24:00:00/;
    my $cmd_fmt = "seqRec -i /project/youngn/zhoup/Data/misc1/seq06/01_id.tbl -o /project/youngn/zhoup/Data/misc1/seq06/13_seq -t 3 -r %s";
    for my $i (1..8) {
        $cmd = sprintf $cmd_fmt, "chr$i";
        my $text = $tmpl->fill_in( HASH=> {
            nodes=>$nodes, ppn=>$ppn, mem=>$mem, walltime=>$walltime, 
            chdir_str=>$chdir_str, module_str=>$module_str, cmd=>$cmd
        } );
        my $f_job = sprintf "job%03d", $i;
        system("echo \"$text\" > $f_job");
        system("qsub $f_job");
        system("rm -rf $f_job");
    }
    get_job_stats();
}
sub run_mpileup {
    ($queue, $nodes, $ppn, $mem, $walltime) = qw/lab-long 1 4 8GB 150:00:00/;
    $chdir_str = "cd /project/youngn/zhoup/Data/repo/mt_35/30_vnt";
    my $cmd_fmt = "samtools mpileup -gD -A -q 1 -f \$data/genome/mt_35/41_genome.fa -b 01_bamlist/acc26.txt -R -r chr%d > 11_bcf/acc26_chr%d.bcf";
    for my $i (1..8) {
        $cmd = sprintf $cmd_fmt, $i, $i;
        my $text = $tmpl->fill_in( HASH=> {
            nodes=>$nodes, ppn=>$ppn, mem=>$mem, walltime=>$walltime, 
            chdir_str=>$chdir_str, module_str=>$module_str, cmd=>$cmd
        } );
        my $f_job = sprintf "job%03d", $i;
        system("echo \"$text\" > $f_job");
        system("qsub $f_job");
        system("rm -rf $f_job");
    }
    get_job_stats();
}
sub run_pipe_bwa {
    ($nodes, $ppn, $mem) = qw/1 4 8GB/;
    my $cmd_fmt = "pipe_bwa.pl --p pipe_bwa --b %d --e %d";

    my $f_task = "/project/youngn/zhoup/Data/misc3/hapmap/11_pipe_bwa/task.tbl";
    my $t = readTable(-in=>$f_task, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($beg, $end, $queue) = $t->row($i);
        my $walltime = $ht->{$queue}.":00:00";
        $cmd = sprintf $cmd_fmt, $beg, $end;
        my $text = $tmpl->fill_in( HASH=> {
            queue=>$queue, nodes=>$nodes, ppn=>$ppn, mem=>$mem, walltime=>$walltime, 
            chdir_str=>$chdir_str, module_str=>$module_str, cmd=>$cmd
        } );
        my $f_job = sprintf "job%03d-%03d", $beg, $end;
        system("echo \"$text\" > $f_job");
        system("qsub $f_job");
        system("rm -rf $f_job");
    }
    get_job_stats();
}
sub run_pipe_bam {
    ($nodes, $ppn, $mem) = qw/1 4 14GB/;
    my $cmd_fmt = "pipe_bam.pl --p pipe_bam --b %d --e %d";

    my $f_task = "/project/youngn/zhoup/Data/misc3/hapmap/15_pipe_bam/task.tbl";
    my $t = readTable(-in=>$f_task, -header=>1);
    for my $i (0..$t->nofRow-1) {
        my ($beg, $end, $queue) = $t->row($i);
        my $walltime = $ht->{$queue}.":00:00";
        $cmd = sprintf $cmd_fmt, $beg, $end;
        my $text = $tmpl->fill_in( HASH=> {
            queue=>$queue, nodes=>$nodes, ppn=>$ppn, mem=>$mem, walltime=>$walltime, 
            chdir_str=>$chdir_str, module_str=>$module_str, cmd=>$cmd
        } );
        my $f_job = sprintf "job%03d-%03d", $beg, $end;
        system("echo \"$text\" > $f_job");
        system("qsub $f_job");
        system("rm -rf $f_job");
    }
    get_job_stats();
}

sub get_job_stats {
    system("qstat -u zhoup > tmp");
    open(FH, "<tmp");
    my @stats;
    while( <FH> ) {
        chomp;
        my @ps = split " ";
        next unless @ps >= 10;
        my ($id, $user, $queue, $name, $id_session, $nodes, $mem, $timeR, $timeE) = @ps;
        if( $name=~/^job/ ) {
            push @stats, [$id, $queue, $name, $nodes, $mem, $timeR, $timeE];
        }
    }
    close FH;
    system("rm tmp");
    for my $i (0..$#stats) {
        print join("\t", @{$stats[$i]})."\n";
    }
    return @stats;
}
sub clear_queue {
    my @stats = get_job_stats();
    print "--------------------------------------------\n";
    for (@stats) {
        my ($id, $user, $queue, $name, $id_session, $nodes, $mem, $timeR, $timeE) = @$_;
        print "qdel $id\n";
        my $id2 = [split(/\./, $id)]->[0];
        print("rm $dir/*$id2\n");
    }
    print "--------------------------------------------\n";
}


