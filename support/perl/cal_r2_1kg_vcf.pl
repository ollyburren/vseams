#!/usr/bin/perl

## Macd location (get from my github)
use lib 'FILL IN PATH/macd/lib';
use Macd;
use FindBin qw($Bin);
use File::Find;
use strict;

## site specific conf file for macd - see http://github.com/ollyburren/macd
## for more details
my $grid_cnf = 'FILL ME IN';
my $log_dir = 'FILL ME IN';
	
my $DRIVER = Macd::GRIDDriverFactory->instantiate('SGE',inifile => $grid_cnf);

my $step = Macd::Step->new(
	logdir=> $log_dir,
	driver=>$DRIVER,
	env_settings=>"R_LIBS=$ENV{R_LIBS}"
);


###############
#CONFIGURATION#
###############
## where to mail to once finished.
my $MAIL_TO = 'nobody@nowhere.com';
my $TEST = 0; #IF set to true allows us to test the script by running only a few jobs
# location of recombination frequency based LD blocks from previous step
my $DATADIR = 'FILL ME IN';
my $OUTDIR = 'FILL ME IN';
# extension of thousand genomes vcf files - if recommended source used extensio should not need modification
my $VCF_EXTENSION='.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.bgz';
# location of 1Kgenome vcf files split by chromosome
my $THOUSAND_GENOME_DATA_DIR='FILL ME IN';
my $FILEPATTERN='RData$';
$FILEPATTERN='chr1_.*\.RData$' if $TEST;
my $RSCRIPT='/usr/bin/Rscript --vanilla $Bin/../R/calc_r2_recomb.R';

find(sub {
		if(/$FILEPATTERN/){
			 my @param = "$RSCRIPT";
      push @param, "region.file=\\'$File::Find::name\\'";
      push @param, "vcfext=\\'$VCF_EXTENSION\\'";
      push @param, "out_dir=\\'$OUTDIR\\'";
      push @param, "thou_gen_data_dir=\\'$THOUSAND_GENOME_DATA_DIR\\'";
      my $cmd = join(" ",@param);
			my $cmd = "$RSCRIPT region.file=\\'$File::Find::name\\'";
			my $job =  Macd::Step::Job->new(command=>$cmd);
			$step->add_job($job);
			print "Added $File::Find::name\n";
			$FILEPATTERN='NOFINDME' if $TEST;
		}
	}, ($DATADIR));


if($step->execute()){
	print "Step submitted successfully\n";
	## we can hold up prog execution as follows
	## in case we have a step below that requires the output
	$step->wait_on_complete();
	`echo "Step completed successfully but please check $log_dir " | mail -s "Step finished $0" $MAIL_TO`;
	print "Step completed successfully\n";
}else{
	#print "Error submitting step\n";
	`echo "Step failed, but please check $log_dir " | mail -s "Step failed $0" $MAIL_TO`;
}
