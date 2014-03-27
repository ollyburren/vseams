#########
#PRAGMA #
#########

use strict;
use FindBin qw($Bin) ;
use lib ("$Bin/perl/lib","$Bin/../macd/lib");
use Getopt::Long;
use Config::IniFiles;
use Macd;
use Data::Dumper;
use VSEAMS::Runnable;
use Fcntl qw/ :flock /; 

#######################################
#OPTION HANDLING AND INI FILE PARSING #
#######################################

my $USAGE=<<EOL;
perl vseams_prioritise.pl --[i]ni_file -[t]arget  {- [n]_perms -[v]erbose -[h]elp)

Prioritise regions identified from VSEAMS enrichment.

	ini_file - path to a suitable inifile (see  http://github.com/ollyburren/ini/default.ini) for a template
	target - Region/Gene set to prioritise.
	n_perms - Number of permutations to run default is number in ini file
	verbose - enable verbose output to STDERR - currently not implemented
	help - print this message.

For help and suggestions please contact Olly Burren (http://github.com/ollyburren)
EOL

my ($ini_file, $target, $verbose,$n_perms,$help,$ERROR_FLAG);
GetOptions (
	'ini_file|i=s' => \$ini_file,
	'target|t=s'	=> \$target,
	'n_perms|n=s' => \$n_perms,
	'help|h' => \$help,
	'verbose|v' => \$verbose
	);

if($help){
	print STDERR $USAGE."\n";
	exit(1);
}
	
if(!-e $ini_file){
	print STDERR "Cannot find ini_file:$ini_file. See http://github.com/ollyburren/ini/default.ini for details\n";
	exit(1);
}

if(! $target){
	print STDERR "No target:$target defined\n";
	exit(1);
}

if($verbose){
	print STDERR "Verbose flag not yet implemented.\n";
}

## configuration


my $cfg = Config::IniFiles->new( -file => $ini_file );

## grab project setting 
## not liking the tied hash technique 
my %PROJECT;
##grab file based parameters 
foreach my $p($cfg->Parameters('PROJECT files')){
	my $set = $cfg->val('PROJECT files',$p);
	## this must be datasets bit
	if($set =~m/\n/){
		my %ds;
		foreach my $d(split("\n",$set)){
			my ($name,$path)=split(":",$d);
			if(! -e $path){
				print STDERR "Cannot access path for dataset $name $path\n";
				$ERROR_FLAG=1;
			}else{
				$ds{$name}=$path;
			}
		}
		$PROJECT{$p}=\%ds;
		next;
	}
	if(! -e $set){
		print STDERR "Cannot access path for $p:$set\n";
		$ERROR_FLAG=1;
	}else{
		$PROJECT{$p}=$cfg->val('PROJECT files',$p);
	}
}

## grab other mandatory settings

foreach  my $p($cfg->Parameters('PROJECT settings')){
	if(my $set = $cfg->val('PROJECT settings',$p)){
		$PROJECT{$p}=$set
	}else{
		print STDERR "No setting for $p\n";
		$ERROR_FLAG=1;
	}
}



exit(1) if $ERROR_FLAG;


##############################################
#END OF OPTION HANDLING AND INI FILE PARSING #
##############################################

#######################
##SET GLOBAL VARIABLES#
#######################
my $RLIB_ENV = "R_LIBS=$ENV{R_LIBS}";
my $BASE_DIR = $PROJECT{base_dir};
my $PROJECT = $PROJECT{name};
my $N_PERMS = $n_perms || $PROJECT{number_perms};
my $SIGMA_CACHE_DIR = $PROJECT{sigma_cache_dir};
my $QUEUE_ENGINE = $PROJECT{queue_engine};
my $MACD_CONF_FILE = $PROJECT{mac_d_conf_file};
my $GENE_SET_FILE = $PROJECT{region_set_file};
my $LD_INDEX_FILE = $PROJECT{ld_index_file};
#my $EXCLUDE_FILE = $PROJECT{exclude_regions};
my $PROJECT_DIR = "$BASE_DIR/$PROJECT";
my %DATASETS = %{$PROJECT{datasets}};
my $DRIVER = Macd::GRIDDriverFactory->instantiate($QUEUE_ENGINE,inifile => $MACD_CONF_FILE);


foreach my $ds(keys %DATASETS){
	my %STEP;
	my $dataset_dir_stub = "$PROJECT_DIR/$ds";
	my $log_dir_stub = "$dataset_dir_stub/log/";
	my $odir = "$dataset_dir_stub/prioritise_genes/";
	$STEP{$ds}{PrioritiseGenes} = VSEAMS::Runnable::PrioritiseRegions->new(
		log_dir=>"$log_dir_stub/prioritise_genes/",
		env_settings=>$RLIB_ENV,                         
		macd_driver=>$DRIVER,
		inputs=>{
			test_name => $target,
			input_file=>$GENE_SET_FILE,		
			support_file=>"$PROJECT_DIR/support/support.RData", ## need to check for this
			snp_manifest => $DATASETS{$ds},
			sigma_cache_dir => $SIGMA_CACHE_DIR,
			n_perms => $N_PERMS,
			index_file => $LD_INDEX_FILE,
		},
		outputs=>{
			out_dir=>$odir,
		}
	)->run_step;
}
