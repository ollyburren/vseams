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
perl vseams.pl --[i]ni_file {-[v]erbose -[h]elp)

Variant set enrichment of GWAS statistics using multivariate sampling.

	ini_file - path to a suitable inifile (see  http://github.com/ollyburren/ini/default.ini) for a template
	verbose - enable verbose output to STDERR - currently not implemented
	help - print this message.

For help and suggestions please contact Olly Burren (http://github.com/ollyburren)
EOL

my ($ini_file,$verbose,$help,$ERROR_FLAG);
GetOptions (
	'ini_file|i=s' => \$ini_file,
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

## grab TEST settings here
my %TESTS;
foreach my $sect(grep{/TEST/}$cfg->Sections()){
	(my $tname = $sect)=~s/TEST (.*)/\1/;
	unless($TESTS{$tname}{test_set}=$cfg->val($sect,'test_set')){
		print STDERR "$tname test set is misconfigured please check\n";
		$ERROR_FLAG=1;
	}
	unless($TESTS{$tname}{control_set}=$cfg->val($sect,'control_set')){
		print STDERR "$tname control set is misconfigured please check\n";
		$ERROR_FLAG=1;
	}
	$TESTS{$tname}{test_set}="SET.".$TESTS{$tname}{test_set};
	$TESTS{$tname}{control_set}="SET.".$TESTS{$tname}{control_set};
}

if((keys %TESTS)==0){
	print STDERR "Cannot find any tests please check ini\n";
	$ERROR_FLAG=1;
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
my $N_PERMS = $PROJECT{number_perms};
my $CHUNK_SIZE = $PROJECT{chunk_size};
my $R2_THRESHOLD = $PROJECT{r2_threshold};
my $MAX_SNPS_PER_SIGMA_BLOCK = $PROJECT{max_snps_sigma_block};
my $SIGMA_CACHE_DIR = $PROJECT{sigma_cache_dir};
my $LD_INDEX_FILE = $PROJECT{ld_index_file};
my $QUEUE_ENGINE = $PROJECT{queue_engine};
my $MACD_CONF_FILE = $PROJECT{mac_d_conf_file};
my $GENE_SET_FILE = $PROJECT{region_set_file};
my $EXCLUDE_FILE = $PROJECT{exclude_regions};
my $PROJECT_DIR = "$BASE_DIR/$PROJECT";
my $TSS_EXTENSION = $PROJECT{tss_extension};
my $MART_HOST = $PROJECT{mart_host};
my %DATASETS = %{$PROJECT{datasets}};
my $DRIVER = Macd::GRIDDriverFactory->instantiate($QUEUE_ENGINE,inifile => $MACD_CONF_FILE);
	
#######################
## STEP 1 GET REGIONS##
#######################

## THIS STEP DOES NOT DEPEND ON DATASET

my $ofile = "$PROJECT_DIR/support/support.RData";   


my $gr = VSEAMS::Runnable::GetRegion->new(
		log_dir=>"$PROJECT_DIR/log/getregion/",
		env_settings=>$RLIB_ENV,                         
		macd_driver=>$DRIVER,
		inputs=>{
			input_file=>$GENE_SET_FILE,
			mart_host=>$MART_HOST,
			tss_extension=>$TSS_EXTENSION
			
		},
		outputs=>{
			output_file=>$ofile,
		}
)->run_step;

sleep(5);
## next setup a dir for writing shared variables to

my $shared_dir = $PROJECT_DIR.'/shared/';
if(-d $shared_dir){
	File::Path::remove_tree($shared_dir,{keep_root => 1,result=> \my $del_dirs});
}else{
	Macd::Utils::mkpath_wrapper($shared_dir);
}         

## THIS MUST BE RUN FOR EACH DATASET
foreach my $ds(keys %DATASETS){
	if (fork() == 0) {
		my %STEP;
		my $dataset_dir_stub = "$PROJECT_DIR/$ds";
		my $log_dir_stub = "$dataset_dir_stub/log/";
		
		###########################
		## STEP 2 PREPARE SUPPORT##
		###########################
		
		my $odir = "$dataset_dir_stub/support/";
		
		$STEP{$ds}{CreateSupport} = VSEAMS::Runnable::CreateSupport->new(
				log_dir=>"$log_dir_stub/createsupport/",
				env_settings=>$RLIB_ENV,
				macd_driver=>$DRIVER,
				previous_step=>$gr,
				inputs=>{
					exclude_bed=>$EXCLUDE_FILE,
					gwas_bed =>$DATASETS{$ds},
					region_file=>$gr->inputs->{output_file},
					ld_index_file=>$LD_INDEX_FILE,
					max_snps_per_sigma_block=>$MAX_SNPS_PER_SIGMA_BLOCK
				},
				outputs=>{
					out_dir=> $odir
				}
		)->run_step();

		################################
		## STEP3 PRUNE AND CACHE SIGMA##
		################################
		
		$odir = "$dataset_dir_stub/pruned/";
		$STEP{$ds}{PruneSNPsCacheSigma} = VSEAMS::Runnable::PruneSNPsCacheSigma->new(
			log_dir=>"$log_dir_stub//prunesnpsigma/",
			env_settings=>$RLIB_ENV,
			macd_driver=>$DRIVER,
			previous_step=>$STEP{$ds}{CreateSupport},
			debug_flag=>1,
			inputs=>{
				in_dir=>$STEP{$ds}{CreateSupport}->inputs->{out_dir},
				sigma_cache_dir=>$SIGMA_CACHE_DIR,
				r2_threshold => $R2_THRESHOLD,
			},
			outputs=>{
				out_dir=>$odir
			}
		)->run_step();
		
		##########################################
		## STEP 3.5 DOWNSAMPLE THE MASSIVE SIGMA##
		##########################################
		## TODO MAKE UPPER LIMIT ON SNP_COUNT A PARAMETER
		
		$odir = "$dataset_dir_stub/ds_sigma/";
		$STEP{$ds}{DownSampleSigma} = VSEAMS::Runnable::DownSampleSigma->new(
			log_dir=>"$log_dir_stub/downsamplesigma/",
			env_settings=>$RLIB_ENV,
			macd_driver=>$DRIVER,
			previous_step=>$STEP{$ds}{PruneSNPsCacheSigma},
			debug_flag=>1,
			chunk_size=>$CHUNK_SIZE,
			inputs=>{
				snp_dir=>$STEP{$ds}{PruneSNPsCacheSigma}->inputs->{out_dir},
				downsample_file=>$STEP{$ds}{CreateSupport}->inputs->{out_dir}."downsample_files.tab",
				max_snps_per_sigma_block=>$MAX_SNPS_PER_SIGMA_BLOCK,
				sigma_cache_dir=>$STEP{$ds}{PruneSNPsCacheSigma}->inputs->{sigma_cache_dir},
				test=>-1
			},
			outputs=>{
				out_dir=>$odir
			}
			)->run_step();
			
		#############################	
		## STEP4 CALCULATE MVS PERM##
		#############################
		
		$odir = "$dataset_dir_stub/perms/";
		
		$STEP{$ds}{CalcMVSPerm} = VSEAMS::Runnable::CalcMVSPerm->new(
			log_dir=>"$log_dir_stub/calcmvsperm/",
			env_settings=>$RLIB_ENV,
			macd_driver=>$DRIVER,
			previous_step=>$STEP{$ds}{PruneSNPsCacheSigma},
			debug_flag=>1,
			chunk_size=>$CHUNK_SIZE,
			inputs=>{
				in_dir=>$STEP{$ds}{PruneSNPsCacheSigma}->inputs->{out_dir},
				sigma_cache_dir=>$STEP{$ds}{PruneSNPsCacheSigma}->inputs->{sigma_cache_dir},
				sample_sigma_dir=>$STEP{$ds}{DownSampleSigma}->inputs->{out_dir},
				test=>-1,
				n_perms=>$N_PERMS
			},
			outputs=>{
				out_dir=>$odir
			}
			)->run_step();
		
		
		## THESE STEPS ARE TEST SPECIFIC
			
		###########################################
		### STEP4B CONSOLIDATE SNPS INTO MANIFEST##
		#############################	##############
		foreach my $T(keys %TESTS){
			$odir = "$dataset_dir_stub/$T/snp_manifest/";
			$STEP{$ds}{$T}{SNPManifest} = VSEAMS::Runnable::SNPManifest->new(
			log_dir=>"$log_dir_stub/$T/snpmanifest/",
					env_settings=>$RLIB_ENV,
					macd_driver=>$DRIVER,
					previous_step=>$STEP{$ds}{CalcMVSPerm},
					debug_flag=>1,
					inputs=>{
						snp_dir=>$STEP{$ds}{PruneSNPsCacheSigma}->inputs->{out_dir},
						test_set => $TESTS{$T}->{test_set},
						control_set=>$TESTS{$T}->{control_set},
						test=>-1,
					},
					outputs=>{
						out_dir=>$odir
					}
				)->run_step();
		}
		## run the below when the above has finished
		## let filesystem catch up
		sleep(5);
		foreach my $T(keys %TESTS){
			my $runObj = $STEP{$ds}{$T}{SNPManifest};
			next if $runObj->skipped;
			$runObj->step->wait_on_complete();
			if($runObj->check_R_jobs){
				print "Done check ".$runObj->outputs->{out_dir}."\n";
			}else{
				die "Encountered an issue with $ds for $T ".ref($runObj)." check ".$runObj->log_dir."\n";
			}
		}
		
		

		##########################
		## STEP5 CALCULATE WSTAR##
		##########################
		foreach my $T(keys %TESTS){
				$odir = "$dataset_dir_stub/$T/wstar/";
				$STEP{$ds}{$T}{CalcWstar} = VSEAMS::Runnable::CalcWstar->new(
					log_dir=>"$log_dir_stub/$T/calcwstar/",
					env_settings=>$RLIB_ENV,
					macd_driver=>$DRIVER,
					previous_step=>$STEP{$ds}{$T}{SNPManifest},
					debug_flag=>1,
					inputs=>{
						perm_dir=>$STEP{$ds}{CalcMVSPerm}->inputs->{out_dir},
						snp_manifest=>$STEP{$ds}{$T}{SNPManifest}->out_file,
						test_set => $TESTS{$T}->{test_set},
						control_set=>$TESTS{$T}->{control_set},
						test=>-1,
						n_perms=>$STEP{$ds}{CalcMVSPerm}->inputs->{n_perms},
						chunk_size=>$STEP{$ds}{CalcMVSPerm}->chunk_size
					},
					outputs=>{
						out_dir=>$odir
					}
				)->run_step();
		}
		
		## run the below when the above has finished
		## let filesystem catch up
		sleep(5);
		foreach my $T(keys %TESTS){
			my $runObj = $STEP{$ds}{$T}{CalcWstar};
			next if $runObj->skipped;
			$runObj->step->wait_on_complete();
			if($runObj->check_R_jobs){
				print "Done check ".$runObj->outputs->{out_dir}."\n";
			}else{
				die "Encountered an issue with $ds for $T ".ref($runObj)." check ".$runObj->log_dir."\n";
			}
		}
		
		
		
		####################################
		## STEP6 CALCULATE WILCOXON PARAMS##
		####################################
		
		foreach my $T(keys %TESTS){
				$odir = "$dataset_dir_stub/$T/wilcoxon_params/";
				$STEP{$ds}{$T}{CalcWilcoxonParams} = VSEAMS::Runnable::CalcWilcoxonParams->new(
					log_dir=>"$log_dir_stub/$T/paramwilcoxon/",
					env_settings=>$RLIB_ENV,
					macd_driver=>$DRIVER,
					previous_step=>$STEP{$ds}{$T}{CalcWstar},
					debug_flag=>1,
					inputs=>{
						wstar_dir=> $STEP{$ds}{$T}{CalcWstar}->inputs->{out_dir},
						snp_manifest=>$STEP{$ds}{$T}{SNPManifest}->out_file,
						test_set=>$STEP{$ds}{$T}{CalcWstar}->inputs->{test_set},
						control_set=>$STEP{$ds}{$T}{CalcWstar}->inputs->{control_set},
						dataset=>$ds
					},
					outputs=>{                                                            
						out_dir=>$odir
					}
				)->run_step();
		}
		sleep(5);
		foreach my $T(keys %TESTS){
			my $runObj = $STEP{$ds}{$T}{CalcWilcoxonParams};
			next if $runObj->skipped;
			$runObj->step->wait_on_complete();
				if($runObj->check_R_jobs){
				print "Done check ".$runObj->outputs->{out_dir}."\n";
			}else{
				die "Encountered an issue with $ds for $T ".ref($runObj)." check ".$runObj->log_dir."\n";
			}
			my @ofile = (
					$ds,
					$runObj->inputs->{test_set},
					$runObj->inputs->{control_set},
					'param',
					'RData'
				);
			test_success(join(".",@ofile[1,2]),$runObj->outputs->{out_dir}.join(".",@ofile),$shared_dir."/success.tab");
		}
		exit; ## required due to fork
	}
}

########################################
## STEP 7: CALCULATE ENRICHMENT SCORES##
########################################

## wait for processes to finish
## read in input files for final step

while (wait() != -1) {}
my %SUCCESSFUL_TESTS;
if(-e $shared_dir."/success.tab"){
	open(SUCCESS,$shared_dir.'/success.tab');
	while(<SUCCESS>){
		chomp();
		my ($t,$f)=split("\t",$_);
		push @{$SUCCESSFUL_TESTS{$t}},$f;
	}
	close(SUCCESS);
}else{
	print "Encountered an error with reading $shared_dir/success.tab - please check !\n";
	exit;
}
		
my $results = VSEAMS::Runnable::CalcZ->new(
	log_dir=>"$PROJECT_DIR/log/calcz/",
	env_settings=>$RLIB_ENV,
	macd_driver=>$DRIVER,
	previous_step=>undef,
		debug_flag=>1,
		inputs=>{    
			tests => \%SUCCESSFUL_TESTS
		},                    
		outputs=>{                       
			out_dir=>$PROJECT_DIR.'/results/'
		}
	)->run_step();

	
sub display_results{
	my $res_dir = shift;
	my $RSCRIPT=<<ERSCRIPT;
rdir<-'$res_dir'
rfiles<-list.files(path=rdir,pattern="*.\\.RData",full.names=TRUE)
test.names<-gsub(".RData","",basename(rfiles))
test.names<-sub("^SET\\.","",test.names)
test.names<-sub("\\.SET\\.","__vs__",test.names)
res<-do.call("rbind",lapply(seq_along(rfiles),function(i){
	r<-get(load(rfiles[i]))
	list(test=test.names[i],p.theoretical=r\$Z.theoretical\$p.value,p.empirical=r\$Z.empirical\$p.value)
	}
))
write.table(res,file='$res_dir/results.tab',row.names=FALSE,quotes=FALSE)
ERSCRIPT
`Rscript --vanilla << $RSCRIPT`;
}
	
sub test_success{
	my ($test,$fname,$filename)=@_;
	open my $fh, ">>", $filename or die  "$0 [$$]: open: $!";
	flock $fh, LOCK_EX or die  "$0 [$$]: flock: $!";
	print $fh join("\t",$test,$fname)."\n" or die  "$0 [$$]: write: $!";
	close $fh or warn "$0 [$$]: close: $!";
}




