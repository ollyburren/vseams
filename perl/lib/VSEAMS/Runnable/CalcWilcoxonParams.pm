package VSEAMS::Runnable::CalcWilcoxonParams;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;
use File::Find;

use constant _inputs_expected => {
	out_dir => 'file',
	wstar_dir => 'file',                      
	snp_manifest => 'file',
	test_set => 'string',
	control_set => 'string',                      
	dataset => 'string',
};

use constant _defaults => {
};                          
	

use constant _script => 'compute_wilcoxon_params.R';

sub run{
	my $self=shift;
  if(-d $self->inputs->{out_dir}){
		File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($self->inputs->{out_dir});
	}         
	my $cmd = $self->build_Rscript_cmd();  
	my $job = Macd::Step::Job->new(command=>$cmd);
	$self->step->add_job($job);
	if($self->step->execute()){
		#$self->step->wait_on_complete();  
		#if($self->check_R_jobs){                    
		#	print "Done check ".$self->outputs->{out_dir}."\n";
			return 1;
		#}
	}
	return 0;
	
}

sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: The output path ".$self->inputs->{out_dir}." already exists and contains files. Would you like to skip this step ?";
}

sub run_conditions{
	my $self=shift;
	my $odir = $self->inputs->{out_dir};
	my $flag=0;
	if(-d $odir){
		opendir (DIR, $odir) or die $!;
		if(scalar(grep{/.*\.RData$/ && -f "$odir/$_"}readdir(DIR))>0){
			$flag=1;
		}
		closedir(DIR);
	}
	return $flag;
}

1;



