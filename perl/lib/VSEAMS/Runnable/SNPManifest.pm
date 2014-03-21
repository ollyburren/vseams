package VSEAMS::Runnable::SNPManifest;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;
use File::Find;

use constant _inputs_expected => {
	out_dir => 'file',                    
	snp_dir => 'file',
	test_set => 'string',
	control_set => 'string'
};

use constant _defaults => {
};                          
	

use constant _script => 'create_snp_manifest.R';

sub run{
	my $self=shift;
  if(-d $self->inputs->{out_dir}){
		File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($self->inputs->{out_dir});
	}
	## set up the output file
	$self->outputs->{out_file}=$self->{out_file};
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

sub out_file{
	my $self = shift;
	my $ofile = $self->outputs->{out_dir}."/".$self->inputs->{test_set}.".".$self->inputs->{control_set}.".RData";
	return $ofile;
}

sub _skip_message{
	my $self=shift;
	return "[".ref($self)."]: The output path ".$self->inputs->{out_dir}." already exists and contains files. Would you like to skip this step ?";
}

sub run_conditions{
	my $self=shift;
	if(-e $self->out_file){
		return 1;
	}else{
		return 0;
	}
}

1;



