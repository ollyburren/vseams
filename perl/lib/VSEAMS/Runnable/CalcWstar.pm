package VSEAMS::Runnable::CalcWstar;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;

use constant _inputs_expected => {
	perm_dir => 'file',
	snp_manifest=> 'file',
	out_dir=>'file',
	test_set => 'string',
	control_set => 'string',
	n_perms => 'number',
	chunk_size=> 'number',
};
	
use constant _defaults => {
};

use constant _script => 'compute_wstar.R';


sub run{
	my $self=shift;
  if(-d $self->inputs->{out_dir}){
		File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($self->inputs->{out_dir});
	}         
	my $max_run_count = $self->inputs->{n_perms}/$self->inputs->{chunk_size};
	for(my $rc=1;$rc <= $max_run_count;$rc++){
		my $cmd = $self->build_Rscript_cmd(['n_perms'],{run_count=>$rc});
		my $job = Macd::Step::Job->new(command=>$cmd);
		$self->step->add_job($job);
	}
	## we want this to be non blocking let the caller handle execution management
	if($self->step->execute()){
		return 1;
		## can get a list of job associated with a step by
		## my @jids = keys %{$step->submitted_jobs};
		## return \@jids;
		#$self->step->wait_on_complete();
		
		#if($self->check_R_jobs){
		#	#print STDERR "Done check ".$self->outputs->{out_dir}."\n";
		#	return 1;
		#}
	}
	return 0;
}

sub max_run_count{
	my $self = shift;
	$self->inputs->{n_perms}/$self->inputs->{chunk_size};
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



