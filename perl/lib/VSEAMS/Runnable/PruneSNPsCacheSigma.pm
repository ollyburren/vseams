package VSEAMS::Runnable::PruneSNPsCacheSigma;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;
use File::Find;

use constant _inputs_expected => {
	in_dir => 'file',
	sigma_cache_dir => 'file',
	r2_threshold => 'number',
	out_dir => 'file',
	test => 'boolean',
	max_no_jobs => 'number'
};

use constant _defaults => {
	test => 1,
	max_no_jobs => 2500 
};
	

use constant _script => 'prune_and_cache_sigma.R';


sub run{
	my $self=shift;
	if(-d $self->inputs->{out_dir}){
		File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($self->inputs->{out_dir});
	}
	if(!-d $self->inputs->{sigma_cache_dir}){
		Macd::Utils::mkpath_wrapper($self->inputs->{sigma_cache_dir});
	}
	#my $test_flag = $self->inputs->{test};
	## if not set we need it to be -1 rather than undef or 0
	#$test_flag = $test_flag?$test_flag:-1;
	my @excl_param=qw/in_dir max_no_jobs/;
	my $pattern='RData$';
	my @include_files;
	find(sub{
			#print $Find::File::name."\n";
			return unless /$pattern/;
			push @include_files,$File::Find::name
	},$self->inputs->{in_dir});
	my $jperbatch = $self->jobs_per_batch(scalar(@include_files));
	for(my $x=0;$x<scalar(@include_files);$x=$x+$jperbatch){
		my %inc_param = (
				support_files=> join(",",@include_files[$x..($x+$jperbatch-1)]),
				#test=>$test_flag?'TRUE':'FALSE'
				test=>'FALSE'
				);
			my $cmd = $self->build_Rscript_cmd(\@excl_param,\%inc_param);
			my $job = Macd::Step::Job->new(command=>$cmd);
			$self->step->add_job($job);
	}
	if($self->step->execute()){
		$self->step->wait_on_complete();
		if($self->check_R_jobs){
			#"Done check ".($self->outputs)->{out_dir}."\n";
			return 1;
		}
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

sub jobs_per_batch{
	my $self = shift;
	my $total_files = shift;
	## not accurate but makes sure that we do
	## the jobs that we need 
	return int(($total_files/$self->inputs->{max_no_jobs})+1);
}

1;



