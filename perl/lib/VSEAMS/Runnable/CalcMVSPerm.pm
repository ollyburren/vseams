package VSEAMS::Runnable::CalcMVSPerm;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;
use File::Find;

use constant _inputs_expected => {
	out_dir => 'file',
	in_dir => 'file',
	sigma_cache_dir => 'file',
	sample_sigma_dir => 'file',
	n_perms => 'number',
	test => 'boolean'
};

use constant _defaults => {
	test => 1
};
	

use constant _script => 'calc_mvs_perm.R';

my $DEFAULT_CHUNK_SIZE = 100;

sub run{
	my $self=shift;
	if(-d $self->inputs->{out_dir}){
		File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($self->inputs->{out_dir});
	}         
	my $test_flag = $self->inputs->{test};
	#die "$test_flag\n";
	## if not set we need it to be -1 rather than undef or 0
	my @excl_param=qw/in_dir/;
	## get a unique list of chromosomes
	my $chr_cmd = '\ls '.$self->inputs->{in_dir}.'| cut -d_ -f1 | sort | uniq';
	my @chr = split("\n",qx/$chr_cmd/);
	foreach my $chr(@chr){
		my $total_perms = 0;
		my $run_count = 1;
		while($total_perms < $self->inputs->{n_perms}){
			my %inc_param = (
					chr_name=>$chr,
					snp_dir=> $self->inputs->{in_dir},
					test=>$test_flag?'TRUE':'FALSE',
					n_perms=>$self->chunk_size,
					run_count=>$run_count
				);
			my $cmd = $self->build_Rscript_cmd(\@excl_param,\%inc_param);
			print STDERR $self->debug($cmd);
			my $job = Macd::Step::Job->new(command=>$cmd);
			$self->step->add_job($job);
			$total_perms += $self->chunk_size;
			$run_count++;
		}
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

sub chunk_size {
	my $self = shift;
	if(@_)	{$self->{chunk_size} = shift}
	return $self->{chunk_size} || $DEFAULT_CHUNK_SIZE;
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



