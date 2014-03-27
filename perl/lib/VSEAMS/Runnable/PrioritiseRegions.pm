package VSEAMS::Runnable::PrioritiseRegions;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;
use File::Find;

use constant _inputs_expected => {
	out_dir => 'file',
	test_name => 'file',
	support_file => 'file',
	input_file => 'file',
	snp_manifest => 'file',
	sigma_cache_dir => 'file',
	#sample_sigma_dir => 'file',
	n_perms => 'number',
	index_file => 'file',
};

#use constant _defaults => {
#	test => 1
#};
	

use constant _script => 'prioritise_regions.R';

my $DEFAULT_CHUNK_SIZE = 100;

sub run{
	my $self=shift;
	if(-d $self->outputs->{out_dir}){
		File::Path::remove_tree($self->inputs->{out_dir},{keep_root => 1,result=> \my $del_dirs});
	}else{
		Macd::Utils::mkpath_wrapper($self->inputs->{out_dir});
	}
	my $test = $self->inputs->{test_name};
	$self->inputs->{test} = "SET.$test";
	
	my @excl_param=qw/input_file test_name/;
	## get a list of target genes
	open(RTAB,$self->inputs->{input_file}) || die "Cannot open input_file".$self->inputs->{input_file};
	my %regions;
	while(<RTAB>){
		chomp;
		my ($id,$set)=split("\t",$_);
		$regions{$id}++ if $set=~/$test$/;
	}
	close(RTAB);
	foreach my $r(keys %regions){
		## we do all perms at the same time
		my %inc_param = (region_id => $r, test=>$self->inputs->{test});
		my $cmd = $self->build_Rscript_cmd(\@excl_param,\%inc_param);
		print $cmd."\n";
		print STDERR $self->debug($cmd);
		my $job = Macd::Step::Job->new(command=>$cmd);
		$self->step->add_job($job);
		last;
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
		if(scalar(readdir(DIR))>0){
			$flag=1;
		}
		closedir(DIR);
	}
	return $flag;
}

1;



