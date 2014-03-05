package VSEAMS::Runnable::DownSampleSigma;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;
use File::Find;
use File::Basename;

use constant _inputs_expected => {
	out_dir => 'file',
	downsample_file => 'file',
	sigma_cache_dir => 'file',
	test => 'boolean',
	max_snps_per_sigma_block => 'number',
	snp_dir=>'file'
};

use constant _defaults => {
	test => 1
};
	

use constant _script => 'downsample_sigma.R';

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
	my @excl_param=qw/downsample_file snp_dir/;
	## get a unique list of chromosomes
	#my $chr_cmd = '\ls '.$self->inputs->{in_dir}.'| cut -d_ -f1 | sort | uniq';
	#my @chr = split("\n",qx/$chr_cmd/);
	#foreach my $chr(@chr){
	open(DS,$self->inputs->{downsample_file}) || die "Cannot open ".$self->inputs->{downsample_file}."\n";
	my @files = <DS>;
	close(DS);
	foreach my $file(@files){
		chomp($file);
		## grab the basename
		my $snp_file = $self->inputs->{snp_dir}.basename($file);
		my %inc_param = (
					snp_file=>$snp_file,
		);
		my $cmd = $self->build_Rscript_cmd(\@excl_param,\%inc_param);
		print STDERR $self->debug($cmd);
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

1;



