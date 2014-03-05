package VSEAMS::Runnable::CalcZ;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;

#tests <= (
#	test.ctrl => [ ds1,ds2] ..
#)

use constant _inputs_expected => {
	out_dir => 'file',
	tests => 'hashref'
};

use constant _defaults => {
	test => 1
};
	

use constant _script => 'compute_Z.R';


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
	my @excl_param=qw/tests out_dir/;
	foreach my $t(keys %{$self->inputs->{tests}}){
		my $wparam = join(",",@{$self->inputs->{tests}->{$t}});
		my %inc_param = (
				w_param_files => $wparam,
				out_file => $self->inputs->{out_dir}."/$t.results.RData"
		);
		my $cmd = $self->build_Rscript_cmd(\@excl_param,\%inc_param);
		print STDERR $cmd."\n";
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



