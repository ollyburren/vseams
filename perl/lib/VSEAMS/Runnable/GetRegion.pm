package VSEAMS::Runnable::GetRegion;

use strict;
use base('VSEAMS::RunnableI');
use Data::Dumper;

use constant _inputs_expected => {
	input_file => 'file',
	tss_extension => 'number',
	mart_host => 'string',
	output_file => 'file'};
	
use constant _defaults => {
		tss_extension => 20000,
		mart_host => 'www.ensembl.org'
};

use constant _script => 'get_regions.gr.R';


sub run{
	my $self=shift;
	## no setup 	die Dumper($step);required just overwrite file
	my $cmd = $self->build_Rscript_cmd;
	my $job = Macd::Step::Job->new(command=>$cmd);
	$self->step->add_job($job);
	if($self->step->execute()){
		$self->step->wait_on_complete();
		if($self->check_R_jobs){
			#"Done check ".($self->outputs)->{output_file}."\n";
			return 1;
		}
	}
	return 0;
}

sub _skip_message{
	my $self=shift;
	return ref($self).":The output file ".$self->inputs->{output_file}." already exists. Would you like to skip this step ?";
}

sub run_conditions{
	my $self=shift;
	if(-e $self->inputs->{output_file}){
		return 1;
	}
	return 0;
}

1;



