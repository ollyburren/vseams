package VSEAMS::RunnableI;


## This is an interface for all things runnable

use strict;
use warnings;
use Macd;
use VSEAMS;
use File::Path qw/make_path remove_tree/;
use Data::Dumper;
use File::Basename qw/dirname/;
use Scalar::Util qw(looks_like_number);

## each runnable needs
## a macd_driver
## a log_dir
## a hash of inputs
## a list of expected outputs
## a run subroutine - this will alter depending

sub new{
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = {@_};
	bless($self,$class);
	$self->init();
	return $self;
}

sub init{
	my $self = shift;
	$self->_set_defaults;
	$self->script;
	$self->inputs({}) if ! $self->inputs;
	## for outputs we check to see if dir is available
	## if it is not then we make it 
	## other areas check to see if files etc. exist.
	$self->_check_outputs;
	$self->_check_inputs;
}

sub _inputs_expected{
	my $self=shift;
	die ref($self)." doesn't seem to have any inputs\n";
}

sub _script{
	my $self = shift;
	if(@_)	{$self->{_script} = shift}
	return $self->{_script}
}

sub log_dir {
	my $self = shift;
	if(@_)	{$self->{log_dir} = shift}
	return $self->{log_dir}
}

sub debug_flag {
	my $self = shift;
	if(@_)	{$self->{debug_flag} = shift}
	return $self->{debug_flag}
}

sub skipped {
	my $self = shift;
	if(@_)	{$self->{skipped} = shift}
	return $self->{skipped}
}

## this is the default we override in subclasses
## to say something more helpful !

sub _skip_message {
	my $self = shift;
	return "Looks as if you have already run step ".ref($self).". Would you like to skip it ?";
}

sub non_interactive{
	return 0;
}

#sub debug{
#	my $self = shift;
#	Macd::Utils::debug(shift,$self->debug_flag);
#	return undef;
#}

sub previous_step {
	my $self = shift;
	if(@_)	{$self->{previous_step} = shift}
	return $self->{previous_step}
}

sub env_settings {
	my $self = shift;
	if(@_)	{$self->{env_settings} = shift}
	return $self->{env_settings}
}

sub macd_driver {
	my $self = shift;
	if(@_)	{ $self->{macd_driver} = shift}
	return $self->{macd_driver}
}


sub inputs {
	my $self = shift;
	if(@_)	{$self->{inputs} = shift}
	return $self->{inputs}
}

sub outputs {
	my $self = shift;
	if(@_)	{$self->{outputs} = shift}
	return $self->{outputs}
}
	

sub step {
	my $self = shift;
	if(!$self->{step}){
		my $step = Macd::Step->new(
			logdir=>$self->log_dir,
			driver=>$self->macd_driver,
			env_settings=>$self->env_settings
		);
		$self->{step} = $step;
	}
	return $self->{step}
}

sub _set_defaults{
	my $self=shift;
	if($self->_defaults){
		my %defaults = %{$self->_defaults};
		my %got = %{$self->inputs};
		foreach my $d(keys %defaults){
			$got{$d} = $defaults{$d} unless $got{$d};
		}
		$self->inputs(\%got);
	}
}

sub _defaults {
	my $self = shift;
	return {};
}

sub _check_inputs{
	my $self = shift;
	my %expected = %{$self->_inputs_expected};
	my %got = %{$self->inputs};
	my $message;
	my $class = ref($self);
	foreach my $i(keys %expected){
		if(!$got{$i}){
			$message.="$class:Expected parameter $i not set\n";
			next;
		}
		#print $expected{$i}." $got{$i}\n";
		if($expected{$i} eq 'file'){
			next if -e $got{$i};
			my $dirname = dirname($got{$i});
			unless($dirname eq '/'){
				next if -d $dirname;
			}
			$message.="$class:Parameter: $i:$got{$i}, is not an accessible file or dir\n";
		}elsif($expected{$i} eq 'number'){
			next if looks_like_number($got{$i});
			$message.="$class:Parameter: $i:$got{$i} doesn't look like a number\n";
		}elsif($expected{$i} eq 'string'){
			next if $got{$i} =~/^\S+$/;
			$message.="$class:Parameter: $i:$got{$i} doesn't look like a string\n";
		}elsif($expected{$i} eq 'boolean'){
			next if $got{$i} =~/^(true|false|[\-0-9]+)$/i;
			$message.="$class:Parameter: $i:$got{$i} doesn't look like a boolean\n";
		}elsif($expected{$i} eq 'hashref'){
			next if ref($got{$i}) eq 'HASH';
			$message.="$class:Parameter: $i:$got{$i} doesn't look like a hashref\n";
		}else{
			$message.="$class:Parameter: $i:$got{$i} $expected{$i} not recognised\n";
		}
	}
	if($message){
		die $message;
	}
	return 1;
}

sub _check_outputs{
	my $self=shift;
	my %outputs = %{$self->outputs};
	foreach my $k(keys %outputs){
		if($k=~m/\_file$/){
			Macd::Utils::mkpath_wrapper(File::Basename::dirname($outputs{$k}));
			$self->{inputs}->{$k}=$outputs{$k};
		}elsif($k=~m/\_dir$/){
			Macd::Utils::mkpath_wrapper($outputs{$k});
			$self->{inputs}->{$k}=$outputs{$k};
		}
	}
	return 1;
}

sub script {
	my $self=shift;
	if(!$self->{script}){
		use FindBin qw/$Bin/;
		my $script = $Bin.'/R/'.$self->_script;
		unless(-e $script){
			die "Cannot find $script please check configuration\n";
		}
		$self->{script} = $script;
	}
	return $self->{script}
}

sub build_Rscript_cmd {
	my $self = shift;
	my ($exclude,$include) =@_;
	#die(ref($include));
	my $cmd =  "Rscript --vanilla ".$self->script;
	my %inputs = %{$self->inputs};
	## add in outputs
	my %outputs = %{$self->outputs};
	$inputs{keys %outputs} = values %outputs;
	my %expected = %{$self->_inputs_expected};
	## here we allow caller to explicitly remove inputs
	## useful when inputs are a directory for which every
	## file needs to be processed 
	if(ref($exclude) =~ 'ARRAY'){
		delete $expected{$_} foreach @$exclude;
	}
	## here we allow caller to explicitly add command line
	## arguments
	if(ref($include) eq 'HASH'){
		foreach my $i(keys %$include){
			$expected{$i}='string' unless  $expected{$i};
		}
		@inputs{keys %$include}=values %$include;
	}
	my @params;
	foreach my $i(keys %expected){
		if($expected{$i} =~ 'number|boolean'){
			push @params,"$i=$inputs{$i}";
		}else{
			push @params, "$i=\\'$inputs{$i}\\'";
		}
	}
	return join(" ",$cmd,@params);
}

## check to see that all jobs have passed 
## assumes that all Rscripts end with print ("Success")

sub check_R_jobs{
	my $self = shift;
	## here we check the contents the log files for any files 
	## for any failures
	my $log_dir=$self->log_dir;
	my $err = `for i in \`\\ls $log_dir*.log\`; do tail -1 \$i; done | grep -v Success`;
	if(length($err)>0){
		return 0;
	}else{
		return 1;
	}
}
		

sub run {
		my $self = shift;
		print STDERR "run needs to be implemented in".ref($self)."\n";
		return 0;
}

sub run_conditions{
	my $self = shift;
	print STDERR "run_conditions should be overridden".ref($self)."\n";
	return 0;
}

sub run_step{
	my $self = shift;
	#if(!$self->isa('VSEAMS::Runnable::CalcZ')){
	#	$self->skipped(1);
	#	return $self;
	#}
	my $force = 0;
	if(my $prev = $self->previous_step){
		$force = $prev->skipped?0:1;
	}
	#die($self->previous_step->skipped);
	if(!$force && $self->run_conditions){
		## ask here if they want to skip (if that is defined)
		## if they do go ahead
		if($self->skip_ask){
			print STDERR "Skipping ".ref($self)."\n";
			$self->skipped(1);
			return $self;
		}
	}
	## here it's nice to remove the logfiles
	if(-d $self->log_dir){
		File::Path::remove_tree($self->log_dir,{keep_root => 1,result=> \my $del_dirs});
	}
	unless($self->run){
		die "Encountered an issue with ".ref($self)." check ".$self->log_dir."\n";
	} 
	return $self;
}

sub skip_ask{
	my $self=shift;
	## here we should get a global setting 
	if($self->non_interactive){
		return 1;
	}else{
		return yesno($self->_skip_message);
	}
}


		

1;

