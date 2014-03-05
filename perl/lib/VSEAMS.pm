package VSEAMS;
require Exporter;
our @ISA = ("Exporter");

use File::Path qw/make_path remove_tree/;
use File::Basename qw/basename/;
use File::Temp qw/tempfile/;
use POSIX qw/strftime ceil/;
use Term::ReadKey;
use Data::Dumper;

my $TEMPLATE_PREFIX = 'vseams';
my $DEBUG=0;

our @EXPORT = qw/debug get_tmp_file mkpath_wrapper yesno basename remove_tree ceil/;

## Base module that we can use to store shared routines and libraries
                                                                          
sub debug{
	my $msg=shift;
	my $debug=shift;
	my $ts = strftime("[%m/%d/%Y %H:%M:%S]", localtime);    
	print join("\t",$ts,$msg)."\n" if $debug;
}


sub get_tmp_file{
	my $dir = shift;
	my $template = "${TEMPLATE_PREFIX}_XXXXX";
	return tempfile($template, DIR => $dir) if $dir;
	return tempfile($template, TMPDIR => 1);
}

sub mkpath_wrapper{
	my @dirs=@_;
	make_path(@dirs,{verbose=>1,error => \my $err});
	if(@$err){
		for my $diag (@$err) {
			my ($file, $message) = %$diag;
			if ($file eq '') {
				debug("general error: $message",$DEBUG);
			}else {
				debug("problem creating $file: $message",$DEBUG);
			}
		}
		return 0;
	}
	return 1;
}

=head2 yesno

 NAME: yesno
 ARGS: STRING - represents question to ask
 FUNCTION: Writes a question to the terminal and waits for a y or n keystroke
 RETURNS : Boolean true or false

=cut

sub yesno{
	my($question)=shift;
	if(!$question){
		warn "[WARNING] No question: returning false\n";
		return 0;
	}
	ReadMode('cbreak');
	print "$question?\n";
	my $dyes = "Please enter y or n. Use Ctrl-C to quit.\n";
	while(1){
		my $char = ReadKey(0);
		if($char!~/[yn]/i){
			print $dyes;
		}else{
			ReadMode('normal');
			return $char =~/y/i?1:0;
		}
	}
}

