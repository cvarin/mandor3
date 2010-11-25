#!/usr/bin/perl -w
use strict;

#
# Design: do not check for many mistakes, only let user read nicely
#         formatted summary and let him eyeball how stuff is recognized
#         so he/she can find where and why mistake occurs.
#

my $debugLevel = 0;		# verbose level (0 - quiet): prints to stderr.
my %pairs = ();			# Hash of all key-value pairs.

my $_INT_ = '[-\+]?\d+';	# Mask for input int value.

# Extract options for filter itself.
my ($makeReport, $makeHelp, $rcFile) = (0, 0);
for my $arg (@ARGV) {			# Gets debug level.
  if ($arg =~ m/\+\+debug:(\d+)/i)  { $debugLevel = $1; };
  if ($arg =~ m/\+\+debug:help/i)   { $makeHelp = 1 };
  if ($arg =~ m/\+\+debug:report/i) { $makeReport = 1 };
  if ($arg =~ m/\+\+file:(\S*)/i)   { $rcFile = $1 };
}

my ($baseDir) = ($0 =~ /^(.*)\/([^\/]+)$/);	# Imports modules for particular options.
do "$baseDir/import_dx.pl";
do "$baseDir/import_fft.pl";
do "$baseDir/import_partitioning.pl";


#########################################################################################

  #
  # Prints message to STDERR if $debugLevel permits.
  #
  sub debugMsg($$) {
    my $msg = $_[1];
    if ($debugLevel >= $_[0]){
      $msg =~ s/(^)(.*$)/$1deb($_[0]): $2/gm;	# Adds deb(lev) prefix to lines.
      print STDERR $msg;
    }
  }

  #
  # Pushes NEW (type, key, value) hash to the hash of pairs.
  #
  sub setPair ($$$) {
    my $pair = {'type' => $_[0], 'key' => $_[1], 'value' => $_[2]};
    if (!exists ($pairs{$_[1]})) {
      debugMsg (1, "NEW PAIR\n  KEY: '$_[1]'\n  VALUE: '$_[2]'.\n");
      $pairs{$_[1]} = $pair;
    }
  }

  #
  # Receives string and converts it to integer.
  # ARGS: key, string.
  #
  sub getInt ($$) {
    my ($key, $s) = ($_[0], $_[1]);
    if ($s =~ /^\s*($_INT_)\s*$/gc) {			# Reads number.
      setPair ('int', $key, $1);
    }
  }
  
  #
  # Receives string and converts it to float.
  # ARGS: key, string.
  #
  sub getDouble ($$) {
    my ($key, $s) = ($_[0], $_[1]);
    if ($s =~ /^\s*([-\+]?\d*\.?\d*(e[-\+]?\d*.?\d*)?)\s*$/gc) {
      setPair ('double', $key, $1);
    }
  }

  #
  # Receives string and converts it to the vector.
  # ARGS: key, string, default.
  #
  sub getVector ($$$) {
    my ($key, $s, $def) = ($_[0], $_[1], $_[2]);
    my %v = ();						# Accumulator hash.
    my @res = ($def->[0], $def->[1], $def->[2]);	# Default vector.
    debugMsg (2, "  VEC: '$key'\n  STRING: '$s'\n");
    $s =~ tr/[A-Z]/[a-z]/;				# Turns low-case.
    if ($s !~ /[^xyz\d,\s\-\+]/) {			# Checks bad symbols.
      while ($s =~ /([xyz])($_INT_)/gc) {		# Reads components.
        (not exists ($v{$1})) and (defined $2 and length ($2)) and
          ($v{$1} = $2);
      }
      (defined $v{'x'}) and $res[0] = $v{'x'};		# Overrides defaults.
      (defined $v{'y'}) and $res[1] = $v{'y'};
      (defined $v{'z'}) and $res[2] = $v{'z'};
      setPair ('vec3i_t', $key, join ("\n", @res));
    }
  }
  
  #
  # Takes options and returns array of the data strings which masks deliver.
  # ARGS: options array, mask
  #
  sub filter($$) {
    my ($arr, $mask) = ($_[0], $_[1]);
    my @strings = grep (/$mask/, @{$arr});		# Filters valid options.
    map { $_ =~ /$mask/; $_ = $1 } @strings;		# Unwraps data.
    return @strings;					# Returns vector strings.
  }

  #
  # Reads rc file and returnes array of options.
  #
  sub prepareFile ($) {
    open (INFILE, $_[0]) or die "Cannot open file '$_[0]': $!\n";
    my $file = do { local $/; <INFILE> };	# Removes line separator and reads file.
    close (INFILE);
  
    debugMsg (3, ">>> INPUT FILE (original '$_[0]'):\n" . $file . "<<< EOM\n");
    $file =~ s/\s*#.*$//gm;				# Strips '#' comments.
    $file =~ s/^\s*//gm;				# Removes empty lines.
    $file =~ s/(.*)%%(.*)/$1\n/gs;			# Removes commented tail.
  
    debugMsg (2, ">>> INPUT FILE (bare data):\n" . $file . "<<< EOM\n");
    return split (/\s*\n\s*/, $file);
  }

  #
  # Processes array of options.
  #
  sub process ($) {
    process_partitioning ($_[0]);
    process_dx ($_[0]);
    process_fft ($_[0]);
  
    my @strings = filter ($_[0], '--VSP:.*weight\(([^\)]*)\)');
    map { getDouble ('VSP:weight', $_)} @strings;
  }

#########################################################################################

  
process (\@ARGV);				# Inputs command line.

my @fileOptions = ();
if (defined $rcFile and -s $rcFile) {		# Inputs rc file.
  debugMsg (1, "Using resource file '$rcFile'\n");
  @fileOptions = prepareFile ($rcFile);
  process (\@fileOptions);
} else {
  undef $rcFile;
}

for my $key (keys % pairs) {			# Exports keys.
  print ">>> DATA BEGIN\n";
  print "$key\n";
  print $pairs{$key}->{'type'}, "\n";
  print $pairs{$key}->{'value'}, "\n";
  print "<<< DATA END\n\n";
}
print "BYE\n";

if ($makeReport == 0 and $makeHelp == 0) { 	# Quits if no report asked for.
  exit (0);
}

select (STDERR);

if ($makeReport) {				# Sends report to stderr.
  print "\nCommand line options:\n  ", join (" ", @ARGV), "\n";
  if (defined $rcFile) {
    print "\nFile:\n  ", $rcFile, "\n";
    print "\nFile options:\n  ", join (" ", @fileOptions), "\n";
  }
  print "\n\n";
  
  report_partitioning (\%pairs);
  report_dx (\%pairs);
  report_fft (\%pairs);
}

if ($makeHelp) {				# Sends help to stderr.
  print "Command line option hints:\n";
  print "  o ++debug:0, .., ++debug:3\n";
  print "  o ++debug:report\n";
  print "  o ++debug:help\n";

  help_partitioning ();
  help_dx ();
  help_fft ();

  print "\nIf output looks unexpected:\n";
  print "  o use ++debug:1, .., ++debug:3 key to see how input is recognized\n";
  print "  o if something overrides your settings you can see why it happends\n";
  print "  o if parser is wrong - report ('dromanov\@phys.ualberta.ca')\n\n";
}

__END__


  Unix is user-friendly - it\'s just choosy about who its friends are.
  --Anonymous

INTRO:

That is input filter to parse all options in the command line. First script
deals with command line and than it opens file and reads options from there.
File supports comments (with '#'), empty lines (for formatting) and tail
comment (after '%%').

Script extracts all data and prints parsed file (one token per line).
Main stages of parsing are:

1. Strips all comments.
2. Removes empty lines and tail.
3. Processes command line.
4. Processes optional file.

NOTES:

1. Use ++debug:level option to see data on different stages (messages go to
   stderr and they are easily separated even when script is used in pipes).
2. Sequence is printed in form of key/type/total length/number of elements/
   and all elements (in low-case to simplify import in C-code).
3. Instead of the bunch of warning and tests I print final summary in human
   readable format using only final data. It is easy to check that script
   understood your wish correctly by simply reading it. If some keys are
   not what is expected use ++debug:1, .., ++debug:3 to see where crap is
   started (double entry in command line or malishous typo or something else).
4. Report any patches and suggestions to dromanov@phys.ualberta.ca.
5. Syntax
   o 'vector' has 3 components given in any order like '{x4, y4}'
     - no spaces between axis char and number
     - comma or space separated
     - no sign
   o 'sequence' has any number of components like '{x10 y3 z4 x34}';
     - vector type item syntax
   o 'int' usually looks like '(10)'.

   USE "" IN SHELL TO AVOID CONFUGION.

  27 / Jan / 2006.
  Dmitry.

----------------------  EXAMPLE OF FILE ----------------------

# -------------------------
# core.out run-time options
# -------------------------

--partition:slice{x1,y1}		# No slicing across z.
--partition:noise{x10, y10, z3}		# Debugging of decomposition algo.
--partition:direct{z4, y2}		# Explicit partitioning sequence.

%%					# End of options' section.

Here comes unstructured comment.
