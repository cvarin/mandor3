#
# OpenDX options processing units.
#

#
# Processes array of options looking for domain decomposition requests.
#
sub process_dx ($) {
  my @strings = filter ($_[0], '--dx:.*res{([^}]*)}');
  map { getVector ('dx:resolution', $_, [1, 1, 1])} @strings;

  @strings = filter ($_[0], '--dx:.*fseek\(([^\)]*)\)');
  map { getInt ('dx:fseek', $_)} @strings;
}

#
# Prints txt report for any key-value pair generated for domain decomposition.
#
sub report_dx ($) {
  my $pairs = $_[0];
  
  if (exists ($pairs->{"dx:resolution"})) {
    if ($pairs->{"dx:resolution"}->{'value'} =~ /(\d+)\n(\d+)\n(\d+)/) {
      print "OpenDX visualizer resolution:\n  ($1, $2, $3)\n\n";
    }
  }

  if (exists ($pairs->{"dx:fseek"})) {
    if ($pairs->{"dx:fseek"}->{'value'} =~ /(\d+)/) {
      print "OpenDX visualizer fseek optimization threshold:\n  $1\n\n";
    }
  }
}

sub help_dx () {
  print "Command line option / OpenDX:\n";
  print "  o --dx:res{x3, z4}\n";
  print "  o --dx:fseek(5)\n";
}
