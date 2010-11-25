#
# Domain-decomposition options processing units.
#

#
# Receives string and converts it to partitioning sequence.
# ARGS: key, string.
#
sub partitioning_getSequence ($$) {
  my ($key, $s, $INT) = ($_[0], $_[1], '\d+');
  my @seq = ();
  if ($s !~ /[^xyz\d,\s]/i) {
    while ($s =~ /([xyz]$INT)/igc) {
      push @seq, "$1";
    }
  }

  if ($#seq >= 0) {
    my $RES = join ("\n", ($#seq + 1), @seq);		# Forms sequence string.
    $RES =~ tr/[A-Z]/[a-z]/;				# Makes string low-case.
    setPair('sequence', $key, join ("\n", length ($RES), $RES));
  }
}

#
# Processes array of options looking for domain decomposition requests.
#
sub process_partitioning ($) {
  my @strings = filter ($_[0], '--partition:.*slice{([^}]*)}');
  map { push @strings, $_ } filter ($_[0], '-p:.*slice{([^}]*)}');
  map { getVector ('decomposition->cut across', $_, [0, 0, 0]) } @strings;

  @strings = filter ($_[0], '--partition:.*noise{([^}]*)}');
  map { push @strings, $_ } filter ($_[0], '-p:.*noise{([^}]*)}');
  map { getVector ('decomposition->rand shifts', $_, [0, 0, 0]) } @strings;

  @strings = filter ($_[0], '--partition:.*direct{([^}]*)}');
  map { push @strings, $_ } filter ($_[0], '-p:.*direct{([^}]*)}');
  map { partitioning_getSequence ('decomposition->direct', $_) } @strings;
}

#
# Prints txt report for any key-value pair generated for domain decomposition.
#
sub report_partitioning ($) {
  my $pairs = $_[0];
  
  if (exists ($pairs->{"decomposition->direct"})) {
    print "User decomposition:\n";
    my $cpuN = 1;
    while ($pairs->{"decomposition->direct"}->{'value'} =~ /([xyz])(\d+)/gc) {
      print "  o across $1-axis on $2 pieces\n";
      $cpuN *= $2;
    }
    print "  NOTE: $cpuN cpus assumed!\n\n";
  }

  if (exists ($pairs->{"decomposition->cut across"})) {
    if ($pairs->{"decomposition->cut across"}->{'value'} =~ /(\d+)\n(\d+)\n(\d+)/) {
      print "Decomposition slicing axises:\n  ",
        ($1 != "0") ? "X" : "", ($2 != "0") ? "Y" : "",  ($3 != "0") ? "Z" : "", "\n\n";
    }
  }

  if (exists ($pairs->{"decomposition->rand shifts"})) {
    if ($pairs->{"decomposition->rand shifts"}->{'value'} =~ /(\d+)\n(\d+)\n(\d+)/) {
      print "Decomposition debugging noise:\n  ($1, $2, $3)\n\n";
    }
  }
}

sub help_partitioning ($) {
  print "Command line option / domain decomposition:\n";
  print "  o [-p:/--partition:]slice{x1, z1}      - vector, set 1 to unlock axis\n";
  print "  o [-p:/--partition:]noise{z4, x0, y5}  - vector, displacement amplitudes\n";
  print "  o [-p:/--partition:]direct{x[Nx], y[Ny]...} - sequence, vec-component like\n";
  print "    elements, space or comma separated\n";
}
