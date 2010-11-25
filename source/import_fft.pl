#
# OpenDX options processing units.
#

# Mask for double input values.
my $DBL = '[-\+]?\d*\.?\d*(e[-\+]?\d+)?';

#
# Receives string and converts it to the vector.
# ARGS: key, string, default.
#
sub getDoubleVector ($$$) {
  my ($key, $s, $def) = ($_[0], $_[1], $_[2]);
  my %v = ();						# Accumulator hash.
  my @res = ($def->[0], $def->[1], $def->[2]);		# Default vector.
  debugMsg (2, "  VEC: '$key'\n  STRING: '$s'\n");
  $s =~ tr/[A-Z]/[a-z]/;				# Turns low-case.
  if ($s !~ /[^xyz\d,\s\-\+e\.]/) {			# Checks bad symbols.
    while ($s =~ /([xyz])($DBL)/gc) {			# Reads components.
      (not exists ($v{$1})) and (defined $2 and length ($2)) and
        ($v{$1} = $2);
    }
    (defined $v{'x'}) and $res[0] = $v{'x'};		# Overrides defaults.
    (defined $v{'y'}) and $res[1] = $v{'y'};
    (defined $v{'z'}) and $res[2] = $v{'z'};
    setPair ('vec3d_t', $key, join ("\n", @res));
  }
}

#
# Processes array of options looking for domain decomposition requests.
#
sub process_fft ($) {
  my @strings = filter ($_[0], '--fft:.*minR{([^}]*)}');
  map { getVector ('fft->min r', $_, [-100000, -100000, -100000]) } @strings;

  @strings = filter ($_[0], '--fft:.*maxR{([^}]*)}');
  map { getVector ('fft->max r', $_, [100000, 100000, 100000]) } @strings;

  @strings = filter ($_[0], '--fft:.*minK{([^}]*)}');
  map { getVector ('fft->min k', $_, [-100000, -100000, -100000]) } @strings;

  @strings = filter ($_[0], '--fft:.*maxK{([^}]*)}');
  map { getVector ('fft->max k', $_, [100000, 100000, 100000]) } @strings;

  @strings = filter ($_[0], '--fft:.*lambda{([^}]*)}');
  map { getDoubleVector ('fft->lambda', $_, [1, 1, 1]) } @strings;

  @strings = filter ($_[0], '--fft:.*vector{([xyz])}');
  map { 
    setPair ('string', 'fft->component', $_);
    setPair ('string', 'fft->type', 'vector');
  } @strings;

  @strings = filter ($_[0], '--fft:.*scalar');
  map { setPair ('string', 'fft->type', 'scalar') } @strings;
}

#
# Prints txt report for any key-value pair generated for domain decomposition.
#
sub report_fft ($) {
  # Alias and masks for input values.
  my ($pairs, $INT) = ($_[0], '[-\+]?\d+');

  if (exists ($pairs->{'fft->min r'})) {
    if ($pairs->{'fft->min r'}->{'value'} =~ /($INT)\n($INT)\n($INT)/) {
      print "Mesh FFT min R:\n  ($1, $2, $3)\n\n";
    }
  }

  if (exists ($pairs->{'fft->max r'})) {
    if ($pairs->{"fft->max r"}->{'value'} =~ /($INT)\n($INT)\n($INT)/) {
      print "Mesh FFT max R:\n  ($1, $2, $3)\n\n";
    }
  }
  
  if (exists ($pairs->{'fft->min k'})) {
    if ($pairs->{'fft->min k'}->{'value'} =~ /($INT)\n($INT)\n($INT)/) {
      print "Mesh FFT min K:\n  ($1, $2, $3)\n\n";
    }
  }

  if (exists ($pairs->{'fft->max k'})) {
    if ($pairs->{'fft->max k'}->{'value'} =~ /($INT)\n($INT)\n($INT)/) {
      print "Mesh FFT max K:\n  ($1, $2, $3)\n\n";
    }
  }
  
  # $1, $3, $5 because of the definition of the $DBL has one more ()-block.
  if (exists ($pairs->{"fft->lambda"})) {
    if ($pairs->{'fft->lambda'}->{'value'} =~ /($DBL)\n($DBL)\n($DBL)/) {
      print "Reference wavelentgth:\n  ($1, $3, $5) [mesh step]\n\n";
    }
  }
  
  if (exists ($pairs->{'fft->type'})) {
    print "Mesh FFT source type:\n   $pairs->{'fft->type'}->{'value'}";
    if ($pairs->{'fft->type'}->{'value'} eq 'vector') {
      print " ($pairs->{'fft->component'}->{'value'}-component)";
    }
    print "\n\n";
  }
}

sub help_fft () {
  print "Command line option / mesh FFT:\n";
  print "  o --fft:minR{<vector>}    subregion in space (low boundary)\n";
  print "  o --fft:maxR{<vector>}    subregion in space (high boundary)\n";
  print "  o --fft:minK{<vector>}    spectral subregion (low boundary)\n";
  print "  o --fft:maxK{<vector>}    spectral subregion (high boundary)\n";
}
