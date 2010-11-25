#!/bin/sh

echo "It is 'total field/scattered field' interface test. Main idea is "
echo "to record field in the region bounded by interface and than to"
echo "play it. To perform test just execute record session, save some"
echo "data to compare (like output of the probe diagnostic), than run"
echo "play session and observe absolutely the same field in the bounded"
echo "region with no field outside the interface. To release field you"
echo "can use cut plane to remove part of the interface - EM field will "
echo "leave the domain through this hole."
echo;
echo "Example output is written in res/ directory."
echo "Use test_rec.sh and test_play.sh to run rec/play session."
echo;

