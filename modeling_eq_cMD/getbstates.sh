#!/bin/bash
#pulls a frame every 1 ns from the last 5 ns of prod.nc

for i in $(seq	1 5); do
for n in $(seq 1 5); do
COMMAND="         parm bound_new.prmtop\n"
COMMAND="$COMMAND trajin $i/prod.nc\n"
COMMAND="$COMMAND trajout $i/$n.ncrst onlyframes $((100-$n))\n"
COMMAND="$COMMAND go\n"
echo -e $COMMAND | cpptraj
done
done
