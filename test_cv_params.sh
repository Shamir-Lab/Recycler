#!/bin/sh
set -e # exit on first error
for ref in 100 200 400 800; do
	for cv in 0.25; do
		make -f map_and_recycle.mk clean
		make -f nucmer.sim.mk clean CV="$cv" REF_CNT="$ref"
		make -f map_and_recycle.mk all
		make -f nucmer.sim.mk all CV="$cv" REF_CNT="$ref"
	done
done
