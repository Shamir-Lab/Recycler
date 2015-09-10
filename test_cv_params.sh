#!/bin/sh
set -e # exit on first error
for ref in 100 200 400 800 1600; do
	for cv in 0.25 0.375 0.5; do
		make -f nucmer.sim.mk clean CV="$cv" REF_CNT="$ref"
		make -f nucmer.sim.mk all CV="$cv" REF_CNT="$ref"
	done
done
