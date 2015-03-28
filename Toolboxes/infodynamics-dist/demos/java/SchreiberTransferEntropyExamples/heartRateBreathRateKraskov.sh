#!/bin/bash
#
# Inputs:
#  $1 - kHistory - destination embedding length. Defaults to 1.
#  $2 - lHistory - source embedding length. Defaults to 1.
#  $3 - knns - a scalar specifying a single, or X:X format specifying a range, or a comma separated list specifying multiple, value(s) of K nearest neighbours to evaluate TE (Kraskov) with. Defaults to 4.
#  $4 - numSurrogates - a scalar specifying the number of surrogates to evaluate TE from null distribution. Defaults to 0 (i.e. don't evaluate surrogates)


# The examples are intended to be compiled and run from the above directory:
cd ..

# Make sure the latest example source file is compiled.
javac -classpath "../../infodynamics.jar" "infodynamics/demos/schreiberTransferEntropyExamples/HeartBreathRateKraskovRunner.java"


# Run the example:
java -classpath ".:../../infodynamics.jar" infodynamics.demos.schreiberTransferEntropyExamples.HeartBreathRateKraskovRunner $1 $2 $3 $4

cd SchreiberTransferEntropyExamples

