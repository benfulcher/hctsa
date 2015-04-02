#!/bin/bash

# Make sure the latest example source file is compiled.
javac -classpath "../../infodynamics.jar" "infodynamics/demos/Example8TeContinuousDataByBinning.java"

# Run the example:
java -classpath ".:../../infodynamics.jar" infodynamics.demos.Example8TeContinuousDataByBinning

