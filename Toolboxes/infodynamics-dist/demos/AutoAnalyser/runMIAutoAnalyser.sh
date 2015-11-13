#!/bin/bash

# Make sure the latest example source file is compiled.
javac -classpath "../java:../../infodynamics.jar" "../java/infodynamics/demos/autoanalysis/AutoAnalyserMI.java"

# Run the example:
java -classpath "../java:../../infodynamics.jar" infodynamics.demos.autoanalysis.AutoAnalyserMI

