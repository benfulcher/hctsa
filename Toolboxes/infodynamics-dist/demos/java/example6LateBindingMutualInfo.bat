@ECHO OFF

REM Make sure the latest example source file is compiled.

REM  Note: this does not nulify our claim that you can compile the code
REM  once and then only change the properties file - we only compile here
REM  every time this script is run so that we can capture any changes the user
REM  made to the demo source code. The demo as written should be compiled once
REM  then one can make dynamic changes to the props file and simply
REM  run the class file without recompiling.
javac -classpath "..\..\infodynamics.jar" "infodynamics\demos\Example6LateBindingMutualInfo.java"

REM Run the example, feeding in the properties file as the command line argument
java -classpath ".;..\..\infodynamics.jar" infodynamics.demos.Example6LateBindingMutualInfo example6LateBindingMutualInfo.props
