Java Information Dynamics Toolkit (JIDT)
Copyright (C) 2012-2014 Joseph T. Lizier
Copyright (C) 2014-2015 Joseph T. Lizier and Ipek Özdemir

Version 1.2.1 (see release notes below)

JIDT provides a standalone, open source code Java implementation (usable in Matlab, Octave and Python) of information-theoretic measures of distributed computation in complex systems: i.e. information storage, transfer and modification.

This includes implementations for:
- both discrete and continuous-valued variables, principally for the measures transfer entropy, mutual information and active information storage;
- using various types of estimators (e.g. Kraskov-Stögbauer-Grassberger estimators, kernel estimation, linear-Gaussian).
    
=============
   License
=============

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=============
   Website
=============

Full information on the JIDT (usage, etc) is provided at the project page and wiki on google code:

http://code.google.com/p/information-dynamics-toolkit/

=============
Installation
=============

"Full" description of any required installation is at: http://code.google.com/p/information-dynamics-toolkit/wiki/Installation

However, if you are reading this file, you've downloaded a distribution and you're halfway there!

There are no dependencies to download; unless:
 a. You don't have java installed - download it from http://www.java.com/
 b. You wish to build the project using the build.xml script - this requires ant: http://ant.apache.org/
 c. You wish to run the JUnit test cases - this requires JUnit: http://www.junit.org/ - for how to run JUnit with our ant script see http://code.google.com/p/information-dynamics-toolkit/wiki/JUnitTestCases

Then just put the jar in a relevant location in your file structure.

That's it.

=============
Documentation
=============

A research paper describing the toolkit is included in the top level directory -- "InfoDynamicsToolkit.pdf".

A tutorial, providing background to the information-theoretic measures, various estimators, and then to the JIDT toolkit itself is included in the tutorial folder (see "JIDT-TutorialSlides.pdf" for the tutorial slides, and "README-TutorialAndExercise.pdf" for further description of the tutorial exercises).

Javadocs for the toolkit are included in the full distribution at javadocs.
They can also be generated using "ant javadocs" (useful if you are on an SVN view).
Further, they will soon be posted on the web.

The project wiki also contains further information on various aspects; see http://code.google.com/p/information-dynamics-toolkit/ to start.

Further documentation is provided by the Usage demo examples below.

You can also join our email discussion group jidt-discuss at http://groups.google.com/d/forum/jidt-discuss

=============
    Usage
=============

Several sets of demonstration code are distributed with the toolkit:

 a. demos/java -- basic examples on easily using the Java toolkit -- run these from the shell scripts in this directory -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/SimpleJavaExamples

 b. Several demo sets mirror the SimpleJavaExamples to demonstrate the use of the toolkit in non-Java environments: 
 
   i. demos/octave -- basic examples on easily using the Java toolkit from Octave or Matlab environments -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/OctaveMatlabExamples
 
   ii. demos/python -- basic examples on easily using the Java toolkit from Python -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/PythonExamples

   iii. demos/r -- basic examples on easily using the Java toolkit from R -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/R_Examples

   iv. demos/julia -- basic examples on easily using the Java toolkit from Julia -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/JuliaExamples

   v. demos/clojure -- basic examples on easily using the Java toolkit from Clojure -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/Clojure_Examples

 c. demos/octave/CellularAutomata -- using the Java toolkit to plot local information dynamics profiles in cellular automata; the toolkit is run under Octave or Matlab -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/CellularAutomataDemos
 
 d. demos/octave/SchreiberTransferEntropyExamples -- recreates the transfer entropy examples in Schreiber's original paper presenting this measure; shows the correct parameter settings to reproduce these results  -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/SchreiberTeDemos
 
 e. demos/octave/DetectingInteractionLags -- demonstration of using the transfer entropy with source-destination lags; the demo is run under Octave or Matlab -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/DetectingInteractionLags

 f. demos/java/InterregionalTransfer -- higher level example using collective transfer entropy to infer effective connections between "regions" of data -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/InterregionalTransfer

 g. demos/octave/NullDistributions --  investigating the correspondence between analytic and bootstrapped distributions for TE and MI under null hypotheses of no relationship; the demo is run under Octave or Matlab -- see description at https://code.google.com/p/information-dynamics-toolkit/wiki/NullDistributions

 h. java/unittests -- the JUnit test cases for the Java toolkit are included in the distribution -- these case also be browsed to see simple use cases for the various calculators in the toolkit -- see description at http://code.google.com/p/information-dynamics-toolkit/wiki/JUnitTestCases

=============
  Citation
=============

Please cite your use of this toolkit as:

Joseph T. Lizier, "JIDT: An information-theoretic toolkit for studying the dynamics of complex systems", Frontiers in Robotics and AI 1:11, 2014; doi:10.3389/frobt.2014.00011

A pre-print of this paper is distributed with this toolkit (InfoDynamicsToolkit.pdf) and is available at arXiv:1408.3270

=============
   Notices
=============

This project includes modified files from the Apache Commons Math library -- http://commons.apache.org/proper/commons-math/
This Apache 2 software is now included as a derivative work in this GPLv3 licensed JIDT project, as per: http://www.apache.org/licenses/GPL-compatibility.html
Notices and license for this software are found in the notices/commons-math directory.

The project includes adapted code from the JAMA project -- http://math.nist.gov/javanumerics/jama/
Notices and license for this software are found in the notices/JAMA directory.

The project includes adapted code from the octave-java package of the Octave-Forge project -- http://octave.sourceforge.net/java/
Notices for this software are found in the notices/JAMA directory.

===============
 Release notes
===============

v1.2.1 12/2/2015 at r621
------------------------
Added tutorial slides, description of exercises and sample exercise solutions;
Made jar target Java 1.6;
Added Schreiber TE heart-breath rate with KSG estimator demo code for Python.

v1.2 28/1/2015 at r601
-----------------------
Dynamic correlation exclusion, or Theiler window, added to all Kraskov estimators;
Added univariate MI calculation to simple demo 6;
Added Java code for Schreiber TE heart-breath rate with KSG estimator, ready for use as a template in Tutorial;
Patch for crashes in KSG conditional MI algorithm 2;

v1.1 14/11/2014 at r576
-----------------------
Implemented Fast Nearest Neighbour Search for Kraskov-Stögbauer-Grassberger (KSG) estimators for MI, conditional MI, TE, conditional TE, AIS, Predictive info, and multi-information. This includes a general (multivariate) k-d tree implementation;
Added multi-threading (using all available processors by default) for the KSG estimators -- code contributed by Ipek Özdemir;
Added Predictive information / Excess entropy implementations for KSG, kernel and Gaussian estimators;
Added R, Julia, and Clojure demos;
Added Windows batch files for the Simple Java Demos;
Added property for adding a small amount of noise to data in all KSG estimators;


v1.0 14/8/2014 at r434
----------------------
Added the draft of the paper on the toolkit to the release;
Javadocs made ready for release;
Switched source->destination arguments for discrete TE calculators to be with source first in line with continuous calculators;
Renamed all discrete calculators to have Discrete suffix -- TE and conditional TE calculators also renamed to remove "Apparent" prefix and change "Complete" to "Conditional";
Kraskov estimators now using 4 nearest neighbours by default;
Unit test for Gaussian TE against ChaLearn Granger causality measurement;
Added Schreiber TE demos; Interregional transfer demos; documentation for Interaction lag demos; added examples 7 and 8 to Simple Java demos;
Added property to add noise to data for Kraskov MI;
Added derivation of Apache Commons Math code for chi square distribution, and included relevant notices in our release;
Inserted translation class for arrays between Octave and Java;
Added analytic statistical significance calculation to Gaussian calculators and discrete TE;
Corrected Kraskov algorithm 2 for conditional MI to follow equation in Wibral et al. 2014.


v0.2.0 20/4/2014 at r284
------------------------
Rearchitected (most) Transfer Entropy and Multivariate TE calculators to use an underlying conditional mutual information calculator, and have arbitrary embedding delay, source-dest delay;
this includes moving Kraskov-Grassberger Transfer Entropy calculator to use a single conditional mutual information estimator instead of two mutual information estimators;
Rearchitected (most) Active Information Storage calculators to use an underlying mutual information calculator;
Added Conditional Transfer Entropy calculators using underlying conditional mutual information calculators;
Moved mixed discrete-continuous calculators to a new "mixed" package;
bug fixes. 

v0.1.4 11/9/2013 at r241
------------------------
added scripts to generate CA figures for 2013 book chapters;
added general Java demo code;
added Python demo code;
made Octave/Matlab demos and CA demos properly compatible for Matlab;
added extra Octave/Matlab general demos;
added more unit tests for MI and conditional MI calculators, including against results from Wibral's TRENTOOL;
bug fixes. 

v0.1.3 13/1/2013 at r151
------------------------
existing Octave/Matlab demo code made compatible with Matlab;
several bug fixes, including using max norm by default in Kraskov calculator (instead of requiring this to be set explicitly);
more unit tests (including against results from Kraskov's own MI implementation) 

v0.1.2 19/11/2012 at r116
-------------------------
Includes demo code for two newly submitted papers

v0.1.1 31/10/2012 at r104
------------------------
No notes

v0.1 24/10/2012 at r65?
------------------------
First distribution

=============

Joseph T. Lizier, 12/02/2015

