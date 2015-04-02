/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.demos;

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.ParsedProperties;

/**
 * = Example 6 - Late binding Mutual info calculator =
 *
 * This class is used to demonstrate the manner in which a user
 *  can code to the interfaces defined in infodynamics.measures.continuous,
 *  and dynamically alter the instantiated class at runtime.
 * We demonstrate this using a multivariate mutual information calculation.
 * 
 * This example also demonstrates how to read simple files of arrays of data
 *  with the toolkit, as well as how to dynamically load properties from a 
 *  java properties file. 
 * 
 * @author Joseph Lizier
 *
 */
public class Example6LateBindingMutualInfo {

	/**
	 * @param args One command line argument taken, specifying location of 
	 *  the properties file. This should be example6LateBindingMutualInfo.props
	 *  in the demos/java directory.
	 */
	public static void main(String[] args) throws Exception {
		
		// 0. Preliminaries (reading in the dynamic properties and the data):
		//  a. Read in the properties file defined as the first
		//     command line argument:
		ParsedProperties props = new ParsedProperties(args[0]);
		//  b. Read in the data file, whose filename is defined in the
		//     property "datafile" in our properties file:
		ArrayFileReader afr = new ArrayFileReader(props.getStringProperty("datafile"));
		double[][] data = afr.getDouble2DMatrix();
		//  c. Pull out the columns from the data set which 
		//     correspond to the univariate and joint variables we will work with:
		//     First the univariate series to compute standard MI between: 
		int univariateSeries1Column = props.getIntProperty("univariateSeries1Column");
		int univariateSeries2Column = props.getIntProperty("univariateSeries2Column");
		double[] univariateSeries1 = MatrixUtils.selectColumn(data, univariateSeries1Column);
		double[] univariateSeries2 = MatrixUtils.selectColumn(data, univariateSeries2Column);
		//     Next the multivariate series to compute joint or multivariate MI between: 
		int[] jointVariable1Columns = props.getIntArrayProperty("jointVariable1Columns");
		int[] jointVariable2Columns = props.getIntArrayProperty("jointVariable2Columns");
		double[][] jointVariable1 = MatrixUtils.selectColumns(data, jointVariable1Columns);
		double[][] jointVariable2 = MatrixUtils.selectColumns(data, jointVariable2Columns);
		
		// 1. Create a reference for our calculator as
		//  an object implementing the interface type:
		MutualInfoCalculatorMultiVariate miCalc;
		
		// 2. Define the name of the class to be instantiated here:
		String implementingClass = props.getStringProperty("implementingClass");
		
		// 3. Dynamically instantiate an object of the given class:
		//  Part 1: Class.forName(implementingClass) grabs a reference to
		//   the class named by implementingClass.
		//  Part 2: .newInstance() creates an object instance of that class.
		//  Part 3: (MutualInfoCalculatorMultiVariate) casts the return
		//   object into an instance of our generic interface type.
		miCalc = (MutualInfoCalculatorMultiVariate)
				Class.forName(implementingClass).newInstance();
		
		// 4. Start using our MI calculator, paying attention to only
		//  call common methods defined in the interface type, not methods
		//  only defined in a given implementation class.
		// a. Initialise the calculator for a univariate calculation:
		miCalc.initialise(1, 1);
		// b. Supply the observations to compute the PDFs from:
		miCalc.setObservations(univariateSeries1, univariateSeries2);
		// c. Make the MI calculation:
		double miUnivariateValue = miCalc.computeAverageLocalOfObservations();

		// 5. Continue onto a multivariate calculation, still only
		//  calling common methods defined in the interface type.
		// a. Initialise the calculator for a multivariate calculation
		//  to use the required number of dimensions for each variable:
		miCalc.initialise(jointVariable1Columns.length, jointVariable2Columns.length);
		// b. Supply the observations to compute the PDFs from:
		miCalc.setObservations(jointVariable1, jointVariable2);
		// c. Make the MI calculation:
		double miJointValue = miCalc.computeAverageLocalOfObservations();
		
		System.out.printf("MI calculator %s computed the univariate MI(%d;%d) as %.5f " +
				" and joint MI as %.5f\n",
				implementingClass, univariateSeries1Column, univariateSeries2Column,
				miUnivariateValue, miJointValue);
	}

}
