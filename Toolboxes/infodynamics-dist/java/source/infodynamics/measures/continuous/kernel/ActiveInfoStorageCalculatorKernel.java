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

package infodynamics.measures.continuous.kernel;

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;

/**
 * An Active Information Storage (AIS) calculator (implementing {@link ActiveInfoStorageCalculator})
 * which is affected using a 
 * box-kernel Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateKernel}) to make the calculations.
 * 
 * <p>
 * That is, this class implements an AIS calculator using box-kernel estimation.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKernel}
 * as the calculator into the parent class {@link ActiveInfoStorageCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #ActiveInfoStorageCalculatorKernel()};</li>
 *  <li>Additional initialisation options {@link #initialise(int, double)}
 *      and {@link #initialise(int, int, double)}; and</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateKernel#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} as outlined
 *      in {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)})</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>J.T. Lizier, M. Prokopenko and A.Y. Zomaya,
 * 		<a href="http://dx.doi.org/10.1016/j.ins.2012.04.016">
 * 		"Local measures of information storage in complex distributed computation"</a>,
 * 		Information Sciences, vol. 208, pp. 39-54, 2012.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ActiveInfoStorageCalculator
 * @see ActiveInfoStorageCalculatorViaMutualInfo
 * @see MutualInfoCalculatorMultiVariateKernel
 */
public class ActiveInfoStorageCalculatorKernel
	extends ActiveInfoStorageCalculatorViaMutualInfo {
	
	public static final String MI_CALCULATOR_KERNEL = MutualInfoCalculatorMultiVariateKernel.class.getName();
		
	/**
	 * Creates a new instance of the box-kernel estimator for AIS
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorKernel() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_KERNEL);
	}

	/**
	 * Initialises the calculator with the existing values for embedding length k,
	 * embedding delay tau and kernel width epsilon
	 * 
	 */
	public void initialise() throws Exception {
		initialise(k, tau, ((MutualInfoCalculatorMultiVariateKernel) miCalc).getKernelWidth());
	}

	/**
	 * Initialises the calculator using existing value for tau
	 * 
	 * @param k history embedding length
	 * @param epsilon kernel width for the box-kernel (same for all dimensions)
	 */
	public void initialise(int k, double epsilon) throws Exception {
		initialise(k, tau, epsilon);
	}

	/**
	 * Initialises the calculator with parameters as supplied here
	 * 
	 * @param k history embedding length
	 * @param tau embedding delay (see {@link ActiveInfoStorageCalculator#initialise(int, int)})
	 * @param epsilon kernel width for the box-kernel (same for all dimensions)
	 */
	public void initialise(int k, int tau, double epsilon) throws Exception {
		// Set the property before the calculator is initialised by the super class
		miCalc.setProperty(MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME, Double.toString(epsilon));
		super.initialise(k, tau);
	}
}
