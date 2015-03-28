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

import infodynamics.measures.continuous.PredictiveInfoCalculator;
import infodynamics.measures.continuous.PredictiveInfoCalculatorViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;

/**
 * A Predictive Information (PI) / Excess Entropy calculator
 * (implementing {@link PredictiveInfoCalculator})
 * which is affected using a 
 * box-kernel Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateKernel}) to make the calculations.
 * 
 * <p>
 * That is, this class implements a PI calculator using box-kernel estimation.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKernel}
 * as the calculator into the parent class {@link PredictiveInfoCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link PredictiveInfoCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #PredictiveInfoCalculatorKernel()};</li>
 *  <li>Additional initialisation options {@link #initialise(int, double)}
 *      and {@link #initialise(int, int, double)}; and</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateKernel#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} as outlined
 *      in {@link PredictiveInfoCalculatorViaMutualInfo#setProperty(String, String)})</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>Bialek, W., Nemenman, I., and Tishby, N.,
 *  <a href="http://dx.doi.org/10.1016/S0378-4371(01)00444-7">
 * 	"Complexity through nonextensivity"</a>,
 *  Physica A, 302, 89-99. (2001).</li>
 * 	<li>J. P. Crutchfield, D. P. Feldman,
 *  <a href="http://dx.doi.org/10.1063/1.1530990">
 * 	"Regularities Unseen, Randomness Observed: Levels of Entropy Convergence"</a>,
 *  Chaos, Vol. 13, No. 1. (2003), pp. 25-54.</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see PredictiveInfoCalculator
 * @see PredictiveInfoCalculatorViaMutualInfo
 * @see MutualInfoCalculatorMultiVariateKernel
 */
public class PredictiveInfoCalculatorKernel
	extends PredictiveInfoCalculatorViaMutualInfo {
	
	public static final String MI_CALCULATOR_KERNEL = MutualInfoCalculatorMultiVariateKernel.class.getName();
		
	/**
	 * Creates a new instance of the box-kernel estimator for PI
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public PredictiveInfoCalculatorKernel() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
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
	 * @param k embedding length
	 * @param epsilon kernel width for the box-kernel (same for all dimensions)
	 */
	public void initialise(int k, double epsilon) throws Exception {
		initialise(k, tau, epsilon);
	}

	/**
	 * Initialises the calculator with parameters as supplied here
	 * 
	 * @param k embedding length
	 * @param tau embedding delay (see {@link PredictiveInfoCalculator#initialise(int, int)})
	 * @param epsilon kernel width for the box-kernel (same for all dimensions)
	 */
	public void initialise(int k, int tau, double epsilon) throws Exception {
		// Set the property before the calculator is initialised by the super class
		miCalc.setProperty(MutualInfoCalculatorMultiVariateKernel.KERNEL_WIDTH_PROP_NAME, Double.toString(epsilon));
		super.initialise(k, tau);
	}
}
