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

package infodynamics.measures.continuous.gaussian;

import infodynamics.measures.continuous.PredictiveInfoCalculator;
import infodynamics.measures.continuous.PredictiveInfoCalculatorViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;

/**
 * A Predictive Information (PI) / Excess Entropy calculator
 * (implementing {@link PredictiveInfoCalculator})
 * which is affected using a 
 * Gaussian Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateGaussian}) to make the calculations.
 * 
 * <p>
 * That is, this class implements a PI calculator using model of 
 * Gaussian variables with linear interactions.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateGaussian}
 * as the calculator into the parent class {@link PredictiveInfoCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link PredictiveInfoCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #PredictiveInfoCalculatorGaussian()}.</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} as outlined
 *      in {@link PredictiveInfoCalculatorViaMutualInfo#setProperty(String, String)})</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
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
 * @see MutualInfoCalculatorMultiVariateGaussian
 */
public class PredictiveInfoCalculatorGaussian
	extends PredictiveInfoCalculatorViaMutualInfo {
	
	public static final String MI_CALCULATOR_GAUSSIAN = MutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Creates a new instance of the Gaussian-estimate style predictive info calculator
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public PredictiveInfoCalculatorGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_GAUSSIAN);
	}
}
