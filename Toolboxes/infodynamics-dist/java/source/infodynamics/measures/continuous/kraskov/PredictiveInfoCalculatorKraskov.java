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

package infodynamics.measures.continuous.kraskov;

import infodynamics.measures.continuous.PredictiveInfoCalculator;
import infodynamics.measures.continuous.PredictiveInfoCalculatorViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;

/**
 * A Predictive Information (PI) / Excess Entropy calculator
 * (implementing {@link PredictiveInfoCalculator})
 * which is affected using a 
 * Kraskov-Stoegbauer-Grassberger (KSG) Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateKraskov}) to make the calculations.
 * 
 * <p>
 * That is, this class implements a PI calculator using the KSG nearest-neighbour approach.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKraskov}
 * as the calculator into the parent class {@link PredictiveInfoCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link PredictiveInfoCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step is either a simple call to {@link #PredictiveInfoCalculatorKraskov()},
 *      or else specifies which KSG algorithm to implement via {@link #PredictiveInfoCalculatorKraskov(int)}
 *      or {@link #PredictiveInfoCalculatorKraskov(String)};</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} as outlined
 *      in {@link PredictiveInfoCalculatorViaMutualInfo#setProperty(String, String)}).</li>
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
 * @see MutualInfoCalculatorMultiVariateKraskov
 *
 */
public class PredictiveInfoCalculatorKraskov
	extends PredictiveInfoCalculatorViaMutualInfo {
	
	/**
	 * Class name for KSG MI estimator via KSG algorithm 1
	 */
	public static final String MI_CALCULATOR_KRASKOV1 = MutualInfoCalculatorMultiVariateKraskov1.class.getName();
	/**
	 * Class name for KSG MI estimator via KSG algorithm 2
	 */
	public static final String MI_CALCULATOR_KRASKOV2 = MutualInfoCalculatorMultiVariateKraskov2.class.getName();
		
	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger style PI calculator.
	 * 
	 * Uses algorithm 2 by default.
	 *
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public PredictiveInfoCalculatorKraskov() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_KRASKOV2);
	}

	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger estimator for PI,
	 *  with the supplied MI calculator name.
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be equal to {@link #MI_CALCULATOR_KRASKOV1} or {@link #MI_CALCULATOR_KRASKOV2}
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public PredictiveInfoCalculatorKraskov(String calculatorName) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(calculatorName);
		// Now check that it was one of our Kraskov-Grassberger calculators:
		if (!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV1) &&
				!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV2)) {
			throw new ClassNotFoundException("Must be an underlying Kraskov-Grassberger calculator");
		}
	}

	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger estimator for PI,
	 *  with the supplied Kraskov-Stoegbauer-Grassberger MI algorithm number
	 * 
	 * @param algorithm must be either 1 or 2 for the first or second KSG algorithm
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public PredictiveInfoCalculatorKraskov(int algorithm) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(algorithm == 1 ? MI_CALCULATOR_KRASKOV1 : MI_CALCULATOR_KRASKOV2);
		if ((algorithm != 1) && (algorithm != 2)) {
			throw new ClassNotFoundException("Algorithm must be 1 or 2");
		}
	}
}
