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

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ActiveInfoStorageCalculatorViaMutualInfo;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;

/**
 * An Active Information Storage (AIS) calculator (implementing {@link ActiveInfoStorageCalculator})
 * which is affected using a 
 * Kraskov-Stoegbauer-Grassberger (KSG) Mutual Information (MI) calculator
 * ({@link MutualInfoCalculatorMultiVariateKraskov}) to make the calculations.
 * 
 * <p>
 * That is, this class implements an AIS calculator using the KSG nearest-neighbour approach.
 * This is achieved by plugging in {@link MutualInfoCalculatorMultiVariateKraskov}
 * as the calculator into the parent class {@link ActiveInfoStorageCalculatorViaMutualInfo}.
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link ActiveInfoStorageCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step is either a simple call to {@link #ActiveInfoStorageCalculatorKraskov()},
 *      or else specifies which KSG algorithm to implement via {@link #ActiveInfoStorageCalculatorKraskov(int)}
 *      or {@link #ActiveInfoStorageCalculatorKraskov(String)};</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties for
 *      {@link MutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}
 *      (except {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF} as outlined
 *      in {@link ActiveInfoStorageCalculatorViaMutualInfo#setProperty(String, String)}).</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
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
 * @see MutualInfoCalculatorMultiVariateKraskov
 *
 */
public class ActiveInfoStorageCalculatorKraskov
	extends ActiveInfoStorageCalculatorViaMutualInfo {
	
	/**
	 * Class name for KSG MI estimator via KSG algorithm 1
	 */
	public static final String MI_CALCULATOR_KRASKOV1 = MutualInfoCalculatorMultiVariateKraskov1.class.getName();
	/**
	 * Class name for KSG MI estimator via KSG algorithm 2
	 */
	public static final String MI_CALCULATOR_KRASKOV2 = MutualInfoCalculatorMultiVariateKraskov2.class.getName();
		
	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger style AIS calculator.
	 * 
	 * Uses algorithm 2 by default.
	 *
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorKraskov() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(MI_CALCULATOR_KRASKOV2);
	}

	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger estimator for AIS,
	 *  with the supplied MI calculator name.
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be equal to {@link #MI_CALCULATOR_KRASKOV1} or {@link #MI_CALCULATOR_KRASKOV2}
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorKraskov(String calculatorName) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(calculatorName);
		// Now check that it was one of our Kraskov-Grassberger calculators:
		if (!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV1) &&
				!calculatorName.equalsIgnoreCase(MI_CALCULATOR_KRASKOV2)) {
			throw new ClassNotFoundException("Must be an underlying Kraskov-Grassberger calculator");
		}
	}

	/**
	 * Creates a new instance of the Kraskov-Stoegbauer-Grassberger estimator for AIS,
	 *  with the supplied Kraskov-Stoegbauer-Grassberger MI algorithm number
	 * 
	 * @param algorithm must be either 1 or 2 for the first or second KSG algorithm
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ActiveInfoStorageCalculatorKraskov(int algorithm) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(algorithm == 1 ? MI_CALCULATOR_KRASKOV1 : MI_CALCULATOR_KRASKOV2);
		if ((algorithm != 1) && (algorithm != 2)) {
			throw new ClassNotFoundException("Algorithm must be 1 or 2");
		}
	}
}
