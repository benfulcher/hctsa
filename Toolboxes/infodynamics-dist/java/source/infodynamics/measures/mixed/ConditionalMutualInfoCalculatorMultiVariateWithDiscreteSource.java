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

package infodynamics.measures.mixed;

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.utils.EmpiricalMeasurementDistribution;

/**
 * A conditional mutual information calculator between a joint set of continuous variables, 
 *  and a discrete variable, conditioned on another continuous variable vector. 
 * 
 * <p>These calculators are <b>EXPERIMENTAL</b> -- not properly tested,
 * and not well documented. The intended calling pattern is similar to
 * {@link ConditionalMutualInfoCalculatorMultiVariate}
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface ConditionalMutualInfoCalculatorMultiVariateWithDiscreteSource {

	/**
	 * 
	 * 
	 * @param dimensions the number of joint continuous variables
	 * @param base the base of the discrete variable
	 * @param dimensionsCond the number of joint continuous variables
	 *     to condition on
	 * @throws Exception
	 */
	public void initialise(int dimensions, int base, int dimensionsCond) throws Exception;
	
	public void setProperty(String propertyName, String propertyValue);
	
	public void startAddObservations();

	public void finaliseAddObservations() throws Exception;

	public void addObservations(double[][] continuousObservations,
			int[] discreteObservations, double[][] conditionedObservations) throws Exception;

	public void setObservations(double[][] continuousObservations,
			int[] discreteObservations, double[][] conditionedObservations) throws Exception;

	public double computeAverageLocalOfObservations() throws Exception;

	public double[] computeLocalOfPreviousObservations() throws Exception;

	public double[] computeLocalUsingPreviousObservations(double[][] contStates,
			int[] discreteStates, double[][] conditionedStates) throws Exception;

	/**
	 * 
	 * @param reorderDiscreteVariable boolean for whether to reorder the discrete variable (true)
	 *    or the continuous variable (false)
	 * @param numPermutationsToCheck
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(boolean reorderDiscreteVariable, int numPermutationsToCheck) throws Exception;

	/**
	 * 
	 * @param reorderDiscreteVariable boolean for whether to reorder the discrete variable (true)
	 *    or the continuous variable (false)
	 * @param newOrderings
	 * @return
	 * @throws Exception
	 */
	public EmpiricalMeasurementDistribution computeSignificance(boolean reorderDiscreteVariable, int[][] newOrderings) throws Exception;

	public void setDebug(boolean debug);

	public double getLastAverage();
	
	public int getNumObservations();

}
