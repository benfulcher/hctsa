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

package infodynamics.measures.discrete;

/**
 * Interface for calculators of information-theoretic measures
 * for single variables (e.g. entropy, active information storage).
 * The interface defines common operations such as
 * adding observations and calculating
 * local and average values, etc. 
 * 
 * <p>Usage is as per {@link InfoMeasureCalculatorDiscrete}, with
 * many methods for supplying observations and making 
 * calculations defined here.</p>
 * 
 * <p>It would ideally be an abstract class to be inherited from, but
 * it's more important for some of our calculators to have inheritance from
 * ContextOfPastCalculator, and since java doesn't allow multiple
 * inheritance, one of them has to miss out.
 * To get around this, we combine the two in
 * {@link SingleAgentMeasureDiscreteInContextOfPastCalculator}.
 * </p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface SingleAgentMeasureDiscrete {

	/**
 	 * Add observations in to our estimates of the pdfs.
	 *
	 * @param states series of samples
	 */
	public void addObservations(int states[]);
	
	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to the PDFs.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 */
	public void addObservations(int states[][]);

	/**
 	 * Add observations for a single variable of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to {@link #addObservations(int[][])}
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param col index of agent
	 */
	public void addObservations(int states[][], int col);

	/**
 	 * Add observations in to our estimates of the pdfs.
 	 * This call suitable only for homogeneous agents, as all
 	 *  agents will contribute to single pdfs.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 */
	public void addObservations(int states[][][]);

	/**
 	 * Add observations for a single agent of the multi-agent system
 	 *  to our estimates of the pdfs.
 	 * This call should be made as opposed to {@link #addObservations(int[][][])}
 	 *  for computing active info for heterogeneous agents.
	 *
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param index1 row index index the variable
	 * @param index2 column index of the variable
	 */
	public void addObservations(int states[][][], int index1, int index2);

	/**
	 * Compute the average value of the measure
	 * from the previously-supplied samples.
	 * 
	 * Must set average, min and max
	 * 
	 * @return the estimate of the measure
	 */
	public double computeAverageLocalOfObservations();

	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 *  
	 * Must set average, min and max
	 * 
	 * @param states time series of samples
	 * @return time-series of local values (indexed as per states)
	 */
	public double[] computeLocalFromPreviousObservations(int states[]);
	
	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method.
	 * This method to be used for homogeneous agents only,
	 *  since the local values will be computed for all variables.
	 *  
	 * Must set average, min and max
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @return 2D time-series of local values (indexed as per states)
	 */
	public double[][] computeLocalFromPreviousObservations(int states[][]);

	/**
	 * Computes local information theoretic measure for the given
	 *  variable in the 2D time-series
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method. 
	 * This method is suitable for heterogeneous agents, since
	 *  the specific variable is identified.
	 *  
	 * Must set average, min and max.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param col index of the given variable
	 * @return time-series of local values for the variable
	 */
	public double[] computeLocalFromPreviousObservations(int states[][], int col);

	/**
	 * Computes local information theoretic measure for the given
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method to be used for homogeneous agents only,
	 *  since the local values will be computed for all variables.
	 *  
	 * Must set average, min and max
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @return 3D time-series of local values (indexed as per states)
	 */
	public double[][][] computeLocalFromPreviousObservations(int states[][][]);

	/**
	 * Computes the local information theoretic measure for the given
	 *  variable in the 3D time-series
	 *  states, using pdfs built up from observations previously
	 *  sent in via the addObservations method 
	 * This method is suitable for heterogeneous agents, since
	 *  the specific variable is identified.
	 *  
	 * Must set average, min and max
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param index1 row index of the given variable
	 * @param index2 column index of the given variable
	 * @return time-series of local values for the variable
	 */
	public double[] computeLocalFromPreviousObservations(int states[][][], int index1, int index2);

	/**
	 * Standalone routine to 
	 * compute the local information-theoretic measure across a
	 *  time-series of states.
	 * Return a time-series array of local values.
	 * First history rows are zeros when the measure must build up
	 *  embedded history of the variable.
	 * 
	 * @param states time series of samples
	 * @return time-series of local values (indexed as per states)
	 */
	public double[] computeLocal(int states[]);
	
	/**
	 * Standalone routine to 
	 * compute the local information-theoretic measure across a 2D spatiotemporal
	 *  array of the states of homogeneous agents,
	 * Return a 2D spatiotemporal array of local values.
	 * First history rows are zeros when the measure must build up
	 *  embedded history of the variable.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @return 2D time-series of local values (indexed as per states)
	 */
	public double[][] computeLocal(int states[][]);

	/**
	 * Standalone routine to 
	 * compute the local information theoretic measure across a 3D spatiotemporal
	 *  array of the states of homogeneous agents
	 * Return a 3D spatiotemporal array of local values.
	 * First history rows are zeros when the measure must build up
	 *  embedded history of the variable.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @return 3D time-series of local values (indexed as per states)
	 */
	public double[][][] computeLocal(int states[][][]);

	/**
	 * Standalone routine to 
	 * compute the average information theoretic measure across a time-series
	 *  of states.
	 * Return the average.
	 * 
	 * @param states time series of samples
	 * @return average of the information-theoretic measure.
	 */
	public double computeAverageLocal(int states[]);
	
	/**
	 * Standalone routine to 
	 * compute the average information theoretic measure across a 2D spatiotemporal
	 *  array of the states of homogeneous agents.
	 * Return the average.
	 * This method to be called for homogeneous agents only,
	 * since all variables are used in the PDFs.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @return average of the information-theoretic measure.
	 */
	public double computeAverageLocal(int states[][]);

	/**
	 * Standalone routine to 
	 * compute the average information theoretic measure across a 3D spatiotemporal
	 *  array of the states of homogeneous agents.
	 * Return the average.
	 * This method to be called for homogeneous agents only,
	 * since all variables are used in the PDFs.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @return average of the information-theoretic measure.
	 */
	public double computeAverageLocal(int states[][][]);
	
	/**
	 * Standalone routine to 
	 * compute local information theoretic measure for one variable
	 *  in a 2D spatiotemporal
	 *  multivariate array.
	 * Return a time-series array of local values.
	 * First history rows are zeros when the measure must build up
	 *  embedded history of the variable.
	 * This method should be used for heterogeneous agents
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param col index of the given variable
	 * @return time-series of local values of the measure for
	 *  the given variable
	 */
	public double[] computeLocal(int states[][], int col);
	
	/**
	 * Standalone routine to 
	 * compute local information theoretic measure for one variable
	 *  in a 3D spatiotemporal
	 *  multivariate array.
	 * Return a time-series array of local values.
	 * First history rows are zeros when the measure must build up
	 *  embedded history of the variable.
	 * This method should be used for heterogeneous agents
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param index1 row index of the given variable
	 * @param index2 column index of the given variable
	 * @return time-series of local values of the measure for
	 *  the given variable
	 */
	public double[] computeLocal(int states[][][], int index1, int index2);
	
	/**
	 * Standalone routine to 
	 * compute the average information theoretic measure 
	 * for a single agent in a multivariate time series.
	 * Returns the average.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable number)
	 * @param col index of the given variable
	 * @return average of the measure for the given variable.
	 */
	public double computeAverageLocal(int states[][], int col);

	/**
	 * Standalone routine to 
	 * compute the average information theoretic measure 
	 * for a single agent in a multivariate time series.
	 * Returns the average.
	 * This method suitable for heterogeneous agents.
	 * 
	 * @param states multivariate time series
	 *  (1st index is time, 2nd index is variable row number,
	 *  3rd is variable column number)
	 * @param index1 row index of the given variable
	 * @param index2 column index of the given variable
	 * @return average of the measure for the given variable.
	 */
	public double computeAverageLocal(int states[][][], int index1, int index2);
}
