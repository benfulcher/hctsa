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

package infodynamics.utils;

/**
 * 
 * Structure to hold a distribution of info-theoretic measurements,
 *  and a significance value for how an original measurement compared 
 *  with these, which is determined empirically by bootstrapping.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EmpiricalMeasurementDistribution extends MeasurementDistribution {

	/**
	 * Distribution of surrogate measurement values
	 */
	public double[] distribution;
	/**
	 * Whether the mean of the surrogate measurement distribution has
	 * been computed
	 */
	protected boolean computedMean = false;
	/**
	 * Computed mean of the surrogate measurement distribution
	 */
	protected double meanOfDist;
	/**
	 * Computed mean of the surrogate measurement distribution
	 */
	protected double stdOfDist;
	
	/**
	 * Construct an instance, ready to fill out the distribution
	 * 
	 * @param size number of surrogates used
	 */
	public EmpiricalMeasurementDistribution(int size) {
		super(); // Creating the super class with mean and pValue 0
		// These value will be filled out by the caller later.
		distribution = new double[size];
	}

	/**
	 * Construct an instance with the given distribution of 
	 * surrogates
	 * 
	 * @param distribution surrogate measurements
	 * @param actualValue actual observed value
	 */
	public EmpiricalMeasurementDistribution(double[] distribution, double actualValue) {
		super(actualValue, 0); // Using pValue = 0 temporarily ...
		this.distribution = distribution;
		int countWhereActualIsNotGreater = 0;
		for (int i = 0; i < distribution.length; i++) {
			if (distribution[i] >= actualValue) {
				countWhereActualIsNotGreater++;
			}
		}
		pValue = (double) countWhereActualIsNotGreater / (double) distribution.length;
	}
	
	// TODO Compute the significance under the assumption of a Gaussian distribution
	/*
	public double computeGaussianSignificance() {
		// Need to conpute the significance based on the assumption of 
		//  an underlying Gaussian distribution.
		// Use the t distribution for analysis, since we have a finite
		//  number of samples to comptue the mean and std from.
		return 0;
	}
	*/
	
	/**
	 * Assuming the distribution is Gaussian, return a t-score
	 * for our observed measurement
	 * 
	 * @return a t-score for our observed measurement
	 */
	public double getTSscore() {
		if (! computedMean) {
			meanOfDist = MatrixUtils.mean(distribution);
			stdOfDist = MatrixUtils.stdDev(distribution, meanOfDist);
			computedMean = true;
		}
		double t = (actualValue - meanOfDist) / stdOfDist;
		return t;
	}
	
	/**
	 * Return the mean of the distribution
	 * 
	 * @return the mean of the distribution
	 */
	public double getMeanOfDistribution() {
		if (! computedMean) {
			meanOfDist = MatrixUtils.mean(distribution);
			stdOfDist = MatrixUtils.stdDev(distribution, meanOfDist);
			computedMean = true;
		}
		return meanOfDist;
	}

	/**
	 * Return the standard deviation of the distribution
	 * 
	 * @return the standard deviation of the distribution
	 */
	public double getStdOfDistribution() {
		if (! computedMean) {
			meanOfDist = MatrixUtils.mean(distribution);
			stdOfDist = MatrixUtils.stdDev(distribution, meanOfDist);
			computedMean = true;
		}
		return stdOfDist;
	}
}
