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

import infodynamics.measures.continuous.EntropyCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential entropy of a given multivariate set of observations
 *  (implementing {@link EntropyCalculatorMultiVariate}, assuming that
 *  the probability distribution function for these observations is Gaussian.</p>
 *
 * <p>Usage is as per the paradigm outlined for {@link EntropyCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #EntropyCalculatorMultiVariateGaussian()}.</li>
 * 	<li>The user can call {@link #setCovariance(double[][])} or
 *     {@link #setCovarianceAndMeans(double[][], double[])}
 *     instead of supplying observations via {@link #setObservations(double[])}.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
    <li>Differential entropy for Gaussian random variables defined at 
 *      <a href="http://mathworld.wolfram.com/DifferentialEntropy.html">MathWorld</a></li>
 *  <li>Multivariate normal distribution on <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Wikipedia</a></li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorMultiVariateGaussian
	implements EntropyCalculatorMultiVariate, Cloneable {

	/**
	 * Cached Cholesky decomposition of the most recently supplied covariance matrix
	 */
	protected double[][] L;

	/**
	 * Means of the most recently supplied observations (source variables
	 *  listed first, destination variables second).
	 */
	protected double[] means;
	
	/**
	 * The set of observations, retained in case the user wants to retrieve the local
	 *  entropy values of these
	 */
	protected double[][] observations;
	
	/**
	 * Number of dimensions for our multivariate data
	 */
	protected int dimensions;
	
	/**
	 * Determinant of the covariance matrix; stored to save computation time
	 */
	protected double detCovariance;

	/**
	 * Last average entropy we computed
	 */
	protected double lastAverage;
	
	/**
	 * Whether we are in debug mode
	 */
	protected boolean debug;
	
	/**
	 * Construct an instance
	 */
	public EntropyCalculatorMultiVariateGaussian() {
		// Nothing to do
	}
	
	public void initialise(int dimensions) {
		means = null;
		L = null;
		observations = null;
		this.dimensions = dimensions;
		detCovariance = 0;
	}

	/**
	 * @throws Exception where the observations do not match the expected number of 
	 *  dimensions, or covariance matrix is not positive definite (reflecting 
	 *  redundant variables in the observations)
	 */
	public void setObservations(double[][] observations) throws Exception {
		// Check that the observations was of the correct number of dimensions:
		if (observations[0].length != dimensions) {
			means = null;
			L = null;
			throw new Exception("Supplied observations does not match initialised number of dimensions");
		}
		means = MatrixUtils.means(observations);
		setCovariance(MatrixUtils.covarianceMatrix(observations, means));
		// And keep a reference to the observations used here (must set this
		//  *after* setCovariance, since setCovariance sets the observations to null
		this.observations = observations;
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  entropy.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}.</p>
	 *  
	 * @param covariance covariance matrix between the variables
	 * @throws Exception if the covariance matrix does not match the dimensions supplied
	 *  in {@link #initialise(int)}, is non-square, is asymmetric or is non-positive
	 *  definite (i.e. there are redundant terms).
	 */
	public void setCovariance(double[][] covariance) throws Exception {
		detCovariance = 0;
		observations = null;
		// Check the dimensions of the covariance matrix:
		int rows = covariance.length;
		if (rows != dimensions) {
			throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
		}
		// Make sure the matrix is symmetric and positive definite, by taking the 
		//  Cholesky decomposition (which we need for the determinant later anyway):
		//  (this will check and throw Exceptions for non-square,
		//   asymmetric, non-positive definite A)
		L = MatrixUtils.CholeskyDecomposition(covariance);
	}
	
	/**
	 * <p>Set the covariance and mean of the distribution for which we will compute the
	 *  entropy.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}.</p>
	 * 
	 * @param covariance covariance matrix of the variables
	 * @param means mean of the variables
	 * @throws Exception where the dimensions of the covariance or means are not correct,
	 *  or the covariance matrix is non-square, is asymmetric or is non-positive
	 *  definite (i.e. there are redundant terms).
	 */
	public void setCovarianceAndMeans(double[][] covariance, double[] means) throws Exception {
		setCovariance(covariance);
		// Only set means after setCovariance has returned ok
		if (means.length != dimensions) {
			throw new Exception("Supplied mean matrix does not match initialised number of dimensions");
		}
		this.means = means;
	}

	/**
 	 * Compute the entropy from the previously supplied observations, or
	 * based on the supplied variance.
	 * 
	 * <p>The joint entropy for a multivariate Gaussian-distribution of dimension n
	 *  with covariance matrix C is 0.5*\log_e{(2*pi*e)^n*|det(C)|},
	 *  where det() is the matrix determinant of C.</p>
	 * 
	 * <p>Here we compute the joint entropy assuming that the recorded estimation of the
	 *  covariance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the joint entropy of the previously provided observations or from the
	 *  supplied covariance matrix. Returned in nats (NOT bits).
	 */
	public double computeAverageLocalOfObservations() {
		// Simple way:
		// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
		// Using cached Cholesky decomposition:
		detCovariance = MatrixUtils.determinantViaCholeskyResult(L);
		lastAverage =  0.5 * (dimensions* (1 + Math.log(2.0*Math.PI)) +
			Math.log(detCovariance));
		return lastAverage;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	/**
	 * No properties are defined here, so this method will have no effect.
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		// No properties to set here
	}

	public double getLastAverage() {
		return lastAverage;
	}

	/**
	 * @return an array of local values in nats (NOT bits).
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states)
			throws Exception {
		if (means == null) {
			throw new Exception("Cannot compute local values without having means either supplied or computed via setObservations()");
		}
		
		// Check that the covariance matrix was positive definite:
		// (this was done earlier in computing the Cholesky decomposition,
		//  we may still need to compute the determinant)
		if (detCovariance == 0) {
			// The determinant has not been computed yet
			// Simple way:
			// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
			// Using cached Cholesky decomposition:
			detCovariance = MatrixUtils.determinantViaCholeskyResult(L);
			if (detCovariance == 0) {
				throw new Exception("Covariance matrix is not positive definite");
			}
		}
		// Now we are clear to take the matrix inverse (via Cholesky decomposition,
		//  since we have a symmetric positive definite matrix):
		double[][] invCovariance = MatrixUtils.solveViaCholeskyResult(L,
				MatrixUtils.identityMatrix(L.length));
		
		// If we have a time delay, slide the local values
		double[] localValues = new double[states.length];
		for (int t = 0; t < states.length; t++) {
			double[] deviationsFromMean =
					MatrixUtils.subtract(states[t], means);
			// Computing PDF
			// (see the PDF defined at the wikipedia page referenced in the method header)
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double pJoint = Math.pow(2.0 * Math.PI, -(double) dimensions / 2.0) *
					Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);
			localValues[t] = - Math.log(pJoint);
		}
		
		// Don't set average if this was the previously supplied observations,
		//  since it won't be the same as what would have been computed 
		//  analytically.
		
		return localValues;
	}

	/**
	 * @throws Exception if {@link #setCovariance(double[][])} or
	 * {@link #setCovarianceAndMeans(double[][], double[])} were used previously instead
	 * of {@link #setObservations(double[][])}
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		if (observations == null) {
			throw new Exception("Cannot compute local values since no observations were supplied");
		}
		return computeLocalUsingPreviousObservations(observations);
	}

	public int getNumObservations() throws Exception {
		if (observations == null) {
			throw new Exception("Cannot return number of observations because either " +
					"this calculator has not had observations supplied or " +
					"the user supplied the covariance matrix instead of observations");
		}
		return observations.length;
	}

	/**
	 * Provide an implementation of the clone() method.
	 * This does not deeply copy all of the underlying data, just providing
	 *  a copy of the references to it all.
	 * This is enough to protect the integrity of the calculator
	 *  however if the clone is supplied different data (though the 
	 *  clone should not alter the data).
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
}
