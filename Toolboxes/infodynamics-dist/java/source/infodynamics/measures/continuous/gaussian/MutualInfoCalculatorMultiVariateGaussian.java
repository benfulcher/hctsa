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

import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoMultiVariateCommon;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential mutual information of two given multivariate
 *  <code>double[][]</code> sets of
 *  observations (implementing {@link MutualInfoCalculatorMultiVariate}),
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link MutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #MutualInfoCalculatorMultiVariateGaussian()}.</li>
 * 	<li>The user can call {@link #setCovariance(double[][], boolean)} or
 *     {@link #setCovariance(double[][], int)} or {@link #setCovarianceAndMeans(double[][], double[], int)}
 *     instead of supplying observations via {@link #setObservations(double[][], double[][])} or
 *     {@link #addObservations(double[][], double[][])} etc.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  <li>Additional method {@link #computeSignificance()} to compute null distribution analytically.</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. M. Cover and J. A. Thomas, 'Elements of Information
Theory' (John Wiley & Sons, New York, 1991).</li>
 * </ul>
 * 
 * @see <a href="http://mathworld.wolfram.com/DifferentialEntropy.html">Differential entropy for Gaussian random variables at Mathworld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Differential_entropy">Differential entropy for Gaussian random variables at Wikipedia</a>
 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MutualInfoCalculatorMultiVariateGaussian 
		extends MutualInfoMultiVariateCommon
		implements MutualInfoCalculatorMultiVariate,
			AnalyticNullDistributionComputer, Cloneable {

	/**
	 * Cached Cholesky decomposition of the covariance matrix
	 * of the most recently supplied observations.
	 * Is a matrix [C_ss, C_sd; C_ds, C_dd], where C_ss is the covariance
	 * matrix of the source observations, C_dd is the covariance matrix
	 * of the destination observations, and C_sd and C_ds are the covariances
	 * of source to destination and destination to source observations.
	 * The covariance matrix is symmetric, and should be positive definite
	 * (otherwise we have linealy dependent variables).
	 */
	protected double[][] L;
	/**
	 * Cached Cholesky decomposition of the source covariance matrix
	 */
	protected double[][] Lsource;
	/**
	 * Cached Cholesky decomposition of the destination covariance matrix
	 */
	protected double[][] Ldest;
	
	/**
	 * Means of the most recently supplied observations (source variables
	 *  listed first, destination variables second).
	 */
	protected double[] means;

	/**
	 * Cached determinant of the joint covariance matrix
	 */
	protected double detCovariance;
	/**
	 * Cached determinant of the source covariance matrix
	 */
	protected double detSourceCovariance;
	/**
	 * Cached determinant of the destination covariance matrix
	 */
	protected double detDestCovariance;
	
	/**
	 * Construct an instance of the Gaussian MI calculator
	 */
	public MutualInfoCalculatorMultiVariateGaussian() {
		// Nothing to do
	}
	
	public void initialise(int sourceDimensions, int destDimensions) {
		super.initialise(sourceDimensions, destDimensions);
		L = null;
		Lsource = null;
		Ldest = null;
		means = null;
		detCovariance = 0;
		detSourceCovariance = 0;
		detDestCovariance = 0;
	}

	/**
	 * @throws Exception if the observation variables are not linearly independent
	 *  (leading to a non-positive definite covariance matrix).
	 */
	public void finaliseAddObservations() throws Exception {

		// Get the observations properly stored in the sourceObservations[][] and
		//  destObservations[][] arrays.
		super.finaliseAddObservations();

		// Store the means of each variable (useful for local values later)
		means = new double[dimensionsSource + dimensionsDest];
		double[] sourceMeans = MatrixUtils.means(sourceObservations);
		double[] destMeans = MatrixUtils.means(destObservations);
		System.arraycopy(sourceMeans, 0, means, 0, dimensionsSource);
		System.arraycopy(destMeans, 0, means, dimensionsSource, dimensionsDest);
		
		// Store the covariances of the variables
		// Generally, this should not throw an exception, since we checked
		//  the observations had the correct number of variables
		//  on receiving them, and in constructing the covariance matrix
		//   ourselves we know it should be symmetric.
		// It could occur however if the covariance matrix was not
		//  positive definite, which would occur if one variable
		//  is linearly redundant.
		setCovariance(MatrixUtils.covarianceMatrix(sourceObservations,
			destObservations), true);
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information.</p>
	 *  
	 * <p>This is an alternative to sequences of calls to {@link #setObservations(double[][], double[][])} or
	 * {@link #addObservations(double[][], double[][])} etc.
	 * Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of the source and destination
	 *  variables, considered together (variable indices start with the source
	 *  and continue into the destination).
	 *  I.e. it is a matrix [C_ss, C_sd; C_ds, C_dd], where C_ss is the covariance
	 * matrix of the source observations, C_dd is the covariance matrix
	 * of the destination observations, and C_sd and C_ds are the covariances
	 * of source to destination and destination to source observations.
	 * @param numObservations the number of observations that the covariance
	 *  was determined from. This is used for later significance calculations
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
	public void setCovariance(double[][] covariance, int numObservations) throws Exception {
		setCovariance(covariance, false);
		totalObservations = numObservations;
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information.</p>
	 *  
	 * <p>This is an alternative to sequences of calls to {@link #setObservations(double[][], double[][])} or
	 * {@link #addObservations(double[][], double[][])} etc.
	 * Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of the source and destination
	 *  variables, considered jointly together (variable indices start with the source
	 *  and continue into the destination).
	 * @param determinedFromObservations whether the covariance matrix
	 *  was determined internally from observations or not
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
	protected void setCovariance(double[][] covariance, boolean determinedFromObservations)
			throws Exception {
		if (!determinedFromObservations) {
			// Make sure we're not keeping any observations
			sourceObservations = null;
			destObservations = null;
		}
		// Make sure the supplied covariance matrix matches the required dimenions:
		int rows = covariance.length;
		if (rows != dimensionsSource + dimensionsDest) {
			throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
		}

		// Make sure the matrix is symmetric and positive definite, by taking the 
		//  Cholesky decomposition (which we need for the determinant later anyway):
		//  (this will check and throw Exceptions for non-square,
		//   asymmetric, non-positive definite A)
		L = MatrixUtils.CholeskyDecomposition(covariance);
		
		// And store the Cholesky decompositions for the source covariance
		//  and dest covariance as well:
		int[] sourceIndicesInCovariance = MatrixUtils.range(0, dimensionsSource - 1);
		double[][] sourceCovariance =
			MatrixUtils.selectRowsAndColumns(covariance, 
					sourceIndicesInCovariance, sourceIndicesInCovariance);
		Lsource = MatrixUtils.CholeskyDecomposition(sourceCovariance);

		int[] destIndicesInCovariance = MatrixUtils.range(dimensionsSource,
				dimensionsSource + dimensionsDest - 1);
		double[][] destCovariance =
			MatrixUtils.selectRowsAndColumns(covariance, 
				destIndicesInCovariance, destIndicesInCovariance);
		Ldest = MatrixUtils.CholeskyDecomposition(destCovariance);
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information.</p>
	 * 
	 * <p>This is an alternative to sequences of calls to {@link #setObservations(double[][], double[][])} or
	 * {@link #addObservations(double[][], double[][])} etc.
	 * Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}.</p>
	 * 
	 * @param covariance covariance matrix of the source and destination
	 *  variables, considered together (variable indices start with the source
	 *  and continue into the destination).
	 * @param means mean of the source and destination variables (as per
	 *  covariance)
	 * @param numObservations the number of observations that the mean and covariance
	 *  were determined from. This is used for later significance calculations
	 */
	public void setCovarianceAndMeans(double[][] covariance, double[] means,
			int numObservations) throws Exception {
		this.means = means;
		setCovariance(covariance, numObservations);
	}

	/**
	 * Compute the MI from the supplied observations or covariances.
	 * 
	 * <p>The joint entropy for a multivariate Gaussian-distribution of dimension n
	 *  with covariance matrix C is -0.5*\log_e{(2*pi*e)^n*|det(C)|},
	 *  where det() is the matrix determinant of C.</p>
	 * 
	 * <p>Here we compute the mutual information from the joint entropies
	 *  of the source variables (H_s), destination variables (H_d), and all variables
	 *  taken together (H_sd), giving MI = H_s + H_d - H_sd.
	 *  We assume that the recorded estimation of the
	 *  covariance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the MI of the previously provided observations or from the
	 *  supplied covariance matrix, in nats (not bits!).
	 *  Returns NaN if any of the determinants are zero
	 *  (because this will make the denominator of the log 0).
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		// Simple way:
		// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
		// Using cached Cholesky decomposition:
		detCovariance = MatrixUtils.determinantViaCholeskyResult(L);
		detSourceCovariance = MatrixUtils.determinantViaCholeskyResult(Lsource);
		detDestCovariance = MatrixUtils.determinantViaCholeskyResult(Ldest);

		lastAverage = 0.5 * Math.log(Math.abs(
				detSourceCovariance * detDestCovariance /
						detCovariance));
		miComputed = true;
		return lastAverage;
	}

	/**
	 * <p>Computes the local values of the MI,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>If the samples were supplied via a single call such as
	 * {@link #setObservations(double[])},
	 * then the return value is a single time-series of local
	 * channel measure values corresponding to these samples.</p>
	 * 
	 * <p>Otherwise where disjoint time-series observations were supplied using several 
	 *  calls such as {@link addObservations(double[])}
	 *  then the local values for each disjoint observation set will be appended here
	 *  to create a single "time-series" return array.</p>
	 * 
	 * <p>If the user supplied covariance matrices rather than observations
	 * then this method cannot be called.</p>
	 * 
	 * @return the "time-series" of local MIs in nats (not bits!)
	 * @throws Exception
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		// Cannot do if destObservations haven't been set
		if (destObservations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"if they have not been set!");
		}
		
		return computeLocalUsingPreviousObservations(sourceObservations,
				destObservations, true);
	}

	/**
	 * Generate an <b>analytic</b> distribution of what the MI would look like,
	 * under a null hypothesis that our variables had no relation.
	 * This is performed without bootstrapping (which is done in
	 * {@link #computeSignificance(int)} and {@link #computeSignificance(int[][])}).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below, and the other papers referenced in
	 * {@link AnalyticNullDistributionComputer#computeSignificance()}
	 * (in particular Brillinger and Geweke),
	 * for a description of how this is done for MI.
	 * Basically, the null distribution is a chi-square distribution 
	 * with degrees of freedom equal to the product of the number of variables
	 * in each joint variable 1 and 2.
	 * </p>
	 * 
	 * @return ChiSquareMeasurementDistribution object which describes
	 * the proportion of MI scores from the null distribution
	 *  which have higher or equal MIs to our actual value.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	public ChiSquareMeasurementDistribution computeSignificance() {
		return new ChiSquareMeasurementDistribution(2*totalObservations*lastAverage,
				dimensionsSource * dimensionsDest);
	}
	
	/**
	 * @throws Exception where the user did not set observations 
	 * but set covariance matrices instead.
	 */
	public int getNumObservations() throws Exception {
		if (destObservations == null) {
			throw new Exception("Cannot return number of observations because either " +
					"this calculator has not had observations supplied or " +
					"the user supplied the covariance matrix instead of observations");
		}
		return super.getNumObservations();
	}

	/**
	 * @return the MI under the new ordering, in nats (not bits!).
	 *  Returns NaN if any of the determinants are zero
	 *  (because this will make the denominator of the log 0).
	 * @throws Exception if the user previously supplied covariance directly rather
	 *  than by setting observations (this means we have no observations
	 *  to reorder).
	 */
	public double computeAverageLocalOfObservations(int[] newOrdering)
			throws Exception {
		// Cannot do if observations haven't been set (i.e. the variances
		//  were directly supplied)
		if (destObservations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"without supplying observations");
		}
		return super.computeAverageLocalOfObservations(newOrdering);
	}

	/**
	 * @return the local values in nats (not bits).
	 *  If the {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 *  property was set to say k, then the local values align with the
	 *  destination value (i.e. after the given delay k). As such, the
	 *  first k values of the array will be zeros.  
	 * @throws Exception if means were not defined by supplying observations
	 *  (eg via {@link #setObservations(double[][], double[][])} etc) 
	 *  or calling {@link #setCovarianceAndMeans(double[][], double[])}
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] newSourceObs,
			double[][] newDestObs) throws Exception {
		return computeLocalUsingPreviousObservations(newSourceObs, newDestObs, false);
	}

	/**
	 * Protected utility function to compute the local MI values for each of the
	 * supplied samples in <code>newSourceObs</code> and <code>newDestObs</code>.
	 * 
	 * <p>PDFs are computed using all of the previously supplied
	 * observations. <code>isPreviousObservations</code> indicates whether
	 * those in <code>states1</code> and <code>states2</code>
	 * were some of the previously supplied samples.</p>
	 * 
	 * @param newSourceObs provided source observations
	 * @param newDestObs provided destination observations
	 * @param isPreviousObservations whether these are our previous
	 *  observations - this determines whether to add zeros for the first
	 *  timeDiff local values, and also 
	 *  whether to set the internal lastAverage field,
	 *  which is returned by later calls to {@link #getLastAverage()}
	 * @return the local values in nats (not bits).
	 *  If the {@link MutualInfoCalculatorMultiVariate#PROP_TIME_DIFF}
	 *  property was set to say k, then the local values align with the
	 *  destination value (i.e. after the given delay k). As such, the
	 *  first k values of the array will be zeros.
	 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
	 * @see <a href="http://en.wikipedia.org/wiki/Positive-definite_matrix>"Positive definite matrix in Wikipedia"</a>
	 * @throws Exception if means were not defined by supplying observations
	 *  (eg via {@link #setObservations(double[][], double[][])} etc) 
	 *  or calling {@link #setCovarianceAndMeans(double[][], double[])}
	 */
	protected double[] computeLocalUsingPreviousObservations(double[][] newSourceObs,
			double[][] newDestObs, boolean isPreviousObservations) throws Exception {
		
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
			detSourceCovariance = MatrixUtils.determinantViaCholeskyResult(Lsource);
			detDestCovariance = MatrixUtils.determinantViaCholeskyResult(Ldest);
		}

		// Now we are clear to take the matrix inverse (via Cholesky decomposition,
		//  since we have a symmetric positive definite matrix):
		double[][] invCovariance = MatrixUtils.solveViaCholeskyResult(L,
				MatrixUtils.identityMatrix(L.length));
		double[][] invSourceCovariance = MatrixUtils.solveViaCholeskyResult(Lsource,
				MatrixUtils.identityMatrix(Lsource.length));
		double[][] invDestCovariance = MatrixUtils.solveViaCholeskyResult(Ldest,
				MatrixUtils.identityMatrix(Ldest.length));
			
		double[] sourceMeans = MatrixUtils.select(means, 0, dimensionsSource);
		double[] destMeans = MatrixUtils.select(means, dimensionsSource, dimensionsDest);
		
		int lengthOfReturnArray, offset;
		if (isPreviousObservations && addedMoreThanOneObservationSet) {
			// We're returning the local values for a set of disjoint
			//  observations. So we don't add timeDiff zeros to the start,
			//  and note that the required timeDiff is already 
			//  built into the supplied observations.
			lengthOfReturnArray = newDestObs.length;
			offset = 0;
		} else {
			lengthOfReturnArray = newDestObs.length + timeDiff;
			offset = timeDiff;
		}
		// If we have a time delay, slide the local values
		double[] localValues = new double[lengthOfReturnArray];
		for (int t = offset; t < newDestObs.length; t++) {
			// Computing local values for:
			//  a. sourceObservations[t - offset]
			//  b. destObservations[t]
			
			double[] sourceDeviationsFromMean =
					MatrixUtils.subtract(newSourceObs[t - offset],
							sourceMeans);
			double[] destDeviationsFromMean =
					MatrixUtils.subtract(newDestObs[t], destMeans);
			double[] deviationsFromMean =
					MatrixUtils.append(sourceDeviationsFromMean,
							destDeviationsFromMean);
			
			// Computing PDFs WITHOUT (2*pi)^dim factor, since these will cancel:
			// (see the PDFs defined at the wikipedia page referenced in the method header)
			double sourceExpArg = MatrixUtils.dotProduct(
						MatrixUtils.matrixProduct(sourceDeviationsFromMean,
								invSourceCovariance),
						sourceDeviationsFromMean);
			double adjustedPSource = Math.exp(-0.5 * sourceExpArg) /
						Math.sqrt(detSourceCovariance);
			double destExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(destDeviationsFromMean,
							invDestCovariance),
					destDeviationsFromMean);
			double adjustedPDest = Math.exp(-0.5 * destExpArg) /
					Math.sqrt(detDestCovariance);
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double adjustedPJoint = Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);
			
			// Returning results in nats:
			double localValue = Math.log(adjustedPJoint /
					(adjustedPSource * adjustedPDest));
			localValues[t] = localValue;			
		}
		
		// if (isPreviousObservations) {
		// Don't store the average value here, since it won't be exactly 
		//  the same as what would have been computed under the analytic expression
		// }
		
		return localValues;
	}

}
