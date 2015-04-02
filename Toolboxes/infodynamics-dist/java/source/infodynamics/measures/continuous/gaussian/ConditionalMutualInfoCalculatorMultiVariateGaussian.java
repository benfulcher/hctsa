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

import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.NonPositiveDefiniteMatrixException;

/**
 * <p>Computes the differential conditional mutual information of two multivariate
 *  <code>double[][]</code> sets of observations, conditioned on another
 *  (implementing {@link ConditionalMutualInfoCalculatorMultiVariate}),
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.
 *  This is achieved by extending the common code base in {@link ConditionalMutualInfoMultiVariateCommon}.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link ConditionalMutualInfoCalculatorMultiVariate},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #ConditionalMutualInfoCalculatorMultiVariateGaussian()}.</li>
 *  <li>The property {@link ConditionalMutualInfoMultiVariateCommon#PROP_NORMALISE}
 *     is set to false by default here (since this makes more sense for
 *     linear-Gaussian analysis), which is different to the parent class.</li>
 * 	<li>The user can call {@link #setCovariance(double[][], boolean)} or
 *     {@link #setCovariance(double[][], int)} or {@link #setCovarianceAndMeans(double[][], double[], int)}
 *     instead of supplying observations via {@link #setObservations(double[][], double[][], double[][])} or
 *     {@link #addObservations(double[][], double[][], double[][])} etc.</li>
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
public class ConditionalMutualInfoCalculatorMultiVariateGaussian 
		extends ConditionalMutualInfoMultiVariateCommon
		implements ConditionalMutualInfoCalculatorMultiVariate,
			AnalyticNullDistributionComputer, Cloneable {

	/**
	 * Cached Cholesky decomposition of the covariance matrix
	 * of the most recently supplied observations.
	 * Is a matrix [C_11, C_12, C_1c; C_21, C_22, C_2c; C_c1, C_c2, C_cc],
	 * where C_xy represents the covariance matrix of variable x to variable y
	 * where x,y are either variable 1, 2 or the conditional.
	 * The covariance matrix is symmetric, and should be positive definite
	 * (otherwise we have linealy dependent variables).
	 */
	protected double[][] L;
	/**
	 * Cached Cholesky decomposition of the (var1, conditional) covariance matrix
	 */
	protected double[][] L_1c;
	/**
	 * Cached Cholesky decomposition of the (var2, conditional) covariance matrix
	 */
	protected double[][] L_2c;
	/**
	 * Cached Cholesky decomposition of the conditional covariance matrix
	 */
	protected double[][] L_cc;
	
	/**
	 * Means of the most recently supplied observations (source variables
	 *  listed first, destination variables second).
	 */
	protected double[] means;

	/**
	 * Cached determinants of the joint covariance matrix
	 */
	protected double detCovariance;
	/**
	 * Cached determinants of the covariance matrix for
	 *  variable 1 and conditional
	 */
	protected double det1cCovariance;
	/**
	 * Cached determinants of the covariance matrix for
	 *  variable 2 and the conditional
	 */
	protected double det2cCovariance;
	/**
	 * Cached determinants of the covariance matrix
	 *  for the conditional
	 */
	protected double detccCovariance;
	
	/**
	 * Cache the sub-variables which are a linearly-independent set
	 *  (and so are used in the covariances) for variable 1
	 */
	protected int[] var1IndicesInCovariance;
	/**
	 * Cache the sub-variables which are a linearly-independent set
	 *  (and so are used in the covariances) for variable 2
	 */
	protected int[] var2IndicesInCovariance;
	/**
	 * Cache the sub-variables which are a linearly-independent set
	 *  (and so are used in the covariances) for the conditional
	 */
	protected int[] condIndicesInCovariance;

	/**
	 * Construct an instance
	 */
	public ConditionalMutualInfoCalculatorMultiVariateGaussian() {
		// Normalising data makes less sense for linear-Gaussian estimation,
		//  so we turn this off by default.
		normalise = false;
	}
	
	@Override
	public void initialise(int var1Dimensions, int var2Dimensions, int condDimensions) {
		super.initialise(var1Dimensions, var2Dimensions, condDimensions);
		L = null;
		L_1c = null;
		L_2c = null;
		L_cc = null;
		means = null;
		detCovariance = 0;
		det1cCovariance = 0;
		det2cCovariance = 0;
		detccCovariance = 0;
		condIndicesInCovariance = null;
		var1IndicesInCovariance = null;
		var2IndicesInCovariance = null;
	}

	/**
	 * @throws Exception if the observation variables are not linearly independent
	 *  (leading to a non-positive definite covariance matrix).
	 */
	@Override
	public void finaliseAddObservations() throws Exception {

		// Get the observations properly stored in the sourceObservations[][] and
		//  destObservations[][] arrays.
		super.finaliseAddObservations();

		// Store the means of each variable (useful for local values later)
		means = new double[dimensionsVar1 + dimensionsVar2 + dimensionsCond];
		double[] var1Means = MatrixUtils.means(var1Observations);
		double[] var2Means = MatrixUtils.means(var2Observations);
		double[] condMeans = MatrixUtils.means(condObservations);
		System.arraycopy(var1Means, 0, means, 0, dimensionsVar1);
		System.arraycopy(var2Means, 0, means, dimensionsVar1, dimensionsVar2);
		System.arraycopy(condMeans, 0, means, dimensionsVar1 + dimensionsVar2,
				dimensionsCond);
		
		// Store the covariances of the variables
		// Generally, this should not throw an exception, since we checked
		//  the observations had the correct number of variables
		//  on receiving them, and in constructing the covariance matrix
		//   ourselves we know it should be symmetric.
		// It could occur however if the covariance matrix was not
		//  positive definite, which would occur if one variable
		//  is linearly redundant.
		setCovariance(
				MatrixUtils.covarianceMatrix(var1Observations, var2Observations, condObservations),
				true);
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  conditional mutual information.</p>
	 *  
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of var1, var2, conditional
	 *  variables, considered jointly together.
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
	 *  conditional mutual information.</p>
	 *  
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}.</p>
	 * 
	 * @param covariance covariance matrix of var1, var2 and the conditional
	 *  variables, considered jointly together.
	 * @param determinedFromObservations whether the covariance matrix
	 *  was determined internally from observations or not
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
	protected void setCovariance(double[][] covariance, boolean determinedFromObservations)
			throws Exception {
		if (!determinedFromObservations) {
			// Make sure we're not keeping any observations
			var1Observations = null;
			var2Observations = null;
			condObservations = null;
		}
		// Make sure the supplied covariance matrix matches the required dimenions:
		int rows = covariance.length;
		if (rows != dimensionsVar1 + dimensionsVar2 + dimensionsCond) {
			throw new Exception("Supplied covariance matrix does not match initialised number of dimensions");
		}
		
		// Now store the Cholesky decompositions for computing the cond MI later.
		// Start with the conditional variable:
		//  In case there are linear redundancies amongst the conditional variable,
		//  remove some of its components until the linear redundancy is gone.
		//  This allows a conditional MI computation to still take place.
		boolean redundanciesRemoved = false;
		condIndicesInCovariance = null;
		for (int v = dimensionsCond; v >= 1; v--) {
			// Select the first v variables in the conditional joint variable:
			condIndicesInCovariance = MatrixUtils.range(dimensionsVar1 + dimensionsVar2,
					dimensionsVar1 + dimensionsVar2 + v - 1);
			double[][] condCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						condIndicesInCovariance, condIndicesInCovariance);
			try {
				L_cc = MatrixUtils.CholeskyDecomposition(condCovariance);
			} catch (NonPositiveDefiniteMatrixException e) {
				// There is a linear redundancy between the variables - so allow
				//  one to be removed in the next loop iteration
				continue;
			}
			// Allow exceptions indicating asymmetric and non-square to be propagated
			// Otherwise, we've found a linearly independent subset of the conditioned variable
			redundanciesRemoved = true;
			break;
		}
		if (!redundanciesRemoved) {
			// This will only happen if the last remaining variable had covariance of 0
			L_cc = null;
			// So we won't condition on anything here:
			condIndicesInCovariance = new int[] {};
		}
		
		// And next store the Cholesky decompositions for var1 with
		//  the conditional variable:
		//  In case there are linear redundancies amongst the conditional variable
		//  with variable 1, remove some of the components of variable 1
		//  until the linear redundancy is gone.
		//  This allows a conditional MI computation to still take place.
		var1IndicesInCovariance = null;
		redundanciesRemoved = false;
		for (int v = dimensionsVar1; v >= 1; v--) {
			var1IndicesInCovariance = MatrixUtils.range(0, v - 1);
			int[] var1AndCondIndicesInCovariance = MatrixUtils.append(var1IndicesInCovariance, condIndicesInCovariance);
			double[][] var1AndCondCovariance =
					MatrixUtils.selectRowsAndColumns(covariance, 
							var1AndCondIndicesInCovariance, var1AndCondIndicesInCovariance);
			try {
				L_1c = MatrixUtils.CholeskyDecomposition(var1AndCondCovariance);
			} catch (NonPositiveDefiniteMatrixException e) {
				// There is a linear redundancy between the variables - so allow
				//  one to be removed in the next loop iteration
				continue;
			}
			// Allow exceptions indicating asymmetric and non-square to be propagated
			// Otherwise, we've found a linearly independent subset of the conditioned variable
			redundanciesRemoved = true;
			break;
		}
		if (!redundanciesRemoved) {
			// This means that variable 1 is *fully* linearly dependent on the conditioned
			//  variable (i.e. not even cutting it down to just one variable helped).
			// Flag this by setting:
			L_1c = null;
		}
		
		// And next store the Cholesky decompositions for var2 with
		//  the conditional variable:
		//  In case there are linear redundancies amongst the conditional variable
		//  with variable 2, remove some of the components of variable 2
		//  until the linear redundancy is gone.
		//  This allows a conditional MI computation to still take place.
		var2IndicesInCovariance = null;
		int[] var2AndCondIndicesInCovariance = null;
		redundanciesRemoved = false;
		for (int v = dimensionsVar2; v >= 1; v--) {
			var2IndicesInCovariance = MatrixUtils.range(dimensionsVar1, dimensionsVar1 + v - 1);
			var2AndCondIndicesInCovariance = MatrixUtils.append(var2IndicesInCovariance, condIndicesInCovariance);
			double[][] var2AndCondCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						var2AndCondIndicesInCovariance, var2AndCondIndicesInCovariance);
			try {
				L_2c = MatrixUtils.CholeskyDecomposition(var2AndCondCovariance);
			} catch (NonPositiveDefiniteMatrixException e) {
				// There is a linear redundancy between the variables - so allow
				//  one to be removed in the next loop iteration
				continue;
			}
			// Allow exceptions indicating asymmetric and non-square to be propagated
			// Otherwise, we've found a linearly independent subset of the conditioned variable
			redundanciesRemoved = true;
			break;
		}
		if (!redundanciesRemoved) {
			// This means that variable 2 is *fully* linearly dependent on the conditioned
			//  variable (i.e. not even cutting it down to just one variable helped).
			// Flag this by setting:
			L_2c = null;
		}

		// Finally, store the Cholesky decomposition for the whole covariance matrix:
		//  first prune the covariance in line with the variable removed
		//  from variable 1, 2 and the conditional:
		int[] prunedIndicesInCovariance = MatrixUtils.append(var1IndicesInCovariance, var2AndCondIndicesInCovariance);
		double[][] prunedCovariance =
				MatrixUtils.selectRowsAndColumns(covariance, 
						prunedIndicesInCovariance, prunedIndicesInCovariance);
		try {
			L = MatrixUtils.CholeskyDecomposition(prunedCovariance);
		} catch (NonPositiveDefiniteMatrixException e) {
			// There is a linear redundancy between the variables - 
			// Flag this by setting:
			L = null;
		}
		// Allow exceptions indicating asymmetric and non-square to be propagated

		// Postcondition: L's contain Cholesky decompositions of covariance
		//  matrices with linearly dependent variables removed (except for 
		//  the whole covariance matrix), using null to flag where this was not possible
	}

	/**
	 * <p>Set the covariance of the distribution for which we will compute the
	 *  mutual information, as well as the means for each variable.</p>
	 * 
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}
	 *  (but you can call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}).</p>
	 * 
	 * @param covariance covariance matrix of var1, var2 and conditional
	 *  variables, considered jointly together.
	 * @param means mean of var1, var2 and conditional variables (as per
	 *  <code>covariance</code>)
	 * @param numObservations the number of observations that the mean and covariance
	 *  were determined from. This is used for later significance calculations
	 */
	public void setCovarianceAndMeans(double[][] covariance, double[] means,
			int numObservations) throws Exception {
		
		this.means = means;
		setCovariance(covariance, numObservations);
	}

	/**
	 * <p>Computes the local values of the conditional mutual information,
	 *  for each valid observation in the previously supplied observations
	 *  (with PDFs computed using all of the previously supplied observation sets).</p>
	 *  
	 * <p>The joint differential entropy for a multivariate Gaussian-distribution of dimension n
	 *  with covariance matrix C is -0.5*\log_e{(2*pi*e)^n*|det(C)|},
	 *  where det() is the matrix determinant of C.</p>
	 * 
	 * <p>Here we compute the conditional mutual information from the joint entropies
	 *  of all variables (H_12c), variable 1 and conditional (H_1c),
	 *  variable 2 and conditional (H_2c) and conditional (H_c),
	 *  giving MI = H_1c + H_2c - H_c - H_12c.
	 *  We assume that the recorded estimation of the
	 *  covariance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the conditional mutual information of the previously provided observations or from the
	 *  supplied covariance matrix, in <b>nats</b> (not bits!).
	 *  Returns NaN if any of the determinants are zero
	 *  (because this will make the denominator of the log zero)
	 */
	public double computeAverageLocalOfObservations() throws Exception {
		// Simple way:
		// detCovariance = MatrixUtils.determinantSymmPosDefMatrix(covariance);
		// Using cached Cholesky decomposition:
		// And also with extended checks for linear redundancies:
		
		// Should always have a valid L_cc (since we can reduce it down to one variable:
		
		if (L_1c == null) {
			// Variable 1 is fully linearly redundant with conditional, so
			//  we will have zero conditional MI:
			lastAverage = 0;
		} else {
			det1cCovariance = MatrixUtils.determinantViaCholeskyResult(L_1c);
			if (L_2c == null) {
				// Variable 2 is fully linearly redundant with conditional, so
				//  we will have zero conditional MI:
				lastAverage = 0;
			} else {
				det2cCovariance = MatrixUtils.determinantViaCholeskyResult(L_2c);
				if (L == null) {
					// There is a linear dependence amongst variables 1 and 2 given the
					//  conditional which did not exist for either with the conditional alone,
					//  so conditional MI diverges:
					lastAverage = Double.POSITIVE_INFINITY;
				} else {
					detCovariance = MatrixUtils.determinantViaCholeskyResult(L);

					// Else all the covariance matrices were ok, except perhaps the
					//  conditional covariance
					if (L_cc == null) {
						// The conditional variables had no covariance, so
						//  just return an MI:
						lastAverage = 0.5 * Math.log(Math.abs(
								det1cCovariance * det2cCovariance /
										detCovariance));
					} else {
						detccCovariance = MatrixUtils.determinantViaCholeskyResult(L_cc);
						// So compute as normal:
						lastAverage = 0.5 * Math.log(Math.abs(
									det1cCovariance * det2cCovariance /
											(detCovariance * detccCovariance)));
					}
				}
			}
		}
		
		condMiComputed = true;
		return lastAverage;
	}

	/**
	 * @return array of the local values in nats (not bits!)
	 * @throws Exception if the user had previously supplied covariances
	 *  directly (ie had not supplied observations).
	 */
	public double[] computeLocalOfPreviousObservations() throws Exception {
		// Cannot do if destObservations haven't been set
		if (var2Observations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"if they have not been set!");
		}
		
		return computeLocalUsingPreviousObservations(var1Observations,
				var2Observations, condObservations, true);
	}

	/**
	 * Generate an <b>analytic</b> distribution of what the conditional MI would look like,
	 * under a null hypothesis that our variables had no relation
	 * (in the context of the conditional value).
	 * This is performed without bootstrapping (which is done in
	 * {@link #computeSignificance(int, int)} and {@link #computeSignificance(int, int[][])}).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below, and the other papers referenced in
	 * {@link AnalyticNullDistributionComputer#computeSignificance()}
	 * (in particular Geweke),
	 * for a description of how this is done for conditional MI.
	 * Basically, the null distribution is a chi-square distribution 
	 * with degrees of freedom equal to the product of the number of variables
	 * in each joint variable 1 and 2.
	 * </p>
	 * 
	 * @return ChiSquareMeasurementDistribution object which describes
	 * the proportion of conditional MI scores from the null distribution
	 *  which have higher or equal conditional MIs to our actual value.
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	@Override
	public ChiSquareMeasurementDistribution computeSignificance() throws Exception {
		if (!condMiComputed) {
			computeAverageLocalOfObservations();
		}
		// Number of extra parameters in the model incorporating the
		//  extra variable is independent of the number of variables
		//  in the conditional:
		// (Assuming that all variables went into the calculation:)
		// return new ChiSquareMeasurementDistribution(2.0*((double)totalObservations)*lastAverage,
		//		dimensionsVar1 * dimensionsVar2);
		// Taking the subsets into account:
		return new ChiSquareMeasurementDistribution(2.0*((double)totalObservations)*lastAverage,
				var1IndicesInCovariance.length * var2IndicesInCovariance.length);
	}
	
	/**
	 * @throws Exception if user passed in covariance matrix rather than observations
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int numPermutationsToCheck) throws Exception {
		if (var2Observations == null) {
			throw new Exception("Cannot compute empirical statistical significance " +
					"if user passed in covariance matrix rather than observations.");
		}

		return super.computeSignificance(variableToReorder, numPermutationsToCheck);
	}

	/**
	 * @throws Exception if user passed in covariance matrix rather than observations
	 */
	@Override
	public EmpiricalMeasurementDistribution computeSignificance(
			int variableToReorder, int[][] newOrderings) throws Exception {
		if (var2Observations == null) {
			throw new Exception("Cannot compute empirical statistical significance " +
					"if user passed in covariance matrix rather than observations.");
		}

		return super.computeSignificance(variableToReorder, newOrderings);
	}

	public int getNumObservations() throws Exception {
		if (var2Observations == null) {
			throw new Exception("Cannot return number of observations because either " +
					"this calculator has not had observations supplied or " +
					"the user supplied the covariance matrix instead of observations");
		}
		return super.getNumObservations();
	}

	/**
	 * @throws Exception if the user previously supplied covariance directly rather
	 *  than by setting observations (this means we have no observations
	 *  to reorder).
	 */
	public double computeAverageLocalOfObservations(int variableToReorder, 
			int[] newOrdering) throws Exception {
		// Cannot do if observations haven't been set (i.e. the variances
		//  were directly supplied)
		if (var1Observations == null) {
			throw new Exception("Cannot compute local values of previous observations " +
					"without supplying observations");
		}
		return super.computeAverageLocalOfObservations(variableToReorder, newOrdering);
	}

	/**
	 * @return the series of local conditional MI values in <b>nats</b> (not bits).
	 * @throws Exception
	 */
	public double[] computeLocalUsingPreviousObservations(double[][] states1,
			double[][] states2, double[][] condStates) throws Exception {
		return computeLocalUsingPreviousObservations(
				states1, states2, condStates, false);
	}

	/**
	 * Utility function to implement either {@link #computeLocalOfPreviousObservations()}
	 * or {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])}
	 * depending on value of <code>isPreviousObservations</code>
	 * 
	 * @param newVar1Obs provided variable 1 observations
	 * @param newVar2Obs provided variable 2 observations
	 * @param newCondObs provided conditional observations
	 * @param isPreviousObservations whether these are our previous
	 *  observations - this determines whether we are implementing
	 *  {@link #computeLocalOfPreviousObservations()} if true 
	 *  or {@link #computeLocalUsingPreviousObservations(double[][], double[][], double[][])} 
	 *  if false. Also indicates whether to 
	 *  set the internal lastAverage field,
	 *  which is returned by later calls to {@link #getLastAverage()} 
	 *  (in case of true).
	 * @return the local conditional MI values in nats (not bits).
	 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
	 * @see <a href="http://en.wikipedia.org/wiki/Positive-definite_matrix>"Positive definite matrix in Wikipedia"</a>
	 * @throws Exception if means were not defined by {@link #setObservations(double[][], double[][])} etc
	 *  or {@link #setCovarianceAndMeans(double[][], double[])}
	 */
	protected double[] computeLocalUsingPreviousObservations(double[][] newVar1Obs,
			double[][] newVar2Obs, double[][] newCondObs, boolean isPreviousObservations) throws Exception {
		
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

			if (L_cc == null) {
				// We will compute local MIs
				detccCovariance = 0;
			} else {
				detccCovariance = MatrixUtils.determinantViaCholeskyResult(L_cc);
			}
			if (L_1c == null) {
				// Variable 1 is fully linearly redundant with conditional, so
				//  we will have zero conditional MI:
				return MatrixUtils.constantArray(newVar2Obs.length, 0);
			} else {
				det1cCovariance = MatrixUtils.determinantViaCholeskyResult(L_1c);
				if (L_2c == null) {
					// Variable 2 is fully linearly redundant with conditional, so
					//  we will have zero conditional MI:
					return MatrixUtils.constantArray(newVar2Obs.length, 0);
				} else {
					det2cCovariance = MatrixUtils.determinantViaCholeskyResult(L_2c);
					if (L == null) {
						// There is a linear dependence amongst variables 1 and 2 given the
						//  conditional which did not exist for either with the conditional alone,
						//  so conditional MI diverges:
						return MatrixUtils.constantArray(newVar2Obs.length, Double.POSITIVE_INFINITY);
					} else {
						detCovariance = MatrixUtils.determinantViaCholeskyResult(L);
					}
				}
			}
		}

		// Now we are clear to take the matrix inverse (via Cholesky decomposition,
		//  since we have a symmetric positive definite matrix):
		double[][] invCovariance = MatrixUtils.solveViaCholeskyResult(L,
				MatrixUtils.identityMatrix(L.length));
		double[][] invVar1CondCovariance = MatrixUtils.solveViaCholeskyResult(L_1c,
				MatrixUtils.identityMatrix(L_1c.length));
		double[][] invVar2CondCovariance = MatrixUtils.solveViaCholeskyResult(L_2c,
				MatrixUtils.identityMatrix(L_2c.length));
		double[][] invCondCovariance = null;
		if (L_cc != null) {
			invCondCovariance = MatrixUtils.solveViaCholeskyResult(L_cc,
					MatrixUtils.identityMatrix(L_cc.length));
		}
			
		// Now, only use the means from the subsets of linearly independent variables:
		// double[] var1Means = MatrixUtils.select(means, 0, dimensionsVar1);
		double[] var1Means = MatrixUtils.select(means, var1IndicesInCovariance);
		// double[] var2Means = MatrixUtils.select(means, dimensionsVar1, dimensionsVar2);
		double[] var2Means = MatrixUtils.select(means, var2IndicesInCovariance);
		// double[] condMeans = MatrixUtils.select(means, dimensionsVar1 + dimensionsVar2, dimensionsCond);
		double[] condMeans = MatrixUtils.select(means, condIndicesInCovariance);
		
		int lengthOfReturnArray;
		lengthOfReturnArray = newVar2Obs.length;

		double[] localValues = new double[lengthOfReturnArray];
		int[] var2IndicesSelected = MatrixUtils.subtract(var2IndicesInCovariance, dimensionsVar1);
		int[] condIndicesSelected = MatrixUtils.subtract(condIndicesInCovariance, dimensionsVar1 + dimensionsVar2);
		for (int t = 0; t < newVar2Obs.length; t++) {
			
			double[] var1DeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newVar1Obs[t], var1IndicesInCovariance),
							var1Means);
			double[] var2DeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newVar2Obs[t], var2IndicesSelected),
							var2Means);
			double[] condDeviationsFromMean =
					MatrixUtils.subtract(
							MatrixUtils.select(newCondObs[t], condIndicesSelected),
							condMeans);
			double[] var1CondDeviationsFromMean =
					MatrixUtils.append(var1DeviationsFromMean,
							condDeviationsFromMean);
			double[] var2CondDeviationsFromMean =
					MatrixUtils.append(var2DeviationsFromMean,
							condDeviationsFromMean);
			double[] tempDeviationsFromMean =
					MatrixUtils.append(var1DeviationsFromMean,
							var2DeviationsFromMean);
			double[] deviationsFromMean =
					MatrixUtils.append(tempDeviationsFromMean,
							condDeviationsFromMean);
			
			// Computing PDFs WITHOUT (2*pi)^dim factor, since these will cancel:
			// (see the PDFs defined at the wikipedia page referenced in the method header)
			double var1CondExpArg = MatrixUtils.dotProduct(
						MatrixUtils.matrixProduct(var1CondDeviationsFromMean,
								invVar1CondCovariance),
						var1CondDeviationsFromMean);
			double adjustedPVar1Cond = Math.exp(-0.5 * var1CondExpArg) /
						Math.sqrt(det1cCovariance);
			double var2CondExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(var2CondDeviationsFromMean,
							invVar2CondCovariance),
					var2CondDeviationsFromMean);
			double adjustedPVar2Cond = Math.exp(-0.5 * var2CondExpArg) /
					Math.sqrt(det2cCovariance);
			double condExpArg = 0;
			double adjustedPCond = 0;
			if (L_cc != null) {
				condExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(condDeviationsFromMean,
							invCondCovariance),
					condDeviationsFromMean);
				adjustedPCond = Math.exp(-0.5 * condExpArg) /
					Math.sqrt(detccCovariance);
			}
			double jointExpArg = MatrixUtils.dotProduct(
					MatrixUtils.matrixProduct(deviationsFromMean,
							invCovariance),
					deviationsFromMean);
			double adjustedPJoint = Math.exp(-0.5 * jointExpArg) /
					Math.sqrt(detCovariance);
			
			if (L_cc != null) {
				// Returning results in nats:
				localValues[t] = Math.log(adjustedPJoint * adjustedPCond /
						(adjustedPVar1Cond * adjustedPVar2Cond));
			} else {
				// Return an MI (no linearly independent, non-zero conditional vars):
				localValues[t] = Math.log(adjustedPJoint /
						(adjustedPVar1Cond * adjustedPVar2Cond));
			}
				
		}
		
		// if (isPreviousObservations) {
		// Don't store the average value here, since it won't be exactly 
		//  the same as what would have been computed under the analytic expression
		// }
		
		return localValues;
	}
}
