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

import infodynamics.measures.continuous.ConditionalTransferEntropyCalculator;
import infodynamics.measures.continuous.ConditionalTransferEntropyCalculatorViaCondMutualInfo;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ChiSquareMeasurementDistribution;

/**
 * <p>Computes the differential conditional transfer entropy (TE) between two univariate
 *  <code>double[]</code> time-series of observations
 *  (implementing {@link ConditionalTransferEntropyCalculator}),
 *  assuming that the probability distribution function for these observations is
 *  a multivariate Gaussian distribution.
 *  TE was defined by Schreiber, conditional TE by Lizier et al.,
 *  and Kaiser and Schreiber showed how to compute
 *  TE via the Gaussian assumption.
 *  This estimator is realised here by plugging in
 *  {@link ConditionalMutualInfoCalculatorMultiVariateGaussian}
 *  as the calculator into the parent class
 *  {@link ConditionalTransferEntropyCalculatorViaCondMutualInfo}.</p>
 *  
 * <p>That is, this class implements a conditional TE calculator using model of 
 * Gaussian variables with linear interactions, making it equivalent
 * (up to a multiplicative constant) to the (partial) Granger causality (see Barnett et al below).
 * </p> 
 * 
 * <p>Usage is as per the paradigm outlined for {@link ConditionalTransferEntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to
 * 		{@link #ConditionalTransferEntropyCalculatorGaussian()}.</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties defined for both
 * 		{@link ConditionalTransferEntropyCalculator#setProperty(String, String)} and
 *      {@link ConditionalMutualInfoCalculatorMultiVariateGaussian#setProperty(String, String)}
 *      as outlined
 *      in {@link ConditionalTransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}).</li>
 * 	<li>The user can call {@link #setCovariance(double[][], int)}
 *     instead of supplying observations via {@link #setObservations(double[], double[], double[][])} or
 *     {@link #addObservations(double[], double[], double[][])} etc.</li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  <li>Additional method {@link #computeSignificance()} to compute null distribution analytically.</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>L. Barnett, A. B. Barrett, A. K. Seth, <a href="http://dx.doi.org/10.1103/physrevlett.103.238701">
 *  "Granger Causality and Transfer Entropy Are Equivalent for Gaussian Variables"</a>,
 *  Physical Review Letters 103 (23) 238701, 2009;</li>
 *  <li>A. Kaiser, T. Schreiber, <a href="http://dx.doi.org/10.1016/s0167-2789(02)00432-3">
 *  "Information transfer in continuous processes"</a>,
 *  Physica D, Vol. 166, No. 1-2., pp. 43-62 (2002).</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href=http://dx.doi.org/10.1063/1.3486801">
 *  "Information modification and particle collisions in distributed computation"</a>
 *  Chaos 20, 3, 037109 (2010).</li>
 * </ul>
 * 
 * @see <a href="http://mathworld.wolfram.com/DifferentialEntropy.html">Differential entropy for Gaussian random variables at Mathworld</a>
 * @see <a href="http://en.wikipedia.org/wiki/Differential_entropy">Differential entropy for Gaussian random variables at Wikipedia</a>
 * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution">Multivariate normal distribution on Wikipedia</a>
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see ConditionalTransferEntropyCalculator
 * @see ConditionalMutualInfoCalculatorMultiVariateGaussian
 */
public class ConditionalTransferEntropyCalculatorGaussian
	extends ConditionalTransferEntropyCalculatorViaCondMutualInfo
	implements AnalyticNullDistributionComputer {
	
	/**
	 * Name of the Gaussian conditional MI calculator we will use 
	 */
	public static final String COND_MI_CALCULATOR_GAUSSIAN = ConditionalMutualInfoCalculatorMultiVariateGaussian.class.getName();
	
	/**
	 * Creates a new instance of the Gaussian-estimate style conditional transfer entropy calculator
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public ConditionalTransferEntropyCalculatorGaussian() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(COND_MI_CALCULATOR_GAUSSIAN);
	}

	/**
	 * <p>Set the joint covariance of the distribution for which we will compute the
	 *  transfer entropy.</p>
	 *  
	 * <p>Note that without setting any observations, you cannot later
	 *  call {@link #computeLocalOfPreviousObservations()}, and without
	 *  providing the means of the variables, you cannot later call
	 *  {@link #computeLocalUsingPreviousObservations(double[][], double[][])}.</p>
	 * 
	 * @param covariance joint covariance matrix of source, dest, dest history
	 *  and conditional variables, considered together.
	 * @param numObservations the number of observations that the covariance
	 *  was determined from. This is used for later significance calculations
	 * @throws Exception for covariance matrix not matching the expected dimensions,
	 *  being non-square, asymmetric or non-positive definite
	 */
	public void setCovariance(double[][] covariance, int numObservations) throws Exception {
		((ConditionalMutualInfoCalculatorMultiVariateGaussian) condMiCalc).
				setCovariance(covariance, numObservations);
	}

	/**
	 * Generate an <b>analytic</b> distribution of what the 
	 * conditional TE would look like,
	 * under a null hypothesis that our source and destination
	 * variables had no relation
	 * (in the context of the conditional value).
	 * This is performed without bootstrapping (which is done in
	 * {@link #computeSignificance(int, int)} and {@link #computeSignificance(int, int[][])}).
	 * The method is implemented using the corresponding method of the
	 *  underlying {@link ConditionalMutualInfoCalculatorMultiVariateGaussian}
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below, and the other papers referenced in
	 * {@link AnalyticNullDistributionComputer#computeSignificance()}
	 * (in particular Geweke),
	 * for a description of how this is done for TE and conditional MI.
	 * Basically, the null distribution is a chi-square distribution.
	 * </p>
	 * 
	 * @return ChiSquareMeasurementDistribution object which describes
	 * the proportion of TE scores from the null distribution
	 *  which have higher or equal conditional MIs to our actual value.
	 * @see {@link ConditionalMutualInfoCalculatorMultiVariateGaussian#computeSignificance()}
	 * @see "J.T. Lizier, 'JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems', 2014."
	 * @throws Exception
	 */
	public ChiSquareMeasurementDistribution computeSignificance()
			throws Exception {
		return ((ConditionalMutualInfoCalculatorMultiVariateGaussian) condMiCalc).computeSignificance();
	}
}
