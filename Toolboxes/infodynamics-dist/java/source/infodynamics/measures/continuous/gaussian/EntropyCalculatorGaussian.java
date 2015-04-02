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

import infodynamics.measures.continuous.EntropyCalculator;
import infodynamics.utils.MatrixUtils;

/**
 * <p>Computes the differential entropy of a given set of observations
 *  (implementing {@link EntropyCalculator}, assuming that
 *  the probability distribution function for these observations is Gaussian.</p>
 *  
 * <p>Usage is as per the paradigm outlined for {@link EntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to {@link #EntropyCalculatorGaussian()}.</li>
 * 	<li>The user can call {@link #setVariance(double)}
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
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class EntropyCalculatorGaussian implements EntropyCalculator {

	/**
	 * Variance of the most recently supplied observations, or set directly
	 */
	protected double variance;
	
	/**
	 * Whether we are in debug mode
	 */
	protected boolean debug;
	
	/**
	 * Construct an instance
	 */
	public EntropyCalculatorGaussian() {
		// Nothing to do
	}
	
	public void initialise() {
		// Nothing to do
	}

	public void setObservations(double[] observations) {
		variance = MatrixUtils.stdDev(observations);
		variance *= variance;
	}

	/**
	 * An alternative to {@link #setObservations(double[])}, allowing user to
	 * set the variance of the distribution for which we will compute the
	 *  entropy.
	 * 
	 * @param variance the variance of the univariate distribution.
	 */
	public void setVariance(double variance) {
		this.variance = variance;
	}
	
	/**
	 * Compute the entropy from the previously supplied observations, or
	 * based on the supplied variance.
	 * 
	 * <p>The entropy for a Gaussian-distribution random variable with
	 *  variance \sigma is 0.5*\log_e{2*pi*e*\sigma}.</p>
	 * 
	 * <p>Here we compute the entropy assuming that the recorded estimation of the
	 *  variance is correct (i.e. we will not make a bias correction for limited
	 *  observations here).</p>
	 * 
	 * @return the entropy of the previously provided observations or from the supplied
	 *   covariance matrix. Entropy returned in <b>nats</b>, not bits!
	 */
	public double computeAverageLocalOfObservations() {
		return 0.5 * Math.log(2.0*Math.PI*Math.E*variance);
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

}
