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
 * Calculators implementing this interface must provide a
 *  {@link #computeSignificance()} method to compute
 *  the statistical significance of their measurement, returning an analytically
 *  determined distribution {@link AnalyticNullDistributionComputer}.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public interface AnalyticNullDistributionComputer {

	/**
	 * Generate an <b>analytic</b> distribution of what a given measure would look like,
	 * under a null hypothesis that the source values of our
	 * samples had no relation to the destination value (possibly
	 * in the context of a conditional value).
	 * 
	 * <p>See Section II.E "Statistical significance testing" of 
	 * the JIDT paper below for a description of how this is done for MI,
	 * conditional MI and TE as per the references below.
	 * </p>
	 * 
	 * <p><b>References:</b><br/>
	 *  <ul>
	 *   <li>J.T. Lizier, "JIDT: An information-theoretic
	 *    toolkit for studying the dynamics of complex systems", 2014.</li>
	 * 	 <li>Brillinger, <a href="http://www.stat.berkeley.edu/~brill/Papers/MIBJPS.pdf">
	 * 		"Some data analyses using mutual information"</a>,
	 * 		Brazilian Journal of Probability and Statistics, <b>18</b>, p. 163, (2004)</li>
	 * 	 <li>Cheng et al., <a href="http://www.jds-online.com/file_download/112/JDS-369.pdf">
	 * 		"Data Information in Contingency Tables: A Fallacy of Hierarchical Loglinear Models"</a>,
	 * 		Journal of Data Science, <b>4</b>, p. 387 (2006).</li>
	 *   <li>Geweke, <a href="http://dx.doi.org/10.1080/01621459.1982.10477803">
	 *   	"Measurement of Linear Dependence and Feedback between Multiple Time Series"</a>,
	 *   	Journal of the American Statistical Association, <b>77</b>, p. 304-313 (1982).</li>
	 * 	 <li>Barnett and Bossomaier, <a href="http://arxiv.org/abs/1205.6339">
	 * 		"Transfer Entropy as a Log-likelihood Ratio"</a>,
	 * 		Physical Review Letters, <b>109</b>, p. 138105+ (2012).</li>
     * </ul>
	 * 
	 * @return the analytic distribution of channel measure scores under this null hypothesis.
	 * @throws Exception
	 */
	public AnalyticMeasurementDistribution computeSignificance() throws Exception;
}
