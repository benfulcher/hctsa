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
 * <p>Structure to hold a distribution of info-theoretic measurements,
 *  and a significance value for how an original measurement compared 
 *  with these.</p>
 * 
 * <p>While in theory this class could be directly used, it is it's
 * children which are intended to be used.</p>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MeasurementDistribution {
	/**
	 * Actual observed value of the measurement
	 */
	public double actualValue;
	/**
	 * Probability that surrogate measurement is greater than
	 *  the observed value.
	 * (Small pValue means that the observed value is highly significant).
	 */
	public double pValue;

	/**
	 * Allow empty constructor for internal use only when actualValue and
	 * pValue will be supplied later
	 */
	protected MeasurementDistribution() {
	}
	
	/**
	 * Construct with supplied actual value and p-value for it.
	 * 
	 * @param actualValue actual observed value
	 * @param pValue p-value that the surrogate measurement is larger
	 *  than the observed value
	 */
	public MeasurementDistribution(double actualValue, double pValue) {
		this.actualValue = actualValue;
		this.pValue = pValue;
	}
}
