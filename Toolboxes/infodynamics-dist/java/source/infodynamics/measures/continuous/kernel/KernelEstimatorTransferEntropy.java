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

package infodynamics.measures.continuous.kernel;

import infodynamics.utils.MatrixUtils;

/**
 * <p>Kernel estimator for computing PDFs ready for the transfer entropy.</p>
 * 
 * <p>Extends KernelEstimatorMultiVariate, using the super class to manage the history of the destination
 * variable, and adds the next state and source on top of this. Any calls to the super class methods will only
 * function on the joint history.
 * </p> 
 * 
 * <p>
 * TODO More thoroughly check the Javadocs here
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 *  <li>H. Kantz and T. Schreiber, "Nonlinear Time Series Analysis"
 *  (Cambridge University Press, Cambridge, MA, 1997).</li>
 * </ul>
 *
 * @see KernelEstimatorMultiVariate
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class KernelEstimatorTransferEntropy extends KernelEstimatorMultiVariate {

	private double suppliedKernelWidthSource;
	private double kernelWidthSourceInUse;
	
	private double[] destNext;
	private double[] source;
	
	// Store the current observations passed in the synchronized method getCount
	//  waiting for callbacks from the underlying kernel estimator
	private double destNextObs;
	private double sourceObs;
	// Counts of destPastNext, destPastSource, destPastNextSource to be filled in 
	//  with the callbacks
	private int countNextPast;
	private int countPastSource;
	private int countNextPastSource;
	
	public KernelEstimatorTransferEntropy() {
		super();
		// Make sure when get a callbacl when correlated points are found
		makeCorrelatedPointAddedCallback = true;
	}

	public void initialise(int dimensions, double epsilon) {
		super.initialise(dimensions, epsilon);
		this.suppliedKernelWidthSource = epsilon;
	}
	
	public void initialise(int dimensions, double epsilonDest, 
			double epsilonSource) {
		super.initialise(dimensions, epsilonDest);
		this.suppliedKernelWidthSource = epsilonSource;
	}
	
	public void setObservations(double[][] destPastVectors,
			double[] destNext, double[] source) {
		setObservations(destPastVectors);
		// epsilonInUse has been computed for the destination.
		// TODO We could compute and set it directly here so
		//  we don't have a mismatch between any of the vector
		//  variables.
		
		if (normalise) {
			double std = MatrixUtils.stdDev(source);
			kernelWidthSourceInUse = suppliedKernelWidthSource * std;
		} else {
			kernelWidthSourceInUse = suppliedKernelWidthSource;
		}

		this.source = source;
		this.destNext = destNext;
	}
	
	/**
	 * Compute the required counts for Transfer Entropy using kernel estimation on the destination's past.
	 * Use callbacks to check if the joint counts need to be incremented.
	 * 
	 * If observationTimeStep &lt; 0, then no dynamic correlation exclusion will be attempted
	 *
	 * @param destPast
	 * @param destNextObs
	 * @param sourceObs
	 * @param observationTimeStep
	 * @return
	 */
	public synchronized TransferEntropyKernelCounts getCount(
			double[] destPast, double destNextObs,
			double sourceObs, int observationTimeStep) {
		
		// Prepare for any callbacks
		countNextPast = 0;
		countPastSource = 0;
		countNextPastSource = 0;
		this.destNextObs = destNextObs;
		this.sourceObs = sourceObs;

		// Get the count, and have the joint counts filled in via callbacks
		int countPast;
		if (observationTimeStep < 0) {
			countPast = super.getCount(destPast);
		} else {
			countPast = super.getCount(destPast, observationTimeStep);
		}
		
		TransferEntropyKernelCounts teKernelCount = 
			new TransferEntropyKernelCounts(countPast, countNextPast, countPastSource, countNextPastSource);
		return teKernelCount;
	}

	public void setEpsSource(double epsilonSource) {
		this.suppliedKernelWidthSource = epsilonSource;
	}

	/**
	 * A callback for where a correlated point is found at
	 * correlatedTimeStep in the destination's past.
	 * Now check whether we need to increment the joint counts.
	 *
	 */
	protected void correlatedPointAddedCallback(int correlatedTimeStep) {
		boolean sourceMatches = false;
		if (Math.abs(sourceObs - source[correlatedTimeStep]) <= kernelWidthSourceInUse) {
			countPastSource++;
			sourceMatches = true;
		}
		// The epsilons across the destination variables should all be approximately 
		//  equal, so just use the first one.
		if (Math.abs(destNextObs - destNext[correlatedTimeStep]) <= kernelWidthsInUse[0]) {
			countNextPast++;
			if (sourceMatches) {
				countNextPastSource++;
			}
		}
	}

	/**
	 * A callback for where a correlated point is removed due to dynamic correlated exclusion.
	 * The removal is for the point at correlatedTimeStep in the destination's past.
	 * Now check whether we need to decrement the joint counts.
	 */
	protected void correlatedPointRemovedCallback(int removedCorrelatedTimeStep) {
		boolean sourceMatches = false;
		if (Math.abs(sourceObs - source[removedCorrelatedTimeStep]) <= kernelWidthSourceInUse) {
			countPastSource--;
			sourceMatches = true;
		}
		// The epsilons across the destination variables should all be approximately 
		//  equal, so just use the first one.
		if (Math.abs(destNextObs - destNext[removedCorrelatedTimeStep]) <= kernelWidthsInUse[0]) {
			countNextPast--;
			if (sourceMatches) {
				countNextPastSource--;
			}
		}
	}
}
