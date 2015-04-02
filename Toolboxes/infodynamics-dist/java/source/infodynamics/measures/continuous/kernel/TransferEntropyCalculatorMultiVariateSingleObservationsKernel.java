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

import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate;
import infodynamics.utils.MatrixUtils;

import java.util.Iterator;
import java.util.Vector;

/**
 * <p>
 * Extends {@link TransferEntropyCalculatorMultiVariateKernel} for
 * computing the differential transfer entropy (TE) between two <b>multivariate</b>
 *  <code>double[][]</code> time-series of observations
 *  using box-kernel estimation.
 * This calculator however only be used to add observation tuples, i.e.
 *  (source, destination next state, destination past)
 *  one at a time. This allows the user to specify the variable that
 *  should be used as the destination past state, for advanced applications
 *  (where the state is somehow captured differently to the past
 *  embedding vector). As an example, see Wang et al. (2012) below.
 * </p> 
 * 
 * <p>Javadocs are somewhat incomplete since this is 
 * a niche class. TODO Finish these properly.</p>
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculatorMultiVariateKernel}
 * (extending {@link TransferEntropyCalculatorMultiVariate}),
 * with:
 * <ul>
 * 	<li>The constructor step being a simple call to
 * 		{@link #TransferEntropyCalculatorMultiVariateSingleObservationsKernel()}.</li>
 *  <li>Adds additional {@link #initialise(double)} and {@link #initialiseAllDimensions(int, int, int)} options</li>
 *  <li>The addition of {@link #addSingleObservation(double[], double[], double[])}
 *      for adding single observations in.</li>
 *  </ul>
 * </p>
 * 
 * <p>
 * TODO Implement dynamic correlation exclusion with multiple observation sets. (see the
 *  way this is done in Plain calculator).
 * TODO Think about added error-trapping code to make sure the user only makes one type of addObservations call.
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J.T. Lizier, J. Heinzle, A. Horstmann, J.-D. Haynes, M. Prokopenko,
 *  <a href="http://dx.doi.org/10.1007/s10827-010-0271-2">
 *  "Multivariate information-theoretic measures reveal directed information
 *  structure and task relevant changes in fMRI connectivity"</a>,
 *  Journal of Computational Neuroscience, vol. 30, pp. 85-107, 2011.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>H. Kantz and T. Schreiber, "Nonlinear Time Series Analysis"
 *  (Cambridge University Press, Cambridge, MA, 1997).</li>
 *  <li>X. R. Wang, J. M. Miller, J. T. Lizier, M. Prokopenko, and L. F. Rossi,
 *  <a href="http://dx.doi.org/10.1371/journal.pone.0040084">
 *  "Quantifying and Tracing Information Cascades in Swarms"</a>,
 *  PLoS ONE 7, e40084+ (2012).</li>
 * </ul>
 *
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class TransferEntropyCalculatorMultiVariateSingleObservationsKernel
		extends TransferEntropyCalculatorMultiVariateKernel {

	/**
	 * Storage for destination history observations for addObservsations
	 */
	protected Vector<double[][]> vectorOfJointDestinationPastObservations;

	protected int destPastDimensions = 1;

	/**
	 * Construct an instance
	 */
	public TransferEntropyCalculatorMultiVariateSingleObservationsKernel() {
		super();
	}

	/**
	 * Initialises the calculator
	 * 
	 * @param epsilon kernel width
	 */
	public void initialise(double epsilon) throws Exception {
		this.kernelWidth = epsilon;
		initialise(1, 1); // assume 1 dimension in source and dest
	}

	@Override
	public void initialise(int sourceDimensions, int destDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		this.destPastDimensions = destDimensions; // assume same
		super.initialise(1); // Feeds k=1 to super and calls initialise();
	}

	/**
	 * Initialise routine where the number of dimensions considered
	 * to be part of the past state may also be supplied.
	 * 
	 * @param sourceDimensions
	 * @param destDimensions
	 * @param destPastDimensions
	 * @throws Exception
	 */
	public void initialiseAllDimensions(int sourceDimensions,
			int destDimensions, int destPastDimensions) throws Exception {
		this.destDimensions = destDimensions;
		this.sourceDimensions = sourceDimensions;
		this.destPastDimensions = destPastDimensions;
		
		// Mimic super.initialise(1) (it would replace dest and source
		//  dimensions if we're not careful)
		addedMoreThanOneObservationSet = false;
		k = 1;
		// Mimic super.initialise() (it would use k * destDimenions in the kernel estimator
		//  for destPast instead of destPastDimensions if we're not careful)
		teKernelEstimator.initialise(destPastDimensions,
				sourceDimensions, kernelWidth, kernelWidth);
		nextStateKernelEstimator.initialise(destDimensions, kernelWidth);
		destPastVectors = null;
		destNextVectors = null;
		sourceVectors = null;
		localProbNextCondPast = null;
	}

	/**
	 * Set the observations to compute the probabilities from 
	 * 
	 * @param source
	 * @param destination
	 */
	public void setObservations(double[][] source, double[][] destination,
			double[][] destinationPast) throws Exception {
		startAddObservations();
		addObservations(source, destination, destinationPast);
		finaliseAddObservations();
	}
	
	@Override
	public void startAddObservations() {
		vectorOfJointDestinationPastObservations = new Vector<double[][]>();
		super.startAddObservations();
	}

	/**
	 * Add observations of the joint source and destinations
	 * 
	 * @param source
	 * @param destination
	 * @param destinationPast
	 * @throws Exception
	 */
	public void addObservations(double[][] source, double[][] destination,
			double[][] destinationPast) throws Exception {
		if (destinationPast.length != destination.length) {
			throw new Exception(String.format("Destination past and destination lengths (%d and %d) must match!",
					destinationPast.length, destination.length));
		}
		int thisDestPastDimensions = destinationPast[0].length;
		if ((thisDestPastDimensions != destPastDimensions)) {
			throw new Exception("Cannot add observsations for destination past variables " +
					" of " + thisDestPastDimensions +
					" dimensions for TE calculator set up for " + destPastDimensions +
					" destination past dimensions");
		}
		if (vectorOfJointDestinationPastObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		vectorOfJointDestinationPastObservations.add(destinationPast);
		super.addObservations(source, destination);
	}

	/**
	 * Add a single observation of the joint source, destinations and
	 * destination past
	 * 
	 * @param source
	 * @param destination
	 * @param destinationPast
	 * @throws Exception
	 */
	public void addSingleObservation(double[] source, double[] destination,
			double[] destinationPast) throws Exception {
		int thisSourceDimensions = source.length;
		int thisDestDimensions = destination.length;
		int thisDestPastDimensions = destinationPast.length;
		if ((thisDestDimensions != destDimensions) ||
				(thisSourceDimensions != sourceDimensions) ||
				(thisDestPastDimensions != destPastDimensions)) {
			throw new Exception("Cannot add observsations for source, destination and destPast variables " +
					" of " + thisSourceDimensions + " and " + thisDestDimensions + " and " +
					thisDestPastDimensions +
					" dimensions respectively for TE calculator set up for " + sourceDimensions + ", " +
					destDimensions + " and " + destPastDimensions +
					" source, destination and destPast dimensions respectively");
		}
		if (vectorOfJointDestinationPastObservations == null) {
			// startAddObservations was not called first
			throw new RuntimeException("User did not call startAddObservations before addObservations");
		}
		// Now make the multidimensional arrays to add in
		double[][] sourceContainer = new double[1][];
		double[][] destContainer = new double[1][];
		double[][] destPastContainer = new double[1][];
		sourceContainer[0] = source;
		destContainer[0] = destination;
		destPastContainer[0] = destinationPast;
		vectorOfJointDestinationPastObservations.add(destPastContainer);
		super.addObservations(sourceContainer, destContainer);
	}

	@Override
	public void finaliseAddObservations() {
		// First work out the size to allocate the joint vectors, and do the allocation:
		totalObservations = 0;
		for (double[][] destination : vectorOfJointDestinationObservations) {
			// No need t jump k values in, since we've got the destination
			//  past values held separately
			totalObservations += destination.length;
		}
		destPastVectors = new double[totalObservations][destPastDimensions];
		destNextVectors = new double[totalObservations][destDimensions];
		sourceVectors = new double[totalObservations][sourceDimensions];
		
		// Construct the joint vectors from the given observations
		int startObservation = 0;
		Iterator<double[][]> iterator = vectorOfJointDestinationObservations.iterator();
		Iterator<double[][]> iteratorDestPast = vectorOfJointDestinationPastObservations.iterator();
		for (double[][] source : vectorOfJointSourceObservations) {
			double[][] destination = iterator.next();
			double[][] destinationPast = iteratorDestPast.next();
			// Add in all observations - no need to offset by k since
			//  we've got the destination past held separately.
			MatrixUtils.arrayCopy(destinationPast, 0, 0,
					destPastVectors, startObservation, 0, destinationPast.length,
					destPastDimensions);
			MatrixUtils.arrayCopy(destination, 0, 0,
					destNextVectors, startObservation, 0,
					destination.length, destDimensions);
			MatrixUtils.arrayCopy(source, 0, 0,
					sourceVectors, startObservation, 0,
					source.length, sourceDimensions);
			startObservation += destination.length;
		}
		
		// Now set the joint vectors in the kernel estimators
		teKernelEstimator.setObservations(destPastVectors, destNextVectors, sourceVectors);

		// Store whether there was more than one observation set:
		addedMoreThanOneObservationSet = vectorOfJointDestinationObservations.size() > 1;
		
		// And clear the vector of observations
		vectorOfJointSourceObservations = null;
		vectorOfJointDestinationObservations = null;
		vectorOfJointDestinationPastObservations = null;
	}

}
