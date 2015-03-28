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

package infodynamics.measures.continuous;

import java.util.Vector;

import infodynamics.utils.MatrixUtils;

/**
 * A Multivariate Transfer Entropy (TE) calculator (implementing
 * {@link TransferEntropyCalculatorMultiVariate})
 * which is affected using a 
 * given Conditional Mutual Information (MI) calculator (implementing
 * {@link ConditionalMutualInfoCalculatorMultiVariate}) to make the calculations.
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculatorMultiVariate},
 * except that in the constructor(s) for this class the implementation for
 * a {@link ConditionalMutualInfoCalculatorMultiVariate} must be supplied.
 * </p>
 * 
 * <p>This class <i>may</i> be used directly, however users are advised that
 * several child classes are available which already plug-in the various
 * conditional MI estimators
 * to provide TE calculators (taking specific caution associated with
 * each type of estimator):</p>
 * <ul>
 * 	<li>{@link infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorMultiVariateGaussian}</li>
 * 	<li>{@link infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov}</li>
 * </ul>
 * 
 * <p>This implementation inherits from the {@link TransferEntropyCalculatorViaCondMutualInfo},
 * so the univariate methods for adding observations etc are still available, however
 * they will throw Exceptions if this calculator was not initialised for single
 * dimensional data sets.
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 *  <li>J.T. Lizier, J. Heinzle, A. Horstmann, J.-D. Haynes, M. Prokopenko,
 *  <a href="http://dx.doi.org/10.1007/s10827-010-0271-2">
 *  "Multivariate information-theoretic measures reveal directed information
 *  structure and task relevant changes in fMRI connectivity"</a>,
 *  Journal of Computational Neuroscience, vol. 30, pp. 85-107, 2011.</li>
 * </ul>
 *  
 * @author Joseph Lizier, <a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>
 * @see TransferEntropyCalculatorViaCondMutualInfo
 * @see TransferEntropyCalculatorMultiVariate
 * @see TransferEntropyCalculator
 */
public class TransferEntropyCalculatorMultiVariateViaCondMutualInfo
		extends TransferEntropyCalculatorViaCondMutualInfo
		// which means we implement TransferEntropyCalculator 
		implements TransferEntropyCalculatorMultiVariate {

	/**
	 * Number of dimensions of the destination
	 */
	protected int destDimensions = 1;
	/**
	 * Number of dimensions of the source
	 */
	protected int sourceDimensions = 1;

	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalculatorClassName as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalculatorClassName fully qualified name of the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorMultiVariateViaCondMutualInfo(String condMiCalculatorClassName)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(condMiCalculatorClassName);
	}

	/**
	 * Construct a transfer entropy calculator using an instance of
	 * condMiCalcClass as the underlying conditional mutual information calculator.
	 * 
	 * @param condMiCalcClass the class which must implement
	 * 	{@link ConditionalMutualInfoCalculatorMultiVariate}
	 * @throws InstantiationException if the given class cannot be instantiated
	 * @throws IllegalAccessException if illegal access occurs while trying to create an instance
	 *   of the class
	 * @throws ClassNotFoundException if the given class is not found
	 */
	public TransferEntropyCalculatorMultiVariateViaCondMutualInfo(Class<ConditionalMutualInfoCalculatorMultiVariate> condMiCalcClass)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(condMiCalcClass);
	}

	/**
	 * Construct this calculator by passing in a constructed but not initialised
	 * underlying Conditional Mutual information calculator.
	 * 
	 * @param condMiCalc An instantiated conditional mutual information calculator.
	 * @throws Exception if the supplied calculator has not yet been instantiated.
	 */
	public TransferEntropyCalculatorMultiVariateViaCondMutualInfo(ConditionalMutualInfoCalculatorMultiVariate condMiCalc) throws Exception {
		super(condMiCalc);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#initialise(int, int)
	 */
	@Override
	public void initialise(int sourceDimensions, int destDimensions)
			throws Exception {
		initialise(sourceDimensions, destDimensions, k, k_tau, l, l_tau, delay);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate#initialise(int, int, int)
	 */
	@Override
	public void initialise(int k, int sourceDimensions, int destDimensions)
			throws Exception {
		initialise(sourceDimensions, destDimensions, k, k_tau, l, l_tau, delay);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo#initialise(int, int, int, int, int)
	 */
	@Override
	public void initialise(int k, int k_tau, int l, int l_tau, int delay)
			throws Exception {
		initialise(sourceDimensions, destDimensions, k, k_tau, l, l_tau, delay);
	}

	/**
	 * Initialise the calculator for re-use with new observations.
	 * New embedding parameters, source and dest dimensions
	 * and source-destination delay
	 * may be supplied here; all other parameters
	 * remain unchanged.
	 * 
	 * @param sourceDimensions dimensionality of the source variable
	 * @param destDimensions dimensionality of the destination variable
	 * @param k Length of destination past history to consider
	 * @param k_tau embedding delay for the destination variable
	 * @param l length of source past history to consider
	 * @param l_tau embedding delay for the source variable
	 * @param delay time lag between last element of source and destination next value
	 */
	public void initialise(int sourceDimensions, int destDimensions, int k, int k_tau, int l, int l_tau, int delay) throws Exception {
		if (delay < 0) {
			throw new Exception("Cannot compute TE with source-destination delay < 0");
		}
		this.sourceDimensions = sourceDimensions;
		this.destDimensions = destDimensions;
		this.k = k;
		this.k_tau = k_tau;
		this.l = l;
		this.l_tau = l_tau;
		this.delay = delay;
		
		// Now check which point we can start taking observations from in any
		//  addObservations call. These two integers represent the last
		//  point of the destination embedding, in the cases where the destination
		//  embedding itself determines where we can start taking observations, or
		//  the case where the source embedding plus delay is longer and so determines
		//  where we can start taking observations.
		int startTimeBasedOnDestPast = (k-1)*k_tau;
		int startTimeBasedOnSourcePast = (l-1)*l_tau + delay - 1;
		startTimeForFirstDestEmbedding = Math.max(startTimeBasedOnDestPast, startTimeBasedOnSourcePast);

		condMiCalc.initialise(l*sourceDimensions, destDimensions, k*destDimensions);
	}

	/**
	 * <p>Sets a single set of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimensions
	 * of both source and destination are 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param source univariate observations for the source variable
	 * @param destination univariate observations for the destination variable
	 * @throws Exception if initialised dimensions were not 1
	 */
	@Override
	public void setObservations(double[] source, double[] destination) throws Exception {
		
		if ((sourceDimensions != 1) || (destDimensions != 1)) {
			throw new Exception("Cannot call the univariate setObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
		}
		
		super.setObservations(source, destination);
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#setObservations(double[][], double[][])
	 */
	@Override
	public void setObservations(double[][] source, double[][] destination)
			throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to add here, the time series is too short
			throw new Exception("Not enough observations to set here given k, k_tau, l, l_tau and delay parameters");
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k, k_tau,
						startTimeForFirstDestEmbedding,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding + 1,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						source.length - startTimeForFirstDestEmbedding - 1);
		condMiCalc.setObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
	}

	/**
	 * <p>Adds a new set of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimensions
	 * of both source and destination are 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param source univariate observations for the source variable
	 * @param destination univariate observations for the destination variable
	 * @throws Exception if initialised dimensions were not 1
	 * @see {@link ChannelCalculator#addObservations(double[], double[])}
	 */
	@Override
	public void addObservations(double[] source, double[] destination) throws Exception {

		if ((sourceDimensions != 1) || (destDimensions != 1)) {
			throw new Exception("Cannot call the univariate addObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
		}
		super.addObservations(source, destination);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#addObservations(double[][], double[][])
	 */
	@Override
	public void addObservations(double[][] source, double[][] destination)
			throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to add here, the time series is too short
			// Don't throw an exception, do nothing since more observations
			//  can be added later.
			return;
		}
		double[][] currentDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(destination, k, k_tau,
						startTimeForFirstDestEmbedding,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(destination, 1,
						startTimeForFirstDestEmbedding + 1,
						destination.length - startTimeForFirstDestEmbedding - 1);
		double[][] currentSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(source, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						source.length - startTimeForFirstDestEmbedding - 1);
		condMiCalc.addObservations(currentSourcePastVectors, currentDestNextVectors, currentDestPastVectors);
	}

	/**
	 * <p>Adds a new sub-series of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimensions
	 * of both source and destination are 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param source univariate observations for the source variable
	 * @param destination univariate observations for the destination variable
	 * @throws Exception if initialised dimensions were not 1
	 * @see {@link ChannelCalculator#addObservations(double[], double[], int, int)}
	 */
	@Override
	public void addObservations(double[] source, double[] destination,
			int startTime, int numTimeSteps) throws Exception {
		if ((sourceDimensions != 1) || (destDimensions != 1)) {
			throw new Exception("Cannot call the univariate addObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
		}
		super.addObservations(source, destination, startTime, numTimeSteps);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#addObservations(double[][], double[][], int, int)
	 */
	@Override
	public void addObservations(double[][] source, double[][] destination,
			int startTime, int numTimeSteps) throws Exception {
		if (source.length != destination.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					source.length, destination.length));
		}
		if (source.length < startTime + numTimeSteps) {
			// There are not enough observations given the arguments here
			throw new Exception("Not enough observations to set here given startTime and numTimeSteps parameters");
		}
		addObservations(MatrixUtils.selectRows(source, startTime, numTimeSteps),
					    MatrixUtils.selectRows(destination, startTime, numTimeSteps));
	}

	/**
	 * <p>Adds a new sub-series of <b>univariate</b> observations to compute the PDFs from.
	 * Can only be called on this multivariate calculator if the dimensions
	 * of both source and destination are 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param source univariate observations for the source variable
	 * @param destination univariate observations for the destination variable
	 * @throws Exception if initialised dimensions were not 1
	 * @see {@link ChannelCalculator#setObservations(double[], double[], boolean[], boolean[])}
	 */
	public void setObservations(double[] source, double[] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		if ((sourceDimensions != 1) || (destDimensions != 1)) {
			throw new Exception("Cannot call the univariate setObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
		}
		super.setObservations(source, destination, sourceValid, destValid);
	}
	
	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#setObservations(double[][], double[][], boolean[], boolean[])
	 */
	@Override
	public void setObservations(double[][] source, double[][] destination,
			boolean[] sourceValid, boolean[] destValid) throws Exception {
		
		Vector<int[]> startAndEndTimePairs = computeStartAndEndTimePairs(sourceValid, destValid);
		
		// We've found the set of start and end times for this pair
		startAddObservations();
		for (int[] timePair : startAndEndTimePairs) {
			int startTime = timePair[0];
			int endTime = timePair[1];
			addObservations(source, destination, startTime, endTime - startTime + 1);
		}
		finaliseAddObservations();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.ChannelCalculatorMultiVariate#setObservations(double[][], double[][], boolean[][], boolean[][])
	 */
	@Override
	public void setObservations(double[][] source, double[][] destination,
			boolean[][] sourceValid, boolean[][] destValid) throws Exception {
		boolean[] jointSourceValid = MatrixUtils.andRows(sourceValid);
		boolean[] jointDestValid = MatrixUtils.andRows(destValid);
		setObservations(source, destination, jointSourceValid, jointDestValid);
	}

	/**
	 * <p>Computes the local values of the transfer entropy
	 *  for each valid observation in the supplied univariate observations
	 * Can only be called on this multivariate calculator if the dimensions
	 * of both source and destination are 1, otherwise throws an exception</p>
	 * 
	 * {@inheritDoc}
	 * 
	 * @param newSourceObservations univariate observations for the source variable
	 * @param newDestObservations univariate observations for the destination variable
	 * @throws Exception if initialised dimensions were not 1
	 */
	public double[] computeLocalUsingPreviousObservations(double[] newSourceObservations,
			double[] newDestObservations) throws Exception {

		if ((sourceDimensions != 1) || (destDimensions != 1)) {
			throw new Exception("Cannot call the univariate computeLocalUsingPreviousObservations if you " +
					"have initialised with dimension > 1 for either source or destination");
		}
		return super.computeLocalUsingPreviousObservations(newSourceObservations, newDestObservations);
	}
	
	@Override
	public double[] computeLocalUsingPreviousObservations(double[][] newSourceObservations,
			double[][] newDestObservations) throws Exception {
		if (newSourceObservations.length != newDestObservations.length) {
			throw new Exception(String.format("Source and destination lengths (%d and %d) must match!",
					newSourceObservations.length, newDestObservations.length));
		}
		if (newDestObservations.length < startTimeForFirstDestEmbedding + 2) {
			// There are no observations to compute for here
			return new double[newDestObservations.length];
		}
		double[][] newDestPastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, k, k_tau,
						startTimeForFirstDestEmbedding,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newDestNextVectors =
				MatrixUtils.makeDelayEmbeddingVector(newDestObservations, 1,
						startTimeForFirstDestEmbedding + 1,
						newDestObservations.length - startTimeForFirstDestEmbedding - 1);
		double[][] newSourcePastVectors = 
				MatrixUtils.makeDelayEmbeddingVector(newSourceObservations, l, l_tau,
						startTimeForFirstDestEmbedding + 1 - delay,
						newSourceObservations.length - startTimeForFirstDestEmbedding - 1);
		double[] local = condMiCalc.computeLocalUsingPreviousObservations(
						newSourcePastVectors, newDestNextVectors, newDestPastVectors);
		// Pad the front of the array with zeros where local TE isn't defined:
		double[] localsToReturn = new double[local.length + startTimeForFirstDestEmbedding + 1];
		System.arraycopy(local, 0, localsToReturn, startTimeForFirstDestEmbedding + 1, local.length);
		return localsToReturn;
	}	
}
