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

package infodynamics.networkinference.interregional;

import infodynamics.measures.continuous.ChannelCalculatorMultiVariate;
import infodynamics.measures.continuous.MutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.TransferEntropyCalculatorMultiVariate;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.ParsedProperties;
import infodynamics.utils.RandomGenerator;

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintStream;
import java.util.Properties;


/**
 * <p>Compute a significance score for an inter-regional channel measure
 * between two sets, based on looking at the channel measure between elements taken 
 * jointVars1 and jointVars2 at a time.
 * </p>
 * <p>Usage:
 * 	<ol>
 * 		<li>newInstance()</li>
 * 		<li>initialise</li>
 * 		<li>setObservations</li>
 * 		<li>computeMean, or computeSignificance, or computeLocals</li>
 *  </ol>
 * </p>
 * 
 * @author Joseph Lizier, joseph.lizier at gmail.com 
 *
 */
public abstract class InterregionalChannelMeasure {

	protected ChannelCalculatorMultiVariate channelCalc;
	/**
	 * The number of joint variables to consider from each region in each calculation
	 */
	protected int jointVars1; // destination
	protected int jointVars2; // source
	/**
	 * The total number of variables in each region
	 */
	protected boolean useAverageRegion1 = false;
	protected int totalVars1;
	protected boolean useAverageRegion2 = false;
	protected int totalVars2;
	protected int maxNumSubsets;
	protected boolean debug;
	protected boolean writeResults;

	protected double[][] region1;
	protected double[][] region2;
	protected boolean allValid;
	protected boolean validityForIndividualElements;
	protected boolean[] jointValidity1;
	protected boolean[] jointValidity2;
	protected boolean[][] individualValidity1;
	protected boolean[][] individualValidity2;
	
	protected double lastAverage;
	
	protected RandomGenerator rg;
	protected long seed;
	
	// Subsets to use
	protected int[][] subsetsForEachRegion1; 
	protected int[][] subsetsForEachRegion2; 
	protected int numOfSets;
	
	// Track whether each subset pair is known to be significant or not
	protected double[] significanceLevelsForEachSubset;
	// Track the t-score for each subset.
	protected double[] tScoreForEachSubject;

	private String calculatorClass;
	protected Properties calculatorProperties;

	public static final String PROP_CALCULATOR_PROPERTIES_PREFIX = "props.interregionalChannel.calculatorProperties.";
	public static final String PROP_CALCULATOR_CLASS = "props.interregionalChannel.calculator";
	public static final String PROP_JOINT_VARS_1 = "props.interregionalChannel.jointVars1";
	public static final String PROP_JOINT_VARS_2 = "props.interregionalChannel.jointVars2";
	public static final String PROP_NUM_SUBSETS = "props.interregionalChannel.maxNumSets";
	public static final String PROP_SEED = "props.interregionalChannel.seed";
	public static final String PROP_GET_SIGNIFICANCE = "props.interregionalChannel.directRun.getSignificance";
	public static final String PROP_SIG_REORDERINGS = "props.interregionalChannel.directRun.numReorderingsForSignificance";
	public static final String PROP_SIGNIFICANCE_COMPARE_POPULATIONS = "props.interregionalChannel.directRun.compareSubsetPopulationsForSignificance";
	public static final String PROP_WRITE_RESULTS = "props.interregionalChannel.writeResultsToStdOut";
	public static final String PROP_DIRECT_FILE1 = "props.interregionalChannel.directRun.file1";
	public static final String PROP_DIRECT_FILE2 = "props.interregionalChannel.directRun.file2";
	public static final String PROP_DIRECT_OUTPUT_PREFIX = "props.interregionalChannel.directRun.outputPrefix";
	public static final String PROP_DEBUG = "props.debug";
	
	public static final String JOINT_VARS_SELECT_AVERAGE = "Av";

	/**
	 * Data structure to store a distribution of channel measurements. 
	 * 
	 * @author Joseph Lizier
	 *
	 */
	public class ChannelMeasurementDistribution {
		public double[] channelMeasurements;
		public double mean;
		public double std;
	}
	
	/**
	 * Data structure to store the means and standard deviations of local values of the 
	 *  channel measurement in time. The averages are taken across a number of subsets
	 *  selected from each region.
	 * 
	 * @author Joseph Lizier
	 *
	 */
	public class LocalChannelMeasurementDistribution {
		public double[] meanLocalChannelMeasurements;
		public double[] stdLocalChannelMeasurements;
		public boolean hasLocalsForSignificantChannelsOnly;
		public int countOfSignificantChannels;
		public double[] meanLocalForSignificantChannelsOnly;
		public double[] stdLocalForSignificantChannelsOnly;
	}

	public static InterregionalChannelMeasure newInstance(ParsedProperties props) throws Exception {
		// Generate a new instance based on the calculator name
		String calcName = props.getStringProperty(PROP_CALCULATOR_CLASS);
		return newInstance(calcName);
	}
	
	public static InterregionalChannelMeasure newInstance(String calcName) throws Exception {
		// Generate a new instance based on the calculator name
		ChannelCalculatorMultiVariate channelCalcInstance = (ChannelCalculatorMultiVariate)
			Class.forName(calcName).newInstance();
		if (TransferEntropyCalculatorMultiVariate.class.isInstance(channelCalcInstance)) {
			return new InterregionalTransferEntropy();
		} else if (MutualInfoCalculatorMultiVariate.class.isInstance(channelCalcInstance)) {
			return new InterregionalMutualInfo();
		} else {
			throw new Exception("Calculator name not recognised");
		}
	}

	protected InterregionalChannelMeasure() {
		super();
	}
	
	public void initialise(ParsedProperties props) throws Exception {
		setCalculatorName(props.getStringProperty(PROP_CALCULATOR_CLASS));
		useAverageRegion1 = props.getProperty(PROP_JOINT_VARS_1).
								equalsIgnoreCase(JOINT_VARS_SELECT_AVERAGE);
		if (useAverageRegion1) {
			jointVars1 = 1;
		} else {
			setJointVars1(props.getIntProperty(PROP_JOINT_VARS_1));
		}
		useAverageRegion2 = props.getProperty(PROP_JOINT_VARS_2).
								equalsIgnoreCase(JOINT_VARS_SELECT_AVERAGE);
		if (useAverageRegion2) {
			jointVars2 = 1;
		} else {
			setJointVars2(props.getIntProperty(PROP_JOINT_VARS_2));
		}
		setDebug(props.getBooleanProperty(PROP_DEBUG));
		if (useAverageRegion1 && useAverageRegion2) {
			// Since we're taking the average of each region, 
			//  there are no subsets to take
			setMaxNumSubsets(1);
		} else {
			setMaxNumSubsets(props.getIntProperty(PROP_NUM_SUBSETS));
		}
		setSeed(props.getLongProperty(PROP_SEED));
		setWriteResults(props.getBooleanProperty(PROP_WRITE_RESULTS));
		readInCalculatorProperties(props);
		initialise();
	}

	public void initialise(int jointVars1, int jointVars2, int maxNumSubsets) throws Exception {
		setJointVars1(jointVars1);
		setJointVars2(jointVars2);
		useAverageRegion1 = false;
		useAverageRegion2 = false;
		setMaxNumSubsets(maxNumSubsets);
		initialise();
	}

	public void initialise() throws Exception {
		// Create our channel calculator
		channelCalc = (ChannelCalculatorMultiVariate)
					Class.forName(calculatorClass).newInstance();
		channelCalc.setDebug(debug);
		setPropertiesOnCalculator(true);
		rg = new RandomGenerator();
		rg.setSeed(seed);
		subsetsForEachRegion1 = null;
		subsetsForEachRegion2= null;
		significanceLevelsForEachSubset = null;
		tScoreForEachSubject = null;
		numOfSets = 0;
	}
	
	private void readInCalculatorProperties(ParsedProperties props) {
		calculatorProperties = new Properties();
		for (Object keyObject : props.getProperties().keySet()) {
			String key = (String) keyObject;
			if (key.startsWith(PROP_CALCULATOR_PROPERTIES_PREFIX)) {
				// Add this key to the calculator properties
				String propertyName = key.replaceFirst(PROP_CALCULATOR_PROPERTIES_PREFIX, "");
				calculatorProperties.setProperty(propertyName, props.getProperty(key));
			}
		}
	}

	public void setObservations(double[][] region1, double[][] region2) {
		setRegion1And2Data(region1, region2);
		allValid = true;
		finaliseSetObservations();
	}
	
	public void setObservations(double[][] region1, double[][] region2,
			boolean[] validity1, boolean[] validity2) {
		setRegion1And2Data(region1, region2);
		jointValidity1 = validity1;
		jointValidity2 = validity2;
		allValid = false;
		validityForIndividualElements = false;
		finaliseSetObservations();
	}

	public void setObservations(double[][] region1, double[][] region2,
			boolean[][] validity1, boolean[][] validity2) {
		setRegion1And2Data(region1, region2);
		individualValidity1 = validity1;
		individualValidity2 = validity2;
		allValid = false;
		validityForIndividualElements = true;
		finaliseSetObservations();
	}
	
	protected void setRegion1And2Data(double[][] region1, double[][] region2) {
		if (useAverageRegion1) {
			// We're taking the average across all of the variables in the region
			this.region1 = new double[region1.length][1];
			try {
				MatrixUtils.copyIntoColumn(this.region1, 0,
						MatrixUtils.meansOfRows(region1));
			} catch (Exception e) {
				// This should never happen
				throw new RuntimeException(e);
			}
		} else {
			// We're using all of the individual variables
			this.region1 = region1;
		}
		if (useAverageRegion2) {
			// We're taking the average across all of the variables in the region
			this.region2 = new double[region2.length][1];
			try {
				MatrixUtils.copyIntoColumn(this.region2, 0,
						MatrixUtils.meansOfRows(region2));
			} catch (Exception e) {
				// This should never happen
				throw new RuntimeException(e);
			}
		} else {
			// We're using all of the individual variables
			this.region2 = region2;
		}
	}
	
	protected void finaliseSetObservations() {
		totalVars1 = region1[0].length;
		totalVars2 = region2[0].length;
	}

	/**
	 * Initialise the calculator used for the calculation between each subset.
	 * Should be overridden by each subclass if if needs a different type of initialisation.
	 *
	 */
	protected void initialiseCalculator() throws Exception {
		channelCalc.initialise(jointVars1, jointVars2);
	}
	
	/**
	 * Set properties for the calculator object to match those we're holding.
	 * 
	 * @throws Exception
	 */
	protected void setPropertiesOnCalculator(boolean forcePrintProperties)
		throws Exception {
		
		boolean oldDebug = debug;
		channelCalc.setDebug(debug || forcePrintProperties);
		for (Object keyObject : calculatorProperties.keySet()) {
			channelCalc.setProperty((String) keyObject,
					calculatorProperties.getProperty((String) keyObject));
		}
		channelCalc.setDebug(oldDebug);
	}
	
	/**
	 * Check that the subsets for each region have been generated
	 *
	 */
	protected void checkSubsetsAreGenerated() {
		// Grab a set of subsets for each region if required:
		if (subsetsForEachRegion1 == null) {
			subsetsForEachRegion1 = rg.generateNRandomSets(totalVars1, jointVars1, maxNumSubsets); 
			subsetsForEachRegion2 = rg.generateNRandomSets(totalVars2, jointVars2, maxNumSubsets); 
			numOfSets = Math.min(subsetsForEachRegion1.length, subsetsForEachRegion2.length);
		}
	}
	
	/**
	 * Compute the mean transfer entropy across maxNumSubsets subsets of 
	 * joint variables from each region
	 * 
	 * @return
	 * @throws Exception
	 */
	public ChannelMeasurementDistribution computeMean() throws Exception {
		// Set up the calculator
		setPropertiesOnCalculator(false);
		
		checkSubsetsAreGenerated();
		
		double[] tes = new double[numOfSets];
		// Store all the results
		for (int s = 0; s < numOfSets; s++) {
			// Select the subset of region 1:
			double[][] r1Subset = MatrixUtils.selectColumns(region1, subsetsForEachRegion1[s]);
			// Select the subset of region 2:
			double[][] r2Subset = MatrixUtils.selectColumns(region2, subsetsForEachRegion2[s]);
			initialiseCalculator();
			// Then add them in as observations
			if (allValid) {
				// all the observations are valid
				channelCalc.setObservations(r1Subset, r2Subset);
			} else if (validityForIndividualElements) {
				// we've been given validity for each individual sub-variable
				// compute the joint validities
				boolean[] computedJointValidity1 = MatrixUtils.andRowsOverSelectedColumns(
						individualValidity1, subsetsForEachRegion1[s]);
				boolean[] computedJointValidity2 = MatrixUtils.andRowsOverSelectedColumns(
						individualValidity2, subsetsForEachRegion2[s]);
				channelCalc.setObservations(r1Subset, r2Subset, computedJointValidity1, computedJointValidity2);
			} else {
				// we've got joint validity for each time series
				channelCalc.setObservations(r1Subset, r2Subset, jointValidity1, jointValidity2);
			}
			
			// Compute the TEs
			tes[s] = channelCalc.computeAverageLocalOfObservations();
			if (writeResults) {
				System.out.print("Done average for subset " + s + " ");
				printSubset(System.out, subsetsForEachRegion1[s]);
				System.out.print(" -> ");
				printSubset(System.out, subsetsForEachRegion2[s]);
				System.out.println(": average = " +
						tes[s]);
			}
		}

		ChannelMeasurementDistribution ted = new ChannelMeasurementDistribution();
		ted.channelMeasurements = tes;
		ted.mean = MatrixUtils.mean(tes);
		ted.std = MatrixUtils.stdDev(tes, ted.mean);
		lastAverage = ted.mean;
		if (writeResults) {
			System.out.println("Average across " + numOfSets + " subsets was " + lastAverage);
		}
		return ted;
	}
	
	/**
	 * <p>Generate reorderingsForSignificance reorderings for the observations supplied
	 * here. Supplied where a user wants to generate the reordering once, and thereafter
	 * supply to a number of InterregionalChannelMeasure objects.
	 * </p> 
	 * 
	 * @param reorderingsForSignificance
	 * @return
	 * @throws Exception
	 */
	public int[][] generateReorderings(int reorderingsForSignificance) throws Exception {
		// Compute the reorderings
		int numObservationsToReorder = computeNumObservationsToReorder();
		return rg.generateDistinctRandomPerturbations(numObservationsToReorder,
				reorderingsForSignificance);
	}
	
	/**
	 * Compute the significance of the average of the channel measure across all subsets, 
	 *  by comparing it to the population of averages that could be obtained where the
	 *  first region is reordered against the second.
	 * Each average in the population of reordering is made with each underlying subset
	 *  reordered the same way. 
	 * 
	 * @param reorderingsForSignificance the number of reorderings to consider
	 * @return the distribution of SxP surrogate measurements for each subset S and permutation P 
	 * @throws Exception
	 */
	public MeasurementDistributionPermutationsOverSubsets computeSignificance(int reorderingsForSignificance) throws Exception {
		return computeSignificance(generateReorderings(reorderingsForSignificance));
	}
	
	/**
	 * Compute the significance of the average of the channel measure across all subsets, 
	 *  by comparing it to the population of averages that could be obtained where the
	 *  first region is reordered against the second.
	 * Each average in the population of reordering is made with each underlying subset
	 *  reordered the same way. 
	 * 
	 * @param reorderings the reorderings to consider
	 * @return the distribution of SxP surrogate measurements for each subset S and permutation P 
	 * @throws Exception
	 */
	public MeasurementDistributionPermutationsOverSubsets computeSignificance(int[][] reorderings) throws Exception {
		int reorderingsForSignificance = reorderings.length;
		
		// Set up the calculator
		setPropertiesOnCalculator(false);
		
		checkSubsetsAreGenerated();

		double[] measureForEachSet = new double[numOfSets];
		double[][] reorderedTesForAllSubsets = new double[numOfSets][];
		
		// Track the significance for each subset
		significanceLevelsForEachSubset = new double[numOfSets];
		tScoreForEachSubject = new double[numOfSets];
		
		// Store all the results
		for (int s = 0; s < numOfSets; s++) {
			// Select the subset of region 1:
			double[][] r1Subset = MatrixUtils.selectColumns(region1, subsetsForEachRegion1[s]);
			// Select the subset of region 2:
			double[][] r2Subset = MatrixUtils.selectColumns(region2, subsetsForEachRegion2[s]);
			initialiseCalculator();
			// Then add them in as observations
			if (allValid) {
				// all the observations are valid
				channelCalc.setObservations(r1Subset, r2Subset);
			} else if (!validityForIndividualElements) {
				// we've got joint validity for each time series
				channelCalc.setObservations(r1Subset, r2Subset, jointValidity1, jointValidity2);
			} else {
				// we've been given validity for each individual sub-variable
				channelCalc.setObservations(r1Subset, r2Subset, individualValidity1, individualValidity2);
			}
			
			// Compute the measure for set s
			measureForEachSet[s] = channelCalc.computeAverageLocalOfObservations();

			EmpiricalMeasurementDistribution measDist;
			if (allValid || !validityForIndividualElements) {
				// Ask the TE calculator to work out the significance for us
				measDist = channelCalc.computeSignificance(reorderings);
			} else {
				// We need to explicitly reorder including the individual validities.
				// Can't ask the calculator to do it, as it will only use the source-dest
				//  pairings it originally used, which won't match across all subsets
				measDist = new EmpiricalMeasurementDistribution(reorderingsForSignificance);
				measDist.actualValue = measureForEachSet[s];
				int countWhereReorderedIsMoreSignificantThanOriginal = 0;
				for (int p = 0; p < reorderingsForSignificance; p++) {
					double[][] reorderedR1Subset =
						MatrixUtils.extractSelectedTimePointsReusingArrays(r1Subset, reorderings[p]);
					boolean[][] reorderedIndividualValidity1 = 
						MatrixUtils.extractSelectedTimePointsReusingArrays(individualValidity1, reorderings[p]);
					initialiseCalculator();
					channelCalc.setObservations(reorderedR1Subset, r2Subset,
							reorderedIndividualValidity1, individualValidity2);
					measDist.distribution[p] = channelCalc.computeAverageLocalOfObservations();
					if (measDist.distribution[p] >= measureForEachSet[s]) {
						countWhereReorderedIsMoreSignificantThanOriginal++;
					}
				}
				measDist.pValue = (double) countWhereReorderedIsMoreSignificantThanOriginal /
										(double) reorderingsForSignificance;
			}
			reorderedTesForAllSubsets[s] = measDist.distribution;
			significanceLevelsForEachSubset[s] = measDist.pValue;
			tScoreForEachSubject[s] = measDist.getTSscore();
			if (writeResults) {
				double meanOfDist = MatrixUtils.mean(measDist.distribution);
				double stdOfDist = MatrixUtils.stdDev(measDist.distribution, meanOfDist);
				System.out.print("Significance for subset " + s + " ");
				printSubset(System.out, subsetsForEachRegion1[s]);
				System.out.print(" -> ");
				printSubset(System.out, subsetsForEachRegion2[s]);
				System.out.printf(" was %.3f (%.4f compared to %.4f +/- %.4f, t=%.4f, factor of %.2f)\n",
						measDist.pValue,
						measDist.actualValue,
						meanOfDist, stdOfDist, tScoreForEachSubject[s], 
						measDist.actualValue / meanOfDist);
			}
		}

		// Compute the averages over subsets for each reordering, and
		//  compute the significance = probability that we could have gotten a 
		//  value >= by chance alone
		MeasurementDistributionPermutationsOverSubsets measDist = 
			new MeasurementDistributionPermutationsOverSubsets(
					reorderedTesForAllSubsets, measureForEachSet);
		
		lastAverage = measDist.actualValue;
		if (writeResults) {
			System.out.println("Significance across " + reorderingsForSignificance + " permutations of " + 
					numOfSets + " subsets averaged was " + measDist.pValue);
		}
		return measDist;
	}

	/**
	 * Compute the average local transfer entropy values
	 *  (local in time, averaged across all subsets)
	 * from the source region to the destination.
	 * 
	 * @return object containing means and std deviations
	 *  of local values in time.
	 */
	public LocalChannelMeasurementDistribution computeLocals() throws Exception {
		return computeLocals(-1.0, true);
	}
	
	/**
	 * Compute the average local transfer entropy values
	 *  (local in time, averaged across all subsets)
	 * from the source region to the destination.
	 * 
	 * @param cutoff level at which to count an individual channel as significant
	 *   for the return distribution. If &lt; 0 don't look at the significance of individual
	 *   channels. If &gt;= 0 one of the computeSignificance methods should have already
	 *   been called.
	 * @param cutoffIsSignificance whether the cutoff is the cutoff significance of the TE value
	 *  or on the cutoff z score
	 * @return object containing means and std deviations
	 *  of local values in time.
	 */
	public LocalChannelMeasurementDistribution computeLocals(double cutoff, boolean cutoffIsSignificance) throws Exception {
		// Set up the calculator
		setPropertiesOnCalculator(false);
		
		checkSubsetsAreGenerated();

		lastAverage = 0.0;
		double[] averageLocalTes = null;
		double[] stdLocalTes = null;
		int countOfSignificantChannels = 0;
		double[] averageLocalTesForSigChannelsOnly = null;
		double[] stdLocalTesForSigChannelsOnly = null;
		// Store all the results
		for (int s = 0; s < numOfSets; s++) {
			// Select the subset of region 1:
			double[][] r1Subset = MatrixUtils.selectColumns(region1, subsetsForEachRegion1[s]);
			// Select the subset of region 2:
			double[][] r2Subset = MatrixUtils.selectColumns(region2, subsetsForEachRegion2[s]);
			initialiseCalculator();
			// Then add them in as observations
			if (allValid) {
				// all the observations are valid
				channelCalc.setObservations(r1Subset, r2Subset);
			} else if (validityForIndividualElements) {
				// we've been given validity for each individual sub-variable
				// compute the joint validities
				boolean[] computedJointValidity1 = MatrixUtils.andRowsOverSelectedColumns(
						individualValidity1, subsetsForEachRegion1[s]);
				boolean[] computedJointValidity2 = MatrixUtils.andRowsOverSelectedColumns(
						individualValidity2, subsetsForEachRegion2[s]);
				channelCalc.setObservations(r1Subset, r2Subset, computedJointValidity1, computedJointValidity2);
				// TODO I don't think adding in place will work for this case
				//  because each subset pair will have a different number of valid observations,
				//  which means averaging over all the local time series doesn't make much sense ??
				// Could average into 0..t-1 but not make a contribution if we don't 
				//  have a valid local value there.
			} else {
				// we've got joint validity for each time series
				channelCalc.setObservations(r1Subset, r2Subset, jointValidity1, jointValidity2);
			}
			
			// Compute the local TEs
			channelCalc.setDebug(true);
			double[] localTes = channelCalc.computeLocalOfPreviousObservations();
			channelCalc.setDebug(debug);
			lastAverage += channelCalc.getLastAverage();
			
			// And add these into the sums for each time step
			if (averageLocalTes == null) {
				// Create the space to store the average 
				//  and std dev of local TEs.
				// Creating here because outside the loop
				//  we don't know the length of arrays to create
				averageLocalTes = new double[localTes.length];
				stdLocalTes = new double[localTes.length];
			}
			MatrixUtils.addInPlace(averageLocalTes, localTes);
			// stdLocalTes holds the sum of squares for now
			MatrixUtils.addSquaresInPlace(stdLocalTes, localTes);
			if ((cutoff >= 0.0) && (significanceLevelsForEachSubset != null)) {
				// We've been asked to give a result for significant channels only,
				//  and the significance of each channel has been calculated.
				if (averageLocalTesForSigChannelsOnly == null) {
					averageLocalTesForSigChannelsOnly = new double[localTes.length];
					stdLocalTesForSigChannelsOnly = new double[localTes.length];
				}
				// Store these locals if the channel was previously found to be significant
				if ((cutoffIsSignificance && significanceLevelsForEachSubset[s] <= cutoff) ||
					(!cutoffIsSignificance && tScoreForEachSubject[s] >= cutoff)) {
					countOfSignificantChannels++;
					MatrixUtils.addInPlace(averageLocalTesForSigChannelsOnly, localTes);
					// stdLocalTes holds the sum of squares for now
					MatrixUtils.addSquaresInPlace(stdLocalTesForSigChannelsOnly, localTes);
				}
			}
			if (writeResults) {
				System.out.print("Done locals for subset " + s + " ");
				printSubset(System.out, subsetsForEachRegion1[s]);
				System.out.print(" -> ");
				printSubset(System.out, subsetsForEachRegion2[s]);
				System.out.printf(": (average = %.5f)\n",
						channelCalc.getLastAverage());
			}
		}

		// Now take the average across all subsets
		//  and compute the standard devs
		for (int t = 0; t < averageLocalTes.length; t++) {
			averageLocalTes[t] /= (double) numOfSets;
			stdLocalTes[t] = Math.sqrt(
					stdLocalTes[t] / (double) numOfSets - 
					averageLocalTes[t] * averageLocalTes[t]);
		}
		lastAverage /= (double) numOfSets;
		if (writeResults) {
			System.out.println("Average (from locals) across " + numOfSets + " subsets was "
					+ lastAverage);
		}
		LocalChannelMeasurementDistribution lcmd = new
			LocalChannelMeasurementDistribution();
		lcmd.meanLocalChannelMeasurements = averageLocalTes;
		lcmd.stdLocalChannelMeasurements = stdLocalTes;
		lcmd.hasLocalsForSignificantChannelsOnly = (cutoff >= 0.0) && (significanceLevelsForEachSubset != null);
		// Include the results for significant channels only if requested
		if (lcmd.hasLocalsForSignificantChannelsOnly) {
			lcmd.countOfSignificantChannels = countOfSignificantChannels;
			if (countOfSignificantChannels > 0) {
				// Find the mean and std for each point in time
				for (int t = 0; t < averageLocalTes.length; t++) {
					averageLocalTesForSigChannelsOnly[t] /= (double) countOfSignificantChannels;
					stdLocalTesForSigChannelsOnly[t] = Math.sqrt(
							stdLocalTesForSigChannelsOnly[t] / (double) countOfSignificantChannels - 
							averageLocalTesForSigChannelsOnly[t] * averageLocalTesForSigChannelsOnly[t]);
				}
			}
			lcmd.meanLocalForSignificantChannelsOnly = averageLocalTesForSigChannelsOnly;
			lcmd.stdLocalForSignificantChannelsOnly = stdLocalTesForSigChannelsOnly;
		}
		return lcmd;
	}
	
	/**
	 * Compute the number of observations that need to be reordered in a significance
	 *  computation, based on how the validities have been supplied here, and
	 *  what type of calculator we are used (this is why it is left abstract).
	 * 
	 * @return
	 * @throws Exception
	 */
	protected abstract int computeNumObservationsToReorder() throws Exception;
	
	/**
	 * Compute the time indices for the local observations generated here
	 * 
	 * @return
	 * @throws Exception
	 */
	public abstract int[] computeTimeIndicesForLocalValues() throws Exception;
	
	/**
	 * Print the subset out to the given stream
	 * 
	 * @param out
	 * @param subset
	 */
	private void printSubset(PrintStream out, int[] subset) {
		out.print("{");
		boolean printedOne = false;
		for (int i = 0; i < subset.length; i++) {
			if (printedOne) {
				out.print(", ");
			}
			out.print(subset[i]);
			printedOne = true;
		}
		out.print("}");
	}

	public void setCalculatorName(String calculatorName) {
		calculatorClass = calculatorName;
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setJointVars1(int jointVars1) {
		this.jointVars1 = jointVars1;
	}

	public void setJointVars2(int jointVars2) {
		this.jointVars2 = jointVars2;
	}
	
	public void setMaxNumSubsets(int maxNumSubsets) {
		this.maxNumSubsets = maxNumSubsets;
	}
	
	public void setSeed(long seed) {
		this.seed = seed;	
	}

	private void setWriteResults(boolean writeResults) {
		this.writeResults = writeResults;
	}

	public double getLastAverage() {
		return lastAverage;
	}
	
	/**
	 * A simple method to evaluate the interregional channel measure 
	 *  for a given pair of files.
	 * 
	 * @param args first arg is a property filename,
	 *  second and third (optional) are the source and target
	 *  data files (if they are not supplied on command line
	 *  they must be named in the properties file)
	 */
	public static void main(String[] args) throws Exception {
		
		if (args.length < 1) {
			System.err.println("Usage: InterregionalChannelMeasure <propsFile> [<dataFile1> <dataFile2>]");
			return;
		}
		String propertiesFilename = args[0];

		Properties properties = new Properties();
		properties.load(new FileInputStream(new File(propertiesFilename)));
		ParsedProperties props = new ParsedProperties(properties);
		InterregionalChannelMeasure interregionalCalculator =
			InterregionalChannelMeasure.newInstance(props);

		// Load in the data
		// for the one region pair
		ArrayFileReader afr1, afr2;
		if (args.length < 3) {
			// Names of *both* data files were not supplied
			//  on command line
			afr1 = new ArrayFileReader(props.getStringProperty(PROP_DIRECT_FILE1));
			afr2 = new ArrayFileReader(props.getStringProperty(PROP_DIRECT_FILE2));
		} else {
			// Names of the data files were suppied on command line: 
			afr1 = new ArrayFileReader(args[1]);
			afr2 = new ArrayFileReader(args[2]);
		}
		double[][] region1 = afr1.getDouble2DMatrix();
		double[][] region2 = afr2.getDouble2DMatrix();

		interregionalCalculator.initialise(props);
		// System.out.println("Using all " + region1.length + " observations");
		interregionalCalculator.setObservations(region1, region2);
		
		ChannelMeasurementDistribution meanDistribution =
			interregionalCalculator.computeMean();
		
		/*
		System.out.print("MeanActuals,StdActuals");
		if (props.getBooleanProperty(PROP_GET_SIGNIFICANCE)) {
			System.out.println(",MeanPsOvS,StdPsOvS,t,degFree");
		}
		*/
		try {
			String outputPrefix = props.getStringProperty(PROP_DIRECT_OUTPUT_PREFIX);
			System.out.printf("%s", outputPrefix);
		} catch (Exception e) {
			// there was no output prefix to print out
		}
		
		System.out.printf("%.4f,%.4f", meanDistribution.mean, meanDistribution.std);
		if (props.getBooleanProperty(PROP_GET_SIGNIFICANCE)) {
			// Generate the reorderings once and use a common set across all subsets
			//  and subject;
			// Unless we are using only pre and post steps, in which case the number
			//  of observations will be different across subjects, so we will need
			//  to regenerate the reorderings subject by subject.
			int permutations = props.getIntProperty(PROP_SIG_REORDERINGS);
			MeasurementDistributionPermutationsOverSubsets measDist =
				interregionalCalculator.computeSignificance(permutations);
			if (props.getBooleanProperty(PROP_SIGNIFICANCE_COMPARE_POPULATIONS)) {
				// Non-standard method that we investigated with Bernstein on Tracking data set.
				//  Compares the *population* of subset measures against the population of surrogate measures.
				//  It is not recommended that you do this.
				int numSubsets = measDist.avDistributionForSubsets.length;
				double meanOfMeanPermsForSubsets = MatrixUtils.mean(measDist.avDistributionForSubsets);
				// TODO I think using the std of mean permutations for subsets makes the standard dev here
				//  artificially small. This is because we're taken a std of averages, rather than of raw values.
				// We could approximately correct for this by multiplying by sqrt(numSubsets).
				// Of course this alters the number of samples used in its calculation to numSubsets * numPermutations,
				//  which alters s2SqrOnN and degFreedon and tVal.
				double stdOfMeanPermsForSubsets = MatrixUtils.stdDev(
						measDist.avDistributionForSubsets, meanOfMeanPermsForSubsets);
				// t-test of populations with unequal variance: see
				//  http://en.wikipedia.org/wiki/Student%27s_t-test
				double s1SqrOnN = (meanDistribution.std * meanDistribution.std) / numSubsets;
				double s2SqrOnN = (stdOfMeanPermsForSubsets * stdOfMeanPermsForSubsets) / numSubsets;
				double s = Math.sqrt(s1SqrOnN + s2SqrOnN);
				double degFreedom = (s1SqrOnN + s2SqrOnN) * (s1SqrOnN + s2SqrOnN) * (numSubsets - 1) /
									( s1SqrOnN*s1SqrOnN + s2SqrOnN*s2SqrOnN);
				double tValForPopsOverSubsets = (meanDistribution.mean - meanOfMeanPermsForSubsets) / s;
				System.out.printf(",%.4f,%.4f,%.4f,%.1f", meanOfMeanPermsForSubsets,
						stdOfMeanPermsForSubsets, tValForPopsOverSubsets, degFreedom);
			} else {
				// Default statistical significance measurement:
				// Comparing the single interregional value to the distribution of interregional values
				//  for each permutation.
				double tValForInterregionalMeasure = (meanDistribution.mean - measDist.getMeanOfDistribution()) / measDist.getStdOfDistribution();
				System.out.printf(",%.4f,%.4f,%.4f,%d", measDist.getMeanOfDistribution(),
						measDist.getStdOfDistribution(), tValForInterregionalMeasure, permutations);
			}
		}
		System.out.println();
	}
}
