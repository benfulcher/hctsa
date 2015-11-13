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

package infodynamics.measures.continuous.kraskov;

import java.util.Hashtable;

import infodynamics.measures.continuous.ActiveInfoStorageCalculator;
import infodynamics.measures.continuous.ConditionalMutualInfoCalculatorMultiVariate;
import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo;

/**
 * <p>Computes the differential transfer entropy (TE) between two univariate
 *  <code>double[]</code> time-series of observations
 *  (implementing {@link TransferEntropyCalculator}),
 *  using Kraskov-Stoegbauer-Grassberger (KSG) estimation (see references below).
 *  This estimator is realised here by plugging in
 *  a {@link ConditionalMutualInfoCalculatorMultiVariateKraskov}
 *  as the calculator into the parent class {@link TransferEntropyCalculatorViaCondMutualInfo}.</p>
 *  
 * <p>Crucially, the calculation is performed by examining
 * neighbours in the full joint space (as specified by Frenzel and Pompe,
 * and Gomez-Herrero et al.)
 * rather than two MI calculators.</p>
 * 
 * <p>Usage is as per the paradigm outlined for {@link TransferEntropyCalculator},
 * with:
 * <ul>
 * 	<li>The constructor step is either a simple call to {@link #TransferEntropyCalculatorKraskov()},
 *      or else specifies which KSG algorithm to implement via
 *      {@link #TransferEntropyCalculatorKraskov(String)};</li>
 * 	<li>{@link #setProperty(String, String)} allowing properties defined for both
 * 		{@link TransferEntropyCalculator#setProperty(String, String)} and
 *      {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)}
 *      as outlined
 *      in {@link TransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)});
 *      as well as for {@link #PROP_KRASKOV_ALG_NUM}.
 *      Embedding parameters may be automatically determined as per the Ragwitz criteria
 *      by setting the property {@link #PROP_AUTO_EMBED_METHOD} to {@link #AUTO_EMBED_METHOD_RAGWITZ}
 *      or {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY}
 *      (plus additional parameter settings for this).</li>
 *      </li>
 *  <li>Computed values are in <b>nats</b>, not bits!</li>
 *  </ul>
 * </p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>T. Schreiber, <a href="http://dx.doi.org/10.1103/PhysRevLett.85.461">
 * "Measuring information transfer"</a>,
 *  Physical Review Letters 85 (2) pp.461-464, 2000.</li>
 * 	<li>Frenzel and Pompe, <a href="http://dx.doi.org/10.1103/physrevlett.99.204101">
 * 	"Partial Mutual Information for Coupling Analysis of Multivariate Time Series"</a>,
 * 	Physical Review Letters, <b>99</b>, p. 204101+ (2007).</li>
 * 	<li>G. Gomez-Herrero, W. Wu, K. Rutanen, M. C. Soriano, G. Pipa, and R. Vicente,
 * 	<a href="http://arxiv.org/abs/1008.0539">
 * 	"Assessing coupling dynamics from an ensemble of time series"</a>,
 * 	arXiv:1008.0539 (2010).</li>
 * 	<li>Kraskov, A., Stoegbauer, H., Grassberger, P., 
 *   <a href="http://dx.doi.org/10.1103/PhysRevE.69.066138">"Estimating mutual information"</a>,
 *   Physical Review E 69, (2004) 066138.</li>
 *  <li>J. T. Lizier, M. Prokopenko and A. Zomaya,
 *  <a href="http://dx.doi.org/10.1103/PhysRevE.77.026110">
 *  "Local information transfer as a spatiotemporal filter for complex systems"</a>
 *  Physical Review E 77, 026110, 2008.</li>
 * 	<li>Ragwitz and Kantz, "Markov models from data by simple nonlinear time series
 *  	predictors in delay embedding spaces", Physical Review E, vol 65, 056201 (2002).</li>
 * </ul>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 * @see TransferEntropyCalculator
 * @see ConditionalMutualInfoCalculatorMultiVariateKraskov
 */
public class TransferEntropyCalculatorKraskov
	extends TransferEntropyCalculatorViaCondMutualInfo {

	/**
	 * Class name for KSG conditional MI estimator via KSG algorithm 1
	 */
	public static final String COND_MI_CALCULATOR_KRASKOV1 = ConditionalMutualInfoCalculatorMultiVariateKraskov1.class.getName();
	/**
	 * Class name for KSG conditional MI estimator via KSG algorithm 2
	 */
	public static final String COND_MI_CALCULATOR_KRASKOV2 = ConditionalMutualInfoCalculatorMultiVariateKraskov2.class.getName();
	
	/**
	 * Property for setting which underlying Kraskov-Grassberger algorithm to use  (1 or 2).
	 * Will only be applied at the next initialisation.
	 */
	public final static String PROP_KRASKOV_ALG_NUM = "ALG_NUM";
	
	/**
	 * Which Kraskov algorithm number we are using
	 */
	protected int kraskovAlgorithmNumber = 1;
	protected boolean algChanged = false;
	/**
	 * Storage for the properties ready to pass onto the underlying conditional MI calculators should they change 
	 */
	protected Hashtable<String,String> props;

	/**
	 * Property name for specifying which (if any) auto-embedding method to use.
	 * Valid values include {@link #AUTO_EMBED_METHOD_RAGWITZ}, {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY}
	 *  and {@link #AUTO_EMBED_METHOD_NONE}.
	 * Defaults to {@link #AUTO_EMBED_METHOD_NONE}
	 */
	public static final String PROP_AUTO_EMBED_METHOD = "AUTO_EMBED_METHOD";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  no auto embedding should be done (i.e. to use manually supplied parameters)
	 */
	public static final String AUTO_EMBED_METHOD_NONE = "NONE";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the Ragwitz optimisation technique should be used for automatic embedding
	 *  for both source and destination time-series
	 */
	public static final String AUTO_EMBED_METHOD_RAGWITZ = "RAGWITZ";
	/**
	 * Valid value for the property {@link #PROP_AUTO_EMBED_METHOD} indicating that
	 *  the Ragwitz optimisation technique should be used for automatic embedding
	 *  for the destination time-series only
	 */
	public static final String AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY = "RAGWITZ_DEST_ONLY";
	/**
	 * Internal variable tracking what type of auto embedding (if any)
	 *  we are using
	 */
	protected String autoEmbeddingMethod = AUTO_EMBED_METHOD_NONE;
	
	/**
	 * Property name for maximum embedding lengths (i.e. k for destination, and l for source if we're auto-embedding
	 *  the source as well) for the auto-embedding search. Defaults to 1
	 */
	public static final String PROP_K_SEARCH_MAX = "AUTO_EMBED_K_SEARCH_MAX";
	/**
	 * Internal variable for storing the maximum embedding length to search up to for
	 *  automating the parameters.
	 */
	protected int k_search_max = 1;

	/**
	 * Property name for maximum embedding delay (i.e. k_tau for destination, and l_tau for source if we're auto-embedding
	 *   the source as well) for the auto-embedding search. Defaults to 1
	 */
	public static final String PROP_TAU_SEARCH_MAX = "AUTO_EMBED_TAU_SEARCH_MAX";
	/**
	 * Internal variable for storing the maximum embedding delay to search up to for
	 *  automating the parameters.
	 */
	protected int tau_search_max = 1;

	/**
	 * Property name for the number of nearest neighbours to use for the auto-embedding search (Ragwitz criteria).
	 * Defaults to match the value in use for {@link MutualInfoCalculatorMultiVariateKraskov#PROP_K}
	 */
	public static final String PROP_RAGWITZ_NUM_NNS = "AUTO_EMBED_RAGWITZ_NUM_NNS";
	/**
	 * Internal variable for storing the number of nearest neighbours to use for the
	 *  auto embedding search (Ragwitz criteria)
	 */
	protected int ragwitz_num_nns = 1;
	/** 
	 * Internal variable to track whether the property {@link #PROP_RAGWITZ_NUM_NNS} has been
	 * set yet
	 */
	protected boolean ragwitz_num_nns_set = false;

	/**
	 * Creates a new instance of the Kraskov-estimate style transfer entropy calculator
	 * 
	 * Uses algorithm 1 by default, as per Gomez-Herro et al.
	 * 
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public TransferEntropyCalculatorKraskov() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(COND_MI_CALCULATOR_KRASKOV1);
		kraskovAlgorithmNumber = 1;
		props = new Hashtable<String,String>();
	}

	/**
	 * Creates a new instance of the Kraskov-Grassberger style transfer entropy calculator,
	 *  with the supplied conditional MI calculator name
	 * 
	 * @param calculatorName fully qualified name of the underlying MI class.
	 *    Must be {@link #COND_MI_CALCULATOR_KRASKOV1} or {@link #COND_MI_CALCULATOR_KRASKOV2}
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws InstantiationException 
	 *
	 */
	public TransferEntropyCalculatorKraskov(String calculatorName) throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		super(calculatorName);
		// Now check that it was one of our Kraskov-Grassberger calculators:
		if (calculatorName.equalsIgnoreCase(COND_MI_CALCULATOR_KRASKOV1)) {
			kraskovAlgorithmNumber = 1;
		} else if (calculatorName.equalsIgnoreCase(COND_MI_CALCULATOR_KRASKOV2)) {
			kraskovAlgorithmNumber = 2;
		} else {
			throw new ClassNotFoundException("Must be an underlying Kraskov-Grassberger conditional MI calculator");
		}
		props = new Hashtable<String,String>();
	}

	/* (non-Javadoc)
	 * @see infodynamics.measures.continuous.TransferEntropyCalculatorViaCondMutualInfo#initialise(int, int, int, int, int)
	 */
	@Override
	public void initialise(int k, int k_tau, int l, int l_tau, int delay)
			throws Exception {
		if (algChanged) {
			// The algorithm number was changed in a setProperties call:
			String newCalcName = COND_MI_CALCULATOR_KRASKOV1;
			if (kraskovAlgorithmNumber == 2) {
				newCalcName = COND_MI_CALCULATOR_KRASKOV2;
			}
			@SuppressWarnings("unchecked")
			Class<ConditionalMutualInfoCalculatorMultiVariate> condMiClass = 
					(Class<ConditionalMutualInfoCalculatorMultiVariate>) Class.forName(newCalcName);
			ConditionalMutualInfoCalculatorMultiVariate newCondMiCalc = condMiClass.newInstance();
			construct(newCondMiCalc);
			// Set the properties for the Kraskov MI calculators (may pass in properties for our super class
			//  as well, but they should be ignored)
			for (String key : props.keySet()) {
				newCondMiCalc.setProperty(key, props.get(key));
			}
			algChanged = false;
		}
		
		super.initialise(k, k_tau, l, l_tau, delay);
	}

	/**
	 * Sets properties for the TE calculator.
	 *  New property values are not guaranteed to take effect until the next call
	 *  to an initialise method. 
	 *  
	 * <p>Valid property names, and what their
	 * values should represent, include:</p>
	 * <ul>
	 * 		<li>{@link #PROP_KRASKOV_ALG_NUM} -- which Kraskov algorithm number to use (1 or 2).</li>
	 * 		<li>{@link #PROP_AUTO_EMBED_METHOD} -- method by which the calculator
	 * 		automatically determines the embedding history length ({@link #K_PROP_NAME})
	 * 		and embedding delay ({@link #TAU_PROP_NAME}) for destination and potentially source.
	 * 		Default is {@link #AUTO_EMBED_METHOD_NONE} meaning
	 * 		values are set manually; other accepted values include: {@link #AUTO_EMBED_METHOD_RAGWITZ} for use
	 * 		of the Ragwitz criteria for both source and destination (searching up to {@link #PROP_K_SEARCH_MAX} and 
	 * 		{@link #PROP_TAU_SEARCH_MAX}), and {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY} for use
	 * 		of the Ragwitz criteria for the destination only. 
	 * 		Use of any value other than {@link #AUTO_EMBED_METHOD_NONE}
	 * 		will lead to previous settings for embedding lengths and delays (via e.g. {@link #initialise(int, int)} or
	 * 		auto-embedding during previous calculations) for the destination and perhaps source to
	 * 		be overwritten after observations are supplied.</li>
	 * 		<li>{@link #PROP_K_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}) for the time-series to be embedded; default is 1</li>
	 * 		<li>{@link #PROP_TAU_SEARCH_MAX} -- maximum embedded history length to search
	 * 		up to if automatically determining the embedding parameters (as set by
	 * 		{@link #PROP_AUTO_EMBED_METHOD}) for the time-series to be embedded; default is 1</li>
	 * 		<li>{@link #PROP_RAGWITZ_NUM_NNS} -- number of nearest neighbours to use
	 * 		in the auto-embedding if the property {@link #PROP_AUTO_EMBED_METHOD}
	 * 		has been set to {@link #AUTO_EMBED_METHOD_RAGWITZ} or {@link #AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY}.
	 * 		Defaults to the property value
	 *      set for {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#PROP_K}</li>
	 * 		<li>Any properties accepted by {@link TransferEntropyCalculatorViaCondMutualInfo#setProperty(String, String)}</li>
	 * 		<li>Or properties accepted by the underlying
	 * 		{@link ConditionalMutualInfoCalculatorMultiVariateKraskov#setProperty(String, String)} implementation.</li>
	 * </ul>
	 * <p>One should set {@link ConditionalMutualInfoCalculatorMultiVariateKraskov#PROP_K} here, the number
	 *  of neighbouring points one should count up to in determining the joint kernel size.</p> 
	 * <p><b>Note:</b> further properties may be defined by child classes.</p>
	 * 
	 * <p>Unknown property values are ignored.</p>
	 * 
	 * @param propertyName name of the property
	 * @param propertyValue value of the property.
	 * @throws Exception if there is a problem with the supplied value).
	 */
	public void setProperty(String propertyName, String propertyValue)
			throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			int previousAlgNumber = kraskovAlgorithmNumber;
			kraskovAlgorithmNumber = Integer.parseInt(propertyValue);
			if ((kraskovAlgorithmNumber != 1) && (kraskovAlgorithmNumber != 2)) {
				throw new Exception("Kraskov algorithm number (" + kraskovAlgorithmNumber
						+ ") must be either 1 or 2");
			}
			if (kraskovAlgorithmNumber != previousAlgNumber) {
				algChanged = true;
			}
			if (debug) {
				System.out.println(this.getClass().getSimpleName() + ": Set property " + propertyName +
						" to " + propertyValue);
			}
		} else if (propertyName.equalsIgnoreCase(PROP_AUTO_EMBED_METHOD)) {
			// New method set for determining the embedding parameters
			autoEmbeddingMethod = propertyValue;
		} else if (propertyName.equalsIgnoreCase(PROP_K_SEARCH_MAX)) {
			// Set max embedding history length for auto determination of embedding
			k_search_max = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_TAU_SEARCH_MAX)) {
			// Set maximum embedding delay for auto determination of embedding
			tau_search_max = Integer.parseInt(propertyValue);
		} else if (propertyName.equalsIgnoreCase(PROP_RAGWITZ_NUM_NNS)) {
			// Set the number of nearest neighbours to use in case of Ragwitz auto embedding:
			ragwitz_num_nns = Integer.parseInt(propertyValue);
			ragwitz_num_nns_set = true;
		} else {
			// Assume it was a property for the parent class or underlying conditional MI calculator
			super.setProperty(propertyName, propertyValue);
			props.put(propertyName, propertyValue); // This will keep properties for the super class as well as the cond MI calculator, but this is ok
		}
	}

	@Override
	public String getProperty(String propertyName) throws Exception {
		if (propertyName.equalsIgnoreCase(PROP_KRASKOV_ALG_NUM)) {
			return Integer.toString(kraskovAlgorithmNumber);
		} else if (propertyName.equalsIgnoreCase(PROP_AUTO_EMBED_METHOD)) {
			return autoEmbeddingMethod;
		} else if (propertyName.equalsIgnoreCase(PROP_K_SEARCH_MAX)) {
			return Integer.toString(k_search_max);
		} else if (propertyName.equalsIgnoreCase(PROP_TAU_SEARCH_MAX)) {
			return Integer.toString(tau_search_max);
		} else if (propertyName.equalsIgnoreCase(PROP_RAGWITZ_NUM_NNS)) {
			if (ragwitz_num_nns_set) {
				return Integer.toString(ragwitz_num_nns);
			} else {
				return condMiCalc.getProperty(ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K);
			}
		} else {
			// Assume it was a property for the parent class or underlying conditional MI calculator
			return super.getProperty(propertyName);
		}
	}

	@Override
	public void preFinaliseAddObservations() throws Exception {
		// Automatically determine the embedding parameters for the given time series
		
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_NONE)) {
			return;
		}
		// Else we need to auto embed
		
		// If we need to check which embedding method later:
		// if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ) ||
		//		autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY)) {
		
		// Use a Kraskov AIS calculator to embed both time-series individually:
		ActiveInfoStorageCalculatorKraskov aisCalc = new ActiveInfoStorageCalculatorKraskov();
		// Set the properties for the underlying MI Kraskov calculator here to match ours:
		for (String key : props.keySet()) {
			aisCalc.setProperty(key, props.get(key));
		}
		// Set the auto-embedding properties as we require:
		aisCalc.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				ActiveInfoStorageCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ);
		aisCalc.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_K_SEARCH_MAX,
				Integer.toString(k_search_max));
		aisCalc.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_TAU_SEARCH_MAX,
				Integer.toString(tau_search_max));
		// In case !ragwitz_num_nns_set and our condMiCalc has a different default number of
		//  kNNs for Kraskov search than miCalc, we had best supply the number directly here:
		aisCalc.setProperty(ActiveInfoStorageCalculatorKraskov.PROP_RAGWITZ_NUM_NNS,
					getProperty(PROP_RAGWITZ_NUM_NNS));
		
		// Embed the destination with the Ragwitz criteria:
		if (debug) {
			System.out.println("Starting embedding of destination:");
		}
		aisCalc.initialise();
		aisCalc.startAddObservations();
		for (double[] destination : vectorOfDestinationTimeSeries) {
			aisCalc.addObservations(destination);
		}
		aisCalc.finaliseAddObservations();
		// Set the auto-embedding parameters for the destination:
		k = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
		k_tau = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
		if (debug) {
			System.out.printf("Embedding parameters for destination set to k=%d,k_tau=%d\n",
				k, k_tau);
		}
	
		if (autoEmbeddingMethod.equalsIgnoreCase(AUTO_EMBED_METHOD_RAGWITZ)) {
			// Embed the source also with the Ragwitz criteria:
			if (debug) {
				System.out.println("Starting embedding of source:");
			}
			aisCalc.initialise();
			aisCalc.startAddObservations();
			for (double[] source : vectorOfSourceTimeSeries) {
				aisCalc.addObservations(source);
			}
			aisCalc.finaliseAddObservations();
			// Set the auto-embedding parameters for the source:
			l = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.K_PROP_NAME));
			l_tau = Integer.parseInt(aisCalc.getProperty(ActiveInfoStorageCalculator.TAU_PROP_NAME));
			if (debug) {
				System.out.printf("Embedding parameters for source set to l=%d,l_tau=%d\n",
					l, l_tau);
			}
		}

		// Now that embedding parameters are finalised:
		setStartTimeForFirstDestEmbedding();
	}
}
