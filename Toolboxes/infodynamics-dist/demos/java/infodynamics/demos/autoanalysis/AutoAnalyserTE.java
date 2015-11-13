/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2015, Joseph T. Lizier
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

package infodynamics.demos.autoanalysis;

import infodynamics.measures.continuous.ChannelCalculatorCommon;
import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian;
import infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;
import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete;

import javax.swing.JOptionPane;
import javax.swing.event.DocumentListener;

import java.awt.event.ActionListener;
import java.awt.event.MouseListener;

/**
 * This class provides a GUI to build a simple transfer entropy calculation,
 *  and supply the code to execute it.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserTE extends AutoAnalyser
	implements ActionListener, DocumentListener, MouseListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	protected static final String DISCRETE_PROPNAME_K = "k_HISTORY";
	
	/**
	 * Constructor to initialise the GUI for TE
	 */
	protected void makeSpecificInitialisations() {
		
		// Set up the properties for TE:
		measureAcronym = "TE";
		appletTitle = "JIDT Transfer Entropy Auto-Analyser"; 
		
		// Discrete:
		discreteClass = TransferEntropyCalculatorDiscrete.class;
		discreteProperties = new String[] {
				DISCRETE_PROPNAME_K,
				DISCRETE_PROPNAME_BASE
		};
		discretePropertyDefaultValues = new String[] {
				"1",
				"2"
		};
		discretePropertyDescriptions = new String[] {
				"Destination history embedding length",
				"Number of discrete states available for each variable (i.e. 2 for binary)"
		};
		
		// Continuous:
		abstractContinuousClass = TransferEntropyCalculator.class;
		// Common properties for all continuous calcs:
		commonContPropertyNames = new String[] {
				TransferEntropyCalculator.K_PROP_NAME
		};
		commonContPropertiesFieldNames = new String[] {
				"K_PROP_NAME"
		};
		commonContPropertyDescriptions = new String[] {
				"Destination history embedding length (k_HISTORY)"
		};
		// Gaussian properties:
		gaussianProperties = new String[] {
				TransferEntropyCalculator.K_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.DELAY_PROP_NAME, // Not common to Kernel
		};
		gaussianPropertiesFieldNames = new String[] {
				"TransferEntropyCalculator.K_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.DELAY_PROP_NAME", // Not common to Kernel
		};
		gaussianPropertyDescriptions = new String[] {
				"Destination history embedding delay (k_TAU)",
				"Source history embedding length (l_HISTORY)",
				"Source history embeding delay (l_TAU)",
				"Delay from source to destination (in time steps)"
		};
		// Kernel:
		kernelProperties = new String[] {
				TransferEntropyCalculatorKernel.KERNEL_WIDTH_PROP_NAME,
				TransferEntropyCalculatorKernel.DYN_CORR_EXCL_TIME_NAME,
				TransferEntropyCalculatorKernel.NORMALISE_PROP_NAME,			
		};
		kernelPropertiesFieldNames = new String[] {
				"KERNEL_WIDTH_PROP_NAME",
				"DYN_CORR_EXCL_TIME_NAME",
				"NORMALISE_PROP_NAME"			
		};
		kernelPropertyDescriptions = new String[] {
				"Kernel width to be used in the calculation. <br/>If the property " +
						TransferEntropyCalculatorKernel.NORMALISE_PROP_NAME +
						" is set, then this is a number of standard deviations; " +
						"otherwise it is an absolute value.",
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
						"0 (default) means no dynamic exclusion window",
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not  (recommended)",
		};
		// KSG (Kraskov):
		kraskovProperties = new String[] {
				TransferEntropyCalculator.K_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.L_TAU_PROP_NAME, // Not common to Kernel
				TransferEntropyCalculator.DELAY_PROP_NAME, // Not common to Kernel
				ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
				ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
				TransferEntropyCalculatorKraskov.PROP_KRASKOV_ALG_NUM,
				TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
				TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX,
				TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX,
				TransferEntropyCalculatorKraskov.PROP_RAGWITZ_NUM_NNS
		};
		kraskovPropertiesFieldNames = new String[] {
				"TransferEntropyCalculator.K_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.L_TAU_PROP_NAME", // Not common to Kernel
				"TransferEntropyCalculator.DELAY_PROP_NAME", // Not common to Kernel
				"ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE",
				"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS",
				"PROP_KRASKOV_ALG_NUM",
				"PROP_AUTO_EMBED_METHOD",
				"PROP_K_SEARCH_MAX",
				"PROP_TAU_SEARCH_MAX",
				"PROP_RAGWITZ_NUM_NNS"
		};
		kraskovPropertyDescriptions = new String[] {
				"Destination history embedding delay (k_TAU)",
				"Source history embedding length (l)",
				"Source history embeding delay (l_TAU)",
				"Delay from source to destination (in time steps)",
				"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not (recommended)",
				"Number of k nearest neighbours to use <br/>in the full joint kernel space in the KSG algorithm",
				"Standard deviation for an amount <br/>of random Gaussian noise to add to each variable, " +
						"to avoid having neighbourhoods with artificially large counts. <br/>" +
						"(\"false\" may be used to indicate \"0\".). The amount is added in after any normalisation.",
				"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
						"0 (default) means no dynamic exclusion window",
				"<br/>Norm type to use in KSG algorithm between the points in each marginal space. <br/>Options are: " +
						"\"MAX_NORM\" (default), otherwise \"EUCLIDEAN\" or \"EUCLIDEAN_SQUARED\" (both equivalent here)",
				"Number of parallel threads to use <br/>in computation: an integer > 0 or \"USE_ALL\" " +
						"(default, to indicate to use all available processors)",
				"Which KSG algorithm to use (1 or 2)",
				"Method to automatically determine embedding lengths (k_HISTORY,l_HISTORY)<br/> and delays (k_TAU, l_TAU) for " +
						"destination and potentially source time-series. Default is \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_NONE +
						"\" meaning values are set manually; other values include: <br/>  -- \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ +
						"\" for use of the Ragwitz criteria for both source and destination (searching up to \"" + TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX +
						"\" and \"" + TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX + "\"); <br/>  -- \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY +
						"\" for use of the Ragwitz criteria for the destination only. <br/>Use of values other than \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_NONE +
						"\" leads to any previous settings for embedding lengths and delays for the destination and perhaps source to be overwritten after observations are supplied",
				"Max. embedding length to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD + ")",
				"Max. embedding delay to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD + ")",
				"Number of k nearest neighbours for <br/>Ragwitz auto embedding (if used; defaults to match property \"k\")"
		};
		
	}

	/**
	 * Method to assign and initialise our continuous calculator class
	 */
	protected ChannelCalculatorCommon assignCalcObjectContinuous(String selectedCalcType) throws Exception {
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
			return new TransferEntropyCalculatorGaussian();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
			return new TransferEntropyCalculatorKraskov();
		} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
			return new TransferEntropyCalculatorKernel();
		} else {
			throw new Exception("No recognised continuous calculator selected: " +
					selectedCalcType);
		}

	}

	/**
	 * Method to assign and initialise our discrete calculator class
	 */
	protected DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception {
		String kPropValueStr, basePropValueStr;
		try {
			kPropValueStr = propertyValues.get(DISCRETE_PROPNAME_K);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot find a value for property " + DISCRETE_PROPNAME_K);
			return null;
		}
		try {
			basePropValueStr = propertyValues.get(DISCRETE_PROPNAME_BASE);
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot find a value for property " + DISCRETE_PROPNAME_BASE);
			return null;
		}
		int k = Integer.parseInt(kPropValueStr);
		int base = Integer.parseInt(basePropValueStr);
		
		return new DiscreteCalcAndArguments(
				new TransferEntropyCalculatorDiscrete(base, k),
				base,
				base + ", " + k);
	}

	protected void setObservations(ChannelCalculatorCommon calc,
			double[] source, double[] dest) throws Exception {
		// We know this is a TransferEntropyCalculator
		TransferEntropyCalculator teCalc = (TransferEntropyCalculator) calc;
		teCalc.setObservations(source, dest);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyserTE();
	}
}
