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

import infodynamics.utils.commonsmath3.special.Gamma;

/**
 * This class implements a number of static 
 * methods for mathematical functions.
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MathsUtils {

	private static final double EULER_MASCHERONI_CONSTANT = 0.5772156;
	
	private static int highestDigammaArgCalced = 0;
	private static final int NUM_STORED_DIGAMMAS = 10000; // commons.math to handle beyond this
	private static double[] storedDigammas = new double[NUM_STORED_DIGAMMAS];
	
	/**
	 * Returns the integer result of base^power
	 * 
	 * @param base base integer of the operation
	 * @param power power that base is raised to
	 * @return base raised to exponent power (rounded by integer operations)
	 */
	public static int power(int base, int power) {
		int result = 1;
		int absPower = Math.abs(power);
		for (int p = 0; p < absPower; p++) {
			result *= base;
		}
		if (power < 0) {
			// This will be zero for any base except 1 or -1
			result = 1 / result;
		}
		return result;
	}
	
	/**
	 * Returns the integer result of base^power
	 * 
	 * Tested - works.
	 * 
	 * @param base
	 * @param power
	 * @return
	 */
	public static long power(long base, long power) {
		long result = 1;
		long absPower = Math.abs(power);
		for (long p = 0; p < absPower; p++) {
			result *= base;
		}
		if (power < 0) {
			// This will be zero for any base except 1 or -1
			result = 1 / result;
		}
		return result;
	}
	
	/**
	 * Compute n! (n factorial), without checking for overflow
	 * 
	 * @param n
	 * @return
	 */
	public static long factorial(int n) {
		long result = 1;
		for (int i = 1; i <= n; i++) {
			result *= (long) i;
		}
		return result;
	}
	
	/**
	 * Compute n! (n factorial), but check if this causes integer overflow
	 * 
	 * @param n
	 * @return the integer (not long value of) n factorial
	 * @throws Exception when the result exceeds {@link Integer#MAX_VALUE}
	 */
	public static int factorialCheckBounds(int n) throws Exception {
		// TODO I'm not convinced this is safe. Should just hard-code the limit.
		long result = 1;
		for (int i = 1; i <= n; i++) {
			result *= (long) i;
			if (result > Integer.MAX_VALUE) {
				throw new Exception("n! causes integer overflow");
			}
		}
		return (int) result;
	}

	public static double factorialAsDouble(int n) {
		double result = 1;
		for (int i = 1; i <= n; i++) {
			result *= (double) i;
		}
		return result;		
	}
	
	/**
	 * Computes n! as a double (to provide extended range over a long).
	 * 
	 * We include the divisor here since n! hits with n at about 340.
	 * So if there is anything the result would have been divided by, we include it here,
	 *  thus extending the range of n the function is suitable for. 
	 * 
	 * @param n
	 * @param divisor
	 * @return
	 */
	public static double factorialAsDoubleIncludeDivisor(int n, double divisor) {
		double result = 1.0 / divisor;
		for (int i = 1; i <= n; i++) {
			result *= (double) i;
		}
		return result;		
	}

	/**
	 * n!!
	 * see http://en.wikipedia.org/wiki/Double_factorial#Double_factorial
	 * 
	 * @param n
	 * @return
	 */
	public static long doubleFactorial(int n) {
		long result = 1;
		int startValue;
		if (n % 2 == 0) {
			// n even
			startValue = 2;
		} else {
			// n odd
			startValue = 3;
		}
		for (int i = startValue; i <= n; i += 2) {
			result *= (long) i;
		}
		return result;
	}

	/**
	 * n!!
	 * see http://en.wikipedia.org/wiki/Double_factorial#Double_factorial
	 * 
	 * @param n
	 * @return
	 */
	public static double doubleFactorialAsDouble(int n) {
		double result = 1.0;
		int startValue;
		if (n % 2 == 0) {
			// n even
			startValue = 2;
		} else {
			// n odd
			startValue = 3;
		}
		for (int i = startValue; i <= n; i += 2) {
			result *= (double) i;
		}
		return result;
	}
	
	/**
	 * n!!
	 * see http://en.wikipedia.org/wiki/Double_factorial#Double_factorial
	 * 
	 * We include the divisor here since the gamma(d/2+1) hits Inf with d just over 300.
	 * So if there is anything the result would have been divided by, we include it here,
	 *  thus extending the range of d the function is suitable for. 
	 * 
	 * @param n
	 * @param divisor
	 * @return
	 */
	public static double doubleFactorialAsDoublewithDivisor(int n, double divisor) {
		double result = 1.0 / divisor;
		int startValue;
		if (n % 2 == 0) {
			// n even
			startValue = 2;
		} else {
			// n odd
			startValue = 3;
		}
		for (int i = startValue; i <= n; i += 2) {
			result *= (double) i;
		}
		return result;
	}

	/**
	 * Computes gamma(d/2 + 1)
	 * See http://en.wikipedia.org/wiki/Gamma_function
	 *  for description of the analytical result for d odd.
	 * For d even, we have gamma of an integer, which is equal to
	 *  (d/2)!
	 * 
	 * @param d
	 * @return
	 */
	public static double gammaOfArgOn2Plus1(int d) {
		if (d % 2 == 0) {
			// d even
			return factorialAsDouble(d/2);
		} else {
			// d odd
			return Math.sqrt(Math.PI) * (double) doubleFactorialAsDouble(d) / 
				(double) Math.pow(2, ((double) (d + 1)) / 2.0);
		}
	}

	/**
	 * Computes gamma(d/2 + 1)/divisor
	 * See http://en.wikipedia.org/wiki/Gamma_function
	 *  for description of the analytical result for d odd.
	 * For d even, we have gamma of an integer, which is equal to
	 *  (d/2)!
	 * 
	 * We include the divisor here since the gamma(d/2+1) hits Inf with d just over 300.
	 * So if there is anything the result would have been divided by, we include it here,
	 *  thus extending the range of d the function is suitable for. 
	 * 
	 * @param d
	 * @param divisor
	 * @return
	 */
	public static double gammaOfArgOn2Plus1IncludeDivisor(int d, double divisor) {
		if (d % 2 == 0) {
			// d even
			return factorialAsDoubleIncludeDivisor(d/2, divisor);
		} else {
			// d odd
			return doubleFactorialAsDoublewithDivisor(d,
					divisor * Math.pow(2, ((double) (d + 1)) / 2.0) / Math.sqrt(Math.PI));
		}
	}

	/**
	 * <p>Compute digamma(d).</p>
	 * 
	 * <p>Stores previous calculations to speed up computation here, though some precision may
	 *  be lost because we're adding in larger numbers first.</p>
	 * 
	 * @param d
	 * @return
	 * @throws Exception
	 */
	public static double digamma(int d) throws Exception {
		if (d < 1) {
			return Double.NaN;
		}
		// Need to take care with multi-threaded applications
		//  here -- another thread could be attempting to update
		//  highestDigammaArgCalced at the same time.
		// Sample the value of highestDigammaArgCalced and
		//  run with it through the rest of the method, to avoid
		//  any race conditions if it is concurrently updated:
		int highestDigammaArgCalcedAtStartForThisThread =
				highestDigammaArgCalced;
		if (d <= highestDigammaArgCalcedAtStartForThisThread) {
			// We've already calculated this one
			return storedDigammas[d];
		}
		if (d >= NUM_STORED_DIGAMMAS) {
			// Don't bother updating our storage,
			//  directly use commons.math:
			return Gamma.digamma(d);
		}
		if (highestDigammaArgCalcedAtStartForThisThread == 0) {
			// We're not initialised yet:
			storedDigammas[0] = Double.NaN;
			storedDigammas[1] = -EULER_MASCHERONI_CONSTANT;
			highestDigammaArgCalcedAtStartForThisThread = 1;
		}
		// Else we'll calculate it and update the storage:
		double result = storedDigammas[highestDigammaArgCalcedAtStartForThisThread];
		for (int n = highestDigammaArgCalcedAtStartForThisThread + 1; n <= d; n++) {
			result += 1.0 / (double) (n-1);
			// n must be < NUM_STORED_DIGAMMAS by earlier if statement on d
			storedDigammas[n] = result;
		}
		// Finally, update highestDigammaArgCalced:
		// We could synchronize on highestDigammaArgCalced, however
		//  this is costly; the following check may be compromised under
		//  a race condition, but the worst outcome is that we have to 
		//  recalculate a few values -- this is probably faster than
		//  enforcing synchronisation
		if (d > highestDigammaArgCalced) {
			// We make the check in case another thread has already set this higher
			highestDigammaArgCalced = d;
		}
		return result;
	}
	
	/**
	 * Compute the digamma function from first principles
	 * 
	 * @param d
	 * @return
	 * @throws Exception
	 */
	public static double digammaByDefinition(int d) throws Exception {
		if (d < 1) {
			return Double.NaN;
		}
		double result = 0;
		for (int n = d; n > 1; n--) {
			result += 1.0 / (double) (n-1);
		}
		// Now add in result for n == 1
		result += -EULER_MASCHERONI_CONSTANT;
		return result;
	}

	/**
	 * <p>Return the value of the cummulative distribution function of the
	 *  chi-square distribution, evaluated at x, for k degrees of freedom.</p>
	 *  
	 * <p>Note that this relies on our approximation of the error function,
	 *  which is the limiting part of the accuracy. Testing against
	 *  values produced by octave indicates this is accurate to 5-6
	 *  decimal places.</p>
	 * 
	 * @param x value at which to evaluate the CDF
	 * @param k degrees of freedom (must have k>0)
	 * @return chi squared CDF evaluated at x given k degrees of freedom
	 * @see {@link http://en.wikipedia.org/wiki/Chi-squared_distribution}
	 */
	public static double chiSquareCdf(double x, int k) {
		if (k <= 0) {
			throw new IllegalArgumentException("k (" + k + ") must be > 0");
		}
		// Old approach: (not numerically stable):
		// return lowerIncompleteGammaFunctionOfArgsOn2(k,x) /
		//		gammaOfArgOn2Plus1(k-2); // denominator is Gamma(k/2)
		// New approach:
		return Gamma.regularizedGammaP(((double)k)/2.0, x/2.0);
	}
	
	/**
	 * Return the value of the lower Incomplete Gamma function,
	 * given arguments s/2 and x/2.
	 * We assume postive integer parameter s (s could be complex in general,
	 * with positive real part, but we restrict it to real and integer for
	 * this method). We make the evaluation using a recurrence relation,
	 * which terminates at s/2 = 1 or 1/2 (i.e. s = 2 or 1)
	 * 
	 * @param s for parameter s/2 to lower incomplete gamma
	 * @param x for value x/2 to lower incomplete gamma
	 * @return value of lower gamma incomplete function
	 * @see {@link http://en.wikipedia.org/wiki/Incomplete_Gamma_function}
	 */
	public static double lowerIncompleteGammaFunctionOfArgsOn2(int s, double x) {
		if (s <= 0) {
			throw new IllegalArgumentException("s must be > 0");
		}
		if (s == 2) {
			// Terminating condition: evaluate lower gamma(1, x/2):
			return 1 - Math.exp(-x/2.0);
		} else if (s == 1) {
			// Terminating condition: evaluate lower gamma(1/2, x/2):
			return Math.sqrt(Math.PI) * erf(Math.sqrt(x/2.0));
		} else {
			// Else evaluate recurrence relation:
			return (s/2.0-1.0)*lowerIncompleteGammaFunctionOfArgsOn2(s-2,x) -
					Math.pow(x/2.0, s/2.0 - 1.0) * Math.exp(-x/2.0);
		}
	}
	
	/**
	 * Return the value of the error function at a given x.
	 * We approximate the error function using elementary functions
	 *  as described at the link below (quoting Abramowitz and Stegun).
	 *  This approximation is quoted to
	 *  have maximum error 1.5e-7 (and indeed this appears to be the
	 *  case in comparison to values produced by octave).
	 * 
	 * @param x value at which to evaluate the error function
	 * @return erf(x)
	 * @see {@link http://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions}
	 * @see Abramowitz, Milton; Stegun, Irene A., eds. (1972),
	 * "Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables",
	 * New York: Dover Publications, ISBN 978-0-486-61272-0
	 */
	public static double erf(double x) {
		// Constants:
		double p = 0.3275911;
		double[] a = {0.254829592, -0.284496736, 1.421413741,
				-1.453152027, 1.061405429};
		boolean negArg = (x < 0);

		if (negArg) {
			// The rest of the method requires x >= 0, but since erf(x)
			//  is an odd function, we just reflect x.
			x = -x;
		}
		
		double t = 1.0 / (1 + p * x);
		double multiplier = 0.0;
		double tToPower = t;
		for (int i = 0; i < 5; i++) {
			multiplier += a[i] * tToPower;
			tToPower *= t;
		}
		double retVal = 1.0 - multiplier * Math.exp(-x*x);
		// Remember that erf(x) was an odd function:
		return negArg ? - retVal: retVal;
	}
	
	/**
	 * Compute the probability density function (PDF) of an observation
	 *  x, given the mean mu, and the standard deviation sigma,
	 *  given that the observations follow a univariate normal distribution.
	 * 
	 * @param x observation
	 * @param mu mean of distribution
	 * @param sigma standard deviation (not the variance)
	 * @return PDF value
	 * @throws Exception when sigma < 0
	 */
	public static double normalPdf(double x, double mu, double sigma)
			throws Exception {
		if (sigma < 0) {
			throw new Exception("Standard deviation cannot be < 0");
		}
		double expArg = (x - mu)/sigma;
		expArg *= expArg;
		double pdf =
				Math.pow(2.0*Math.PI, -0.5) /
					sigma *
						Math.exp(-0.5 * expArg);
		return pdf;
	}
	
	/**
	 * Compute the cumulative density function (CDF) of an observation
	 *  x, given the mean mu, and the standard deviation sigma,
	 *  given that the observations follow a univariate normal distribution.
	 * 
	 * @param x observation
	 * @param mu mean of distribution
	 * @param sigma standard deviation (not the variance)
	 * @return CDF value
	 * @throws Exception when sigma < 0
	 */
	public static double normalCdf(double x, double mu, double sigma)
			throws Exception {
		if (sigma < 0) {
			throw new Exception("Standard deviation cannot be < 0");
		}
		double erfArg = (x - mu)/Math.sqrt(2.0*sigma*sigma);
		double cdf = 0.5 * (1 +
				erf(erfArg));
		return cdf;
	}
	

	/**
	 * Compute the probability density function (PDF) of a vector of observations
	 *  x, given the means, and the
	 *  covariance of the variables, given that the observations
	 *  follow a multivariate normal distribution
	 * 
	 * @param x observations
	 * @param means means of each variable
	 * @param covariance covariance matrix
	 * @return PDF value
	 * @throws Exception when the lengths of x and covariance
	 *  do not match, or if supplied a non-square matrix, or if the covariance
	 *  is not symmetric or not positive-definite.
	 */
	public static double normalPdf(double[] x, double[] means,
			double[][] covariance) throws Exception {
		if (x.length != means.length) {
			throw new Exception("Length of observations must match means");
		}
		return normalPdf(MatrixUtils.subtract(x, means), covariance);
	}
	
	/**
	 * Compute the probability density function (PDF) of a vector of observations
	 *  x, given the deviationsFromMean of x, i.e. (x - \mu), and the
	 *  covariance of the variables, given that the observations
	 *  follow a multivariate normal distribution
	 * 
	 * @param deviationsFromMean x - \mu
	 * @param covariance
	 * @return PDF value
	 * @throws Exception when the lengths of deviationsFromMean and covariance
	 *  do not match, or if supplied a non-square matrix, or if the covariance
	 *  is not symmetric or not positive-definite.
	 */
	public static double normalPdf(double[] deviationsFromMean,
			double[][] covariance) throws Exception {
		if (deviationsFromMean.length != covariance.length) {
			throw new Exception("Vector length of deviations does not " +
					"match the size of the covariance matrix");
		}
		double det = MatrixUtils.determinant(covariance);
		double[][] invCovariance = MatrixUtils.invertSymmPosDefMatrix(covariance);
		double expArg = MatrixUtils.dotProduct(
				MatrixUtils.matrixProduct(deviationsFromMean,
						invCovariance),
				deviationsFromMean);
		double pdf =
				Math.pow(2.0*Math.PI, -deviationsFromMean.length / 2.0) /
					Math.sqrt(det) *
						Math.exp(-0.5 * expArg);
		return pdf;
	}
	
	/**
	 * Return the number of possible combinations of p from n (i.e. n choose p)
	 * 
	 * @param n
	 * @param p
	 * @return
	 * @throws Exception if the number would be greater than Integer.MAX_INT
	 */
	public static int numOfSets(int n, int p) throws Exception {
		// Compute how many sets there will be
		long counter = n;
		long numSets = 1;
		for (int x = 1; x <= p; x++) {
			numSets *= counter;
			numSets /= x;
			if (numSets > Integer.MAX_VALUE) {
				throw new Exception("nCp causes integer overflow");
			}
			counter--;
		}
		// numSets counts the number of permutations of n.
		// Need to get rid of repeats to make is combinations:
		return (int) numSets;
	}

	/**
	 * Return an array of all possible combinations of p from n 
	 * 
	 * @param n
	 * @param p
	 * @return
	 * @throws Exception when the number of sets is greaterr than Integer.MAX_INT
	 */
	public static int[][] generateAllSets(int n, int p) throws Exception {
		int numOfSets = numOfSets(n,p);
		int[][] sets = new int[numOfSets][p];
		int[] currentSet = new int[p];
		writeSetsIn(n, p, 0, 0, currentSet, sets, 0);
		return sets;
	}

	/**
	 * Recursive call used by generateAllSets.
	 * 
	 * @param n
	 * @param p
	 * @param currentIndexInSet current index in currentSet that we are writing into
	 * @param currentSet current set containing indices already written into the upper parts
	 * @param sets array to write generated sets into
	 * @param upToSetNum
	 * @return new value of upToSetNum
	 */
	private static int writeSetsIn(int n, int p, int currentIndexInSet,
			int firstCandidate, int[] currentSet, int[][] sets, int upToSetNum) {
		/*
		String indent = "";
		for (int i = 0; i < currentIndexInSet; i++) {
			indent += " ";
		}
		System.out.println(indent + String.format("currentIndex=%d", currentIndexInSet));
		*/
		// Put every candidate into this position:
		for (int candidate = firstCandidate; candidate < n - p + currentIndexInSet + 1; candidate++) {
			// System.out.println(indent + candidate);
			currentSet[currentIndexInSet] = candidate;
			if (currentIndexInSet == p - 1) {
				// We just wrote the last index, so copy this one in and return
				// System.out.println(indent + "writing into line " + upToSetNum);
				System.arraycopy(currentSet, 0, sets[upToSetNum++], 0, p);
			} else {
				// There are more indices to be written in, so make a recursive call to write the 
				//  next ones in
				upToSetNum = writeSetsIn(n, p, currentIndexInSet + 1, candidate + 1, currentSet, sets, upToSetNum);
			}
		}
		return upToSetNum;
	}

	/**
	 * Perform some testing:
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String args[]) throws Exception {
		/*
		System.out.println(numOfSets(158,4));
		System.out.println(numOfSets(158,3));
		System.out.println(numOfSets(158,2));
		*/
		// int[][] sets = generateAllSets(6,4);
		// MatrixUtils.printMatrix(System.out, sets);
		
		/*
		System.out.printf("digamma()  digammaOld()\n");
		for (int n = 0; n < 100; n++) {
			System.out.printf("%d  %.3f  %.3f\n", n, MathsUtils.digamma(n), MathsUtils.digammaByDefinition(n));
		}
		for (int n = 0; n < 101; n++) {
			System.out.printf("%d  %.3f  %.3f\n", n, MathsUtils.digamma(n), MathsUtils.digammaByDefinition(n));
		}
		*/
		
		/*
		System.out.println("erf(" + 1 + ")= " + MathsUtils.erf(1));
		System.out.println("erf(" + 2 + ")= " + MathsUtils.erf(2));
		System.out.println("erf(" + 0 + ")= " + MathsUtils.erf(0));
		for (int n=0; n<100; n++) {
			System.out.println("erf(" + n*0.1 + ")= " + MathsUtils.erf(n*0.1));
		}*/
		
		int degFree = 10;
		System.out.println("chi2cdf(1," + degFree +")= " + MathsUtils.chiSquareCdf(1, degFree));
		System.out.println("chi2cdf(2," + degFree +")= " + MathsUtils.chiSquareCdf(2, degFree));
		System.out.println("chi2cdf(3," + degFree +")= " + MathsUtils.chiSquareCdf(3, degFree));
		for (int n=0; n<100; n++) {
			System.out.println("chi2cdf(" + n*0.1 + "," + degFree +")= " + MathsUtils.chiSquareCdf(n*0.1, degFree));
		}
		
	}
}
