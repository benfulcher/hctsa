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

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Vector;

/**
 * Utilities for computations on arrays and matrices of data.
 * Matrices are represented as either one-dimensional
 *  arrays of int[] or double[], or two-dimensional
 *  arrays of double[][] or int[][] - it is assumed that all
 *  multidimensional matrices have consistent lengths in each dimension
 *  matrix[i].
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class MatrixUtils {

	/**
	 * Generate the identity matrix of the given size
	 * 
	 * @param size (size along one dimension)
	 * @return two dimensional double array representing the identity matrix
	 */
	public static double[][] identityMatrix(int size) {
		double[][] I = new double[size][size];
		for (int r = 0; r < size; r++) {
			I[r][r] = 1.0;
		}
		return I;
	}
	
	/**
	 * Return an array with values enumerated through the given range
	 * 
	 * @param startValue first value for the array
	 * @param endValue last value for the array
	 * @return
	 */
	public static int[] range(int startValue, int endValue) {
		int[] array = new int[endValue - startValue + 1];
		for (int i = 0; i < endValue - startValue + 1; i++) {
			array[i] = startValue + i;
		}
		return array;
	}
	
	/**
	 * Return an array with the given value at every index
	 * 
	 * @param length length of array
	 * @param value value for every element of the arry
	 * @return
	 */
	public static double[] constantArray(int length, double value) {
		double[] array = new double[length];
		Arrays.fill(array, value);
		return array;
	}
	
	public static double sum(double[] input) {
		double total = 0;
		for (int i = 0; i < input.length; i++) {
			total += input[i];
		}
		return total;
	}
	
	public static double sum(double[] input, int startIndex, int length) {
		double total = 0;
		for (int i = startIndex; i < startIndex + length; i++) {
			total += input[i];
		}
		return total;
	}
	
	public static double sumSpecificIndices(double[] input, int[] indices) {
		double total = 0;
		for (int i = 0; i < indices.length; i++) {
			total += input[indices[i]];
		}
		return total;
	}

	public static double sumSpecificIndices(double[] input, int[][] indices, int columnInIndices) {
		double total = 0;
		for (int i = 0; i < indices.length; i++) {
			total += input[indices[i][columnInIndices]];
		}
		return total;
	}

	public static double sumSpecificIndices(double[] input, int[][] indices, int columnInIndices,
			int indicesOffset) {
		double total = 0;
		for (int i = 0; i < indices.length; i++) {
			total += input[indices[i][columnInIndices] + indicesOffset];
		}
		return total;
	}

	public static double sum(double[][] input) {
		double total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[i].length; j++) {
				total += input[i][j];
			}
		}
		return total;
	}

	public static double sum(double[][] input, int column) {
		double total = 0;
		for (int i = 0; i < input.length; i++) {
			total += input[i][column];
		}
		return total;
	}

	public static int sum(int[] input) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			total += input[i];
		}
		return total;
	}

	public static int sum(int[][] input) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[i].length; j++) {
				total += input[i][j];
			}
		}
		return total;
	}

	/**
	 * Return an array of the sums for each column in the 2D input
	 * 
	 * @param input
	 * @return
	 */
	public static double[] sums(double[][] input) {
		double[] theSums = new double[input[0].length];
		for (int r = 0; r < input.length; r++) {
			for (int c = 0; c < input[r].length; c++) {
				theSums[c] += input[r][c];
			}
		}
		return theSums;
	}
	
	/**
	 * Return an array of the sums for each column in the 2D input
	 * 
	 * @param input
	 * @param startRow which row to start from
	 * @param length how many rows to take the sum over
	 * @return
	 */
	public static double[] sums(double[][] input, int startRow, int length) {
		double[] theSums = new double[input[0].length];
		for (int r = startRow; r < startRow + length; r++) {
			for (int c = 0; c < input[r].length; c++) {
				theSums[c] += input[r][c];
			}
		}
		return theSums;
	}
	
	public static int countIf(int[] input, int condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			if (input[i] == condition)
				total++;
		}
		return total;
	}

	public static int countIf(int[][] input, int condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[0].length; j++) {
				if (input[i][j] == condition)
					total++;				
			}
		}
		return total;
	}

	public static int countIf(long[][] input, long condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[0].length; j++) {
				if (input[i][j] == condition)
					total++;				
			}
		}
		return total;
	}

	public static int countIf(int[] input1, int condition1, int[] input2, int condition2) 
		throws Exception {
		
		if (input1.length != input2.length)
			throw new Exception("MatrixUtils.sumIf() - arguments are not of equal length (" +
					input1.length + " != " + input2.length + ")");
		int total = 0;
		for (int i = 0; i < input1.length; i++) {
			if ((input1[i] == condition1) && (input2[i] == condition2))
				total++;
		}
		return total;
	}

	public static int countIf(int[] input1, int condition1, int[] input2, int condition2,
			int[] input3, int condition3) throws Exception {
	
		if ((input1.length != input2.length) || (input1.length != input3.length))
			throw new Exception("MatrixUtils.sumIf() - arguments are not of equal length (" +
					input1.length + " != " + input2.length + " != " + input3.length + ")");
		int total = 0;
		for (int i = 0; i < input1.length; i++) {
			if ((input1[i] == condition1) && (input2[i] == condition2) && (input3[i] == condition3))
				total++;
		}
		return total;
	}
	
	public static int countIf(boolean[] input, boolean condition) {
		int total = 0;
		for (int i = 0; i < input.length; i++) {
			if (input[i] == condition)
				total++;
		}
		return total;
	}

	public static double mean(int[] input) {
		return sum(input) / (double) input.length;
	}
	
	public static double mean(double[] input) {
		return sum(input) / (double) input.length;
	}
	
	public static double mean(double[] input, int startIndex, int length) {
		return sum(input, startIndex, length) / (double) length;
	}

	public static double mean(double[][] input) {
		return sum(input) / (double) (input.length * input[0].length);
	}

	/**
	 * Compute the mean along the given column 
	 * 
	 * @param input
	 * @param column
	 * @return
	 */
	public static double mean(double[][] input, int column) {
		return sum(input, column) / (double) input.length;
	}

	/**
	 * Return an array of the means of each column in the 2D input
	 * 
	 * @param input
	 * @return
	 */
	public static double[] means(double[][] input) {
		double[] theMeans = sums(input);
		for (int i = 0; i < theMeans.length; i++) {
			theMeans[i] = theMeans[i] / input.length;
		}
		return theMeans;
	}

	/**
	 * Return an array of the means of each column in the 2D input
	 * 
	 * @param input
	 * @param startRow which row to start from
	 * @param length how many rows to take the mean over
	 * @return
	 */
	public static double[] means(double[][] input, int startRow, int length) {
		double[] theMeans = sums(input, startRow, length);
		for (int i = 0; i < theMeans.length; i++) {
			theMeans[i] = theMeans[i] / length;
		}
		return theMeans;
	}

	/**
	 * Return an array of the means of each row in the 2D input matrix
	 * 
	 * @param input
	 * @return
	 */
	public static double[] meansOfRows(double[][] input) {
		double[] theMeans = new double[input.length];
		for (int i = 0; i < input.length; i++) {
			theMeans[i] = mean(input[i]);
		}
		return theMeans;
	}

	public static int[][] columnShift(int[][] input, int shiftBy){
		if (shiftBy == 0) {
			return input;
		}
		
		int rows = input.length;
		int columns = input[0].length;
		for ( ; shiftBy < 0; shiftBy += columns) {
			// Using % mod operator to come back to a +ve column value wont work.
			// So we're shifting the shiftBy value (above) until it's in the appropriate
			// range 0 .. columns-1
		}
		int[][] output = new int[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				output[r][(c + shiftBy) % columns] = input[r][c];
			}
		}
		return output;
	}
	
	public static double[][] columnShift(double[][] input, int shiftBy){
		if (shiftBy == 0) {
			return input;
		}
		
		int rows = input.length;
		int columns = input[0].length;
		for ( ; shiftBy < 0; shiftBy += columns) {
			// Using % mod operator to come back to a +ve column value wont work.
			// So we're shifting the shiftBy value (above) until it's in the appropriate
			// range 0 .. columns-1
		}
		double[][] output = new double[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				output[r][(c + shiftBy) % columns] = input[r][c];
			}
		}
		return output;
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows or columns. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @param fromColumn
	 * @param colums
	 * @return Single dimensional array containing the required data
	 */
	public static int[] matrixToArray(int[][] input, int fromRow, int rows, int fromColumn, int columns) {
		int[] output = new int[rows * columns];
		for (int c = 0; c < columns; c++) {
			for (int r = 0; r < rows; r++) {
				output[c * rows + r] = input[r + fromRow][c + fromColumn];
			}
		}
		return output;
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows only. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @return Single dimensional array containing the required data
	 */
	public static int[] matrixToArray(int[][] input, int fromRow, int rows) {
		return matrixToArray(input, fromRow, rows, 0, input[0].length);
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * 
	 * @param input
	 * @return Single dimensional array containing the required data
	 */
	public static int[] matrixToArray(int[][] input) {
		return matrixToArray(input, 0, input.length, 0, input[0].length);
	}
	
	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows or columns. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @param fromColumn
	 * @param colums
	 * @return Single dimensional array containing the required data
	 */
	public static double[] matrixToArray(double[][] input, int fromRow, int rows, int fromColumn, int columns) {
		double[] output = new double[rows * columns];
		for (int c = 0; c < columns; c++) {
			for (int r = 0; r < rows; r++) {
				output[c * rows + r] = input[r + fromRow][c + fromColumn];
			}
		}
		return output;
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * Provides controllers for selecting a subset of rows only. 
	 * 
	 * @param input
	 * @param fromRow
	 * @param rows
	 * @return Single dimensional array containing the required data
	 */
	public static double[] matrixToArray(double[][] input, int fromRow, int rows) {
		return matrixToArray(input, fromRow, rows, 0, input[0].length);
	}

	/**
	 * Converts a 2 dimensional array into a single dimension.
	 * Places each column on top of each other.
	 * 
	 * @param input
	 * @return Single dimensional array containing the required data
	 */
	public static double[] matrixToArray(double[][] input) {
		return matrixToArray(input, 0, input.length, 0, input[0].length);
	}

	/**
	 * Adds two arrays together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 */
	public static int[] add(int[] input1, int[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		int[] returnValues = new int[input1.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = input1[i] + input2[i];
		}
		return returnValues;
	}

	/**
	 * Adds two arrays together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 */
	public static double[] add(double[] input1, double[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		double[] returnValues = new double[input1.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = input1[i] + input2[i];
		}
		return returnValues;
	}

	/**
	 * Adds two arrays together, returning the result in input1
	 * 
	 * @param input1
	 * @param input2
	 */
	public static void addInPlace(double[] input1, double[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		for (int i = 0; i < input1.length; i++) {
			input1[i] = input1[i] + input2[i];
		}
	}

	/**
	 * Adds the squares of the second array to the first,
	 *  returning the result in input1
	 * 
	 * @param input1
	 * @param input2
	 */
	public static void addSquaresInPlace(double[] input1, double[] input2) throws Exception {
		if (input1.length != input2.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		for (int i = 0; i < input1.length; i++) {
			input1[i] = input1[i] + input2[i] * input2[i];
		}
	}

	/**
	 * Adds two matrices together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 * @throws Exception
	 */
	public static int[][] add(int[][] input1, int[][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		int[][] returnValues = new int[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				returnValues[r][c] = input1[r][c] + input2[r][c];
			}
		}
		return returnValues;
	}
	
	/**
	 * Adds two matrices together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 * @throws Exception
	 */
	public static double[][] add(double[][] input1, double[][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		double[][] returnValues = new double[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				returnValues[r][c] = input1[r][c] + input2[r][c];
			}
		}
		return returnValues;
	}

	/**
	 * Adds two matrices together
	 * 
	 * @param input1
	 * @param input2
	 * @return
	 * @throws Exception
	 */
	public static double[][][] add(double[][][] input1, double[][][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		int height = input1[0][0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		if (input2[0][0].length != height) {
			throw new Exception("Heights (3rd dim) of arrays are not equal");
		}
		double[][][] returnValues = new double[rows][columns][height];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				for (int h = 0; h < height; h++) {
					returnValues[r][c][h] = input1[r][c][h] + input2[r][c][h];
				}
			}
		}
		return returnValues;
	}

	/**
	 * Subtracts second array from the first
	 * 
	 * @param first
	 * @param second
	 * @return first - second
	 */
	public static double[] subtract(double[] first, double[] second) throws Exception {
		if (first.length != second.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		double[] returnValues = new double[first.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = first[i] - second[i];
		}
		return returnValues;
	}

	/**
	 * Subtracts a constant value from all items in an array
	 * 
	 * @param array
	 * @param value
	 * @return array - constant value
	 */
	public static double[] subtract(double[] array, double value) throws Exception {
		double[] returnValues = new double[array.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = array[i] - value;
		}
		return returnValues;
	}

	/**
	 * Subtracts second array from the first, overwriting the
	 *  values in first
	 * 
	 * @param first
	 * @param second
	 */
	public static void subtractInPlace(double[] first, double[] second) throws Exception {
		if (first.length != second.length) {
			throw new Exception("Lengths of arrays are not equal");
		}
		for (int i = 0; i < first.length; i++) {
			first[i] = first[i] - second[i];
		}
	}

	/**
	 * Subtract one matrix from another
	 * 
	 * @param input1
	 * @param input2
	 * @return input1 - input2
	 * @throws Exception
	 */
	public static int[][] subtract(int[][] input1, int[][] input2) throws Exception {
		int rows = input1.length;
		int columns = input1[0].length;
		if (input2.length != rows) {
			throw new Exception("Row length of arrays are not equal");
		}
		if (input2[0].length != columns) {
			throw new Exception("Column length of arrays are not equal");
		}
		int[][] returnValues = new int[rows][columns];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				returnValues[r][c] = input1[r][c] - input2[r][c];
			}
		}
		return returnValues;
	}

	/**
	 * Subtracts a constant value from all items in an array
	 * 
	 * @param array
	 * @param value
	 * @return array - constant value
	 */
	public static int[] subtract(int[] array, int value) throws Exception {
		int[] returnValues = new int[array.length];
		for (int i = 0; i < returnValues.length; i++) {
			returnValues[i] = array[i] - value;
		}
		return returnValues;
	}

	/**
	 * Return the matrix product A x B
	 * 
	 * @param A mxn matrix
	 * @param B nxq matrix
	 * @return mxq matrix product of A and B
	 */
	public static double[][] matrixProduct(double[][] A, double[][] B) throws Exception {
		if (A[0].length != B.length) {
			throw new Exception("Number of columns of a " + A[0].length +
					" does not match the number of rows of b " + B.length);
		}
		double[][] result = new double[A.length][B[0].length];
		for (int r = 0; r < result.length; r++) {
			for (int c = 0; c < result[r].length; c++) {
				result[r][c] = 0;
				for (int k = 0; k < A[r].length; k++) {
					result[r][c] += A[r][k] * B[k][c];
				}
			}
		}
		return result;
	}
	
	/**
	 * Return the matrix product v A
	 *  (i.e. a left multiplication of the 1xn vector and the nxm matrix A)
	 * 
	 * @param v a 1xn vector
	 * @param A an nxm matrix
	 * @return a 1xm vector output
	 */
	public static double[] matrixProduct(double[] v, double[][] A) throws Exception {
		if (v.length != A.length) {
			throw new Exception("Number of entries of v " + v.length +
					" does not match the number of rows of A " + A.length);
		}
		// Result length is the number of columns of A
		double[] result = new double[A[0].length];
		for (int c = 0; c < result.length; c++) {
			result[c] = 0;
			for (int r = 0; r < v.length; r++) {
				result[c] += v[r]*A[r][c];
			}
		}
		return result;
	}

	/**
	 * Return the matrix product A v
	 *  (i.e. a right multiplication of the nxm matrix A and the 1xn vector)
	 * 
	 * @param A an mxn matrix
	 * @param v a nx1 vector
	 * @return a mx1 vector output
	 */
	public static double[] matrixProduct(double[][] A, double[] v) throws Exception {
		if (v.length != A[0].length) {
			throw new Exception("Number of entries of v " + v.length +
					" does not match the number of columns of A " + A[0].length);
		}
		// Result length is the number of rows of A
		double[] result = new double[A.length];
		for (int r = 0; r < result.length; r++) {
			result[r] = 0;
			for (int c = 0; c < v.length; c++) {
				result[r] += A[r][c] * v[c];
			}
		}
		return result;
	}

	/**
	 * Return the dot product of two vectors v u
	 * 
	 * @param v a nx1 vector
	 * @param u a nx1 vector
	 * @return the scalar dot product
	 */
	public static double dotProduct(double[] v, double[] u) throws Exception {
		if (v.length != u.length) {
			throw new Exception("Number of entries of v " + v.length +
					" does not match the number of entries of u " + u.length);
		}
		double result = 0;
		for (int r = 0; r < v.length; r++) {
			result += v[r] * u[r];
		}
		return result;
	}

	/**
	 * Duplicates a matrix; handles different number of columns
	 *  for each row
	 * 
	 * @param src
	 * @return
	 */
	public static int[][] duplicateMatrix(int[][] src) {
		int[][] dest = new int[src.length][];
		for (int r = 0; r < src.length; r++) {
			dest[r] = new int[src[r].length];
			System.arraycopy(src[r], 0, dest[r], 0, src[r].length);
		}
		return dest;
	}
	
	/**
	 * Copies all rows and columns between two double arrays
	 * 
	 * @param src
	 * @param dest
	 */
	public static void arrayCopy(double[][] src, double[][] dest) {
		for (int r = 0; r < src.length; r++) {
			System.arraycopy(src[r], 0,	dest[r], 0,	src[r].length);
		}
	}

	/**
	 * Copies all rows and columns between two double arrays
	 * 
	 * @param src
	 * @param dest
	 */
	public static double[][] arrayCopy(double[][] src) {
		double[][] dest = new double[src.length][];
		for (int r = 0; r < src.length; r++) {
			dest[r] = new double[src[r].length];
			System.arraycopy(src[r], 0,	dest[r], 0,	src[r].length);
		}
		return dest;
	}

	/**
	 * Copies the required rows and columns between two 
	 * double arrays
	 * 
	 * @param src
	 * @param srcStartRow
	 * @param srcStartCol
	 * @param dest
	 * @param destStartRow
	 * @param destStartCol
	 * @param rows
	 * @param cols
	 */
	public static void arrayCopy(double[][] src, int srcStartRow, int srcStartCol,
								double[][] dest, int destStartRow, int destStartCol,
								int rows, int cols) {
		
		for (int r = 0; r < rows; r++) {
			System.arraycopy(src[srcStartRow + r], srcStartCol,
					dest[destStartRow + r], destStartCol,
					cols);
		}
	}

	/**
	 * Copies all rows and columns between two int arrays
	 * 
	 * @param src
	 * @param dest
	 */
	public static void arrayCopy(int[][] src, int[][] dest) {
		for (int r = 0; r < src.length; r++) {
			System.arraycopy(src[r], 0,	dest[r], 0,	src[r].length);
		}
	}

	/**
	 * Copies the required rows and columns between two 
	 * double arrays
	 * 
	 * @param src
	 * @param srcStartRow
	 * @param srcStartCol
	 * @param dest
	 * @param destStartRow
	 * @param destStartCol
	 * @param rows
	 * @param cols
	 */
	public static void arrayCopy(int[][] src, int srcStartRow, int srcStartCol,
								int[][] dest, int destStartRow, int destStartCol,
								int rows, int cols) {
		
		for (int r = 0; r < rows; r++) {
			System.arraycopy(src[srcStartRow + r], srcStartCol,
					dest[destStartRow + r], destStartCol,
					cols);
		}
	}
	
	/**
	 * Copies the given source array into the required column number of the destination
	 * @param destination
	 * @param column
	 * @param source
	 */
	public static void copyIntoColumn(int[][] destination, int column, int[] source) throws Exception {
		if (source.length != destination.length) {
			throw new Exception("Destination column is not of the same length as the source (" +
					destination.length + " vs " + source.length + ")");
		}
		for (int r = 0; r < destination.length; r++) {
			destination[r][column] = source[r];
		}
	}

	/**
	 * Copies the given source array into the required column number of the destination
	 * @param destination
	 * @param column
	 * @param source
	 */
	public static void copyIntoColumn(double[][] destination, int column, 
			int destFromRowNumber, double[] source, int sourceFromRowNumber,
			int rows) throws Exception {
		if (sourceFromRowNumber + rows > source.length) {
			throw new Exception("Attempting to copy too many rows " + rows +
					" after the start row " + sourceFromRowNumber +
					" from the source of length " + source.length);
		}
		if (destFromRowNumber + rows > destination.length) {
			throw new Exception("Attempting to copy too many rows " + rows +
					" after the start row " + destFromRowNumber +
					" from the destination of length " + destination.length);
		}
		for (int r = 0; r < rows; r++) {
			destination[r + destFromRowNumber][column] = source[r + sourceFromRowNumber];
		}
	}

	/**
	 * Copies the given source array into the required column number of the destination
	 * @param destination
	 * @param column
	 * @param source
	 */
	public static void copyIntoColumn(int[][] destination, int column, 
			int destFromRowNumber, int[] source, int sourceFromRowNumber,
			int rows) throws Exception {
		if (sourceFromRowNumber + rows > source.length) {
			throw new Exception("Attempting to copy too many rows " + rows +
					" after the start row " + sourceFromRowNumber +
					" from the source of length " + source.length);
		}
		if (destFromRowNumber + rows > destination.length) {
			throw new Exception("Attempting to copy too many rows " + rows +
					" after the start row " + destFromRowNumber +
					" from the destination of length " + destination.length);
		}
		for (int r = 0; r < rows; r++) {
			destination[r + destFromRowNumber][column] = source[r + sourceFromRowNumber];
		}
	}

	/**
	 * Copies the given source array into the required column number of the destination
	 * @param destination
	 * @param column
	 * @param source
	 */
	public static void copyIntoColumn(double[][] destination, int column, double[] source) throws Exception {
		if (source.length != destination.length) {
			throw new Exception("Destination column is not of the same length as the source (" +
					destination.length + " vs " + source.length + ")");
		}
		for (int r = 0; r < destination.length; r++) {
			destination[r][column] = source[r];
		}
	}

	/**
	 * Return a new matrix with the columns of matrix1 joined on the back of matrix2
	 * 
	 * @param matrix1
	 * @param matrix2
	 * @return
	 * @throws Exception 
	 */
	public static double[][] appendColumns(double[][] matrix1, double[][] matrix2) throws Exception {
		double[][] data = new double[matrix1.length][];
		
		if (matrix1.length != matrix2.length) {
			throw new Exception("matrix1 and matrix2 have different lengths");
		}
		if (matrix1.length == 0) {
			return data;
		}
		for (int r = 0; r < matrix1.length; r++) {
			data[r] = append(matrix1[r], matrix2[r]);
		}
		
		return data;
	}
	
	/**
	 * Append the vector u to the vector v and return the result
	 * 
	 * @param v vector 1
	 * @param u vector 2
	 * @return [v, u] appended result
	 */
	public static double[] append(double[] v, double[] u) {
		double[] result = new double[v.length + u.length];
		System.arraycopy(v, 0, result, 0, v.length);
		System.arraycopy(u, 0, result, v.length, u.length);
		return result;
	}
	
	/**
	 * Append the vector u to the vector v and return the result
	 * 
	 * @param v vector 1
	 * @param u vector 2
	 * @return [v, u] appended result
	 */
	public static int[] append(int[] v, int[] u) {
		int[] result = new int[v.length + u.length];
		System.arraycopy(v, 0, result, 0, v.length);
		System.arraycopy(u, 0, result, v.length, u.length);
		return result;
	}

	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply each column
	 *  by a different power of the base.
	 */
	public static int[] computeCombinedValues(int separateValues[][], int base) throws Exception {
		// Number of columns (second index) is sizeof first element
		int columns = separateValues[0].length;
		
		return computeCombinedValues(separateValues, columns, base);
		/*
		// rows = first index specifies rows
		int rows = separateValues.length;
		int[] combinedValues = new int[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			int combinedRowValue = 0;
			int multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += separateValues[r][c] * multiplier;
				multiplier *= base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
		*/
	}
	
	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply the first "columbs" columns
	 *  by a different power of the base.
	 */
	public static int[] computeCombinedValues(int separateValues[][], int columns, int base) throws Exception {
		if (columns > separateValues[0].length) {
			throw new Exception("computeCombinedValues: computation request across more columns " +
					columns + " than are available " + separateValues[0].length);
		}
		// Make sure we won't get any overflow here
		if (combinedValuesOverflow(columns, base)) {
			// multiplier has overflown
			throw new Exception("Too many columns " + columns + " for the given base " + base +
					" for this call to computeCombinedValues");
		}
		
		// rows = first index specifies rows
		int rows = separateValues.length;
		int[] combinedValues = new int[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			int combinedRowValue = 0;
			int multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += separateValues[r][c] * multiplier;
				multiplier *= base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
	}

	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply each column
	 *  by a different power of the base.
	 */
	public static long[] computeCombinedValuesLong(int separateValues[][], int base) throws Exception {
		// Number of columns (second index) is sizeof first element
		int columns = separateValues[0].length;
		
		return computeCombinedValuesLong(separateValues, columns, base);
		/*
		// rows = first index specifies rows
		int rows = separateValues.length;
		long[] combinedValues = new long[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			long combinedRowValue = 0;
			long multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += ((long) separateValues[r][c]) * multiplier;
				multiplier *= (long) base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
		*/
	}

	/**
	 * 
	 * @param separateValues
	 * @return Single dimensional matrix where each row
	 *  has been combined into a single output value, unique
	 *  to the input row. We basically multiply the first "columns" columns
	 *  by a different power of the base.
	 */
	public static long[] computeCombinedValuesLong(int separateValues[][], int columns, int base) throws Exception {
		if (columns > separateValues[0].length) {
			throw new Exception("computeCombinedValuesLong: computation request across more columns " +
					columns + " than are available " + separateValues[0].length);
		}
		
		// Make sure we won't get any overflow here
		if (combinedValuesOverflowLong(columns, base)) {
			// multiplier has overflown
			throw new Exception("Too many columns " + columns + " for the given base " + base +
					" for this call to computeCombinedValuesLong");
		}
		
		// rows = first index specifies rows
		int rows = separateValues.length;
		long[] combinedValues = new long[rows];
		for (int r = 0; r < rows; r++) {
			// For each row in vec1
			long combinedRowValue = 0;
			long multiplier = 1;
			for (int c = columns - 1; c >= 0; c--) {
				// Add in the contribution from each column
				combinedRowValue += ((long) separateValues[r][c]) * multiplier;
				multiplier *= (long) base;
			}
			combinedValues[r] = combinedRowValue;
		}
		return combinedValues;
	}

	public static boolean combinedValuesOverflow(int columns, int base) {
		// Make sure we won't get any overflow here
		int multiplier = 1;
		for (int c = columns - 1; c >= 0; c--) {
			if (multiplier < 0) {
				// multiplier has overflown.
				// Technically, it's possible to use one negative value if we were using base-2,
				//  but realistically it's safer if we just call it off now.
				return true;
			}
			multiplier *= (long) base;
		}
		return false;
	}

	public static boolean combinedValuesOverflowLong(int columns, int base) {
		// Make sure we won't get any overflow here
		long multiplier = 1;
		for (int c = columns - 1; c >= 0; c--) {
			if (multiplier < 0) {
				// multiplier has overflown.
				// Technically, it's possible to use one negative value if we were using base-2,
				//  but realistically it's safer if we just call it off now.
				return true;
			}
			multiplier *= (long) base;
		}
		return false;
	}
	
	/**
	 * Select out part of an array.
	 * 
	 * @param data
	 * @param fromIndex
	 * @param length
	 * @return
	 */
	public static double[] select(double[] data, int fromIndex, int length) {
		double[] returnData = new double[length];
		System.arraycopy(data, fromIndex, returnData, 0, length);
		return returnData;
	}
	
	/**
	 * Select out part of an array.
	 * 
	 * @param data
	 * @param indices which array indices to pull out
	 * @return
	 */
	public static double[] select(double[] data, int[] indices) {
		double[] returnData = new double[indices.length];
		for (int i = 0; i < indices.length; i++) {
			returnData[i] = data[indices[i]];
		}
		return returnData;
	}
	
	/**
	 * Select out part of an array.
	 * 
	 * @param data
	 * @param fromIndex
	 * @param length
	 * @return
	 */
	public static int[] select(int[] data, int fromIndex, int length) {
		int[] returnData = new int[length];
		System.arraycopy(data, fromIndex, returnData, 0, length);
		return returnData;
	}
	
	/**
	 * Select out part of an array.
	 * 
	 * @param data
	 * @param indices which array indices to pull out
	 * @return
	 */
	public static int[] select(int[] data, int[] indices) {
		int[] returnData = new int[indices.length];
		for (int i = 0; i < indices.length; i++) {
			returnData[i] = data[indices[i]];
		}
		return returnData;
	}

	public static int[] selectColumn(int matrix[][], int columnNo) {
		int[] column = new int[matrix.length];
		for (int r = 0; r < matrix.length; r++) {
			column[r] = matrix[r][columnNo];
		}
		return column;
	}
	
	public static double[] selectColumn(double matrix[][], int columnNo) {
		double[] column = new double[matrix.length];
		for (int r = 0; r < matrix.length; r++) {
			column[r] = matrix[r][columnNo];
		}
		return column;
	}

	public static double[] selectColumn(double matrix[][], int columnNo,
			int startRow, int rows) {
		double[] column = new double[rows];
		for (int r = 0; r < rows; r++) {
			column[r] = matrix[startRow + r][columnNo];
		}
		return column;
	}

	public static int[] selectColumn(int matrix[][], int columnNo,
			int startRow, int rows) {
		int[] column = new int[rows];
		for (int r = 0; r < rows; r++) {
			column[r] = matrix[startRow + r][columnNo];
		}
		return column;
	}

	public static byte[] selectColumn(byte matrix[][], int columnNo,
			int startRow, int rows) {
		byte[] column = new byte[rows];
		for (int r = 0; r < rows; r++) {
			column[r] = matrix[startRow + r][columnNo];
		}
		return column;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param columns
	 * @return
	 */
	public static double[][] selectColumns(double matrix[][], int columns[]) {
		double[][] data = new double[matrix.length][columns.length];
		for (int r = 0; r < matrix.length; r++) {
			for (int cIndex = 0; cIndex < columns.length; cIndex++) {
				data[r][cIndex] = matrix[r][columns[cIndex]];
			}
		}
		return data;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param includeColumnFlags
	 * @return
	 */
	public static double[][] selectColumns(double matrix[][], boolean includeColumnFlags[]) {
		Vector<Integer> v = new Vector<Integer>();
		
		for (int i = 0; i < includeColumnFlags.length; i++) {
			if (includeColumnFlags[i]) {
				v.add(new Integer(i));
			}
		}
		double[][] data = new double[matrix.length][v.size()];
		for (int r = 0; r < matrix.length; r++) {
			for (int outputColumnIndex = 0; outputColumnIndex < v.size(); outputColumnIndex++) {
				int outputColumn = v.get(outputColumnIndex);
				data[r][outputColumnIndex] = matrix[r][outputColumn];
			}
		}
		return data;
	}

	/**
	 * Extract the required columns from the matrix
	 * 
	 * @param matrix
	 * @param columns
	 * @return
	 */
	public static double[][] selectColumns(double matrix[][], Vector<Integer> columns) {
		double[][] data = new double[matrix.length][columns.size()];
		for (int r = 0; r < matrix.length; r++) {
			for (int cIndex = 0; cIndex < columns.size(); cIndex++) {
				data[r][cIndex] = matrix[r][columns.elementAt(cIndex).intValue()];
			}
		}
		return data;
	}

	/**
	 * Extract the required rows from the matrix
	 * 
	 * @param matrix 2D data array
	 * @param fromRow index of the first row to return
	 * @param rows number of rows (including the first) to return
	 * @return a 2D data array of the selected rows
	 */
	public static double[][] selectRows(double matrix[][], int fromRow, int rows) {
		double[][] data = new double[rows][];
		for (int rIndex = 0; rIndex < rows; rIndex++) {
			data[rIndex] = matrix[rIndex + fromRow];
		}
		return data;
	}

	/**
	 * Extract the required rows and columns from the matrix
	 * 
	 * @param matrix 2D data array
	 * @param rows indices of the rows to select
	 * @param columns indices of the columns to select
	 * @return a 2D data array of the selected rows and columns
	 */
	public static double[][] selectRowsAndColumns(double matrix[][], int rows[], int columns[]) {
		double[][] data = new double[rows.length][columns.length];
		for (int rIndex = 0; rIndex < rows.length; rIndex++) {
			for (int cIndex = 0; cIndex < columns.length; cIndex++) {
				data[rIndex][cIndex] = matrix[rows[rIndex]][columns[cIndex]];
			}
		}
		return data;
	}

	public static double[][] selectFirstTwoDimenions(double[][][][] matrix, int d2, int d3) {
		double[][] newMatrix = new double[matrix.length][];
		for (int i = 0; i < matrix.length; i++) {
			newMatrix[i] = new double[matrix[i].length];
			for (int j = 0; j < matrix[i].length; j++) {
				newMatrix[i][j] = matrix[i][j][d2][d3];
			}
		}
		return newMatrix;
	}
	
	public static double[][] copyMatrixEliminateRowAndColumn(double[][] matrix,
			int rowToEliminate, int colToEliminate) {
		double[][] newMatrix = new double[matrix.length - 1][matrix[0].length - 1];
		for (int r = 0; r < matrix.length; r++) {
			if (r == rowToEliminate) {
				continue;
			}
			for (int c = 0; c < matrix.length; c++) {
				if (c == colToEliminate) {
					continue;
				}
				int newRow = r;
				int newCol = c;
				if (newRow > rowToEliminate) {
					newRow--;
				}
				if (newCol > colToEliminate) {
					newCol--;
				}
				newMatrix[newRow][newCol] = matrix[r][c];
			}
		}
		return newMatrix;
	}
	
	public static double[] extractSelectedTimePoints(double[] data, int[] timePoints) {
		double[] extracted = new double[timePoints.length];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	public static int[] extractSelectedTimePoints(int[] data, int[] timePoints) {
		int[] extracted = new int[timePoints.length];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	public static double[][] extractSelectedTimePoints(double[][] data, int[] timePoints) {
		int columns = data[0].length;
		double[][] extracted = new double[timePoints.length][columns];
		for (int t = 0; t < timePoints.length; t++) {
			System.arraycopy(data[timePoints[t]], 0, extracted[t], 0, columns);
		}
		return extracted;
	}

	/**
	 * Extraxts the double[] vectors at each of the selected time points.
	 * The return double[][] array is an array of points to the existing
	 * double[] vectors.
	 * 
	 * @param data
	 * @param timePoints
	 * @return
	 */
	public static double[][] extractSelectedTimePointsReusingArrays(double[][] data, int[] timePoints) {
		double[][] extracted = new double[timePoints.length][];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	/**
	 * Extracts the boolean[] vectors at each of the selected time points.
	 * The return boolean[][] array is an array of points to the existing
	 * boolean[] vectors.
	 * 
	 * @param data
	 * @param timePoints
	 * @return
	 */
	public static boolean[][] extractSelectedTimePointsReusingArrays(boolean[][] data, int[] timePoints) {
		boolean[][] extracted = new boolean[timePoints.length][];
		for (int t = 0; t < timePoints.length; t++) {
			extracted[t] = data[timePoints[t]];
		}
		return extracted;
	}

	
	public static double[][] extractSelectedTimePoints(double[][] data, int[][] timePoints,
			int columnInTimePoints) {
		int columns = data[0].length;
		double[][] extracted = new double[timePoints.length][columns];
		for (int t = 0; t < timePoints.length; t++) {
			System.arraycopy(data[timePoints[t][columnInTimePoints]], 0, extracted[t], 0, columns);
		}
		return extracted;
	}

	/**
	 * Extract from data the vectors for rows corresponding to the time values in
	 *  column columnInTimePoints of each row of timePoints.
	 * 
	 * @param data
	 * @param timePoints
	 * @param columnInTimePoints
	 * @param timeOffset
	 * @return a 2D array of doubles, with timePoints.length rows and data[0].length columns
	 */
	public static double[][] extractSelectedTimePoints(double[][] data, int[][] timePoints,
			int columnInTimePoints, int timeOffset) {
		int columns = data[0].length;
		double[][] extracted = new double[timePoints.length][columns];
		for (int t = 0; t < timePoints.length; t++) {
			System.arraycopy(data[timePoints[t][columnInTimePoints] + timeOffset], 0, extracted[t], 0, columns);
		}
		return extracted;
	}
	
	/**
	 * Return the rows of data, where the conditionalData matched the
	 *  conditionalValue for that given row.
	 * Assumes data.length == conditionalData.length.
	 * 
	 * @param data
	 * @param conditionalData
	 * @param conditionalValue
	 * @return a 2D array of doubles where the conditionalData matched the
	 *  conditionalValue for those rows.
	 */
	public static double[][] extractSelectedPointsMatchingCondition(
			double[][] data, int[] conditionalData, int conditionalValue) {
		
		// Count the number of matching points first.
		int numNewRows = 0;
		for (int t = 0; t < data.length; t++) {
			if (conditionalData[t] == conditionalValue) {
				numNewRows++;
			}
		}
		// Create the new extracted data
		return extractSelectedPointsMatchingCondition(data, conditionalData,
				conditionalValue, numNewRows);
	}

	/**
	 * Return the rows of data, where the conditionalData matched the
	 *  conditionalValue for that given row.
	 * Assumes data.length == conditionalData.length.
	 * Here, the caller knows that there will be at minimum knownNumExtractedValues
	 *  values to be extracted, and only wants those values.
	 * 
	 * @param data
	 * @param conditionalData
	 * @param conditionalValue
	 * @param knownNumExtractedValues the known number of matching values
	 * @return a 2D array of doubles where the conditionalData matched the
	 *  conditionalValue for those rows.
	 */
	public static double[][] extractSelectedPointsMatchingCondition(
			double[][] data, int[] conditionalData, int conditionalValue,
			int knownNumExtractedValues) {
		
		// Create the new extracted data
		int columns = data[0].length;
		double[][] extracted = new double[knownNumExtractedValues][columns];
		int rowsCopied = 0;
		if (knownNumExtractedValues == 0) {
			return extracted;
		}
		for (int t = 0; t < data.length; t++) {
			if (conditionalData[t] == conditionalValue) {
				System.arraycopy(data[t], 0, extracted[rowsCopied++],
						0, columns);
			}
			if (rowsCopied == knownNumExtractedValues) {
				// We've extracted enough values
				break;
			}
		}
		return extracted;
	}

	/**
	 * Inserts the given time points (in the order prescribed in timePoints) 
	 * from the vector originalSourceValuesInJoint into the given column in matrix
	 * 
	 * @param inputValues
	 * @param timePoints
	 * @param matrix
	 * @param column
	 */
	public static void reorderVectorIntoMatrix(double[] inputValues, int[] timePoints,
			double[][] matrix, int column) {
		for (int i = 0; i < timePoints.length; i++) {
			int t = timePoints[i];
			matrix[i][column] = inputValues[t];
		}
	}
	
	/**
	 * Return data[x][y]:
	 *  - y==0: inputValues[x][0]
	 *  - y>0:  inputValues[reordering[y-1][x]][y]
	 * 
	 * @param inputValues holds the raw data values
	 * @param reordering outlines how to rearrange the raw data values for each variable or column.
	 *  First index is variable
	 *  or column number. Reorderings may be supplied for all of the columns of the inputValues,
	 *  or for one less than all of the columns, in which case the first column is not
	 *  reordered. Second index is for the row number or time step. The value at that
	 *  point states which row number to pull the data from.
	 * @return
	 */
	public static double[][] reorderDataForVariables(double[][] inputValues, int[][] reordering) {
		int rows = inputValues.length;
		int columns = inputValues[0].length;
		boolean reorderingFirstColumn = (reordering.length == columns);
		double[][] data = new double[rows][columns];
		for (int r = 0; r < rows; r++) {
			int reorderIndex = 0;
			if (reorderingFirstColumn) {
				data[r][0] = inputValues[reordering[reorderIndex++][r]][0];
			} else {
				data[r][0] = inputValues[r][0];
			}
			for (int c = 1; c < columns; c++) {
				data[r][c] = inputValues[reordering[reorderIndex++][r]][c];
			}
		}
		return data;
	}
	
	/**
	 * Reshapes the given single dimensional array into a 2D array of the given
	 *  size
	 * 
	 * @param data
	 * @param rows
	 * @param columns
	 * @return
	 */
	public static double[][] reshape(double[] data, int rows, int columns) {
		double[][] matrix = new double[rows][columns];
		int i = 0;
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < columns; c++) {
				matrix[r][c] = data[i++];
			}
		}
		return matrix;
	}
	
	/**
	 * Constructs all embedding vectors of size k for the data.
	 * There will be (data.length - k + 1) of these vectors returned.
	 * 
	 * @param data time series data
	 * @param k embedding length
	 * @return An array of k-length embedding vectors
	 */
	public static double[][] makeDelayEmbeddingVector(double[] data, int k) {
		try {
			return makeDelayEmbeddingVector(data, k, k - 1, data.length - k + 1);
		} catch (Exception e) {
			// The above call should not throw an Exception, handle here 
			//  in a RuntimeException so this method doesn't throw one
			throw new RuntimeException(e);
		}
	}

	/**
	 * Constructs numEmbeddingVectors embedding vectors of size k for the data,
	 * with the first embedding vector having it's last time point at t=startKthPoint
	 * 
	 * @param data time series data
	 * @param k embedding length
	 * @param startKthPoint last time point of the first embedding vector
	 *   (i.e. use k-1 if you want to go from the start)
	 * @param numEmbeddingVectors the number of embedding vectors to return
	 *   (i.e. use data.length-k+1 if you go from the start and want all
	 *   of them extracted)
	 * @return a 2D array of k-length embedding vectors.
	 */
	public static double[][] makeDelayEmbeddingVector(double[] data, int k,
			int startKthPoint, int numEmbeddingVectors) throws Exception {
		if (startKthPoint < k - 1) {
			throw new Exception("Start point t=" + startKthPoint + " is too early for a " +
					k + " length embedding vector");
		}
		if (numEmbeddingVectors + startKthPoint > data.length) {
			throw new Exception("Too many embedding vectors " + numEmbeddingVectors +
					" requested for the given startPoint " + startKthPoint +
					" and time series length " + data.length);
		}
		double[][] embeddingVectors = new double[numEmbeddingVectors][k];
		for (int t = startKthPoint; t < numEmbeddingVectors + startKthPoint; t++) {
			for (int i = 0; i < k; i++) {
				embeddingVectors[t - startKthPoint][i] = data[t - i];
			}
		}
		return embeddingVectors;
	}

	/**
	 * Constructs numEmbeddingVectors embedding vectors of size k for the data,
	 *  with embedding delay tau between each time sample for the vectors,
	 *  with the first embedding vector having it's last time point at t=startKthPoint
	 * 
	 * @param data time series data
	 * @param k embedding length
	 * @param tau embedding delay between each point in the original time series
	 * 		selected into each embedding vector
	 * @param startKthPoint last time point of the first embedding vector
	 *   (i.e. use (k-1)*tau if you want to go from the start)
	 * @param numEmbeddingVectors the number of embedding vectors to return
	 *   (i.e. use data.length-(k-1)*tau if you go from the start and want all
	 *   of them extracted)
	 * @return  a 2D array of k-length embedding vectors.
	 */
	public static double[][] makeDelayEmbeddingVector(double[] data, int k, int tau,
			int startKthPoint, int numEmbeddingVectors) throws Exception {
		if (startKthPoint < (k - 1)*tau) {
			throw new Exception("Start point t=" + startKthPoint + " is too early for a " +
					k + " length embedding vector with delay " + tau);
		}
		if (numEmbeddingVectors + startKthPoint > data.length) {
			throw new Exception("Too many embedding vectors " + numEmbeddingVectors +
					" requested for the given startPoint " + startKthPoint +
					" and time series length " + data.length);
		}
		double[][] embeddingVectors = new double[numEmbeddingVectors][k];
		for (int t = startKthPoint; t < numEmbeddingVectors + startKthPoint; t++) {
			for (int i = 0; i < k; i++) {
				embeddingVectors[t - startKthPoint][i] = data[t - i*tau];
			}
		}
		return embeddingVectors;
	}

	/**
	 * Constructs all embedding vectors of k time points for the data, including
	 *  all multivariate values at each time point.
	 * Will be data.length - k + 1 of these vectors returned
	 * 
	 * @param data 2D time series data (time is first second, second is variable number),
	 *   all of which is embedded
	 * @param k embedding length (i.e. number of time extractions for each vector)
	 * @return a 2D array of embedding vectors, which are of length
	 *   k x data[0].length.
	 */
	public static double[][] makeDelayEmbeddingVector(double[][] data, int k) {
		try {
			return makeDelayEmbeddingVector(data, k, k - 1, data.length - k + 1);
		} catch (Exception e) {
			// The above call should not throw an Exception, handle here 
			//  in a RuntimeException so this method doesn't throw one
			throw new RuntimeException(e);
		}
	}

	/**
	 * Constructs numEmbeddingVectors embedding vectors of k time points for the data, including
	 *  all multivariate values at each time point.
	 * Return only a subset, with the first embedding vector having it's last time point at t=startKthPoint
	 * 
	 * @param data 2D time series data (time is first second, second is variable number),
	 *   all of which is embedded
	 * @param k embedding length (i.e. number of time extractions for each vector)
	 * @param startKthPoint last time point of the first embedding vector
	 *   (i.e. use k-1 if you want to go from the start)
	 * @param numEmbeddingVectors the number of embedding vectors to return
	 *   (i.e. use data.length-k+1 if you go from the start and want all
	 *   of them extracted)
	 * @return a 2D array of embedding vectors, which are each of length
	 *   k x data[0].length.
	 */
	public static double[][] makeDelayEmbeddingVector(double[][] data, int k,
			int startKthPoint, int numEmbeddingVectors) throws Exception {
		if (startKthPoint < k - 1) {
			throw new Exception("Start point t=" + startKthPoint + " is too early for a " +
					k + " length embedding vector");
		}
		if (numEmbeddingVectors + startKthPoint > data.length) {
			throw new Exception("Too many embedding vectors " + numEmbeddingVectors +
					" requested for the given startPoint " + startKthPoint +
					" and time series length " + data.length);
		}
		int columns = data[0].length;
		double[][] embeddingVectors = new double[numEmbeddingVectors][k * columns];
		for (int t = startKthPoint; t < numEmbeddingVectors + startKthPoint; t++) {
			for (int i = 0; i < k; i++) {
				for (int c = 0; c < columns; c++) {
					embeddingVectors[t - startKthPoint][i*columns + c] = data[t - i][c];
				}
			}
		}
		return embeddingVectors;
	}

	/**
	 * Constructs numEmbeddingVectors embedding vectors of k time points for the data, including
	 *  all multivariate values at each time point,
	 *  with embedding delay tau between each time sample for the vectors,
	 *  with the first embedding vector having it's last time point at t=startKthPoint
	 * 
	 * @param data 2D time series data (time is first second, second is variable number),
	 *   all of which is embedded
	 * @param k embedding length (i.e. number of time extractions for each vector)
	 * @param tau embedding delay between each point in the original time series
	 * 		selected into each embedding vector
	 * @param startKthPoint last time point of the first embedding vector
	 *   (i.e. use k-1 if you want to go from the start)
	 * @param numEmbeddingVectors the number of embedding vectors to return
	 *   (i.e. use data.length-k+1 if you go from the start and want all
	 *   of them extracted)
	 * @return a 2D array of numEmbeddingVectors embedding vectors, which are each of length
	 *   k x data[0].length.
	 */
	public static double[][] makeDelayEmbeddingVector(double[][] data, int k, int tau,
			int startKthPoint, int numEmbeddingVectors) throws Exception {
		if (startKthPoint < (k - 1)*tau) {
			throw new Exception("Start point t=" + startKthPoint + " is too early for a " +
					k + " length embedding vector with delay " + tau);
		}
		if (numEmbeddingVectors + startKthPoint > data.length) {
			throw new Exception("Too many embedding vectors " + numEmbeddingVectors +
					" requested for the given startPoint " + startKthPoint +
					" and time series length " + data.length);
		}
		int columns = data[0].length;
		double[][] embeddingVectors = new double[numEmbeddingVectors][k * columns];
		for (int t = startKthPoint; t < numEmbeddingVectors + startKthPoint; t++) {
			for (int i = 0; i < k; i++) {
				for (int c = 0; c < columns; c++) {
					embeddingVectors[t - startKthPoint][i*columns + c] = data[t - i*tau][c];
				}
			}
		}
		return embeddingVectors;
	}

	/**
	 * Constructs numEmbeddingVectors embedding vectors of k time points for a single column of
	 *  the data,
	 *  with embedding delay tau between each time sample for the vectors,
	 *  with the first embedding vector having it's last time point at t=startKthPoint
	 * 
	 * @param data 2D time series data (time is first second, second is variable number),
	 *   only one particular column of which is embedded
	 * @param column the column index to embed
	 * @param k embedding length (i.e. number of time extractions for each vector)
	 * @param tau embedding delay between each point in the original time series
	 * 		selected into each embedding vector
	 * @param startKthPoint last time point of the first embedding vector
	 *   (i.e. use k-1 if you want to go from the start)
	 * @param numEmbeddingVectors the number of embedding vectors to return
	 *   (i.e. use data.length-k+1 if you go from the start and want all
	 *   of them extracted)
	 * @return a 2D array of numEmbeddingVectors embedding vectors, which are each of length k.
	 */
	public static double[][] makeDelayEmbeddingVector(double[][] data, int column, int k, int tau,
			int startKthPoint, int numEmbeddingVectors) throws Exception {
		if (startKthPoint < (k - 1)*tau) {
			throw new Exception("Start point t=" + startKthPoint + " is too early for a " +
					k + " length embedding vector with delay " + tau);
		}
		if (numEmbeddingVectors + startKthPoint > data.length) {
			throw new Exception("Too many embedding vectors " + numEmbeddingVectors +
					" requested for the given startPoint " + startKthPoint +
					" and time series length " + data.length);
		}
		double[][] embeddingVectors = new double[numEmbeddingVectors][k];
		for (int t = startKthPoint; t < numEmbeddingVectors + startKthPoint; t++) {
			for (int i = 0; i < k; i++) {
				embeddingVectors[t - startKthPoint][i] = data[t - i*tau][column];
			}
		}
		return embeddingVectors;
	}

	public static double stdDev(double[] array) {
		double mean = 0.0;
		double total = 0.0;
		for (int m = 0; m < array.length; m++) {
			total += array[m];
		}
		mean = total / (double) array.length;
		
		return stdDev(array, mean);
	}
	
	public static double stdDev(double[][] matrix, int column) {
		double mean = 0.0;
		double total = 0.0;
		for (int m = 0; m < matrix.length; m++) {
			total += matrix[m][column];
		}
		mean = total / (double) matrix.length;
		
		return stdDev(matrix, column, mean);
	}
	
	/**
	 * Return the standard deviation of all the elements in array
	 * 
	 * @param array
	 * @param mean
	 * @return
	 */
	public static double stdDev(double[] array, double mean) {
		return stdDev(array, mean, array.length);
	}
	
	/**
	 * Standard deviation for the first arrayLength terms of array
	 * 
	 * @param array
	 * @param mean
	 * @param arrayLength
	 * @return
	 */
	public static double stdDev(double[] array, double mean, int arrayLength) {
		if (arrayLength == 0) {
			return 0.0;
		}
		double sumSqs = 0.0;
		for (int m = 0; m < arrayLength; m++) {
			sumSqs += (array[m] - mean) * (array[m] - mean);
		}
		double std = sumSqs / (double) (arrayLength - 1);
		std = Math.sqrt(std);
		return std;
	}
	
	/**
	 * Compute the standard deviation along the given column, with the known
	 *  given mean.
	 * 
	 * @param matrix
	 * @param column
	 * @param mean
	 * @return
	 */
	public static double stdDev(double[][] matrix, int column, double mean) {
		if (matrix.length == 0) {
			return 0.0;
		}
		double sumSqs = 0.0;
		for (int m = 0; m < matrix.length; m++) {
			sumSqs += (matrix[m][column] - mean) * (matrix[m][column] - mean);
		}
		double std = sumSqs / (double) (matrix.length - 1);
		std = Math.sqrt(std);
		return std;
	}
	
	/**
	 * Compute the standard deviation across all values in the 2D matrix
	 * 
	 * @param matrix
	 * @return
	 */
	public static double stdDev(double[][] matrix) {
		double mean = mean(matrix);
		return stdDev(matrix, mean);
	}
	
	/**
	 * Compute the standard deviation across all values in the 2D matrix
	 * 
	 * @param matrix
	 * @param mean
	 * @return
	 */
	public static double stdDev(double[][] matrix, double mean) {
		if (matrix.length == 0) {
			return 0.0;
		}
		double sumSqs = 0.0;
		for (int m = 0; m < matrix.length; m++) {
			for (int c = 0; c < matrix[m].length; c++) {
				sumSqs += (matrix[m][c] - mean) * (matrix[m][c] - mean);
			}
		}
		double std = sumSqs / (double) ((matrix.length * matrix[0].length) - 1);
		std = Math.sqrt(std);
		return std;
	}


	/**
	 * Compute the standard deviations along each column
	 * 
	 * @param matrix
	 * @param means
	 * @return
	 */
	public static double[] stdDevs(double[][] matrix, double[] means) {
		double[] sumSqs = new double[means.length];
		for (int m = 0; m < matrix.length; m++) {
			for (int c = 0; c < matrix[m].length; c++) {
				sumSqs[c] += (matrix[m][c] - means[c]) * (matrix[m][c] - means[c]);
			}
		}
		double[] stds = new double[means.length];
		for (int c = 0; c < stds.length; c++) {
			stds[c] = sumSqs[c] / (double) (matrix.length - 1);
			stds[c] = Math.sqrt(stds[c]);
		}
		return stds;
	}

	/**
	 * Compute the standard deviations along each row
	 * 
	 * @param matrix
	 * @param means
	 * @return
	 */
	public static double[] stdDevsOfRows(double[][] matrix, double[] means) {
		double[] stds = new double[matrix.length];
		for (int r = 0; r < matrix.length; r++) {
			double sumSqs = 0.0;
			for (int c = 0; c < matrix[r].length; c++) {
				sumSqs += (matrix[r][c] - means[r]) * (matrix[r][c] - means[r]);
			}
			stds[r] = sumSqs / (double) (matrix[r].length - 1);
			stds[r] = Math.sqrt(stds[r]);
		}
		return stds;
	}

	public static double max(double[][][] matrix) {
		// double max = 0.0;
		double max = matrix[0][0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				for (int k = 0; k < matrix[i][j].length; k++) {
					if (matrix[i][j][k] > max) {
						max = matrix[i][j][k];
					}
				}
			}
		}
		return max;
	}

	/**
	 * Normalises the elements in the given array
	 * 
	 * @param array
	 */
	public static void normalise(double[] array) {
		double mean = MatrixUtils.mean(array);
		double stdDev = MatrixUtils.stdDev(array, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < array.length; t++) {
				array[t] = (array[t] - mean);
			}			
		} else {
			// stdDev is non zero
			for (int t = 0; t < array.length; t++) {
				array[t] = (array[t] - mean) / stdDev;
			}
		}
	}
	
	/**
	 * Returns a normalised array of the elements in the given array 
	 * 
	 * @param array
	 */
	public static double[] normaliseIntoNewArray(double[] array) {
		double[] newArray = new double[array.length];
		double mean = MatrixUtils.mean(array);
		double stdDev = MatrixUtils.stdDev(array, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < array.length; t++) {
				newArray[t] = (array[t] - mean);
			}
		} else {
			// stdDev is non zero
			for (int t = 0; t < array.length; t++) {
				newArray[t] = (array[t] - mean) / stdDev;
			}
		}
		return newArray;
	}

	/**
	 * Normalises the elements in the given column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 * @param column column number to be normalised
	 */
	public static void normalise(double[][] matrix, int column) {
		double mean = MatrixUtils.mean(matrix, column);
		double stdDev = MatrixUtils.stdDev(matrix, column, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < matrix.length; t++) {
				matrix[t][column] = (matrix[t][column] - mean);
			}
		} else {
			// stdDev is non zero
			for (int t = 0; t < matrix.length; t++) {
				matrix[t][column] = (matrix[t][column] - mean) / stdDev;
			}
		}
	}

	/**
	 * Normalises the elements in the given column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 * @param column column number to be normalised
	 */
	public static double[] normaliseIntoNewArray(double[][] matrix, int column) {
		double[] newArray = new double[matrix.length];
		double mean = MatrixUtils.mean(matrix, column);
		double stdDev = MatrixUtils.stdDev(matrix, column, mean);
		if (Double.isInfinite(1.0 / stdDev)) {
			// The stdDev is 0, just subtract off mean
			for (int t = 0; t < matrix.length; t++) {
				newArray[t] = (matrix[t][column] - mean);
			}
		} else {
			// stdDev is non zero
			for (int t = 0; t < matrix.length; t++) {
				newArray[t] = (matrix[t][column] - mean) / stdDev;
			}
		}
		return newArray;
	}

	/**
	 * Normalises the elements along each column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 */
	public static void normalise(double[][] matrix) {
		double[] means = means(matrix);
		double[] stds = stdDevs(matrix, means);
		
		boolean[] nonZeroStds = new boolean[stds.length];
		for (int c = 0; c < matrix[0].length; c++) {
			nonZeroStds[c] = !Double.isInfinite(1.0 /  stds[c]);
		}
		
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				matrix[r][c] = matrix[r][c] - means[c];
				if (nonZeroStds[c]) {
					matrix[r][c] /= stds[c];
				} // else we just subtract off the mean
			}
		}
	}

	/**
	 * Normalises the elements along each column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 */
	public static double[][] normaliseIntoNewArray(double[][] matrix) {
		double[] means = means(matrix);
		double[] stds = stdDevs(matrix, means);
		return normaliseIntoNewArray(matrix, means, stds);
	}

	/**
	 * Normalises the elements along each column of the matrix
	 * 
	 * @param matrix 2D matrix of doubles
	 */
	public static double[][] normaliseIntoNewArray(double[][] matrix, double[] means, double[] stds) {
		double[][] newMatrix = new double[matrix.length][matrix[0].length];
		for (int r = 0; r < newMatrix.length; r++) {
			for (int c = 0; c < newMatrix[r].length; c++) {
				newMatrix[r][c] = matrix[r][c] - means[c];
				if (!Double.isInfinite(1.0 /  stds[c])) {
					newMatrix[r][c] /= stds[c];
				} // else we just subtract off the mean
			}
		}
		return newMatrix;		
	}

	public static double max(double[][] matrix) {
		// double max = 0.0;
		double max = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (Double.isNaN(max) || (matrix[i][j] > max)) {
					max = matrix[i][j];
				}
			}
		}
		return max;
	}

	public static int max(int[][] matrix) {
		// int max = 0;
		int max = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (matrix[i][j] > max) {
					max = matrix[i][j];
				}
			}
		}
		return max;
	}

	public static double max(double[] array) {
		return maxStartFromIndex(array, 0);
	}

	public static double maxStartFromIndex(double[] array, int startFromIndex) {
		// double max = 0.0;
		double max = array[startFromIndex];
		for (int i = startFromIndex; i < array.length; i++) {
			// TODO Check where we used this and if it's still
			//  the approach we want to take
			if (Double.isNaN(max) || (array[i] > max)) {
				max = array[i];
			}
		}
		return max;
	}

	public static int maxIndex(double[] array) {
		// double max = 0.0;
		double max = array[0];
		int maxIndex = 0;
		for (int i = 1; i < array.length; i++) {
			if (array[i] > max) {
				max = array[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}

	public static int max(int[] array) {
		// int max = 0;
		int max = array[0];
		for (int i = 0; i < array.length; i++) {
			if (array[i] > max) {
				max = array[i];
			}
		}
		return max;
	}

	/**
	 * Works out the maximum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static double max(double[][] matrix, int column) {
		// double max = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		double max = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (Double.isNaN(max) || (matrix[i][column] > max)) {
				max = matrix[i][column];
			}
		}
		return max;
	}

	/**
	 * Works out the maximum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static int max(int[][] matrix, int column) {
		// double max = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		int max = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (matrix[i][column] > max) {
				max = matrix[i][column];
			}
		}
		return max;
	}

	public static double min(double[][] matrix) {
		// double min = 0.0;
		double min = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (Double.isNaN(min) || (matrix[i][j] < min)) {
					min = matrix[i][j];
				}
			}
		}
		return min;
	}

	public static int min(int[][] matrix) {
		// int min = 0;
		int min = matrix[0][0];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (matrix[i][j] < min) {
					min = matrix[i][j];
				}
			}
		}
		return min;
	}

	public static double min(double[] array) {
		return minStartFromIndex(array, 0);
	}

	/**
	 * Find the kth minimum value in the array.
	 * 
	 * @param array
	 * @param k
	 * @return
	 * @throws Exception 
	 */
	public static double kthMin(double[] array, int k) throws Exception {
		if (k == 1) {
			return min(array); 
		}
		if (array.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					array.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
		}
		for (int t = 0; t < array.length; t++) {
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (array[t] < mins[k - 1]) {
				mins[k - 1] = array[t];
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (array[t] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						mins[i] = array[t];
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the kth min
		return mins[k-1];
	}

	/**
	 * Find the kth minimum value in the array subject to
	 *  a given condition.
	 * Assumes that the condition is satisfied at least k
	 *  times in the array (this is not checked in here)
	 * 
	 * @param array
	 * @param k
	 * @param extraData
	 * @param extraCondition
	 * @return
	 * @throws Exception 
	 */
	public static double kthMinSubjectTo(double[] array, int k, int[] extraData, int condition) throws Exception {
		// Can't do a quickie for k==1 here since we're subject to 
		//  checking the condition
		if (array.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					array.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
		}
		for (int t = 0; t < array.length; t++) {
			if (extraData[t] != condition) {
				continue;
			}
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (array[t] < mins[k - 1]) {
				mins[k - 1] = array[t];
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (array[t] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						mins[i] = array[t];
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the kth min
		return mins[k-1];
	}

	public static double minIgnoreIndex(double[] array, int indexToIgnore) {
		// double min = 0.0;
		double min;
		if (indexToIgnore != 0) {
			min = array[0];
		} else {
			min = array[1];
		}
		for (int i = 0; i < array.length; i++) {
			if (indexToIgnore == i) {
				continue;
			}
			if (Double.isNaN(min) || (array[i] < min)) {
				min = array[i];
			}
		}
		return min;
	}

	public static double minStartFromIndex(double[] array, int startFromIndex) {
		// double min = 0.0;
		double min = array[startFromIndex];
		for (int i = startFromIndex; i < array.length; i++) {
			if (Double.isNaN(min) || (array[i] < min)) {
				min = array[i];
			}
		}
		return min;
	}

	public static int min(int[] array) {
		// int min = 0;
		int min = array[0];
		for (int i = 0; i < array.length; i++) {
			if (array[i] < min) {
				min = array[i];
			}
		}
		return min;
	}
	
	/**
	 * Works out the minimum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static double min(double[][] matrix, int column) {
		// double min = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		double min = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (Double.isNaN(min) || (matrix[i][column] < min)) {
				min = matrix[i][column];
			}
		}
		return min;
	}

	/**
	 * Works out the minimum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static int min(int[][] matrix, int column) {
		// double min = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		int min = matrix[0][column];
		for (int i = 1; i < matrix.length; i++) {
			if (matrix[i][column] < min) {
				min = matrix[i][column];
			}
		}
		return min;
	}

	/**
	 * Works out the index of the minimum value in the matrix in a given column 
	 * 
	 * @param matrix
	 * @param column
	 * @return
	 */
	public static int minIndex(double[][] matrix, int column) {
		// double min = 0.0;
		// Allow ArrayIndexOutOfBoundsException if matrix is size 0
		double min = matrix[0][column];
		int minIndex = 0;
		for (int i = 1; i < matrix.length; i++) {
			if (Double.isNaN(min) || (matrix[i][column] < min)) {
				min = matrix[i][column];
				minIndex = i;
			}
		}
		return minIndex;
	}

	/**
	 * Works out the index of the k minimum values in the matrix in a given column 
	 * 
	 * @param matrix data
	 * @param column which column of the data to find the min values from
	 * @param k how many min values to return
	 * @return an array of the (row) indices in the array with the k min values,
	 *  with closest match first.
	 * @throws Exception 
	 */
	public static int[] kMinIndices(double[][] matrix, int column, int k) throws Exception {
		if (matrix.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					matrix.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		int[] minIndices = new int[k];
		if (k == 1) {
			minIndices[0] = minIndex(matrix, column);
			return minIndices;
		}
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
			minIndices[i] = -1;
		}
		for (int t = 0; t < matrix.length; t++) {
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (matrix[t][column] < mins[k - 1]) {
				mins[k - 1] = matrix[t][column];
				minIndices[k-1] = t;
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (matrix[t][column] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						minIndices[i+1] = minIndices[i];
						mins[i] = matrix[t][column];
						minIndices[i] = t;
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the indices of the k mins
		return minIndices;
	}

	/**
	 * Works out the index of the k minimum values in the matrix in a given column 
	 *  subject to the extraData matching a given condition.
	 * We do not check whether there are k matches for the extraData to
	 *  the condition here - the caller should check this themselves.
	 * 
	 * @param matrix
	 * @param column
	 * @param k
	 * @param extraData
	 * @param condition
	 * @return
	 * @throws Exception 
	 */
	public static int[] kMinIndicesSubjectTo(double[][] matrix, int column,
			int k, int[] extraData, int condition) throws Exception {
		if (matrix.length < k) {
			throw new Exception(String.format("Length of array (%d) is less than k (%d)",
					matrix.length, k));
		}
		// Hold the k minimum elements in strictly increasing order from 0 .. k-1
		double[] mins = new double[k];
		int[] minIndices = new int[k];
		// no quick check for k==1 since we need to check the extra condition
		for (int i = 0; i < k; i++) {
			mins[i] = Double.POSITIVE_INFINITY;
			minIndices[i] = -1;
		}
		for (int t = 0; t < matrix.length; t++) {
			if (extraData[t] != condition) {
				continue;
			}
			// Assume that k is small enough that there is no point doing binary
			//  searches to find the best place to insert this element in the minimums (if required).
			// First check if it's smaller than the current kth min:
			if (matrix[t][column] < mins[k - 1]) {
				mins[k - 1] = matrix[t][column];
				minIndices[k-1] = t;
				// Now check if we need to reorder the array of minimums, keeping it sorted
				for (int i = k - 2; i >= 0; i--) {
					if (matrix[t][column] < mins[i]) {
						// Swap array[t] along from mins[i+1]:
						mins[i+1] = mins[i];
						minIndices[i+1] = minIndices[i];
						mins[i] = matrix[t][column];
						minIndices[i] = t;
						continue;
					}
					// else no need to keep checking the array is sorted correctly
					break;
				}
			}
		}
		// Return the index of the kth min
		return minIndices;
	}

	/**
	 * Sort array and return the original indices of each item in the
	 * sorted list, such that array[returnValue[k]] is the kth item in the
	 * sorted list.
	 * Sorting is done from smallest to largest.
	 * 
	 * @param array array of doubles to sort
	 * @return list of original indices, in the sorted order of the array
	 */
	public static int[] sortIndices(double[] array) {
		// Need an instance of MatrixUtils to get to member classes:
		MatrixUtils mUtils = new MatrixUtils();
		
		// First create the array of DoubleWithIndexForSort objects:
		DoubleWithIndexForSort[] objectArray = new DoubleWithIndexForSort[array.length];
		for (int i = 0; i < array.length; i++) {
			objectArray[i] = mUtils.new DoubleWithIndexForSort(array[i], i);
		}
		// Sort the array with original indices in place:
		Arrays.sort(objectArray, mUtils.new DoubleWithIndexForSortComparator());
		// Pull out the original indices:
		int[] arrayOfOriginalIndices = new int[array.length];
		for (int i = 0; i < array.length; i++) {
			arrayOfOriginalIndices[i] = objectArray[i].originalIndex;
		}
		return arrayOfOriginalIndices;
	}
	
	/**
	 * Structure used by {@link MatrixUtils#sortIndices(double[])}
	 * 
	 * @author Joseph Lizier
	 *
	 */
	private class DoubleWithIndexForSort {
		double value;
		int originalIndex;
		
		public DoubleWithIndexForSort(double value, int originalIndex) {
			this.value = value;
			this.originalIndex = originalIndex;
		}
	}
	
	/**
	 * Comparator for {@link DoubleWithIndexForSort}
	 * 
	 * @author Joseph Lizier
	 *
	 */
	private class DoubleWithIndexForSortComparator implements Comparator<DoubleWithIndexForSort> {

		public int compare(DoubleWithIndexForSort arg0,
				DoubleWithIndexForSort arg1) {
			if (arg0.value < arg1.value) {
				return -1;
			}
			if (arg0.value > arg1.value) {
				return 1;
			}
			// values are equal:
			return 0;
		}
		
	}
	
	/**
	 * Mirrors the matrix in both coordinates
	 * 
	 * @param matrix
	 * @return
	 */
	public static int[][] mirrorMatrixBothCoords(int[][] matrix) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		
		int[][] mirrored = new int[rows][cols];
		
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				mirrored[(rows - 1) - r][(cols - 1) - c] = matrix[r][c];
			}
		}
		return mirrored;
	}
	
	/**
	 * Mirrors the matrix in both coordinates
	 * 
	 * @param matrix
	 * @return
	 */
	public static double[][] mirrorMatrixBothCoords(double[][] matrix) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		
		double[][] mirrored = new double[rows][cols];
		
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				mirrored[(rows - 1) - r][(cols - 1) - c] = matrix[r][c];
			}
		}
		return mirrored;
	}

	/**
	 * Moves the rows of the array up by upBy.
	 * Inserts zeros at the bottom
	 * 
	 * @param matrix
	 * @param upBy
	 */
	public static void moveRowsUp(double[][] matrix, int upBy) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		
		for (int r = 0; r < rows - upBy; r++) {
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = matrix[r + upBy][c];
			}
		}
		for (int r = rows - upBy; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = 0;
			}
		}
	}
	
	/**
	 * <p>Returns the covariance between the two arrays of data.</p>
	 * <p>See - <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @return the covariance
	 */
	public static double covariance(double[] x, double[] y) {
		double c = 0;
		double meanX = mean(x);
		double meanY = mean(y);
		for (int t = 0; t < x.length; t++) {
			c += (x[t] - meanX)*(y[t]-meanY);
		}
		return c / (double) (x.length - 1); // -1 for sample covariance
	}

	/**
	 * <p>Returns the covariance between the two arrays of data, with
	 * a given lag between the first and second.</p>
	 * <p>See - <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 * </p>
	 * 
	 * @param x time series 1
	 * @param y time series 2
	 * @param delay delay >= 0 to compute the covariance across (from first to second time series)
	 * @return the covariance
	 */
	public static double covariance(double[] x, double[] y, int delay) {
		double meanX = 0, meanY = 0;
		// No error checking if y is same length as x
		for (int n = 0; n < x.length - delay; n++) {
			meanX += x[n];
			meanY += y[n + delay];
		}
		meanX /= (double) (x.length - delay);
		meanY /= (double) (x.length - delay);
		double c = 0;
		for (int t = 0; t < x.length - delay; t++) {
			c += (x[t] - meanX)*(y[t + delay]-meanY);
		}
		return c / (double) (x.length - delay - 1); // -1 for sample covariance
	}

	/**
	 * <p>Returns the covariance between the first two columns of data.</p>
	 * 
	 * @param data
	 * @return the covariance
	 * @see <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 */
	public static double covarianceFirstTwoColumns(double[][] data) {
		return covarianceTwoColumns(data, 0, 1);
	}
	
	/**
	 * <p>Returns the covariance between two columns of data in
	 *  a multivariate array.</p>
	 * <p>See - <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 * </p>
	 * 
	 * @param data multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param col1 variable number 1 to compute the covariance to
	 * @param col2 variable number 2 to compute the covariance to
	 * @return the covariance
	 */
	public static double covarianceTwoColumns(double[][] data, int col1, int col2) {
		double mean1 = mean(data, col1);
		double mean2 = mean(data, col2);
		return covarianceTwoColumns(data, col1, col2, mean1, mean2);
	}
	
	/**
	 * <p>Returns the covariance between two columns of data in
	 *  a multivariate array.</p>
	 * <p>See - <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 * </p>
	 * 
	 * @param data multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param col1 variable number 1 to compute the covariance to
	 * @param col2 variable number 2 to compute the covariance to
	 * @param mean1 mean of variable 1
	 * @param mean2 mean of variable 2
	 * @return the covariance
	 */
	public static double covarianceTwoColumns(double[][] data, int col1, int col2,
			double mean1, double mean2) {
		double c = 0;
		for (int t = 0; t < data.length; t++) {
			c += (data[t][col1] - mean1)*(data[t][col2]-mean2);
		}
		return c / (double) (data.length - 1);
	}
	
	/**
	 * <p>Returns the covariance between two columns of data in
	 *  two multivariate arrays.</p>
	 * <p>See - <a href="http://mathworld.wolfram.com/Covariance.html">Mathworld</a>
	 * </p>
	 * 
	 * @param data1 first multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param data2 second multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param col1 variable number 1 to compute the covariance to
	 * @param col2 variable number 2 to compute the covariance to
	 * @param mean1 mean of variable 1
	 * @param mean2 mean of variable 2
	 * @return the covariance
	 */
	public static double covarianceTwoColumns(
			double[][] data1, double[][] data2, int col1, int col2,
			double mean1, double mean2) {
		double c = 0;
		for (int t = 0; t < data1.length; t++) {
			c += (data1[t][col1] - mean1)*(data2[t][col2]-mean2);
		}
		return c / (double) (data1.length - 1);
	}
	
	/**
	 * Compute the covariance matrix between all column pairs (variables) in the
	 *  multivariate data set
	 * 
	 * @param data multivariate array of data; first index is time, second is 
	 *    variable number
	 * @return covariance matrix
	 */
	public static double[][] covarianceMatrix(double[][] data) {
		return covarianceMatrix(data, means(data));
	}
	
	/**
	 * Compute the covariance matrix between all column pairs (variables) in the
	 *  multivariate data set
	 * 
	 * @param data multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param means the mean of each variable (column) in the data
	 * @return covariance matrix
	 */
	public static double[][] covarianceMatrix(double[][] data, double[] means) {
		int numVariables = data[0].length;
		double[][] covariances = new double[numVariables][numVariables];
		for (int r = 0; r < numVariables; r++) {
			for (int c = r; c < numVariables; c++) {
				// Compute the covariance between variable r and c:
				covariances[r][c] = covarianceTwoColumns(data, r, c,
						means[r], means[c]);
				// And of course this is symmetric between c and r:
				covariances[c][r] = covariances[r][c];
			}
		}
		return covariances;
	}
	
	/**
	 * Compute the covariance matrix between all column pairs (variables) in the
	 *  multivariate data set, which consists of two separate
	 *  multivariate vectors.
	 * 
	 * @param data1 multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param data2 a second multivariate array of data, which can be though
	 *    of as extensions of rows of the first.
	 * @return covariance matrix, where the columns of dat1 are numbered
	 *   first, and the columns of data2 after that.
	 */
	public static double[][] covarianceMatrix(
			double[][] data1, double[][] data2) {
		return covarianceMatrix(data1, data2, 0);
	}

	/**
	 * Compute the covariance matrix between all column pairs (variables) in the
	 *  multivariate data set, which consists of two separate
	 *  multivariate vectors.
	 * 
	 * @param data1 multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param data2 a second multivariate array of data, which can be thought
	 *    of as extensions of rows of the first.
	 * @param delay compute the lagged covariance of the given delay from
	 *    data1 to data2 (assumes delay >= 0); i.e. compute correlation
	 *    between data1[x] and data2[x+delay]. 
	 * @return covariance matrix, where the columns of data1 are numbered
	 *   first, and the columns of data2 after that.
	 */
	public static double[][] covarianceMatrix(
			double[][] data1, double[][] data2, int delay) {
		if (delay > 0) {
			// Trim out the last delay rows of data1, and the
			//  first delay rows of data2:
			double[][] data1Trimmed = new double[data1.length - delay][];
			double[][] data2Trimmed = new double[data2.length - delay][];
			for (int x = 0; x < data1.length - delay; x++) {
				data1Trimmed[x] = data1[x];
				data2Trimmed[x] = data2[x + delay];
			}
			// Just overwrite our local copy of the pointers to the
			//  original data
			data1 = data1Trimmed;
			data2 = data2Trimmed;
		}
		
		int numVariables1 = data1[0].length;
		int numVariables2 = data2[0].length;
		int numVariables = numVariables1 + numVariables2;
		double[][] covariances = new double[numVariables][numVariables];
		// Compute means of each variable once up front to save time
		double[] means1 = new double[numVariables1];
		double[] means2 = new double[numVariables2];
		for (int r = 0; r < numVariables1; r++) {
			means1[r] = mean(data1, r);
		}
		for (int r = 0; r < numVariables2; r++) {
			means2[r] = mean(data2, r);
		}
		// Now compute the covariances:
		for (int r = 0; r < numVariables1; r++) {
			// Compute the covariances internal to data1:
			for (int c = r; c < numVariables1; c++) {
				// Compute the covariance between variable r and c:
				covariances[r][c] = covarianceTwoColumns(data1, r, c,
						means1[r], means1[c]);
				// And of course this is symmetric between c and r:
				covariances[c][r] = covariances[r][c];
			}
			// Compute the covariances between data1 and data2
			for (int c = 0; c < numVariables2; c++) {
				// Compute the covariance between variable r and c:
				covariances[r][numVariables1 + c] =
						covarianceTwoColumns(data1, data2,
								r, c, means1[r], means2[c]);
				// And of course this is symmetric between c and r:
				covariances[numVariables1 + c][r] =
						covariances[r][numVariables1 + c];
			}
		}
		// Now compute the covariances internal to data2:
		for (int r = 0; r < numVariables2; r++) {
			for (int c = r; c < numVariables2; c++) {
				// Compute the covariance between variable r and c:
				covariances[numVariables1 + r][numVariables1 + c] =
						covarianceTwoColumns(data2, r, c,
						means2[r], means2[c]);
				// And of course this is symmetric between c and r:
				covariances[numVariables1 + c][numVariables1 + r] =
						covariances[numVariables1 + r][numVariables1 + c];
			}
		}		
		return covariances;
	}

	/**
	 * Compute the covariance matrix between all column pairs (variables) in the
	 *  multivariate data set, which consists of three separate
	 *  multivariate vectors.
	 * 
	 * @param data1 multivariate array of data; first index is time, second is 
	 *    variable number
	 * @param data2 a second multivariate array of data, which can be thought
	 *    of as extensions of rows of the first.
	 * @param data2 a third multivariate array of data, which can be thought
	 *    of as extensions of rows of the first and second.
	 * @return covariance matrix, where the columns of data1 are numbered
	 *   first, the columns of data2 after that, and finally the columns
	 *   of data3.
	 */
	public static double[][] covarianceMatrix(
			double[][] data1, double[][] data2, double[][] data3) {
		int numVariables1 = data1[0].length;
		int numVariables2 = data2[0].length;
		int numVariables3 = data3[0].length;
		int numVariables = numVariables1 + numVariables2 + numVariables3;
		double[][] covariances = new double[numVariables][numVariables];
		// Compute means of each variable once up front to save time
		double[] means1 = new double[numVariables1];
		double[] means2 = new double[numVariables2];
		double[] means3 = new double[numVariables3];
		for (int r = 0; r < numVariables1; r++) {
			means1[r] = mean(data1, r);
		}
		for (int r = 0; r < numVariables2; r++) {
			means2[r] = mean(data2, r);
		}
		for (int r = 0; r < numVariables3; r++) {
			means3[r] = mean(data3, r);
		}
		// Now compute the covariances:
		for (int r = 0; r < numVariables1; r++) {
			// Compute the covariances internal to data1:
			for (int c = r; c < numVariables1; c++) {
				// Compute the covariance between variable r and c:
				covariances[r][c] = covarianceTwoColumns(data1, r, c,
						means1[r], means1[c]);
				// And of course this is symmetric between c and r:
				covariances[c][r] = covariances[r][c];
			}
			// Compute the covariances between data1 and data2
			for (int c = 0; c < numVariables2; c++) {
				// Compute the covariance between variable r and c:
				covariances[r][numVariables1 + c] =
						covarianceTwoColumns(data1, data2,
								r, c, means1[r], means2[c]);
				// And of course this is symmetric between c and r:
				covariances[numVariables1 + c][r] =
						covariances[r][numVariables1 + c];
			}
			// Compute the covariances between data1 and data3
			for (int c = 0; c < numVariables3; c++) {
				// Compute the covariance between variable r and c:
				covariances[r][numVariables1 + numVariables2 + c] =
						covarianceTwoColumns(data1, data3,
								r, c, means1[r], means3[c]);
				// And of course this is symmetric between c and r:
				covariances[numVariables1 + numVariables2 + c][r] =
						covariances[r][numVariables1 + numVariables2 + c];
			}
		}
		// Compute the other covariances for data2
		for (int r = 0; r < numVariables2; r++) {
			// Compute the covariances internal to data2:
			for (int c = r; c < numVariables2; c++) {
				// Compute the covariance between variable r and c:
				covariances[numVariables1 + r][numVariables1 + c] =
						covarianceTwoColumns(data2, r, c,
						means2[r], means2[c]);
				// And of course this is symmetric between c and r:
				covariances[numVariables1 + c][numVariables1 + r] =
						covariances[numVariables1 + r][numVariables1 + c];
			}
			// Compute the covariances between data2 and data3
			for (int c = 0; c < numVariables3; c++) {
				// Compute the covariance between variable r and c:
				covariances[numVariables1 + r][numVariables1 + numVariables2 + c] =
						covarianceTwoColumns(data2, data3,
								r, c, means2[r], means3[c]);
				// And of course this is symmetric between c and r:
				covariances[numVariables1 + numVariables2 + c][numVariables1 + r] =
						covariances[numVariables1 + r][numVariables1 + numVariables2 + c];
			}
		}		
		// Compute the internal covariances for data3
		for (int r = 0; r < numVariables3; r++) {
			for (int c = r; c < numVariables3; c++) {
				// Compute the covariance between variable r and c:
				covariances[numVariables1 + numVariables2 + r][numVariables1 + numVariables2 + c] =
						covarianceTwoColumns(data3, r, c,
						means3[r], means3[c]);
				// And of course this is symmetric between c and r:
				covariances[numVariables1 + numVariables2 + c][numVariables1 + numVariables2 + r] =
						covariances[numVariables1 + numVariables2 + r][numVariables1 + numVariables2 + c];
			}
		}
		return covariances;
	}

	/**
	 * <p>Returns the correlation between the two arrays of data.</p>
	 * <p>The arrays are asssumed to have the same lengths</p>
	 * <p>See - <a href="http://en.wikipedia.org/wiki/Correlation">Wikipedia</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @return the correlation
	 */
	public static double correlation(double[] x, double[] y) {
		return correlation(x, y, x.length);
	}
	
	/**
	 * <p>Returns the correlation between the two arrays of data.</p>
	 * <p>See - <a href="http://en.wikipedia.org/wiki/Correlation">Wikipedia</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @param dataLength - number of terms in each vector to consider (we look at the first dataLength terms).
	 * 	Precondition: dataLength is less than min(x.length, y.length)
	 * @return the correlation
	 */
	public static double correlation(double[] x, double[] y, int dataLength) {
		// return covariance(x, y) / stdDev(x) / stdDev(y);
		// Save some code time by reusing the code from inside covariance:
		double c = 0;
		double meanX = mean(x, 0, dataLength);
		double meanY = mean(y, 0, dataLength);
		for (int t = 0; t < dataLength; t++) {
			c += (x[t] - meanX)*(y[t]-meanY);
		}
		double covariance = c / (double) (dataLength - 1);
		return covariance / stdDev(x, meanX, dataLength) / stdDev(y, meanY, dataLength);
	}
	
	/**
	 * <p>Returns the correlation between the two arrays of data,
	 *    ignoring any Nan values</p>
	 * <p>See - <a href="http://en.wikipedia.org/wiki/Correlation">Wikipedia</a>
	 * </p>
	 * 
	 * @param x
	 * @param y
	 * @param dataLength - number of terms in each vector to consider (we look at the first dataLength terms).
	 * 	Precondition: dataLength is less than min(x.length, y.length)
	 * @return the covariance
	 */
	public static double correlationIgnoreNans(double[] x, double[] y, int dataLength) {
		// return covariance(x, y) / stdDev(x) / stdDev(y);
		// Save some code time by reusing the code from inside covariance:
		double c = 0;
		double meanX = 0;
		double meanY = 0;
		int count = 0;
		for (int i = 0; i < dataLength; i++) {
			if ((!Double.isNaN(x[i])) && (!Double.isNaN(y[i]))) {
				// Only add the values in if they are not NaN
				meanX += x[i];
				meanY += y[i];
				count++;
			}
		}
		// Adjust for the values we've skipped:
		meanX = meanX / count;
		meanY = meanY / count;
		
		for (int t = 0; t < dataLength; t++) {
			if ((!Double.isNaN(x[t])) && (!Double.isNaN(y[t]))) {
				// Only add the product in if it is not NaN
				c += (x[t] - meanX)*(y[t]-meanY);
			}
		}
		double covariance = c / (double) (count - 1);
		
		// Now work out the std devs of each:
		double sumSqsX = 0.0;
		double sumSqsY = 0.0;
		for (int m = 0; m < dataLength; m++) {
			if ((!Double.isNaN(x[m])) && (!Double.isNaN(y[m]))) {
				// Ignore if one is NaN
				sumSqsX += (x[m] - meanX) * (x[m] - meanX);
				sumSqsY += (y[m] - meanY) * (y[m] - meanY);
			}
		}
		double stdX = sumSqsX / (double) (count - 1);
		stdX = Math.sqrt(stdX);
		double stdY = sumSqsY / (double) (count - 1);
		stdY = Math.sqrt(stdY);
		
		return covariance / stdX / stdY;
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			matrix[r] = value;
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 * @param offset where in the array to start from
	 * @param length length in the array to fill
	 */
	public static void fill(int[] matrix, int value, int offset, int length) {
		for (int r = offset; r < offset + length; r++) {
			matrix[r] = value;
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[][] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = value;
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[][][] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					matrix[r][c][h] = value;
				}
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(int[][][][] matrix, int value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					int depth = matrix[r][c][h].length;
					for (int d = 0; d < depth; d++) {
						matrix[r][c][h][d] = value;
					}
				}
			}
		}
	}
	
	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			matrix[r] = value;
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[][] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				matrix[r][c] = value;
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[][][] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					matrix[r][c][h] = value;
				}
			}
		}
	}

	/**
	 * Initialises all values in the matrix to the given value
	 * 
	 * @param matrix
	 * @param value
	 */
	public static void fill(long[][][][] matrix, long value) {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			int cols = matrix[r].length;
			for (int c = 0; c < cols; c++) {
				int height = matrix[r][c].length;
				for (int h = 0; h < height; h++) {
					int depth = matrix[r][c][h].length;
					for (int d = 0; d < depth; d++) {
						matrix[r][c][h][d] = value;
					}
				}
			}
		}
	}
	
	public static double[][] transpose(double[][] matrix) {
		double[][] newMatrix = new double[matrix[0].length][matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				newMatrix[j][i] = matrix[i][j];
			}
		}
		return newMatrix;
	}
	
	public static int[][] transpose(int[][] matrix) {
		int[][] newMatrix = new int[matrix[0].length][matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				newMatrix[j][i] = matrix[i][j];
			}
		}
		return newMatrix;
	}

	/**
	 * Converts an int array to a double array
	 * 
	 * @param input
	 * @return
	 */
	public static double[][] convertMatrix(int[][] input) {
		double[][] outputArray = new double[input.length][];
		for (int i = 0; i < input.length; i++) {
			outputArray[i] = new double[input[i].length];
			for (int j = 0; j < input[i].length; j++) {
				outputArray[i][j] = input[i][j];
			}
		}
		return outputArray;
	}
	
	/**
	 * Converts a double array to an int array
	 * 
	 * @param input
	 * @param valueOffset value to be subtracted from each value
	 * @return
	 */
	public static int[] convertMatrix(double[] input, int valueOffset) {
		int[] outputArray = new int[input.length];
		for (int i = 0; i < input.length; i++) {
			outputArray[i] = (int) (input[i]) - valueOffset;
		}
		return outputArray;
	}

	/**
	 * Converts a double array to an int array
	 * 
	 * @param input
	 * @return
	 */
	public static int[] convertMatrix(double[] input) {
		return convertMatrix(input, 0);
	}

	/**
	 * <p>Returns the determinant of the input matrix.
	 * </p>
	 * 
	 * <p>This uses a fairly naive calculation - it will work for small sized
	 * matrices but will not be efficient enough for larger sizes.</p>
	 * 
	 * @param matrix
	 * @return determinant of matrix
	 * @throws Exception if supplied a non-square matrix
	 */
	public static double determinant(double[][] matrix) throws Exception {
		int rows = matrix.length;
		for (int r = 0; r < rows; r++) {
			if (matrix[r].length != rows) {
				throw new Exception("Cannot compute the determinant of a non-square matrix");
			}
		}
		return recursiveDeterminant(matrix);
	}
	
	/**
	 * <p>Private method to compute the determinant recursively.
	 * {@link determinant()} calls this after checking the matrix dimensions. <br/>
	 * @see {@link http://mathworld.wolfram.com/Determinant.html}
	 * </p>
	 * 
	 * @param matrix
	 * @return
	 */
	private static double recursiveDeterminant(double[][] matrix) {
		int rows = matrix.length;
		double result = 0;
		// Base cases:
		if (rows == 1) {
			return matrix[0][0];
		}
		if (rows == 2) {
			return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
		}
		// Recursive case
		int multiplier = 1;
		for(int col = 0; col < rows; col++) {
			// Construct the next sub-matrix to compute the determinant of
			double minor[][] = copyMatrixEliminateRowAndColumn(matrix, 0, col);
			result += (double) multiplier * matrix[0][col] * recursiveDeterminant(minor);
			multiplier *= -1;
		}

		return result; 
	}
	
	/*
	 * The method CholeskyDecomposition() was adapted from the 
	 *  JAMA project -- http://math.nist.gov/javanumerics/jama/
	 *  
	 * This code was distributed with the following original license:
	 * 
	 * Copyright Notice
	 * 
	 * This software is a cooperative product of The MathWorks and the
	 * National Institute of Standards and Technology (NIST) which has
	 * been released to the public domain. Neither The MathWorks nor
	 * NIST assumes any responsibility whatsoever for its use by other
	 * parties, and makes no guarantees, expressed or implied, about
	 * its quality, reliability, or any other characteristic.
	 * 
	 * As Jama is in the public domain other developers are free to
	 * adopt and adapt this code to other styles of programming or to
	 * extend or modernize the API.
	 * Make note, however, that NIST makes no endorsement of these projects.
	 */
	/**
	 * <p>Make the Cholesky decomposition L of a given input matrix A,
	 * where:
	 * 	<ol>
	 * 		<li>A is symmetric and positive definite (has full rank)</li>
	 * 		<li>A = L L^T (L^T is the transpose of L - here A has real
	 * 		entries only, though a Cholesky decomposition is possible
	 * 		with complex entries)</li>
	 * 		<li>L is a lower triangular matrix</li>
	 * 	</ol>
	 * We perform the decomposition using the Cholesky–Banachiewicz
	 *  algorithm, computing L from the top left, row by row (see wikipedia)
	 * </p>
	 * 
	 * <p>This method has been adapted from the JAMA project (public domain software)
	 * </p>
	 * 
	 * @param A input matrix
	 * @return L 
	 * @throws Exception when the matrix A is not square, is asymmetric, or
	 * 		not positive definite
	 * @see {@link en.wikipedia.org/wiki/Cholesky_decomposition}
	 * @see {@link http://mathworld.wolfram.com/CholeskyDecomposition.html}
	 * @see {@link http://en.wikipedia.org/wiki/Positive-definite_matrix}
	 * @see {@link http://math.nist.gov/javanumerics/jama/}
	 * @see {@link http://www2.gsu.edu/~mkteer/npdmatri.html}
	 */
	public static double[][] CholeskyDecomposition(double[][] A) throws Exception {
		int n = A.length;
		double[][] L = new double[n][n];
		// Loop over all rows:
		for (int j = 0; j < n; j++) {
			// Check length of row keeps this a square matrix:
			if (A[j].length != n) {
				throw new Exception("CholeskyDecomposition is only performed on square matrices");
			}
			double d = 0.0;
			for (int k = 0; k < j; k++) {
				double s = 0.0;
				for (int i = 0; i < k; i++) {
					s += L[k][i]*L[j][i];
				}
				L[j][k] = s = (A[j][k] - s)/L[k][k];
				d = d + s*s;
				// Check that these matrix entries remain symmetric:
				if (A[k][j] != A[j][k]) {
					throw new Exception("CholeskyDecomposition is only performed on symmetric matrices");
				}
			}
			d = A[j][j] - d;
			// Check the positive definite condition:
			if (d <= 0.0) {
				// Throw an error with some suggestions. The last suggestion is from my observations
				//  from a simple test with Matlab - I should find a reference for this ...
				throw new NonPositiveDefiniteMatrixException("CholeskyDecomposition is only performed on positive-definite matrices. " + 
						"Some reasons for non-positive-definite matrix are listed at http://www2.gsu.edu/~mkteer/npdmatri.html - " +
						"note: a correlation matrix is non-positive-definite if you have more variables than observations");
			}
			L[j][j] = Math.sqrt(d);
			// Set the upper triangular part to all zeros:
			for (int k = j+1; k < n; k++) {
				L[j][k] = 0.0;
			}
		}
		return L;
	}
	
	/**
	 * Compute matrix inversion of a symmetric, positive definite matrix
	 *  by using the Cholesky Decomposition L of the matrix A.
	 * Since A = L L^T, then A^-1 = (L^T)^-1 L^-1, and the inverses
	 *  of 
	 * 
	 * @param A matrix to be inverted
	 * @return the inverse of A
	 * @throws Exception when the matrix is not symmetric or positive definite
	 * @see #CholeskyDecomposition(double[][])
	 */
	public static double[][] invertSymmPosDefMatrix(double[][] A) throws Exception {
		// First do the Cholesky Decomposition:
		
		double[][] L = CholeskyDecomposition(A);
		
		return solveViaCholeskyResult(L, identityMatrix(A.length));
	}
	// TODO implement solve for identity matrix
	
	/*
	 * The method solveViaCholeskyResult() was adapted from the 
	 *  JAMA project -- http://math.nist.gov/javanumerics/jama/
	 *  
	 * See license above the CholeskyDecomposition() method
	 */
	/**
	 * <p>Solve A*X = B, where A = L*L^T via Cholesky decomposition.
	 * </p>
	 * 
	 * <p>This method has been adapted from the JAMA project (public domain software)
	 * </p>
	 *
	 * @param L Cholesky decomposition of the matrix A
	 * @param B matrix with as many rows as A and any number of columns
	 * @return X so that A*X = B
	 * @see {@link http://math.nist.gov/javanumerics/jama/}
	 * @see #CholeskyDecomposition(double[][])
	 */
	public static double[][] solveViaCholeskyResult(double[][] L, double[][] B) {
		int aRows = L.length;
		if (aRows != B.length) {
			throw new IllegalArgumentException("Matrix row dimensions must agree.");
		}

		// Copy B matrix 
		double[][] X = MatrixUtils.arrayCopy(B);
		int bCols = B[0].length;
		
		// Solve L*Y = B;
		for (int k = 0; k < aRows; k++) {
			for (int j = 0; j < bCols; j++) {
				for (int i = 0; i < k ; i++) {
					X[k][j] -= X[i][j]*L[k][i];
				}
				X[k][j] /= L[k][k];
			}
		}

		// Solve L'*X = Y;
		for (int k = aRows-1; k >= 0; k--) {
			for (int j = 0; j < bCols; j++) {
				for (int i = k+1; i < aRows; i++) {
					X[k][j] -= X[i][j]*L[i][k];
				}
				X[k][j] /= L[k][k];
			}
		}
		return X;
	}
	
	/**
	 * <p>Compute determinant(A), where A = L*L^T via Cholesky decomposition
	 * is a symmetric, positive definite matrix.
	 * </p>
	 * 
	 * <p>This is an efficient computation, since det(A) = det(L*L^T) = det(L)*det(L^T)
	 * and L is triangular so det(L) is just the product along the diagonal
	 * </p>
	 * 
	 * @param A symmetric, positive definite matrix to take the determinant of
	 * @return the determinant of A
	 * @throws Exception when the matrix is not symmetric or positive definite
	 */
	public static double determinantSymmPosDefMatrix(double[][] A) throws Exception {
		double[][] L = CholeskyDecomposition(A);
		
		return determinantViaCholeskyResult(L);
	}

	/**
	 * <p>Compute determinant(A), where A = L*L^T via Cholesky decomposition
	 * is a symmetric, positive definite matrix.
	 * </p>
	 * 
	 * <p>This is an efficient computation, since det(A) = det(L*L^T) = det(L)*det(L^T)
	 * and L is triangular so det(L) is just the product along the diagonal
	 * </p>
	 * 
	 * @param L Cholesky decomposition L of the symmetric positive definite matrix A
	 * @return det(A)
	 * @see #CholeskyDecomposition(double[][])
	 */
	public static double determinantViaCholeskyResult(double[][] L) {
		double detL = 1.0;
		int n = L.length;
		for (int i = 0; i < n; i++) {
			detL *= L[i][i];
		}
		return detL * detL;
	}

	public static void printMatrix(PrintStream out, double[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				out.print(matrix[r][c] + " ");
			}
			out.println();
		}
	}
	
	public static void printMatrix(PrintStream out, int[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				out.print(matrix[r][c] + " ");
			}
			out.println();
		}
	}

	public static void printArray(PrintStream out, double[] array) {
		for (int r = 0; r < array.length; r++) {
				out.print(array[r] + " ");
		}
		out.println();
	}

	public static void printArray(PrintStream out, int[] array) {
		for (int r = 0; r < array.length; r++) {
				out.print(array[r] + " ");
		}
		out.println();
	}
	
	/**
	 * Discretizes using even bin sizes
	 * 
	 * @param data
	 * @param numBins
	 * @return
	 */
	public static int[] discretise(double data[], int numBins) {
		int[] discretised = new int[data.length];
		double min = min(data);
		double max = max(data);
		double binInterval = (max - min) / numBins;
		
		for (int t = 0; t < data.length; t++) {
			discretised[t] = (int) ((data[t] - min) / binInterval);
			if (discretised[t] == numBins) {
				// This occurs for the maximum value; put it in the largest bin (base - 1)
				discretised[t]--;
			}
		}
		return discretised;
	}
	
	/**
	 * Discretizes using a maximum entropy partitioning
	 * 
	 * @param data
	 * @param numBins
	 * @return
	 */
	public static int[] discretiseMaxEntropy(double data[], int numBins){
		int[] newData = new int[data.length];
		
		double[] tempData =  new double[data.length];
		System.arraycopy(data, 0, tempData, 0, data.length);
		Arrays.sort(tempData);
		int compartmentSize;
		double[] cutOffValues = new double[numBins]; 
		for(int i=0;i<numBins;i++){
			compartmentSize = (int)((double)(i+1)*(double)(data.length)/(double)numBins)-1;
			// System.out.println(compartmentSize);
			cutOffValues[i]=tempData[compartmentSize];
		}
		
		for (int i=0;i<data.length;i++){				
			for(int m=0;m<numBins;m++){
				if (data[i] <= cutOffValues[m]){
					newData[i] = m;
					break;
				}
			}
		}
		return newData;
	}

	/**
	 * Discretizes each column of the data independently,
	 * using a maximum entropy partitioning
	 * 
	 * @param data
	 * @param numBins
	 * @return
	 */
	public static int[][] discretiseMaxEntropy(double data[][], int numBins){
		int lastCol = data[0].length;	
		int lastRow = data.length;
		int[][] newData = new int[lastRow][lastCol];
		for(int j=0;j<lastCol;j++){
			double[] tempData = new double[lastRow];
			for (int i=0;i<lastRow;i++){
				tempData[i] = data[i][j]; 			 
			}
		
			Arrays.sort(tempData);

			int compartmentSize;
			double[] cutOffValues = new double[numBins]; 
			for(int i=0;i<numBins;i++){
				compartmentSize = (int)((double)(i+1)*(double)(lastRow)/(double)numBins)-1;
//				System.out.println(compartmentSize);
				cutOffValues[i]=tempData[compartmentSize];
			}
			
			for (int i=0;i<lastRow;i++){				
				for(int m=0;m<numBins;m++){
					if (data[i][j] <= cutOffValues[m]){
						newData[i][j] = m;
						m = numBins;
					}
				}
			}
		}
		return newData;
	}
	
	/**
	 * Take the logical AND of all variables in each row 
	 * 
	 * @param data
	 * @return
	 */
	public static boolean[] andRows(boolean[][] data) {
		boolean[] result = new boolean[data.length];
		
		for (int i = 0; i < data.length; i++) {
			result[i] = true;
			for (int j = 0; j < data[i].length; j++) {
				result[i] &= data[i][j];
			}
		}
		return result;
	}
	
	/**
	 * Take the logical AND of selected variables in each row 
	 * 
	 * @param data
	 * @param columns which variables to take the AND over
	 * @return
	 */
	public static boolean[] andRowsOverSelectedColumns(boolean[][] data, int[] columns) {
		boolean[] result = new boolean[data.length];
		
		for (int i = 0; i < data.length; i++) {
			result[i] = true;
			for (int c = 0; c < columns.length; c++) {
				result[i] &= data[i][columns[c]];
			}
		}
		return result;
	}
	
	/**
	 * Convert a single dimensional double array to a 2D array,
	 *  where the first dimension matches that of the original array
	 *  and second is always 0 (with our convention for first index being
	 *  time and second variable number, this means we have a 2D array
	 *  where the original array becomes the first variable).
	 * 
	 * @param array
	 * @return 2D array with the original array as the only column.
	 */
	public static double[][] doubleTo2DArray(double[] array) {
		double[][] twoDArray = new double[array.length][1];
		for (int i = 0; i < array.length; i++) {
			twoDArray[i][0] = array[i];
		}
		return twoDArray;
	}
	
	/**
	 * Convert a double array to an int array.
	 * This is designed specifically for use of the toolkit in Octave
	 *  where all native arrays are considered as doubles for Java.
	 * To use an integer 1D array (to supply to a java method), one must
	 *  first create the native Octave 1D array, then create a java 1D Double
	 *  array, then use this method to convert that java 1D Double array
	 *  to a java 1D integer array.
	 * 
	 * @param array
	 * @return
	 */
	public static int[] doubleToIntArray(double[] array) {
		if (array == null) {
			return null;
		}
		int[] intArray = new int[array.length];
		for (int i = 0; i < array.length; i++) {
			intArray[i] = (int) array[i];
		}
		return intArray;
	}

	/**
	 * Convert a 2D double array to an int array.
	 * This is designed specifically for use of the toolkit in Octave
	 *  where all native arrays are considered as doubles for Java.
	 * To use an integer 2D array (to supply to a java method), one must
	 *  first create the native Octave 2D array, then create a java 2D Double
	 *  array, then use this method to convert that java 2D Double array
	 *  to a java 2D integer array.
	 * 
	 * @param array
	 * @return
	 */
	public static int[][] doubleToIntArray(double[][] array) {
		if (array == null) {
			return null;
		}
		int[][] intArray = new int[array.length][];
		for (int i = 0; i < array.length; i++) {
			if (array[i] == null) {
				intArray[i] = null;
			} else {
				intArray[i] = new int[array[i].length];
				for (int j = 0; j < array[i].length; j++) {
					intArray[i][j] = (int) array[i][j];
				}
			}
		}
		return intArray;
	}
}
