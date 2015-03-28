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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Writes basic files with one or two dimensional arrays of integers
 *  and doubles and booleans,
 *  delimited by tab.
 *  
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class ArrayFileWriter {
	
	/**
	 * Outputs boolean[][] to a file
	 * 
	 * @param matrixthe data to write
	 * @param outputFileName filename to write to
	 */
	public static void makeBooleanMatrixFile(boolean matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;
		
		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				if (matrix[i][j]==false){
					out.write("0\t");
				}else{
					out.write("1\t");
				}
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}
	
	/**
	 * Outputs double[][] to a file
	 * 
	 * @param matrixthe data to write
	 * @param outputFileName filename to write to
	 */
	public static void makeDoubleMatrixFile(double matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.valueOf(matrix[i][j])+"\t");
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	/**
	 * Outputs double[][] to a file
	 * 
	 * @param matrix the data to write
	 * @param outputFileName filename to write to
	 * @param decimalPlaces number of decimal places to write
	 */
	public static void makeDoubleMatrixFile(double matrix[][],String outputFileName, int decimalPlaces) throws IOException {
		String template = String.format("%%.%df\t", decimalPlaces);
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.format(template, matrix[i][j]));
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	/**
	 * Write the matrix to a file.
	 * Assumes the 2D array is regular - i.e. all rows have the same number of columns as the first
	 * 
	 * @param matrix data to write
	 * @param outputFileName output filename
	 * @throws IOException
	 */
	public static void makeIntMatrixFile(int matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.valueOf(matrix[i][j])+"\t");
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	/**
	 * Outputs double[] to a file
	 * 
	 * @param matrix the data to write
	 * @param outputFileName filename to write to
	 */
	public static void makeMatrixFile(double matrix[],String outputFileName) throws IOException {
    	int rowSize = matrix.length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			out.write(String.valueOf(matrix[i]) + "\n");
		}
		out.close();
	}

	/**
	 * Outputs Number[][] to a file
	 * 
	 * @param matrix the data to write
	 * @param outputFileName filename to write to
	 */
	public static void makeMatrixFile(Number matrix[][],String outputFileName) throws IOException {
    	int rowSize = matrix.length;
		int colSize = matrix[0].length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			for(int j=0;j<colSize;j++){
				out.write(String.valueOf(matrix[i][j])+"\t");
				if (j==colSize-1){
					out.write("\n");
				}
			}
		}
		out.close();
	}

	/**
	 * Outputs Number[] to a file
	 * 
	 * @param matrix the data to write
	 * @param outputFileName filename to write to
	 */
	public static void makeMatrixFile(Number matrix[],String outputFileName) throws IOException {
    	int rowSize = matrix.length;

		createDirectories(outputFileName);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
		for(int i=0;i<rowSize;i++){
			out.write(String.valueOf(matrix[i]) + "\n");
		}
		out.close();
	}
	
	private static void createDirectories(String filename) {
		File file = new File(filename);
		File parentDir = file.getParentFile();
		if ((parentDir != null) && !parentDir.isDirectory()) {
			parentDir.mkdirs();
		}
	}
}
