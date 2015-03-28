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
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.net.InetAddress;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;

/**
 * Octave text file format writer.
 * Usage:
 * <ol>
 *  <li>call constructor</li>
 *  <li>call put for each variable to be stored.</li>
 *  <li>call writeFile(outputFilename) or setFilename(outputFilename)
 *    then writeFile()</li>
 * 
 * @author Joseph Lizier (<a href="joseph.lizier at gmail.com">email</a>,
 * <a href="http://lizier.me/joseph/">www</a>)
 */
public class OctaveFileWriter extends HashMap<String, Object> {

	public static final long serialVersionUID = 1;
	
	private boolean writeLFOnly = true;
	private String octaveFilename = null;
	private static final String LINE_SEPARATOR_PROPERTY = "line.separator";
	private static final String LINE_SEPARATOR = "\n";
	private Hashtable<String, Integer> roundingHT = new Hashtable<String, Integer>();
		
	public OctaveFileWriter() {
	}
	
	public void setFilename(String outputFilename) {
		octaveFilename = outputFilename;
	}
	
	public void writeFile(String outputFilename) throws IOException {
		setFilename(outputFilename);
		writeFile();
	}

	/**
	 * Writes the contained variables to a file formatted for octave to read
	 * 
	 * @param octaveFilename
	 */
	public void writeFile() throws IOException {
		int maxLength = 0; // required for arrays of strings
		
		if (octaveFilename == null) {
			throw new IOException("No filename has been set");
		}
		
		String originalLSValue = null;
		if (writeLFOnly) {
			// We only want to write \n. Save the current line separator for later restoration
			originalLSValue = System.getProperty(LINE_SEPARATOR_PROPERTY);
			System.setProperty(LINE_SEPARATOR_PROPERTY, LINE_SEPARATOR);
		}
		
		// Create the directory if required
    	createDirectories(octaveFilename);
		
		// Open the file for writing
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(octaveFilename)));
		
		// Write header
		String hostname = "null";
		try {
			InetAddress localMachine = InetAddress.getLocalHost();
			hostname = localMachine.getHostName();
		} catch (Exception e) {
			// do nothing
		}
		pw.println("# Created on " + hostname + " by infodynamics.utils.OctaveFileWriter, " + (new Date()));
		
		// Have a decimal format object ready in case we have to do any rounding
		DecimalFormat decFormat = new DecimalFormat();
		int rounding = -1;
		
		Set keySet = this.keySet();
		Iterator iterator = keySet.iterator();
		for (; iterator.hasNext();) {
			String key = (String) iterator.next();
			Object value = this.get(key);
			try {
				// Check if we need to round off this value
				Integer roundObject = roundingHT.get(key);
				if (roundObject != null) {
					rounding = roundObject.intValue();
					decFormat.setMaximumFractionDigits(rounding);
				} else {
					rounding = -1;
				}
				// Start header for this variable
				pw.println("# name: " + key);
				// Check what the type of the object is
				if (value.getClass().isArray()) {
					// we have an array item
					// is it a 1D or 2D matrix? Check if the first item is an array
					if (Array.getLength(value) > 0) {
						Object item1 = Array.get(value, 0);
						// Check whether this is an array of strings
						boolean isStrings = String.class.isInstance(item1); 
						if (isStrings) {
							pw.println("# type: string");
							// Now need to find the maximum length of all strings in the array
							maxLength = 0;
							for (int i = 0; i < Array.getLength(value); i++) {
								if (((String)Array.get(value, i)).length() > maxLength) {
									maxLength = ((String)Array.get(value, i)).length();
								}
							}
							pw.println("# elements: " + Array.getLength(value));
						} else {
							// Ordinary array
							pw.println("# type: matrix");
						}
						if (item1.getClass().isArray()) {
							// We have a 2D array - (cannot have 2D arrays of strings in Octave,
							//  so assume they haven't been passed in)
							pw.println("# rows: " + Array.getLength(value));
							pw.println("# columns: " + Array.getLength(item1));
							for (int i = 0; i < Array.getLength(value); i++) {
								// Grab the next row
								Object row = Array.get(value, i);
								// If the first element is a float or double we'll
								//  round off every element.
								for (int j = 0; j < Array.getLength(row); j++) {
									// Print each index item
									Object thisValue = Array.get(row, j);
									if (Boolean.class.isInstance(thisValue)) {
										pw.print(" " + (((Boolean)thisValue) ? "1" : "0"));
									} else {
										if (rounding >= 0) {
											// print with rounding
											pw.print(" " + decFormat.format(thisValue));
										} else {
											pw.print(" " + thisValue);
										}
									}
								}
								// Close off the row
								pw.println();
							}
						} else {
							// We have a 1D array
							if (isStrings) {
								for (int i = 0; i < Array.getLength(value); i++) {
									pw.println("# length: " + maxLength);
									String thisString = (String) Array.get(value, i); 
									pw.print(thisString);
									for (int l = thisString.length(); l < maxLength; l++) {
										pw.print(' ');
									}
									pw.println();
								}
							} else {
								// ordinary array
								pw.println("# rows: " + Array.getLength(value));
								pw.println("# columns: 1");
								for (int i = 0; i < Array.getLength(value); i++) {
									// Print each index item on it's own row
									Object thisValue = Array.get(value, i);
									if (Boolean.class.isInstance(thisValue)) {
										pw.print(" " + (((Boolean)thisValue) ? "1" : "0"));
									} else {
										if (rounding >= 0) {
											// print with rounding
											pw.print(" " + decFormat.format(thisValue));
										} else {
											pw.println(" " + thisValue);
										}
									}
								}
							}
						}
					} else {
						// Empty array ...
						// just write null values. Not 100% sure this is valid ...
						pw.println("# type: matrix");
						pw.println("# rows: 0");
						pw.println("# columns: 0");
					}
				} else {
					if (String.class.isInstance(value)) {
						pw.println("# type: string");
						pw.println("# elements: 1");
						pw.println("# length: " + ((String)value).length());
						pw.println(value);
					} else if (Boolean.class.isInstance(value)) {
						// we have a boolean value
						pw.println("# type: bool");
						pw.println(((Boolean)value) ? "1" : "0");
					} else {
						// we have a general scalar value
						pw.println("# type: scalar");
						if (rounding >= 0) {
							pw.println(decFormat.format(value));
						} else {
							pw.println(value);
						}
					}
				}
			} catch (Exception e) {
				System.out.println("Problem writing variable " +
						key + " to the output file (value = " +
						value + "):");
				e.printStackTrace();
				System.out.println("Continuing with file.");
			}
		}
		// Add one extra newline at the end
		pw.println();
		pw.close();
		
		if (writeLFOnly) {
			// Restore the saved line separator
			System.setProperty(LINE_SEPARATOR_PROPERTY, originalLSValue);
		}
	}
	
	/**
	 * JL: I have no idea why i put this method here. Nothing appears
	 *  to be using it
	 * 
	 * @param key
	 * @param value
	 * @param roundDPs
	 */
	public void put(String key, Object value, int roundDPs) {
		roundingHT.put(key, new Integer(roundDPs));
		put(key, value);
	}
	
	/** 
	 * Create the relevant directories if required
	 * 
	 * @param filename
	 */
	private static void createDirectories(String filename) {
		File file = new File(filename);
		File parentDir = file.getParentFile();
		if ((parentDir != null) && !parentDir.isDirectory()) {
			parentDir.mkdirs();
		}
	}
	
	public void putAll(Properties props) {
		for (Object keyObject : props.keySet()) {
			String key = (String) keyObject;
			// Alter the key name to remvoe "." characters
			key = key.replaceAll("\\.", "__");
			// Add this pair in
			put(key, props.get(keyObject));
		}
	}
}
