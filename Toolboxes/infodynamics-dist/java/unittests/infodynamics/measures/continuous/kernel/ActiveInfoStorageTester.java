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

import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;
import junit.framework.TestCase;

public class ActiveInfoStorageTester extends TestCase {

	public void testEmbedding() throws Exception {
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();

		int[] ks = {1,2,3,4};
		int[] taus = {1,2,3,4};
		
		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int tauIndex = 0; tauIndex < taus.length; tauIndex++) {
				int tau = taus[tauIndex];
		
				ActiveInfoStorageCalculatorKernel ais = new ActiveInfoStorageCalculatorKernel();
				ais.initialise(k, tau, 0.5);
				ais.setObservations(MatrixUtils.selectColumn(data, 0));

				assertEquals(data.length -  (k-1)*tau - 1,
						ais.getNumObservations());
			}
		}
	}
	
	
	public void testAgainstDirectAISCalculator() throws Exception {
		
		int[] ks = {1,2,3,4};
		double[] kernelwidths = {0.2, 0.5, 1};
		
		ArrayFileReader afr = new ArrayFileReader("demos/data/2coupledRandomCols-1.txt");
		double[][] data = afr.getDouble2DMatrix();

		for (int kIndex = 0; kIndex < ks.length; kIndex++) {
			int k = ks[kIndex];
			for (int kwIndex = 0; kwIndex < kernelwidths.length; kwIndex++) {
				double kernelwidth = kernelwidths[kwIndex];
				
				ActiveInfoStorageCalculatorKernelDirect aisOld = new ActiveInfoStorageCalculatorKernelDirect();
				aisOld.initialise(k, kernelwidth);
				aisOld.setObservations(MatrixUtils.selectColumn(data, 0));
				double aisOldValue = aisOld.computeAverageLocalOfObservations();
				
				ActiveInfoStorageCalculatorKernel aisNew = new ActiveInfoStorageCalculatorKernel();
				aisNew.initialise(k, kernelwidth);
				aisNew.setObservations(MatrixUtils.selectColumn(data, 0));
				double aisNewValue = aisNew.computeAverageLocalOfObservations();
				
				
				System.out.printf("k=%d, kw=%.2f, Kernel calculator: %.5f bits; Direct kernel calculator: %.5f bits\n",
						k, kernelwidth, aisNewValue, aisOldValue);
				assertEquals(aisOldValue, aisNewValue, 0.00000001);
			}
		}
		
	}
}
