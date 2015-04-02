;
;  Java Information Dynamics Toolkit (JIDT)
;  Copyright (C) 2012, Joseph T. Lizier
;  
;  This program is free software: you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation, either version 3 of the License, or
;  (at your option) any later version.
;  
;  This program is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;  
;  You should have received a copy of the GNU General Public License
;  along with this program.  If not, see <http://www.gnu.org/licenses/>.
;

; = Example 3 - Transfer entropy on continuous data using kernel estimators =

; Simple transfer entropy (TE) calculation on continuous-valued data using the (box) kernel-estimator TE calculator.

; Import relevant classes:
(import infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel)
(import java.util.Random)
(def rg (Random.))

(let 
    [numObservations 1000
     covariance 0.4
     ; Generate some random normalised data.
     sourceArray (double-array (take numObservations (repeatedly #(.nextGaussian rg))))
     destArray (double-array 
	(cons 0 
	    (map + 
		(map (partial * covariance) (butlast sourceArray)) 
		(map (partial * (- covariance 1)) (double-array (take (- numObservations 1) (repeatedly #(.nextGaussian rg))))) )))
     sourceArray2 (double-array (take numObservations (repeatedly #(.nextGaussian rg))))
     teCalc (TransferEntropyCalculatorKernel. )
    ]

; Set up the calculator
(.setProperty teCalc "NORMALISE" "true")
(.initialise teCalc 1 0.5) ; Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units

(.setObservations teCalc sourceArray destArray)
; For copied source, should give something close to expected value for correlated Gaussians:
(println "TE result " (.computeAverageLocalOfObservations teCalc)
	" expected to be close to " (/ (Math/log (/ 1 (- 1 (* covariance covariance)))) (Math/log 2))
	" for these correlated Gaussians but biased upward")

(.initialise teCalc ) ; Initialise leaving the parameters the same
(.setObservations teCalc sourceArray2 destArray)
; For random source, it should give something close to 0 bits
(println "TE result " (.computeAverageLocalOfObservations teCalc)
	" expected to be close to 0 bits for these uncorrelated Gaussians but will be biased upward")

; We can get insight into the bias by examining the null distribution:
(def nullDist (.computeSignificance teCalc 100))
(println "Null distribution for unrelated source and destination "
	"(i.e. the bias) has mean " (.getMeanOfDistribution nullDist)
	" and standard deviation " (.getStdOfDistribution nullDist))

)

