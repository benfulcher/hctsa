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

; = Example 4 - Transfer entropy on continuous data using Kraskov estimators =

; Simple transfer entropy (TE) calculation on continuous-valued data using the Kraskov-estimator TE calculator.

; Import relevant classes:
(import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov)
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
     teCalc (TransferEntropyCalculatorKraskov. )
    ]

; Set up the calculator
(.setProperty teCalc "k" "4") ; Use Kraskov parameter K=4 for 4 nearest points
(.initialise teCalc 1) ; Use history length 1 (Schreiber k=1)

; Perform calculation with correlated source:
(.setObservations teCalc sourceArray destArray)
; Note that the calculation is a random variable (because the generated
;  data is a set of random variables) - the result will be of the order
;  of what we expect, but not exactly equal to it; in fact, there will
;  be a large variance around it.
(println "TE result " (.computeAverageLocalOfObservations teCalc)
	" nats expected to be close to " (Math/log (/ 1 (- 1 (* covariance covariance))))
	" nats for these correlated Gaussians")

; Perform calculation with uncorrelated source:
(.initialise teCalc ) ; Initialise leaving the parameters the same
(.setObservations teCalc sourceArray2 destArray)
; For random source, it should give something close to 0 bits
(println "TE result " (.computeAverageLocalOfObservations teCalc)
	" nats expected to be close to 0 nats for these uncorrelated Gaussians")

; We can also compute the local TE values for the time-series samples here:
;  (See more about utility of local TE in the CA demos)
(def localTE (.computeLocalOfPreviousObservations teCalc))

(println "Notice that the mean of locals, "
	(/ (reduce + localTE) (- numObservations 1))
	" nats, equals the previous result")

)
