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

;; = Example 1 - Transfer entropy on binary data =

; Simple transfer entropy (TE) calculation on binary data using the discrete TE calculator:

; Import relevant classes:
(import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete)

(let
    ; Generate some random binary data.
    [sourceArray (int-array (take 100 (repeatedly #(rand-int 2))))
     destArray (int-array (cons 0 (butlast sourceArray))) ; shifts sourceArray by 1
     sourceArray2 (int-array (take 100 (repeatedly #(rand-int 2))))
     ; Create TE calculator
     teCalc (TransferEntropyCalculatorDiscrete. 2 1)
    ]

; Initialise the TE calculator and run it:
(.initialise teCalc)
(.addObservations teCalc sourceArray destArray)
(println "For copied source, result should be close to 1 bit : "
	(.computeAverageLocalOfObservations teCalc))

(.initialise teCalc)
(.addObservations teCalc sourceArray2 destArray)
(println "For random source, result should be close to 0 bits : "
	(.computeAverageLocalOfObservations teCalc))

)
