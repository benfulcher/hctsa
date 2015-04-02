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

; = Example 2 - Transfer entropy on multidimensional binary data =

; Simple transfer entropy (TE) calculation on multidimensional binary data using the discrete TE calculator.

; Import relevant classes:
(import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete)

(let
    ; Create many columns in a multidimensional array (2 rows by 100 columns),
    ;  where the next time step (row 2) copies the value of the column on the left
    ;  from the previous time step (row 1):
    [row1 (int-array (take 100 (repeatedly #(rand-int 2))))
     row2 (int-array (cons (aget row1 99) (butlast row1))) ; shifts row1 by 1
     twoDTimeSeriesClojure (into-array (map int-array [row1 row2]))
     ; Create TE calculator
     teCalc (TransferEntropyCalculatorDiscrete. 2 1)
    ]


; Initialise the TE calculator and run it:
(.initialise teCalc)
; Add observations of transfer across one cell to the right per time step:
(.addObservations teCalc twoDTimeSeriesClojure 1)
(println "The result should be close to 1 bit here, since we are executing copy operations of what is effectively a random bit to each cell here:"
	(.computeAverageLocalOfObservations teCalc))

)
