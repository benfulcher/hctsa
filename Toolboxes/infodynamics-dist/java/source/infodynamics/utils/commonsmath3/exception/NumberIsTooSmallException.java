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

/*
 * This class was originally distributed as part of the Apache Commons
 *  Math3 library, under the Apache License Version 2.0, which is 
 *  copied below. This Apache 2 software is now included as a derivative
 *  work in the GPLv3 licensed JIDT project, as per:
 *  http://www.apache.org/licenses/GPL-compatibility.html
 *  
 * The original Apache source code has been modified as follows:
 * -- We have modified package names to sit inside the JIDT structure.
 */

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package infodynamics.utils.commonsmath3.exception;

import infodynamics.utils.commonsmath3.exception.util.Localizable;
import infodynamics.utils.commonsmath3.exception.util.LocalizedFormats;

/**
 * Exception to be thrown when a number is too small.
 *
 * @since 2.2
 * @version $Id$
 */
public class NumberIsTooSmallException extends MathIllegalNumberException {
    /** Serializable version Id. */
    private static final long serialVersionUID = -6100997100383932834L;
    /**
     * Higher bound.
     */
    private final Number min;
    /**
     * Whether the maximum is included in the allowed range.
     */
    private final boolean boundIsAllowed;

    /**
     * Construct the exception.
     *
     * @param wrong Value that is smaller than the minimum.
     * @param min Minimum.
     * @param boundIsAllowed Whether {@code min} is included in the allowed range.
     */
    public NumberIsTooSmallException(Number wrong,
                                     Number min,
                                     boolean boundIsAllowed) {
        this(boundIsAllowed ?
             LocalizedFormats.NUMBER_TOO_SMALL :
             LocalizedFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED,
             wrong, min, boundIsAllowed);
    }

    /**
     * Construct the exception with a specific context.
     *
     * @param specific Specific context pattern.
     * @param wrong Value that is smaller than the minimum.
     * @param min Minimum.
     * @param boundIsAllowed Whether {@code min} is included in the allowed range.
     */
    public NumberIsTooSmallException(Localizable specific,
                                     Number wrong,
                                     Number min,
                                     boolean boundIsAllowed) {
        super(specific, wrong, min);

        this.min = min;
        this.boundIsAllowed = boundIsAllowed;
    }

    /**
     * @return {@code true} if the minimum is included in the allowed range.
     */
    public boolean getBoundIsAllowed() {
        return boundIsAllowed;
    }

    /**
     * @return the minimum.
     */
    public Number getMin() {
        return min;
    }
}
