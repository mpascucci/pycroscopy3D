
/*
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
 *  
 *   Copyright (C) 2018-2020
 *    Laboratory of Systems Biology, Department of Cybernetics,
 *    School of Science, Tallinn University of Technology
 *   This file is part of project: IOCBIO Deconvolve
 */


#ifndef IOCBIO_CONSTANTS_HPP
#define IOCBIO_CONSTANTS_HPP

#define NDIMS 3 // limited to 3D case

#define IOCBIO__S(x) #x
#define IOCBIO__S_(x) IOCBIO__S(x)
#define IOCBIO__SLINE__ IOCBIO__S_(__LINE__)

#define EXCPT_START "DeconvolveLib [" __FILE__ ":" IOCBIO__SLINE__ "]"
#define EXCPT_INTERNAL EXCPT_START " InternalError: "
#define EXCPT_MEMORY EXCPT_START ": "
#define EXCPT_NOBODYS_FAULT EXCPT_START ": "
#define EXCPT_USER EXCPT_START " UserError: "

#endif
