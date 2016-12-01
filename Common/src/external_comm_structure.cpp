/*!
 * \file geometry_structure.cpp
 * \brief Main subroutines for creating the primal grid and multigrid structure.
 * \author A. Jirasek
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/internal_comm_structure.hpp"

CExtIntrf::CExtIntrf(void) {

  BCnames = NULL;    
  Type = NULL;
  Iname = NULL; 
  PN = 0;
  IP = NULL;   
  Displacement = NULL;
  PartNr = NULL;
  communicator = -1;
  
}

CExtIntrf::~CGeometry(void) {
  
  unsigned long iBC;
  
  if (BCnames != NULL) {
    for (iBC = 0; iBC < nCB; iBC++)
      if (BCnames[iBC] != NULL) delete BCnames[iBC];
    delete[] BCnames;
  }
  
  if (Type != NULL) {
    delete[] Type;
  }

  if (Name != NULL) {
    delete[] Name;
  }

  if (IP != NULL) {
    delete[] IP;
  }
  
  if (Displacement != NULL) {
    delete[] Displacement;
  }
  
  if (PartNr != NULL) {
    delete[] PartNr;
  }
  
}

su2double CGeometry::Point2Plane_Distance(su2double *Coord, su2double *iCoord, su2double *jCoord, su2double *kCoord) {
  su2double CrossProduct[3], iVector[3], jVector[3], distance, modulus;
  unsigned short iDim;
  
  for (iDim = 0; iDim < 3; iDim ++) {
    iVector[iDim] = jCoord[iDim] - iCoord[iDim];
    jVector[iDim] = kCoord[iDim] - iCoord[iDim];
  }
  
  CrossProduct[0] = iVector[1]*jVector[2] - iVector[2]*jVector[1];
  CrossProduct[1] = iVector[2]*jVector[0] - iVector[0]*jVector[2];
  CrossProduct[2] = iVector[0]*jVector[1] - iVector[1]*jVector[0];
  
  modulus = sqrt(CrossProduct[0]*CrossProduct[0]+CrossProduct[1]*CrossProduct[1]+CrossProduct[2]*CrossProduct[2]);
  
  distance = 0.0;
  for (iDim = 0; iDim < 3; iDim ++)
    distance += CrossProduct[iDim]*(Coord[iDim]-iCoord[iDim]);
  distance /= modulus;
  
  return distance;
  
}
  
}
