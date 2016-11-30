/*!
 * \file iteration_structure.cpp
 * \brief Main subroutine managing communication to external solvers
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

#include "../include/iteration_structure.hpp"

#include "../include/ext_man_header.hpp"


void External_Comm(CGeometry ***geometry_container, CSurfaceMovement **surface_movement,
         CVolumetricMovement **grid_movement, CFreeFormDefBox ***FFDBox,CSolver ****solver_container,
	 CConfig **config_container, unsigned short val_iZone, unsigned long IntIter,
	 unsigned long ExtIter)   {

  

  unsigned short iDim, iMGlevel, nMGlevels = config_container[val_iZone]->GetnMGLevels();
  unsigned short Kind_Grid_Movement = config_container[val_iZone]->GetKind_GridMovement(val_iZone);
  unsigned long nIterMesh;
  unsigned long iPoint;

  d6dof_t d6dofdata, d6dofdata_old, *p_6DOFdata, *p_6DOFdata_old;
  conn_t conn, *pconn;
  struct timespec now, tmstart;
  double seconds;

  bool stat_mesh = true;
  bool adjoint = config_container[val_iZone]->GetContinuous_Adjoint();
  bool harmonic_balance = (config_container[val_iZone]->GetUnsteady_Simulation() == HARMONIC_BALANCE);

  p_6DOFdata = &d6dofdata;
  p_6DOFdata_old = &d6dofdata_old;
  pconn = &conn;

        cout << " Outside:   --------------------------  " << endl;
 
}






void CEulerSolver::EC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                         CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, iVar, iMarker;
  
  bool implicit 	 = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
  su2double Density, sq_vel, StaticEnergy;
  
  su2double *Normal = new su2double[nDim];
  su2double *PrimVar_i = new su2double[nPrimVar];
  su2double *PrimVar_j = new su2double[nPrimVar];
  
  su2double P_static, rho_static;
  
  int rank = MASTER_NODE, irank = 1;
  
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE)	{

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->node[iPoint]->GetDomain()) {

          for (iVar = 0; iVar < nPrimVar; iVar++) {
            PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
            PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar);
          }
          
          /*--- Set primitive variables ---*/

          numerics->SetPrimitive( PrimVar_i, PrimVar_j );
          
          if( !( config->GetKind_FluidModel() == STANDARD_AIR || config->GetKind_FluidModel() == IDEAL_GAS ) ){
	        Secondary_i = node[iPoint]->GetSecondary();

		P_static   = PrimVar_j[nDim+1];
		rho_static = PrimVar_j[nDim+2];			  
                FluidModel->SetTDState_Prho(P_static, rho_static);
			
                Secondary_j[0] = FluidModel->GetdPdrho_e();
                Secondary_j[1] = FluidModel->GetdPde_rho();  
                   
                numerics->SetSecondary(Secondary_i, Secondary_j);
	   }

          /*--- Set the normal vector ---*/

          geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
          for (iDim = 0; iDim < nDim; iDim++) 
            Normal[iDim] = -Normal[iDim];

          numerics->SetNormal(Normal);

          if (grid_movement)
            numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

          /*--- Compute the convective residual using an upwind scheme ---*/

          numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

          /*--- Add Residuals and Jacobians ---*/

          LinSysRes.AddBlock(iPoint, Residual);
          if (implicit) 
            Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        }
      }
    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;
  delete [] PrimVar_i;
  delete [] PrimVar_j;
}
