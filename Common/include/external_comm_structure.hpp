/*!
 * \file primal_grid_structure.hpp
 * \brief Headers of the main subroutines for storing data sets
 *   for external solver communication
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

#pragma once

#include "./mpi_structure.hpp"

#include <iostream>
#include <vector>
#include <cstdlib>

#include "dual_grid_structure.hpp"
#include "config_structure.hpp"

using namespace std;

/*!
 * \class CPrimalGrid
 * \brief Class to define the numerical primal grid.
 * \author F. Palacios
 * \version 4.3.0 "Cardinal"
 */
class CExtIntrf {
protected:
	int  nBC;         /*!< \brief Number of BCs */
	char **BCnames;         /*!< \brief Names of BCs included in interface. */
	char *Type;             /*!< \brief Solver type */
	char *Iname;             /*!< \brief Name of interface */
	static int PN;		/*!< \brief Port number */
	char *IP;             /*!< \brief hostname or ip address */
	unsigned long *Displacement;	/*!< \brief Displacement vector for MPI operation*/
	unsigned int *PartNr;	/*!< \brief Vector containing number of partitions in interface */
        int           communicator  /*!< \brief MPI communicator */

public:
	/*!
	 * \brief Constructor of the class.
	 */
	CExtIntrf(void);
	/*!
	 * \overload
	 * \param[in] val_nNodes - Number of nodes of the element.
	 */
	CExtIntrf(char *name);
	
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CExtIntrf(void);
	/*!
	 * \brief Get the elements that surround an element.
	 */
	char **GetBCnames() {return BCnames};
	/*!
	 * \param[in] val_face - Local index of the face.
	 */
	void SetBCnames(char **names);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	void SetType(char *type);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	char *GetType() {return Type};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	void SetPN(int *pn);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	int GetPN() {return PN};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	void SetIP(char *ipaddr);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	char *GetIP() {return IP};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	void SetDisplacement(int *Displacement);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	int *GetDisplacement() {return Displacement};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	void SetPartNr(int *PartNr);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	int *GetPartNr() {return PartNr};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	void Setcommunicator(int communicator);
	
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	int Getcommunicator() {return communicator};



	/*!
	 * \param[in] val_face - Local index of the face.
	 */
	virtual void SetBCnames(char **names);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	virtual void SetPN(int *pn);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	virtual void SetIP(char *ipaddr);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	virtual void SetDisplacement(int *Displacement);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	virtual void SetPartNr(int *PartNr);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	virtual void Setcommunicator(int communicator);
};


#include "external_comm_structure.inl"
