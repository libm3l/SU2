/*!
 * \file primal_grid_structure.inl
 * \brief In-Line subroutines of the <i>external_comm_structure</i> file.
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
	/*!
	 * \brief Get the elements that surround an element.
	 */
	inline char **GetBCnames() {return BCnames};
	/*!
	 * \param[in] val_face - Local index of the face.
	 */
	inline void SetBCnames(char **names);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline void SetType(char *type);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline char *GetType() {return Type};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline void SetPN(int *pn);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline int GetPN() {return PN};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline void SetIP(char *ipaddr);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline char *GetIP() {return IP};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline void SetDisplacement(int *Displacement);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline int *GetDisplacement() {return Displacement};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline void SetPartNr(int *PartNr);
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline int *GetPartNr() {return PartNr};
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline void Setcommunicator(int communicator);
	
	/*!
	 * \param[in] val_coord - Coordinates of the element.
	 */
	inline int Getcommunicator() {return communicator};


