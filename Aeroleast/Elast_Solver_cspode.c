

/*
 *     Copyright (C) 2012  Adam Jirasek
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *     
 *     contact: libm3l@gmail.com
 * 
 */



/*
 *     Function Linear Modal Solver specific to BSCW wing test case
 *
 *     Date: 2017-07-10
 * 
 * 
 *     Description:
 * 
 *
 *     Input parameters:
 * 
 *
 *     Return value:
 * 
 * 
 *
 *     Modifications:
 *     Date		Version		Patch number		CLA 
 *     2017-09-23 Modified by Cleber Spode
 *     1.Modified from central difference to Newmark Method
 *     2.Change the pointer dt to deta_t and variable t to dt
 *    
 *
 *     Description
 * 
 */


//#include "/home/jiraseka/Cprograms/Sources/libm3l/Source/data_util//libm3l.h"
//#include "/home/jiraseka/Cprograms/Sources/lsipdx/Source/lsipdx.h"

#include "Elast_Solver.h"


int main(int argc, char *argv[]) {
  node_t *Gnode = NULL, *Snode = NULL, *FoundNode = NULL, *TmpNode = NULL;
  size_t i, niter, dim[1];

  lmint_t sockfd, portno, restart, *comfreq;

  socklen_t clilen;
  struct sockaddr_in cli_addr;
  lmchar_t *name = "CFD2SIM";
  lmchar_t *name1 = "SIM2CFD";
  lmchar_t *action;

  lmdouble_t *tmpfloat, *time, sign, *ModalForce, f1, f2, w1, w2, md1, md2, mo1, mo2;
  lmdouble_t A11, A12, A13, A21, A22, A23, mf1_n, mf1_n1, mf1_n2, mf2_n, mf2_n1, mf2_n2;
  lmdouble_t psi, theta, phi, *delta_t, dt, Ytranslation;
  lmdouble_t Q1n3, Q1n2, Q1n1, Q2n3, Q2n2, Q2n1;

  // Modal accel, velocity and displacements for mode 1 - Plunge
  lmdouble_t q1_ddot_n, q1_ddot_np1, q1_dot_n, q1_dot_np1, q1_n, q1_np1, p1_np1;
  // Newmark variables for mode 1 - Plunge
  lmdouble_t A1, b11, b12, b13, b14;

  // Modal accel, velocity and displacements for mode 2 - Pitch
  lmdouble_t q2_ddot_n, q2_ddot_np1, q2_dot_n, q2_dot_np1, q2_n, q2_np1, p2_np1;
  // Newmark variables for mode 2 - Pitch
  lmdouble_t A2, b21, b22, b23, b24;
  
  // Newmark method parameters
  lmdouble_t plunge;
  lmdouble_t gamma;
  lmdouble_t betha;
  
  // Structural parameters
  lmdouble_t mode_1, mode_2;
          
  lmdouble_t a, b, c, d, timef, count;

  find_t *SFounds;

  opts_t opts, *Popts_1;

  FILE *fp;

  client_fce_struct_t InpPar, *PInpPar;

  PInpPar = &InpPar;
  
/*
 * get port number
 */
	if (argc < 3) {
		fprintf(stderr,"ERROR, no IPaddress and port number provided\n");
		exit(1);
	}
 	portno = atoi(argv[2]);


        mf1_n2  = 0;
        mf1_n1  = 0;
        mf2_n2  = 0;
        mf2_n1  = 0;
        Q1n3 = 0.0;
        Q1n2 = 0.0;
        Q2n3 = 0;
        Q2n2 = 0;
        
        
        /*
         * Initialization
         */
        q1_n = 0.0;                      // Generalized coordinate mode 1 - Plunge
        q1_dot_n = 0.0;                  // Generalized velocity mode 1 - Plunge
        q1_ddot_n = 0.0;                 // Generalized acceleration mode 1 - Plunge
        
        q2_n = 0.0;                      // Generalized coordinate mode 2 - Pitch
        q2_dot_n = 0.0;                  // Generalized velocity mode 2 - Pitch
        q2_ddot_n = 0.0;                 // Generalized acceleration mode 2 - Pitch

        /*
         * Modal matrices ( mass and damping)
         */        
        
//        printf("Restart [1] or not [0]\n");
//        scanf("%d", &restart);
//        printf("Modal mass\n");
//        scanf("%lf %lf", &mo1, &mo2);
//        printf("Modal damping coefficients\n");
//        scanf("%lf %lf", &md1, &md2);
        
        //Manual setting of modal matrices
        restart = 0;
        mo1 = 1.0;
        mo2 = 1.0;
        md1 = 0.0;
        md2 = 0.0;
        
        f1 = 3.32991158567781;  // frequency
        w1 = 2*3.1415926*f1;
        
        f2 = 5.19961338983582;   // frequency
        w2 = 2*3.1415926*f2;
        
        mode_1 = 0.096153176967948;                    // Eigenvector for plunge   
        mode_2 = -0.464625688784754*180.0/3.141592654; // Eigenvector for pitch
       
        /*
         * Newmark Method constants
         */
        gamma = 0.5;
        betha = 1.0/12.0;
        
        /*
         * If restart is true, start reading the COORDS file to get positions
         * and velocities
         */

        if(restart == 1){

		if ( (fp = fopen("COORDS","r")) == NULL)
		   Perror("fopen");


		while (fscanf(fp,"%lf %lf %lf %lf %lf", &timef, &a, &b, &c, &d) != EOF) {

			mf1_n2  = mf1_n1;
			mf1_n1  = mf1_n;
			mf2_n2  = mf2_n1;
			mf2_n1  = mf2_n;

			Q1n3 = Q1n2;
			Q1n2 = Q1n1;
			Q2n3 = Q2n2;
			Q2n2 = Q2n1;

			mf1_n = c;
			mf2_n = d;
			Q2n1 = a/mode_2;
			Q1n1 = b/mode_1;
		}

                psi    = Q2n1*mode_2;
                theta  = 0;
                phi    = 0;
                Ytranslation = Q1n1*mode_1;
                printf("Initial pitching angle and plunge is %lf  %lf \n", psi, Ytranslation);

		if( fclose (fp) != 0)
			Perror("fclose");
	}

        if ( (fp = fopen("COORDS","aw")) == NULL)
           Perror("fopen");
/*
 * open socket - because we use more then just send - receive scenario
 * we need to open socket manualy and used Send_receive function with hostname = NULL, ie. as server
 * portno is then replaced by socket number
 */
	niter = 0;
 	while(1){

 		printf("\n\n--------------------------------    i = %ld\n\n", niter);
/*
 * open socket
 */
		PInpPar->channel_name = name;
		PInpPar->SR_MODE = 'R';
		if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
			Error("wrong client mode");
		Popts_1 = &opts;
		m3l_set_Send_receive_tcpipsocket(&Popts_1);
	
		if( (sockfd = open_connection_to_server(argv[1], portno, PInpPar, Popts_1)) < 1)
			Error("client_sender: Error when opening socket");

		Gnode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL);
/*
 * get modal forces
 */
	if( (SFounds = m3l_Locate(Gnode, "/CFD_2_STR/ModalF", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("Elast_str:: More then one ModalF data set found");
      
        /* 
       * pointer to list of found nodes
       */
      if ((FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
        Error("Elast_str:: Did not find 1st data pointer");

      if ((ModalForce = (lmdouble_t *) m3l_get_data_pointer(FoundNode)) == NULL)
        Error("Elast_str:: Did not find ModalF data pointer");

      printf("Modal forces are %lf  %lf \n", ModalForce[0], ModalForce[1]);
      
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("Elast_str:: ModalF not found\n");
	}
/*
 * find current time of simulations
 */
	SFounds = m3l_Locate(Gnode, "/CFD_2_STR/Time", "/*/*",  (lmchar_t *)NULL);
	FoundNode = m3l_get_Found_node(SFounds, 0);
	time = (lmdouble_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);

	printf("Time is %lf \n", *time);
/*
 * find time step of simulations
 */
	SFounds = m3l_Locate(Gnode, "/CFD_2_STR/DT", "/*/*",  (lmchar_t *)NULL);
	FoundNode = m3l_get_Found_node(SFounds, 0);
	delta_t = (lmdouble_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);
/*
 * find position in loop
 */
	SFounds = m3l_Locate(Gnode, "/CFD_2_STR/Action", "/*/*",  (lmchar_t *)NULL);
	FoundNode = m3l_get_Found_node(SFounds, 0);
	action = (lmchar_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);
/*
 * find frequency of communiction
 */
	SFounds = m3l_Locate(Gnode, "/CFD_2_STR/COMFREQ", "/*/*",  (lmchar_t *)NULL);
	FoundNode = m3l_get_Found_node(SFounds, 0);
	comfreq = (lmint_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);
/*
 *  close socket
 */
	if( close(sockfd) == -1)
		Perror("close");	
/*
 * If frequency of communication > 0, multiple time step 
 */
        dt = *delta_t;
        if(*comfreq > 0){
          dt = *comfreq*dt;
        }
/*
 * modal forces - they were recevied from SU2
 */
	mf1_n = ModalForce[0];
	mf2_n = ModalForce[1];
    
    if(niter ==0  && restart == 0){
    
       if( fclose (fp) != 0)
		Perror("fclose");
       if ( (fp = fopen("COORDS","w")) == NULL)
           Perror("fopen");
       printf(" Nullifying \n");
       
//       mf1_n=0;
//       mf2_n=0;
       
    }
    
/*
     * Newmark mode 1 - Plunge
     * cspode - 2017-09-30
     */
    p1_np1 = mf1_n;

    //Initialize acceleration if is not restart
    if (niter == 0 && restart == 0){
      q1_ddot_n = (p1_np1-2.0*md1*w1*q1_dot_n - w1*w1*q1_n)/mo1;
    }

    printf("q_ddot_n = %f \n", q1_ddot_n);

    b11 = q1_dot_n + (1 - gamma) * dt*q1_ddot_n;
    b13 = q1_n + dt * q1_dot_n + (0.5 - betha) * dt * dt*q1_ddot_n;
    A1 = mo1 + gamma * dt * 2.0 * md1 * w1 + betha * dt * dt * w1*w1;
    q1_ddot_np1 = (p1_np1 - 2.0 * md1 * w1 * b11 - w1 * w1 * b13) / A1;

    b22 = gamma * dt*q1_ddot_np1;
    b24 = dt * dt * betha*q1_ddot_np1;

    q1_dot_np1 = b11 + b22;
    q1_np1 = b13 + b24;
    
     /*
     * Newmark mode 2 - Pitch
     * cspode - 2017-09-30
     */
    p2_np1 = mf2_n;

    //Initialize acceleration if is not restart
    if (niter == 0 && restart == 0){
      q2_ddot_n = (p2_np1-2.0*md2*w2*q2_dot_n - w2*w2*q2_n)/mo2;
    }
    
    b21 = q2_dot_n + (1-gamma)*dt*q2_ddot_n;
    b23 = q2_n + dt*q2_dot_n + (0.5 - betha)*dt*dt*q2_ddot_n;
    A2 = mo2 + gamma*dt*2.0*md2*w2 + betha*dt*dt*w2*w2;
    q2_ddot_np1 = (p2_np1 - 2.0*md2*w2*b21 - w2*w2*b23)/A2;
    
    b22 = gamma*dt*q2_ddot_np1;
    b24 = dt*dt*betha*q2_ddot_np1;
    
    q2_dot_np1 = b21 + b22;
    q2_np1     = b23 + b24;
    
    printf("Rodando Versao CSPODE mod - somente pitch \n");

    /*
     * shift solution
     */
    if (restart == 1 && niter == 0) {
      printf("Skipping\n");
    } else {
      if (strncmp(action, "shift", 5) == 0) {

//        double analitic_plunge, wd;
//
//        wd = w1 * pow((1.0 - pow(md1, 2.0)), 0.5);
//
//        analitic_plunge = exp(-md1 * w1 * (*time))*((0.0 * cos(wd * (*time))+((10.0 * 3.1415926 / 180.0 + md1 * w1 * 0.0) / wd) * sin(wd * (*time))));
//
//        fprintf(fp, "%lf %lf %lf %lf %lf %lf  %lf\n", *time, Q2n2 * mode_2, Q1n2, mf1_n, mf2_n, q2_n, analitic_plunge);
        fprintf(fp, "%lf %lf %lf %lf %lf \n", *time, p2_np1, p1_np1, q2_np1* mode_2, q1_np1*mode_1);

        fflush(fp);

        ++niter;

        mf1_n2 = mf1_n1;
        mf1_n1 = mf1_n;

        mf2_n2 = mf2_n1;
        mf2_n1 = mf2_n;

        q1_n = q1_np1;
        q1_dot_n = q1_dot_np1;
        q1_ddot_n = q1_ddot_np1;

        q2_n = q2_np1;
        q2_dot_n = q2_dot_np1;
        q2_ddot_n = q2_ddot_np1;

      }
      /*
       * get pitching angle and translation
       */
      psi = q2_np1 * mode_2;
      theta = 0.0;
      phi = 0.0;
      Ytranslation = q1_np1 * mode_1;
    }

     printf("Cl and Cm are %lf  %lf \n", p1_np1/(mode_1*8082.32803*0.4064), p2_np1/(0.464625688784754*8082.32803*0.4064*0.4064));
     printf("Pitching angle and plunge is %lf  %lf \n", psi, Ytranslation);
/*
 * send modal coordinates (translation and angles) back to SU2
 * open socket for sending data back
 */
	PInpPar->channel_name = name1;
	PInpPar->SR_MODE = 'S';
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("wrong client mode");
	Popts_1 = &opts;
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	
	if( (sockfd = open_connection_to_server(argv[1], portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");

	if(  (Snode = m3l_Mklist("STR_2_CFD", "DIR", 0, 0, (node_t **)NULL, (const char *)NULL, (const char *)NULL, (char *)NULL)) == 0)
		Perror("m3l_Mklist");

	dim[0] = 3;
/*
 * add Angles 
 */
	if(  (TmpNode = m3l_Mklist("Angles", "D", 1, dim, &Snode, "/STR_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);

	tmpfloat[0] = psi;
	tmpfloat[1] = theta;
	tmpfloat[2] = phi;
/*
 * add RotCenter
 */	
	if(  (TmpNode = m3l_Mklist("RotCenter", "D", 1, dim, &Snode, "/STR_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
	tmpfloat[0] = 0.2032;
	tmpfloat[1] = 0.;
	tmpfloat[2] = 0.;
/*
 * add translation 
 */	
	if(  (TmpNode = m3l_Mklist("TransVec", "D", 1, dim, &Snode, "/STR_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
	tmpfloat[0] = 0;
	tmpfloat[1] = 0;
	tmpfloat[2] = Ytranslation;

	client_sender(Snode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL);

	
	if(m3l_Umount(&Gnode) != 1)
		Perror("m3l_Umount");
	if(m3l_Umount(&Snode) != 1)
		Perror("m3l_Umount");
/* 
 * close socket
 */
	if( close(sockfd) == -1)
		Perror("close");

 	}

	if( fclose (fp) != 0)
		Perror("fclose");

     return 0; 
}



