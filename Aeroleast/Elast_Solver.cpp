

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
 *     Function Linear Modal Solver specific ti BSCW wing test case
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
 *
 *
 *     Description
 * 
 */


#include "/home/jka/Cprograms/Sources/libm3l/Source/data_util/libm3l.h"
#include "/home/jka/Cprograms/Sources/lsipdx/Source/lsipdx.h"

int main(int argc, char *argv[])
{
	node_t *Gnode=NULL, *Snode=NULL, *FoundNode=NULL, *TmpNode=NULL;
	size_t i, niter, dim[1];

	lmint_t sockfd, portno, restart;

        socklen_t clilen;
        struct sockaddr_in cli_addr;
	lmchar_t *name ="CFD2SIM";
	lmchar_t *name1="SIM2CFD";
        lmchar_t *action;

	lmdouble_t *tmpfloat, *time, sign, *ModalForce, f1, f2, w1, w2, md1, md2, mo1, mo2;
        lmdouble_t A11,A12,A13,A21,A22,A23,mf1_n,mf1_n1,mf1_n2,mf2_n,mf2_n1,mf2_n2;
	lmdouble_t psi,theta,phi,*dt,t,Ytranslation;
        lmdouble_t Q1n3, Q1n2, Q1n1, Q2n3, Q2n2, Q2n1;
	lmdouble_t a,b,c,d,timef, count;
	
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
        Q1n3 = 0;
        Q1n2 = 0;
        Q2n3 = 0;
        Q2n2 = 0;

        printf("Restart [1] or not [0]\n");
        scanf("%d", &restart);

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

			mf1_n = a;
			mf2_n = b;
			Q1n1 = c;
			Q2n1 = d;
		}


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
	
// 		if(m3l_Cat(Gnode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
// 			Error("CatData");
// 		if(m3l_Cat(Gnode, "--detailed", "-P", "-L",  "*",   (char *)NULL) != 0)
// 		Error("CatData");

/*
 * get modal forces
 */
	if( (SFounds = m3l_Locate(Gnode, "/CFD_2_STR/ModalF", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("Elast_str:: More then one ModalF data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("Elast_str:: Did not find 1st data pointer");
		
		if( (ModalForce = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("Elast_str:: Did not find ModalF data pointer");

                printf("Modal forces are %lf  %lf \n",ModalForce[0], ModalForce[1]);
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
	dt = (lmdouble_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);
/*
 * find position in loop
 */
	SFounds = m3l_Locate(Gnode, "/CFD_2_STR/Action", "/*/*",  (lmchar_t *)NULL);
	FoundNode = m3l_get_Found_node(SFounds, 0);
	action = (lmchar_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);
/*
 *  close socket
 */
	if( close(sockfd) == -1)
		Perror("close");	
/*
 * get modal coordinates, set specific parameters for BSCW test case
 */

/*
 * first mode - plunging mode
 */
        t = *dt;
	printf("Time step is %lf \n", t);
	f1 = 3.33000000000000;  // frequency
	w1 = 2*3.1415926*f1;
	md1 = 0;                // damping

	A11 = 0.25*w1*w1 + md1*w1/t+1./(t*t);
	A12 = 0.25*w1*w1 - 2./(t*t);
	A13 = 0.25*w1*w1 - md1*w1/t+1./(t*t);
/*
 * first mode - pitching mode
 */
	f2 = 5.2000000000000;   // frequency
	w2 = 2*3.1415926*f2;
	md2 = 0;                // damping

	A21 = 0.25*w2*w2 + md2*w2/t+1./(t*t);
	A22 = 0.25*w2*w2 - 2./(t*t);
	A23 = 0.25*w2*w2 - md2*w2/t+1./(t*t);
/*
 * modal forces - they were recevied from SU2
 */
	mf1_n = ModalForce[0];
	mf2_n = ModalForce[1];
    
    if(niter ==0  && restart == 0){
   printf(" Nullifying \n");
    mf1_n=0;
    mf2_n=0;
    }

/*
 * get new modal coordinates
 */
	mo1 = 1;     // modal mass
	Q1n1 = ((mf1_n + 2*mf1_n1 +mf1_n2)/(4*mo1) - A12*Q1n2 -A13*Q1n3)/A11;

	mo2 = 1;     // modal mass
	Q2n1 = ((mf2_n + 2*mf2_n1 +mf2_n2)/(4*mo2) - A22*Q2n2 -A23*Q2n3)/A21;
/*
 * shift solution
 */
     if( strncmp(action, "shift", 5) == 0){

        fprintf(fp, "%lf %lf %lf %lf %lf\n", *time, Q2n1*27.264, Q1n1*0.1066528,mf1_n,mf2_n);
        fflush(fp);

        ++niter;

        mf1_n2  = mf1_n1;
        mf1_n1  = mf1_n;

        mf2_n2  = mf2_n1;
        mf2_n1  = mf2_n;

        Q1n3 = Q1n2;
        Q1n2 = Q1n1;
        Q2n3 = Q2n2;
        Q2n2 = Q2n1; }
/*
 * get pitching angle and translation
 */
	psi    = Q2n1*27.264;
	theta  = 0;
	phi    = 0;

        Ytranslation = Q1n1*0.1066528;


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
 * add Angles 
 */	
	
	if(  (TmpNode = m3l_Mklist("RotCenter", "D", 1, dim, &Snode, "/STR_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
	tmpfloat[0] = 0.2023;
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



