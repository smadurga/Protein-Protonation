#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define MAXRES 6000

#define ASP 1
#define GLU 2
#define ARG 3
#define LYS 4
#define CYS 5
#define HIS 6
#define TYR 7
#define Cterm 8
#define Nterm 9

//Para oplsaa.f

#define PROTONAT_Nterm 0
#define DESPROTONAT_Nterm 2
#define PROTONAT_Cterm 2
#define DESPROTONAT_Cterm 0

#define PROTONAT_HIS 2
#define PROTONAT 1
#define DESPROTONAT 0

#define PROTONAR 10
#define DESPROTONAR 20

#define HISD 0
#define HISE 1
#define HISH 2

int main(int narg, char **arg);
float ran2(long *idum);

float alea;
long idum;

int main(int narg, char **arg)
{
FILE *fpIN,*fpOUT;
char nameIN[100],nameOUT[100],nameOUText[100];
char ext[11][4]={".0",".1",".2",".3",".4",".5",".6",".7",".8",".9",".10"};
int n,niter,nconfs,numres,totalres,nResults;
char comentari[500];
float pkares[MAXRES],pkaini[MAXRES];
int indexres[MAXRES],resprotonacio[MAXRES];
int restype[MAXRES];
char nomres[MAXRES][10],nomcad[MAXRES][2];
float pH;
int Charge,nMC;
int numASP,numGLU,numARG,numLYS,numCYS,numHIS,numTYR,numCterm,numNterm;
double difE,RT,bf,SGComega;
int procesSGC;
int qTotal,qASP,qGLU,qARG,qLYS,qHIS,qCterm,qNterm;
char chain[30][3]={"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y"};
char chain1[3];
int first_chain,last_chain,num_chain;
//int HIS1=486;   /* Protonada en ND1 . Con gromacs código: 0*/
//int HIS2=374, HIS3=376, HIS4=512;  /* Protonada en NE2. Con gromacs código: 1 */
                            /* Histidina cargada. Doble protonada. Con gromacs código 2: */
                                  




  if(narg==6)
    {
          sscanf(arg[1],"%99s",nameIN);
          sscanf(arg[2],"%99s",nameOUT);
          sscanf(arg[3],"%f",&pH);
        //  sscanf(arg[4],"%d",&Charge);
          sscanf(arg[4],"%d",&nMC);
          sscanf(arg[5],"%ld",&idum);
    }
    else
    {
 	printf("EXEC.x  fileIN fileOUT pH nMC idum\n");
	printf("FileIN:  Summary of Propka\n");
	printf(" Line1: NumDades\n");
	printf(" Line2: Comment\n");
	printf(" .... Dades \n");
	return(1);
    }



printf("microQ - per IM\n");
printf("IN: %s, OUT %s, ph %f, Charge %d, nMC %d, idum %ld\n",nameIN,nameOUT,pH,Charge,nMC,idum);

    fpIN = fopen(nameIN,"r");
    if (fpIN == NULL)
    {    printf("Error al abrir el archivo %s\n",nameIN);
         exit (1);  }

numASP=0;
numGLU=0;
numARG=0;
numLYS=0;
numHIS=0;
numCYS=0;
numCterm=0;
numNterm=0;

    if(fscanf(fpIN,"%d\n",&numres)!=1)
         {    printf("Error al leer el numero de residuos en el archivo %s\n",nameIN);
         exit (1);  }
    if(fgets(comentari,200,fpIN)==NULL)
		{    printf("Error al leer el comentari en el archivo %s\n",nameIN);
         exit (1);  }
    
    printf("Numatoms: %d, Comentari: %s",numres,comentari);


    for(n=1;n<=numres;n++)
      {  if(fscanf(fpIN, "%3s%4d%2s%f%f\n", nomres[n],&indexres[n],nomcad[n],&pkares[n],&pkaini[n])!=5)
		     {    printf("Error al leer los 5 parametros en una linea en el archivo %s\n",nameIN);
             exit (1);  }


       printf("%s %d %s %f %f\n",
         nomres[n],indexres[n],nomcad[n],pkares[n],pkaini[n]);

      if(strcmp(nomres[n],"ASP")==0)
        {  // printf("ES ASP** %d\n",n);
           numASP++;
	   restype[n]=ASP;  
	   resprotonacio[n]=PROTONAT;     
        }

      if(strcmp(nomres[n],"GLU")==0)
        { // printf("ES GLU**\n");       
           numGLU++;
	   restype[n]=GLU;       
	   resprotonacio[n]=PROTONAT;     
        }

      if(strcmp(nomres[n],"ARG")==0)
        { // printf("ES ARG**\n");       
           numARG++;
	   restype[n]=ARG;       
	   resprotonacio[n]=PROTONAT;     
        }

      if(strcmp(nomres[n],"LYS")==0)
        { // printf("ES LYS**\n");       
           numLYS++;
	   restype[n]=LYS;       
	   resprotonacio[n]=PROTONAT;     
        }

      if(strcmp(nomres[n],"CYS")==0)
        { // printf("ES CYS**\n");       
           numCYS++;
	   restype[n]=CYS;       
	   resprotonacio[n]=PROTONAT;     
        }

      if(strcmp(nomres[n],"HIS")==0)
        { // printf("ES HIS**\n");       
           numHIS++;
	   restype[n]=HIS;       
	   resprotonacio[n]=PROTONAT_HIS;     
        }

      if(strcmp(nomres[n],"TYR")==0)
        { // printf("ES TYR**\n");       
           numTYR++;
	   restype[n]=TYR;       
	   resprotonacio[n]=PROTONAT;     
        }


      if(strcmp(nomres[n],"C-")==0)
        { // printf("ES C-**\n");       
           numCterm++;
	   restype[n]=Cterm;       
	   resprotonacio[n]=DESPROTONAT_Cterm;    /* PER AMBER COO- */ 
        }


      if(strcmp(nomres[n],"N+")==0)
        { // printf("ES N+**\n");       
           numNterm++;
	   restype[n]=Nterm;       
	   resprotonacio[n]=PROTONAT_Nterm;    /* PER AMBER NH3+ */ 
        }

      }  /* for n per dades*/

fclose(fpIN);

printf("numASP: %d\n",numASP);
printf("numGLU: %d\n",numGLU);
printf("numARG: %d\n",numARG);
printf("numLYS: %d\n",numLYS);
printf("numCYS: %d\n",numCYS);
printf("numHIS: %d\n",numHIS);
printf("numTYR: %d\n",numTYR);
printf("numC-:  %d\n",numCterm);
printf("numN+:  %d\n",numNterm);

totalres=numASP+numGLU+numARG+numLYS+numCYS+numHIS+numTYR+numCterm+numNterm;
printf("Total:  %d\n\n",totalres);

if(totalres!=numres)
 {
	printf("WARNING: totalres != numres \n");
	exit(1);
 }

printf("CYS and TYR are not considered for protonation\n");


/* CAMBIO ESTADO PROTONACION: HISTIDINAS CENTRO ACTIVO */
//for(n=1;n<=numres;n++)
//{   if(indexres[n]%808==HIS1)   /***********************/
//     {  resprotonacio[n]=HISD; }
//        
//   if(indexres[n]%808==HIS2 || indexres[n]%808==HIS3 || indexres[n]%808== HIS4)
//     {   		 resprotonacio[n]=HISE; 
//		 printf("fix HIS prot: %d ; indexres: %d\n",resprotonacio[n],indexres[n]); 
//		 }

//}



/***** CALCUL  ****/
nconfs=0;
nResults=0;
do
 {
  for(niter=1;niter<=nMC;niter++)
  {
   nconfs++;
    for(n=1;n<=numres;n++)
    {

     // El número de veces que se prueba el cambio que sea aleatorio.
     // Importante para casos en que pH=Pka, para que no se cambie en cada iteración.
     alea=ran2(&idum);
     if(alea<=0.01)              // Un 1% dels casos no testejat. Per aleatorietat del nombre de casos per residu.
       { continue; }
     

//    if(indexres[n]%808==HIS1 || indexres[n]%808==HIS2 || indexres[n]%808==HIS3 || indexres[n]%808== HIS4)
//      { continue; }  /* NO CAMBIAR ESTADO PROTONACION DE HISTIDINAS CENTRO ACTIVO */

//////////  Acceptar PROTONACIO/DESPROTONACIO? //////////

difE=0;
RT=1; /* No efecto si difE=0 */

/* Para residuos standard */
  if(restype[n]==ASP || restype[n]==GLU || restype[n]==LYS || restype[n]==ARG ) 
     { 
           if(resprotonacio[n]==PROTONAT || ( resprotonacio[n]==PROTONAT_HIS) )
            { procesSGC=DESPROTONAR; }
            else
            { procesSGC=PROTONAR; }
      }
      
/* Para Histidina */
  if( restype[n]==HIS) 
     { 
           if( resprotonacio[n]==PROTONAT_HIS )
            { procesSGC=DESPROTONAR; }
            else
            { procesSGC=PROTONAR; }
      }
      
/* Para Cterm */
  if( restype[n]==Cterm) 
     { 
           if( resprotonacio[n]==PROTONAT_Cterm )
            { procesSGC=DESPROTONAR; }
            else
            { procesSGC=PROTONAR; }
      }

/* Para Nterm */
  if( restype[n]==Nterm) 
     { 
           if( resprotonacio[n]==PROTONAT_Nterm )
            { procesSGC=DESPROTONAR; }
            else
            { procesSGC=PROTONAR; }
      }

      
      



	     if(procesSGC==PROTONAR)
	      { bf=-difE/RT+log(10)*(pkares[n]-pH);}
	     else
	      { bf=-difE/RT-log(10)*(pkares[n]-pH);} 

//printf("%d %f %f\n",n,pkares[n],bf);

     if(bf<-709)
       { /* No es acceptat el SGC: No fer res */ } /* REBUTJAT PER INCREMENT MASSA GRAN DE L'ENERGIA */
   else
  {

      SGComega=exp(bf);

      alea=ran2(&idum);
      if(SGComega>=alea)
        { /* Canvi acceptat */
		/* Per residu standard */
	      if(restype[n]==ASP || restype[n]==GLU || restype[n]==LYS || restype[n]==ARG ) 
            { if(procesSGC==PROTONAR)
               {  resprotonacio[n]=PROTONAT; }   
			        else
			    { resprotonacio[n]=DESPROTONAT; }
			        
            }
          /* Residu HIS */  
          if( restype[n]==HIS) 
            { if(procesSGC==PROTONAR)
               {  resprotonacio[n]=PROTONAT_HIS; }   
			        else
			    { resprotonacio[n]=DESPROTONAT; }
			        
            }
           /* Para Cterm */
          if( restype[n]==Cterm) 
           { if(procesSGC==PROTONAR)
               {  resprotonacio[n]=PROTONAT_Cterm; }   
			        else
			    { resprotonacio[n]=DESPROTONAT_Cterm; }
			        
            }

           /* Para Nterm */
          if( restype[n]==Nterm) 
           { if(procesSGC==PROTONAR)
               {  resprotonacio[n]=PROTONAT_Nterm; }   
			        else
			    { resprotonacio[n]=DESPROTONAT_Nterm; }
			        
            }

  
    }
       

   } 

    //  } /* Fi if per restype que poden canviar destat de protonacio */               
    }  /* for n per calcul*/
  } /* for niter de nMC */


/*************************************/
/* Calcul carrega total */

qASP=0;
qGLU=0;
qARG=0;
qLYS=0;
qHIS=0;
qCterm=0;
qNterm=0;

    for(n=1;n<=numres;n++)
    {
      if(restype[n]==ASP)
       { if(resprotonacio[n]==DESPROTONAT) { qASP--; }  }

      if(restype[n]==GLU)
       { if(resprotonacio[n]==DESPROTONAT) { qGLU--; }  }

      if(restype[n]==ARG)
       { if(resprotonacio[n]==PROTONAT) { qARG++; }  }

      if(restype[n]==LYS)
       { if(resprotonacio[n]==PROTONAT) { qLYS++; }  }

      if(restype[n]==HIS)
       { if(resprotonacio[n]==PROTONAT_HIS) { qHIS++; }  }

      if(restype[n]==Cterm)
       { if(resprotonacio[n]==DESPROTONAT) { qCterm--; }  }

      if(restype[n]==Nterm)
       { if(resprotonacio[n]==PROTONAT) { qNterm++; }  }

    }
qTotal=qASP+qGLU+qARG+qLYS+qHIS+qCterm+qNterm;

if(nconfs<1000 || nconfs%1000==0)
{
printf("Protein Charge: %d   Nconf: %d \n",qTotal,nconfs);
printf("Charge ASP %d, GLU %d, ARG %d, LYS %d, HIS %d, Cterm %d, Nterm %d\n",qASP,qGLU,qARG,qLYS,qHIS,qCterm,qNterm);
}

/***********************************/
/*GENERA OUTPUT PER GROMACS */


nResults++;

printf("**CALCUL: %d \n",nResults);
printf("**Protein Charge: %d   Nconf: %d \n",qTotal,nconfs);
printf("**Charge ASP %d, GLU %d, ARG %d, LYS %d, HIS %d, Cterm %d, Nterm %d\n",qASP,qGLU,qARG,qLYS,qHIS,qCterm,qNterm);

strcpy(nameOUText,nameOUT);
strcat(nameOUText,ext[nResults]);
    fpOUT = fopen(nameOUText,"w");
    if (fpOUT == NULL)
    {    printf("Error al abrir el archivo %s\n",nameOUText);
         exit (1);  }

/* Per force-field AMBER */
/*   fprintf(fpOUT,"1\n");  /* force-filed AMBER03 */
/*   fprintf(fpOUT,"1\n");  /* recommended type of Water */




first_chain=0;
last_chain=1;   //NUMERO DE CADENAS -1

printf("NUMERO DE CADENAS: %d\n",last_chain+1); 




/* Per LYS */
for(n=1;n<=numres;n++)
{
  for(num_chain=first_chain;num_chain<=last_chain;num_chain++)
  {
   strcpy(chain1,chain[num_chain]);

   //printf("cadena: %s ",chain1);	
	
    if(restype[n]==LYS and strcmp(nomcad[n],chain1)==0)
     { printf("LYS %d ; %d\n",resprotonacio[n],indexres[n]);  
       fprintf(fpOUT,"%d \n",resprotonacio[n]); }
   }
}

/* Per ARG */
for(n=1;n<=numres;n++)
{   
   for(num_chain=first_chain;num_chain<=last_chain;num_chain++)
   {
   strcpy(chain1,chain[num_chain]);

   //printf("cadena: %s ",chain1);
    if(restype[n]==ARG and strcmp(nomcad[n],chain1)==0)
     { printf("ARG %d ; %d\n",resprotonacio[n],indexres[n]);  
       fprintf(fpOUT,"%d \n",resprotonacio[n]); }
   }
}     

/* Per ASP */
for(n=1;n<=numres;n++)
{  
   for(num_chain=first_chain;num_chain<=last_chain;num_chain++)
   {
   strcpy(chain1,chain[num_chain]);

   //printf("cadena: %s ",chain1);
	    
    if(restype[n]==ASP and strcmp(nomcad[n],chain1)==0)
     { printf("ASP %d  ; %d\n",resprotonacio[n],indexres[n]);  
       fprintf(fpOUT,"%d \n",resprotonacio[n]); }
   }
}    

/* Per GLU */
for(n=1;n<=numres;n++)
{  
   for(num_chain=first_chain;num_chain<=last_chain;num_chain++)
   {
   strcpy(chain1,chain[num_chain]);

   //printf("cadena: %s ",chain1);     
    if(restype[n]==GLU and strcmp(nomcad[n],chain1)==0)
     { printf("GLU %d  ; %d\n",resprotonacio[n],indexres[n]);  
       fprintf(fpOUT,"%d \n",resprotonacio[n]); }
   }
}     


/* Per HIS */
for(n=1;n<=numres;n++)
{   
   for(num_chain=first_chain;num_chain<=last_chain;num_chain++)
   {
   strcpy(chain1,chain[num_chain]);

   //printf("cadena: %s ",chain1);
//	if(indexres[n]==HIS1 and strcmp(nomcad[n],chain1)==0)
//     { printf("HIS %d  ; %d\n",HISD,indexres[n]);   /* Protonado en ND1 */ 
//       fprintf(fpOUT,"%d \n",HISD); 
//
//       continue; } 
//       
//   if((indexres[n]==HIS2 || indexres[n]==HIS3 || indexres[n]== HIS4) and strcmp(nomcad[n],chain1)==0)
//     { printf("HIS %d  ; %d\n",HISE,indexres[n]);  /* Protonado en NE2 */
//       fprintf(fpOUT,"%d \n",HISE); 
//
//       continue; } 

    if(restype[n]==HIS and strcmp(nomcad[n],chain1)==0)
     { printf("HIS %d  ; %d\n",resprotonacio[n],indexres[n]);    /* 0,1 es DESPROTONAT i 2 es PROTONAT */
       fprintf(fpOUT,"%d \n",resprotonacio[n]); 
       }   /* 0,1 es DESPROTONAT i 2 es PROTONAT */
   }
}     

/* Per Nterm */
for(n=1;n<=numres;n++)
{  
   for(num_chain=first_chain;num_chain<=last_chain;num_chain++)
   {
   strcpy(chain1,chain[num_chain]);

   //printf("cadena: %s ",chain1);
	     
    if(restype[n]==Nterm and strcmp(nomcad[n],chain1)==0)
     { printf("N+ %d ; %d\n",resprotonacio[n],indexres[n]);
       fprintf(fpOUT,"%d \n",resprotonacio[n]); }
   }
}

/* Per Cterm */
for(n=1;n<=numres;n++)
{   
   for(num_chain=first_chain;num_chain<=last_chain;num_chain++)
   {
   strcpy(chain1,chain[num_chain]);

   //printf("cadena: %s ",chain1);
	    
    if(restype[n]==Cterm and strcmp(nomcad[n],chain1)==0)
     {   printf("C- %d ; %d\n",resprotonacio[n],indexres[n]); 
		 fprintf(fpOUT,"%d \n",resprotonacio[n]); }     
   }
} 




fclose(fpOUT);
 /* Fin del output per gromacs */

}while(nResults<10);
return(0);
}

/*********************************************************************
  3. This random number generator is from William H. Press, et al.,
     _Numerical Recipes in C_, Second Ed. with corrections (1994),
     p. 282.  This excellent book is available through the
     WWW at http://nr.harvard.edu/nr/bookc.html.
     The specific section concerning ran2, Section 7.1, is in
     http://cfatab.harvard.edu/nr/bookc/c7-1.ps
*********************************************************************/

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)

/* ran2() - Return a random floating point value between 0.0 and
   1.0 exclusive.  If idum is negative, a new series starts (and
   idum is made positive so that subsequent calls using an unchanged
   idum will continue in the same sequence). */

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {                             /* initialize */
    if (-(*idum) < 1)                           /* prevent idum == 0 */
      *idum = 1;
    else
      *idum = -(*idum);                         /* make idum positive */
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0)
        *idum += IM1;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;                                /* avoid endpoint */
  else
    return temp;
}
/**************************************************************************/
