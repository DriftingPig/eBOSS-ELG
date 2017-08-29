#include<Python.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
#include <math.h>
#define PI 3.1415926
#define RESu 100
#define RESs 40

#include<time.h>
void headreader (char *filename)

/**********************************************************************/
/* Print out all the header keywords in all extensions of a FITS file */
/**********************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    
    int status=0, nkeys, keypos, hdutype, ii, jj, column;
    char *template="dec";
    //    char filename[]  = "atestfil.fit";     /* name of existing FITS file   */
    char card[FLEN_CARD];   /* standard string lengths defined in fitsioc.h */
    
    status = 0;
    
    if ( fits_open_image(&fptr, filename, READONLY, &status) )
        fits_report_error(stderr, status);
    
    /*attempt to move to next HDU, until we get an EOF error */
    for (ii = 1; !(fits_movabs_hdu(fptr, ii, &hdutype, &status) ); ii++)
    {
        //       get no. of keywords
        if (fits_get_hdrpos(fptr, &nkeys, &keypos, &status) )
            fits_report_error(stderr, status);
        
        printf("Header listing for HDU #%d:\n", ii);
        for (jj = 1; jj <= nkeys; jj++)  {
            if ( fits_read_record(fptr, jj, card, &status) )
                fits_report_error(stderr, status);
            printf("%s\n", card);
            //print the keyword card
        }
        printf("END\n\n");  // terminate listing with END
    }
    //fits_get_num_cols(fptr,&column, &status);
    
    
    if (status == END_OF_FILE)   /* status values are defined in fitsioc.h */
        status = 0;              /* got the expected EOF error; reset = 0  */
    else
        fits_report_error(stderr, status);
    if ( fits_close_file(fptr, &status) )
        fits_report_error(stderr, status);
    return;
}
double ra2phi(double ra){return ra*PI/180;}
double dec2theta(double dec){return dec*PI/180.;}
double dc(double z){
    double distance;
    PyObject *pModule,*pFunc;
    PyObject *pArgs, *pValue;
    /*import*/
    pModule = PyImport_Import(PyString_FromString("romberg"));
    /* great_module.great_function */
    pFunc = PyObject_GetAttrString(pModule, "dc");
    /* build args */
    pArgs = PyTuple_New(1);
    PyTuple_SetItem(pArgs,0, PyFloat_FromDouble(z));
    
    /* call */
    pValue = PyObject_CallObject(pFunc, pArgs);
    
    distance = PyFloat_AsDouble(pValue);
    return distance;
    
}
double Separation(double thi1,double phi1,double r1,double thi2, double phi2,double r2){
    double temp;
    //	printf("\nLINE80@@: %lf %lf %lf %lf %lf %lf\n",thi1,phi1,r1,thi2,phi2,r2);
    //	printf("******LINE81*****\n");
    //	temp=cos(thi1);
    //	printf("temp: %lf\n",temp);
    //  temp=thi2;
    // printf("temp: %lf\n",temp);
    return sqrt(r1*r1+r2*r2-2*r1*r2*(cos(thi1)*cos(thi2)*cos(phi1-phi2)+sin(thi1)*sin(thi2)));
}
double Orientation(double r1,double r2,double r){
    return fabs((r1-r2)/r);
}
double Separation2(double COS1,double SIN1,double R1_2,double R1,double COS2,double SIN2,double R2_2,double R2,double Phi1_Phi2){
    return sqrt(R1_2+R2_2-2*R1*R2*(COS1*COS2*cos(Phi1_Phi2)+SIN1*SIN2));
}


int main(int argc, char* argv[]){
   int status=0, nkeys, keypos, hdutype,nfound, ii, jj, column[4],column2[4],status2=0,flag1,flag2;
   long naxes[2],naxes2[2];
   static double *buffer[4],*buffer2[4],buff,buff2;//DD-2DR+RR & RR
   static double DD_total=0,DR_total=0,RR_total=0;
   static double Distribution[RESs][RESu];
   int frow = 1, felem = 1, nullval = -99.,anynulls,ibuff,ibuff2,outside=0;
   double *Cos1,*Sin1,*Cos2,*Sin2,*r1_2,*r2_2,min_dc,max_dc;
   double buff_phi,buff_theta;
   char temp[80]="";
   fitsfile *fptr,*fptr2;
   FILE *fpXX,*fptotal;
   clock_t t;
   int opt=atoi(argv[3]);
   int count1=0,count2=0;
  // printf("*Running...\n");
  //  headreader(argv[1]);
  //  headreader(argv[2]);
    /*open file1*/
    if ( fits_open_image(&fptr,argv[1], READONLY, &status) )
        fits_report_error(stderr, status);
    /*open file2*/
    if ( fits_open_image(&fptr2,argv[2], READONLY, &status2) )
        fits_report_error(stderr, status2);
    /*move hdu to the data area for file1*/
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    if(hdutype !=BINARY_TBL){
        printf("Error: expected to find a binary table in this HDU\n");
    }
    /*move hdu to the data area for file2*/
    fits_movabs_hdu(fptr2, 2, &hdutype, &status2);
    if(hdutype !=BINARY_TBL){
        printf("Error: expected to find a binary table in this HDU\n");
    }
    /*get the row number of file1*/
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
        fits_report_error(stderr, status);
    /*get the row number of file2*/
    if (fits_read_keys_lng(fptr2, "NAXIS", 1, 2, naxes2, &nfound, &status2) )
        fits_report_error(stderr, status2);
    /*print the total column number for data file and random file*/
    printf("total event1: %ld total event2: %ld\n",naxes[1],naxes2[1]);
    if(opt==1 || opt==5) {//DD self,cross
        /*open an arrary to store the columns for data file*/
        for(ii=0;ii<4;ii++)
            buffer[ii] = (double *) malloc(naxes[1]*(sizeof(double)));
        /*get the column number for desired name of columns for data file*/
        if (fits_get_colnum(fptr, CASEINSEN,"ra", &column[0], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr, CASEINSEN,"dec", &column[1], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr, CASEINSEN,"z", &column[2], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr, CASEINSEN,"w", &column[3], &status))
            fits_report_error(stderr, status);
        /*transfer the data to the column arrary for data file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr, TDOUBLE, column[ii], frow, felem, naxes[1],&nullval, buffer[ii], &anynulls, &status))
                fits_report_error(stderr, status);
        
        /*open an arrary to store the columns for data file*/
        for(ii=0;ii<4;ii++)
            buffer2[ii] = (double *) malloc(naxes2[1]*(sizeof(double)));
        /*get the column number for desired name of columns for data file*/
        if (fits_get_colnum(fptr2, CASEINSEN,"ra", &column2[0], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr2, CASEINSEN,"dec", &column2[1], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr2, CASEINSEN,"z", &column2[2], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr2, CASEINSEN,"w", &column2[3], &status))
            fits_report_error(stderr, status);
        /*transfer the data to the column arrary for data file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr2, TDOUBLE, column2[ii], frow, felem, naxes2[1],&nullval, buffer2[ii], &anynulls, &status))
                fits_report_error(stderr, status);
    }

    if(opt==2) {//RR,self
        /*open an arrary to store the columns for random file*/
        for(ii=0;ii<4;ii++)
            buffer[ii] = (double *) malloc(naxes[1]*(sizeof(double)));
        /*get the column number for desired name of columns for random file*/
        if (fits_get_colnum(fptr, CASEINSEN,"ra", &column[0], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr, CASEINSEN,"dec", &column[1], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr, CASEINSEN,"z", &column[2], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr, CASEINSEN,"w", &column[3], &status))
           fits_report_error(stderr, status);
        /*transfer the data to the column arrary for random file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr, TDOUBLE, column[ii], frow, felem, naxes[1],&nullval, buffer[ii], &anynulls, &status2))
                fits_report_error(stderr, status2);
        
        /*open an arrary to store the columns for random file*/
        for(ii=0;ii<4;ii++)
            buffer2[ii] = (double *) malloc(naxes[1]*(sizeof(double)));
        /*get the column number for desired name of columns for random file*/
        if (fits_get_colnum(fptr2, CASEINSEN,"ra", &column2[0], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"dec", &column2[1], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"z", &column2[2], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"w", &column2[3], &status))
            fits_report_error(stderr, status);
        /*transfer the data to the column arrary for random file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr2, TDOUBLE, column2[ii], frow, felem, naxes2[1],&nullval, buffer2[ii], &anynulls, &status2))
                fits_report_error(stderr, status2);
    }

    if(opt==3) {//DR
        /*open an arrary to store the columns for data file*/
        for(ii=0;ii<4;ii++)
            buffer[ii] = (double *) malloc(naxes[1]*(sizeof(double)));
        /*get the column number for desired name of columns for data file*/
        if (fits_get_colnum(fptr, CASEINSEN,"ra", &column[0], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr, CASEINSEN,"dec", &column[1], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr, CASEINSEN,"z", &column[2], &status))
            fits_report_error(stderr, status);
        if (fits_get_colnum(fptr, CASEINSEN,"w", &column[3], &status))
            fits_report_error(stderr, status);
        /*transfer the data to the column arrary for data file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr, TDOUBLE, column[ii], frow, felem, naxes[1],&nullval, buffer[ii], &anynulls, &status))
                fits_report_error(stderr, status);
        /*open an arrary to store the columns for random file*/
        for(ii=0;ii<4;ii++)
            buffer2[ii]= (double *) malloc(naxes2[1]*(sizeof(double)));
        /*get the column number for desired name of columns for random file*/
        if (fits_get_colnum(fptr2, CASEINSEN,"ra", &column2[0], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"dec", &column2[1], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"z", &column2[2], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"w", &column2[3], &status))
            fits_report_error(stderr, status);
        /*transfer the data to the column arrary for random file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr2, TDOUBLE, column2[ii], frow, felem, naxes2[1],&nullval, buffer2[ii], &anynulls, &status2))
                fits_report_error(stderr, status2);
    }
    
    if(opt==4) {//RR,cross
        /*open an arrary to store the columns for random file*/
        for(ii=0;ii<4;ii++)
            buffer[ii]= (double *) malloc(naxes[1]*(sizeof(double)));
        /*get the column number for desired name of columns for random file*/
        if (fits_get_colnum(fptr, CASEINSEN,"ra", &column[0], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr, CASEINSEN,"dec", &column[1], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr, CASEINSEN,"z", &column[2], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr, CASEINSEN,"w", &column[3], &status))
           fits_report_error(stderr, status2);
        /*transfer the data to the column arrary for random file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr, TDOUBLE, column[ii], frow, felem, naxes[1],&nullval, buffer[ii], &anynulls, &status2))
                fits_report_error(stderr, status2);
        
        /*open an arrary to store the columns for random file*/
        for(ii=0;ii<4;ii++)
            buffer2[ii]= (double *) malloc(naxes2[1]*(sizeof(double)));
        /*get the column number for desired name of columns for random file*/
        if (fits_get_colnum(fptr2, CASEINSEN,"ra", &column2[0], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"dec", &column2[1], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"z", &column2[2], &status2))
            fits_report_error(stderr, status2);
        if (fits_get_colnum(fptr2, CASEINSEN,"w", &column2[3], &status2))
             fits_report_error(stderr, status2);
        /*transfer the data to the column arrary for random file*/
        for(ii=0;ii<4;ii++)
            if (fits_read_col(fptr2, TDOUBLE, column2[ii], frow, felem, naxes2[1],&nullval, buffer2[ii], &anynulls, &status2))
                fits_report_error(stderr, status2);
    }
//TODO
//for(ii=0;ii<naxes[1];ii++){
//if(buffer[2][ii]>0.7 && buffer[2][ii]<1.0){
//buffer[2][count1]=buffer[2][ii];
//count1++;
//}
//}
//for(ii=0;ii<naxes2[1];ii++){
//if(buffer2[2][ii]>0.7 && buffer2[2][ii]<1.0){
//buffer[2][count2]=buffer[2][ii];
//count2++;
//}
//}
//printf("***count1: %d count2: %d\n",count1,count2);
//naxes[1]=count1;naxes2[1]=count2;
//TODO

    /*transport (ra,dec,z) to (theta,phi,r) for data file*/
    Py_Initialize();
    for(ii=0;ii<naxes[1];ii++){
        buff_phi = ra2phi(buffer[0][ii]);
        buff_theta = buffer[0][ii]=dec2theta(buffer[1][ii]);
        buffer[2][ii] = dc(buffer[2][ii]);
        buffer[1][ii] = buff_phi;
        buffer[0][ii] = buff_theta;
    }

    /*transport (ra,dec,z) to (theta,phi,r) for random file*/
    for(ii=0;ii<naxes2[1];ii++){
        buff_phi = ra2phi(buffer2[0][ii]);
        buff_theta = dec2theta(buffer2[1][ii]);
        buffer2[2][ii]=dc(buffer2[2][ii]);
        buffer2[1][ii] = buff_phi;
        buffer2[0][ii] = buff_theta;    
}
    min_dc=dc(0.7);
    max_dc=dc(1.0);
    Py_Finalize();
    /*store cos theta, sin theta ...*/
    Cos1 = (double *) malloc(naxes[1]*(sizeof(double)));
    Sin1 = (double *) malloc(naxes[1]*(sizeof(double)));
    r1_2= (double *) malloc(naxes[1]*(sizeof(double)));
    Cos2 = (double *) malloc(naxes2[1]*(sizeof(double)));
    Sin2 = (double *) malloc(naxes2[1]*(sizeof(double)));
    r2_2= (double *) malloc(naxes2[1]*(sizeof(double)));
    
    for(ii=0;ii<naxes[1];ii++){
        Cos1[ii]=cos(buffer[0][ii]);
        Sin1[ii]=sin(buffer[0][ii]);
        r1_2[ii]=buffer[2][ii]*buffer[2][ii];
    }
    
    for(ii=0;ii<naxes2[1];ii++){
        Cos2[ii]=cos(buffer2[0][ii]);
        Sin2[ii]=sin(buffer2[0][ii]);
        r2_2[ii]=buffer2[2][ii]*buffer2[2][ii];
    }
    if(opt==1)//DD,self
    {
        //calculate pairs for DD,DR,RR
        flag1=naxes[1];flag2=naxes2[1];
        for(ii=0;ii<flag1;ii++)//DD
            for(jj=ii+1;jj<flag2;jj++){
                if(buffer[2][ii]>min_dc && buffer[2][ii]<max_dc && buffer2[2][jj]>min_dc && buffer2[2][jj]<max_dc){
                    DD_total+=(buffer[3][ii])*(buffer2[3][jj]);
                    buff=Separation2(Cos1[ii],Sin1[ii],r1_2[ii],buffer[2][ii],Cos2[jj],Sin2[jj],r2_2[jj],buffer2[2][jj],buffer[1][ii]-buffer2[1][jj]);
                    buff2=Orientation(buffer[2][ii],buffer2[2][jj],buff);
                    //if(ii<10 && buff<10){ 
 //                   printf("ii=%d theta=%f phi=%f dc=%f Sep=%f Ori=%f\n",ii,buffer[0][ii],buffer[1][ii],buffer[2][ii],buff,buff2);
//                    printf("jj=%d theta2=%f phi2=%f dc=%f\n",jj,buffer2[0][jj],buffer2[1][jj],buffer2[2][jj]);
//}
                    ibuff=floor(buff/5);
                    ibuff2=floor(buff2/0.01);
                    if(ibuff>=0 && ibuff<RESs && ibuff2>=0 && ibuff2<RESu){
                        Distribution[ibuff][ibuff2]+=(buffer[3][ii])*(buffer2[3][jj]);
                }
            }}
        printf("DD total:%lf\n",DD_total);
        strcat(temp,"./NewData_subfiles/D");strcat(temp,argv[4]);strcat(temp,"D");strcat(temp,argv[5]);strcat(temp,".txt");
        fpXX=fopen(temp,"w+");
        
        for(ii=0;ii<RESs;ii++){
            for(jj=0;jj<RESu;jj++){
                fprintf(fpXX,"%f ",Distribution[ii][jj]);
            }
            fprintf(fpXX,"\n");
        }
        fclose(fpXX);
        
    }
    if(opt==2)//RR,self
    {
        flag1=naxes[1];flag2=naxes2[1];
        for(ii=0;ii<flag2;ii++)//RR
            for(jj=ii+1;jj<flag2;jj++){
                if(buffer[2][ii]>min_dc && buffer[2][ii]<max_dc && buffer2[2][jj]>min_dc && buffer2[2][jj]<max_dc){
                    RR_total+=(buffer[3][ii])*(buffer2[3][jj]);
                    buff=Separation(buffer[0][ii],buffer[1][ii],buffer[2][ii],buffer2[0][jj],buffer2[1][jj],buffer2[2][jj]);
                    buff2=Orientation(buffer[2][ii],buffer2[2][jj],buff);
                    ibuff=floor(buff/5);ibuff2=floor(buff2/0.01);
                    if(ibuff>=0 && ibuff<RESs && ibuff2>=0 && ibuff2<RESu){
                        Distribution[ibuff][ibuff2]+=(buffer[3][ii])*(buffer2[3][jj]);
                    }
                }
            }
        printf("RR self total:%lf\n",RR_total);
        strcat(temp,"./NewData_subfiles/R");strcat(temp,argv[4]);strcat(temp,"R");strcat(temp,argv[5]);strcat(temp,".txt");
        fpXX=fopen(temp,"w+");
        
        for(ii=0;ii<RESs;ii++){
            for(jj=0;jj<RESu;jj++){
                fprintf(fpXX,"%f ",Distribution[ii][jj]);
            }
            fprintf(fpXX,"\n");
        }
        fclose(fpXX);
        
    }
    if(opt==3)//DR
    {
        flag1=naxes[1];flag2=naxes2[1];
        for(ii=0;ii<flag1;ii++)//DR
            for(jj=0;jj<flag2;jj++){
                if(buffer[2][ii]>min_dc && buffer[2][ii]<max_dc && buffer2[2][jj]>min_dc && buffer2[2][jj]<max_dc){
                    DR_total+=(buffer[3][ii])*(buffer2[3][jj]);
                    buff=Separation2(Cos1[ii],Sin1[ii],r1_2[ii],buffer[2][ii],Cos2[jj],Sin2[jj],r2_2[jj],buffer2[2][jj],buffer[1][ii]-buffer2[1][jj]);
                    buff2=Orientation(buffer[2][ii],buffer2[2][jj],buff);
                    ibuff=floor(buff/5);ibuff2=floor(buff2/0.01);
                    if(ibuff>=0 && ibuff<RESs && ibuff2>=0 && ibuff2<RESu){
                        Distribution[ibuff][ibuff2]+=(buffer[3][ii])*(buffer2[3][jj]);
                       // if(ibuff==0 && ibuff2==0) printf("DR_1: D:%d R:%d Dz:%f Rz: %f Sep: %f Orientation: %f\n",ii,jj,buffer[2][ii],buffer2[2][jj],buff,buff2);
                        
                    }
                }
            }
        printf("DR total:%lf\n",DR_total);
        strcat(temp,"./NewData_subfiles/D");strcat(temp,argv[4]);strcat(temp,"R");strcat(temp,argv[5]);strcat(temp,".txt");
        fpXX=fopen(temp,"w+");
        
        for(ii=0;ii<RESs;ii++){
            for(jj=0;jj<RESu;jj++){
                fprintf(fpXX,"%f ",Distribution[ii][jj]);
            }
            fprintf(fpXX,"\n");
        }
        fclose(fpXX);
    }
    if(opt==4)//RR,cross
    {
        flag1=naxes[1];flag2=naxes2[1];
        for(ii=0;ii<flag1;ii++)
            for(jj=0;jj<flag2;jj++){
                if(buffer[2][ii]>min_dc && buffer[2][ii]<max_dc && buffer2[2][jj]>min_dc && buffer2[2][jj]<max_dc){
                    RR_total+=(buffer[3][ii])*(buffer2[3][jj]);
                    buff=Separation2(Cos1[ii],Sin1[ii],r1_2[ii],buffer[2][ii],Cos2[jj],Sin2[jj],r2_2[jj],buffer2[2][jj],buffer[1][ii]-buffer2[1][jj]);
                    buff2=Orientation(buffer[2][ii],buffer2[2][jj],buff);
                    ibuff=floor(buff/5);ibuff2=floor(buff2/0.01);
                    if(ibuff>=0 && ibuff<RESs && ibuff2>=0 && ibuff2<RESu){
                        Distribution[ibuff][ibuff2]+=(buffer[3][ii])*(buffer2[3][jj]);
                    }
                }
            }
        
        printf("RR cross total:%lf\n",RR_total);
        strcat(temp,"./NewData_subfiles/R");strcat(temp,argv[4]);strcat(temp,"R"),strcat(temp,argv[5]);strcat(temp,".txt");
        fpXX=fopen(temp,"w+");
        
        for(ii=0;ii<RESs;ii++){
            for(jj=0;jj<RESu;jj++){
                fprintf(fpXX,"%f ",Distribution[ii][jj]);
            }
            fprintf(fpXX,"\n");
        }
        fclose(fpXX);
    }
        if(opt==5)//DD,cross
    {
        flag1=naxes[1];flag2=naxes2[1];
        for(ii=0;ii<flag1;ii++)
            for(jj=0;jj<flag2;jj++){
                if(buffer[2][ii]>min_dc && buffer[2][ii]<max_dc && buffer2[2][jj]>min_dc && buffer2[2][jj]<max_dc){
                    DD_total+=(buffer[3][ii])*(buffer2[3][jj]);
                    buff=Separation2(Cos1[ii],Sin1[ii],r1_2[ii],buffer[2][ii],Cos2[jj],Sin2[jj],r2_2[jj],buffer2[2][jj],buffer[1][ii]-buffer2[1][jj]);
                    buff2=Orientation(buffer[2][ii],buffer2[2][jj],buff);
                    ibuff=floor(buff/5);ibuff2=floor(buff2/0.01);
                    if(ibuff>=0 && ibuff<RESs && ibuff2>=0 && ibuff2<RESu){
                        Distribution[ibuff][ibuff2]+=(buffer[3][ii])*(buffer2[3][jj]);

                    }
                }
            }

        printf("DD cross total:%lf\n",DD_total);
        strcat(temp,"./NewData_subfiles/D");strcat(temp,argv[4]);strcat(temp,"D"),strcat(temp,argv[5]);strcat(temp,".txt");
        fpXX=fopen(temp,"w+");

        for(ii=0;ii<RESs;ii++){
            for(jj=0;jj<RESu;jj++){
                fprintf(fpXX,"%f ",Distribution[ii][jj]);
            }
            fprintf(fpXX,"\n");
        }
        fclose(fpXX);
    }
    if(atoi(argv[3])==1 && atoi(argv[4])==0 && atoi(argv[5])==0){
        fptotal=fopen("./NewData_subfiles/totalpoints.txt","w");
        printf("YES!\n");
    }
    else{
        fptotal=fopen("./NewData_subfiles/totalpoints.txt","a");
    }
    fprintf(fptotal,"%lf %lf %lf %d %d\n",DD_total,DR_total,RR_total,atoi(argv[4]),atoi(argv[5]));
    fclose(fptotal);
}

