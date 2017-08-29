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
double ra2phi(double ra){return ra*PI/180;}
double dec2theta(double dec){return dec*PI/180.;}
double Separation(double thi1,double phi1,double r1,double thi2, double phi2,double r2){
    double temp;
    return sqrt(r1*r1+r2*r2-2*r1*r2*(cos(thi1)*cos(thi2)*cos(phi1-phi2)+sin(thi1)*sin(thi2)));
}
double Orientation(double r1,double r2,double r){
    return fabs((r1-r2)/r);
}
double Separation2(double COS1,double SIN1,double R1_2,double R1,double COS2,double SIN2,double R2_2,double R2,double Phi1_Phi2){
    return sqrt(R1_2+R2_2-2*R1*R2*(COS1*COS2*cos(Phi1_Phi2)+SIN1*SIN2));
}
int linenum(char* filename){
    FILE *pf = fopen(filename,"r");
    char buf[1000];
    int lineCnt = 0;
    if(!pf){
       printf("fail to open the file\n")
       return -1;}
    while (fgets(buf,1000,pf)) lineCnt++;
    fclose(pf);
    return lineCnt;
}

int main(int argc, char* argv[]){
   int status=0, nkeys, keypos, hdutype,nfound, ii, jj, column[4],column2[4],status2=0,flag1,flag2;
   long naxes,naxes2;
   static double *buffer[4],*buffer2[4],buff,buff2;//DD-2DR+RR & RR
   static double DD_total=0,DR_total=0,RR_total=0;
   static double Distribution[RESs][RESu];
   double *Cos1,*Sin1,*Cos2,*Sin2,*r1_2,*r2_2,min_dc,max_dc;
   double buff_phi,buff_theta;
   char temp[80]="";
   fitsfile *fptr,*fptr2;
   FILE *fpXX,*fptotal,*fp1,*fp2;
   clock_t t;
   int opt=atoi(argv[3]);
   int count1=0,count2=0;
   double *phi,*theta,*r,*w,*phi2,*theta2,*r2,*z2;
    /*open file1*/
    fp1 = fopen(argv[1],"r")
    /*open file2*/
    fp2 = fopen(argv[2],"r")
    /*get the row number of file1*/
    naxes = linenum(argv[1])
    /*get the row number of file2*/
    naxes2 = linenum(argv[2])
    /*print the total row number for file1 and file2*/
    printf("total event1: %ld total event2: %ld\n",naxes,naxes2);
    /*make space to save ra,dec,z,w*/
    phi = (double *)malloc(naxes*sizeof(double));
    theta = (double *)malloc(naxes*sizeof(double));
    r = (double *)malloc(naxes*sizeof(double));
    w = (double *)malloc(naxes*sizeof(double));
    phi = (double *)malloc(naxes2*sizeof(double));
    theta = (double *)malloc(naxes2*sizeof(double));
    r2 = (double *)malloc(naxes2*sizeof(double));
    w2 = (double *)malloc(naxes2*sizeof(double));
    /*save the data in the files*/
    for(ii=0;ii<naxes;ii++){
    fscanf(fp, "%lf",&phi[i]);phi[i] = phi[i]*PI/180;
    fscanf(fp,"%lf",&theta[i]);theta[i] = theta[i]*PI/180.
    fscanf(fp,"%lf",&r[i]);
    fscanf(fp,"%lf",&w);
    }
    for(ii=0;ii<naxes;ii++){
    fscanf(fp, "%lf",&phi2[i]);phi2[i] = phi2[i]*PI/180;
    fscanf(fp,"%lf",&theta2[i]);theta2[i] = theta2[i]*PI/180.
    fscanf(fp,"%lf",&r2[i]);
    fscanf(fp,"%lf",&w2);
}
    /*store cos theta, sin theta ...*/
    Cos1 = (double *) malloc(naxes*(sizeof(double)));
    Sin1 = (double *) malloc(naxes*(sizeof(double)));
    r1_2= (double *) malloc(naxes*(sizeof(double)));
    Cos2 = (double *) malloc(naxes2*(sizeof(double)));
    Sin2 = (double *) malloc(naxes2*(sizeof(double)));
    r2_2= (double *) malloc(naxes2*(sizeof(double)));
    
    for(ii=0;ii<naxes[1];ii++){
        Cos1[ii]=cos(theta[ii]);
        Sin1[ii]=sin(theta[ii]);
        r1_2[ii]=r[ii]*r[ii];
    }
    
    for(ii=0;ii<naxes2[1];ii++){
        Cos2[ii]=cos(theta2[ii]);
        Sin2[ii]=sin(theta2[ii]);
        r2_2[ii]=r2[ii]*r2[ii];
    }
//TODO
    if(opt==1)//DD,self
    {
        //calculate pairs for DD,DR,RR
        flag1=naxes;flag2=naxes2;
        for(ii=0;ii<flag1;ii++)//DD
            for(jj=ii+1;jj<flag2;jj++){
                    DD_total+=w[ii]*w2[jj];
                    buff=Separation2(Cos1[ii],Sin1[ii],r1_2[ii],r[ii],Cos2[jj],Sin2[jj],r2_2[jj],r2[ii],phi[ii]-phi2[jj]);
                    buff2=Orientation(r[ii],r2[jj],buff);
                    ibuff=floor(buff/5);
                    ibuff2=floor(buff2/0.01);
                    if(ibuff>=0 && ibuff<RESs && ibuff2>=0 && ibuff2<RESu){
                        Distribution[ibuff][ibuff2]+=(w[ii])*(w2[jj]);
                }
            }
        printf("DD total:%lf\n",DD_total);
        strcat(temp,"./NewData_subfiles/D");strcat(temp,argv[4]);strcat(temp,"D");strcat(temp,argv[5]);strcat(temp,".txt");
        fpXX=fopen(temp,"w+");
        //TODO
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

