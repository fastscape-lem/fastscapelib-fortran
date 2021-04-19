#include <stdio.h>
#include <stdlib.h>
#include <time.h>
void fastscape_init_(int* ierr);
void fastscape_set_nx_ny_(int* nx,int* ny);
void fastscape_set_xl_yl_(double* xl,double* yl);
void fastscape_set_dt_(double* dt);
void fastscape_set_erosional_parameters_(double* kf1_arr,double* kf2,double* m,double* n,double* kd1_arr,double* kd2,double* g1,double* g2, double* p_flow_dir_exp);
//void fastscape_set_precipitation_rate_(double* preci_rate);
void fastscape_set_marine_parameters_(double* sealevel,double* poro1,double* poro2,double* z1,double* z2, double* ratio,double* L,double* kds1,double* kds2);
void fastscape_set_bc_(int* bc);
void fastscape_init_h_(double *h);
void fastscape_setup_();
void fastscape_copy_h_(double *h);
void fastscape_set_u_(double *u);
void fastscape_execute_step_();
void fastscape_get_step_(int* istep);
void fastscape_debug_();
void fastscape_copy_total_erosion_(double* etot);
void fastscape_copy_erosion_rate_(double* erate);
void fastscape_copy_basement_(double* b);
void fastscape_vtk_(double* etot, double* vexp);
void fastscape_view_();

int main(void){
    srand(time(NULL));
    int  nx,ny,istep,nstep,nn,nfreq,bc,nr_fields;
    // double* u[],ux[],uy[],h[],b[],etot[],erate[],a[],chi[],catchment[],sedflux[],sedflux_shore[];
    double kf1,kf2,m,n,kd1,kd2,g1,g2,preci_rate,p_flow_dir_exp;
    double xl,yl,dx,dy,dt;
    double sealevel, poro1, poro2, ratio, L, kds1, kds2, z1, z2;
    double vexp;
    int i,j,ij;
    int ierr;

    bc = 1010;
    nx=201;
    ny=201;
    nn=nx*ny;
    xl=200.e3;
    yl=200.e3;
    dx=xl/(nx-1);
    dy=yl/(ny-1);
    dt=1.e3;
    kf1=1.e-5;
    kf2=2.e-5;
    // !kf2=kf1;
    m=0.4e0;
    n=1.e0;
    kd1=1.e-2;
    kd2=5.e-2;
    kd1 = 0.e0;
    kd2 = 0.e0;
    // !kd2=kd1;
    g1=1.e0;
    g2=1.e0;
    p_flow_dir_exp = -2.e0;


    preci_rate = 1.e0; // precipitation rate
    sealevel = 0.e0;
    poro1 = 0.0e0;
    poro2 = 0.0e0;
    ratio = 0.5e0;
    L = 0.5e2;
    kds1 = 2.e2;
    kds2 = 1.e2;
    z1 = 1.e3;
    z2 = 1.e3;

    // allocate memory for working arrays
    //double *u creates pointer to start of array. (double*) typecasts the allocated array
    // malloc returns void pointer. Thus we need typecasting. malloc does not assign a value
    // using calloc(nn,sizeof(double)) would initialise array with 0 values.
    double *u = (double*)malloc( nn * sizeof( double ) );
    double *ux= (double*)malloc( nn * sizeof( double ) );
    double *uy= (double*)malloc( nn * sizeof( double ) );
    double *h= (double*)malloc( nn * sizeof( double ) );
    double *b= (double*)malloc( nn * sizeof( double ) );
    double *etot= (double*)malloc( nn * sizeof( double ) );
    double *erate= (double*)malloc( nn * sizeof( double ) );
    double *a= (double*)malloc( nn * sizeof( double ) );
    double *chi= (double*)malloc( nn * sizeof( double ) );
    double *catchment= (double*)malloc( nn * sizeof( double ) );
    double *sedflux= (double*)malloc( nn * sizeof( double ) );
    double *sedflux_shore= (double*)malloc( nn * sizeof( double ) );
    double *kf1_arr = (double*)malloc( nn * sizeof( double ) );
    double *kd1_arr =(double*)malloc( nn * sizeof( double ) );
    nr_fields=2;
    // double field[nr_fields][nn];
    // double **field = (double**)malloc(sizeof(double*)*nr_fields + sizeof(double*)*nr_fields*nn);
    // double **field = (double **)malloc(nr_fields * sizeof(double*));
    // for (i=0; i<nr_fields; i++){
    //      field[i] = (double *)malloc(nn * sizeof(double));
    //  }
     // printf("%d\n",field[0] );
     // printf("%d\n",field[1] );
     fastscape_init_(&ierr);
     fastscape_setup_();
     fastscape_set_nx_ny_(&nx,&ny);
     fastscape_set_xl_yl_(&xl,&yl);
     fastscape_set_dt_(&dt);
     fastscape_set_marine_parameters_(&sealevel, &poro1, &poro2, &z1, &z2, &ratio, &L, &kds1, &kds2);
     fastscape_set_bc_(&bc);
    for (i=0;i <nn; i++){
        h[i] = rand() % 10 + 0;
        u[i] = 0.0;
        kf1_arr[i] = kf1;
        kd1_arr[i] = kd1;
    }
    for(j=0;j<ny;j++){
        for(i=0;i<ny;i++){
            ij=(j)*nx+i;
            if (j<ny/2){
                h[ij]=h[ij]-200.e0;
            }
            if (j>ny/2){
                h[ij]=h[ij]+1000.e0;
            }
        }
    }
    fastscape_set_erosional_parameters_(kf1_arr,&kf2,&m,&n,kd1_arr,&kd2,&g1,&g2,&p_flow_dir_exp);
    // for (i=0;i <nn; i++){
    //     printf("Element #%d: %f\n",i,h[i]);
    // }

    for(j=0;j<ny;j++){
        for(i=0;i<ny;i++){
            ij=(j)*nx+i;
            u[ij]=0.0;
        }
    }

    fastscape_init_h_(h);
    fastscape_set_u_(u);
    fastscape_copy_h_(h);
    // for (i=0;i <nn; i++){
    //     printf("Element #%d: %f\n",i,h[i]);
    // }
    nstep=500;
    nfreq=100;
    istep=nstep;
    fastscape_view_();

    while (istep<=nstep){
        fastscape_execute_step_();
        fastscape_get_step_(&istep);
        fastscape_debug_();
        if (istep%nfreq==0){
            printf("%d\n",istep );
            fastscape_copy_total_erosion_(etot);
            fastscape_copy_erosion_rate_(erate);
            fastscape_copy_h_(h);
            fastscape_copy_basement_(b);
            // for (i=0;i <nn; i++){
            //     printf("Element #%d: %f\n",i,h[i]);
            // }
            vexp = 2.0;
            fastscape_vtk_(etot,&vexp);
            vexp = -2.0;
            fastscape_vtk_(erate,&vexp);
        }
    }
    // fastscape_copy_h_(h);
    // for (i=0;i <nn; i++){
    //     printf("Element #%d: %f\n",i,h[i]);
    // }
    //  free memory again
    free(u);
    free(ux);
    free(uy);
    free(h);
    free(b);
    free(etot);
    free(erate);
    free(a);
    free(chi);
    free(catchment);
    free(sedflux);
    free(sedflux_shore);
    free(kf1_arr);
    free(kd1_arr);

}
