/**
 * @file    example_spectral_unmixing.c
 * @brief   Subpixel land cover mapping via spectral unmixing
 * L6: Subpixel target detection | L7: Urban impervious surface mapping
 * L8: Fan bilinear mixing model
 */
#include "../include/hyperspectral_core.h"
#include "../include/hyperspectral_unmixing.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(void) {
    printf("=== Spectral Unmixing: Subpixel Land Cover ===\n\n");
    size_t nrows=8,ncols=10,nbands=30;
    hspec_datacube_t dc=hspec_datacube_alloc(nrows,ncols,nbands);
    double wl[30],fw[30];
    hspec_create_uniform_bands(400.0,900.0,nbands,wl,fw);
    memcpy(dc.wavelengths,wl,nbands*sizeof(double));
    memcpy(dc.fwhm,fw,nbands*sizeof(double));
    strcpy(dc.sensor_name,"HyMap-Sim");
    dc.spatial_resolution=4.0;
    printf("Sensor: %s | %zux%zu @ %.1fm GSD\n\n",dc.sensor_name,nrows,ncols,dc.spatial_resolution);
    /* Land cover endmembers: Vegetation, BareSoil, Asphalt, Concrete, Water */
    double em[5*30];
    for(size_t b=0;b<nbands;b++){
        double lam=wl[b];
        if(lam<700)em[0*30+b]=0.03+(lam-400)*0.001;else em[0*30+b]=0.03+0.40*(1-exp(-(lam-700)/70));
        em[1*30+b]=0.10+(lam-400)/900*0.25;
        em[2*30+b]=0.04+(lam-400)/900*0.02;
        em[3*30+b]=0.30;
        em[4*30+b]=0.06-(lam-400)/900*0.04;
    }
    printf("Scene Layout:\n");
    for(size_t r=0;r<nrows;r++){
        printf("  ");
        for(size_t c=0;c<ncols;c++){
            double a[5]={0,0,0,0,0};char ch;
            if(r<2&&c<4){a[0]=0.6;a[1]=0.3;a[3]=0.1;ch='V';}
            else if(r<2){a[3]=0.5;a[2]=0.3;a[1]=0.2;ch='C';}
            else if(r<5&&c<3){a[4]=0.8;a[1]=0.1;a[0]=0.1;ch='W';}
            else if(r<5&&c<7){a[2]=0.5;a[3]=0.3;a[0]=0.2;ch='U';}
            else if(r<5){a[1]=0.7;a[0]=0.2;a[3]=0.1;ch='S';}
            else{a[0]=0.4;a[1]=0.3;a[3]=0.2;a[2]=0.1;ch='M';}
            printf("%c ",ch);
            size_t base=c*nrows*nbands+r*nbands;
            for(size_t b=0;b<nbands;b++){
                double val=0;
                for(size_t k=0;k<5;k++)val+=a[k]*em[k*30+b];
                val+=((double)rand()/RAND_MAX-0.5)*0.01;
                dc.data[base+b]=val;
            }
        }
        printf("\n");
    }
    printf("\nV=Vegetation C=Concrete W=Water U=Urban S=Suburban M=Mixed\n");
    printf("\nN-FINDR endmember extraction...\n");
    hspec_nfindr_params_t np;np.n_endmembers=5;np.max_iterations=10;np.n_random_starts=3;
    hspec_unmixing_result_t nr=hspec_nfindr_extract(&dc,np);
    printf("Extracted %zu endmembers\n",nr.n_endmembers);
    hspec_fcls_params_t fp;fp.max_iterations=100;fp.tolerance=1e-5;fp.regularization=0;fp.enforce_sum_to_one=1;
    hspec_unmixing_result_t ur=hspec_unmix_scene(&dc,nr.endmembers,nr.n_endmembers,fp);
    printf("RMSE: %.5f  SAD: %.3f deg\n",ur.rmse,ur.mean_sad*180/M_PI);
    printf("Sum-to-one err: %.5f  Non-neg err: %.5f\n",ur.sum_to_one_error,ur.nonneg_error);
    printf("\nBilinear mixing model test:\n");
    double ta[5]={0.3,0.3,0.1,0.2,0.1};
    double tg[10];for(int i=0;i<10;i++)tg[i]=0.5;
    double *bs=malloc(nbands*sizeof(double));
    hspec_fan_mixing_model(em,ta,tg,5,nbands,bs);
    printf("  Fan bilinear model applied (gamma=0.5)\n");
    free(bs);
    hspec_unmixing_result_free(&nr);hspec_unmixing_result_free(&ur);
    hspec_datacube_free(&dc);
    printf("\nSpectral unmixing completed.\n");
    return 0;
}
