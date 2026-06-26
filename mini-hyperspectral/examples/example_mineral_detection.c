/**
 * @file    example_mineral_detection.c
 * @brief   Example: Mineral detection from hyperspectral imagery
 * @details AVIRIS-Cuprite mineral mapping simulation. Uses VCA endmember
 *          extraction, FCLS abundance estimation, and SAM classification.
 * L6: Canonical Problem ? Mineral mapping from AVIRIS data
 * L7: Application ? Mineral exploration (USGS/NASA AVIRIS)
 */
#include "../include/hyperspectral_core.h"
#include "../include/hyperspectral_unmixing.h"
#include "../include/hyperspectral_classification.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(void) {
    printf("=== Mineral Detection from Hyperspectral Imagery ===\n");
    printf("    AVIRIS-Cuprite (Nevada) Simulation\n\n");
    size_t nrows=6,ncols=8,nbands=40;
    hspec_datacube_t dc=hspec_datacube_alloc(nrows,ncols,nbands);
    double wl[40],fw[40];
    hspec_create_uniform_bands(2000.0,2500.0,nbands,wl,fw);
    memcpy(dc.wavelengths,wl,nbands*sizeof(double));
    memcpy(dc.fwhm,fw,nbands*sizeof(double));
    strcpy(dc.sensor_name,"AVIRIS-Classic-Sim");
    dc.spatial_resolution=17.0;
    const char *minerals[]={"Kaolinite","Alunite","Calcite","Muscovite"};
    double em[4*40];
    for(size_t i=0;i<nbands;i++){
        double lam=wl[i];
        em[0*nbands+i]=0.6-0.3*exp(-pow(lam-2160,2)/(2*25*25))-0.4*exp(-pow(lam-2210,2)/(2*20*20));
        em[1*nbands+i]=0.7-0.5*exp(-pow(lam-2180,2)/(2*30*30));
        em[2*nbands+i]=0.5-0.6*exp(-pow(lam-2340,2)/(2*20*20));
        em[3*nbands+i]=0.55-0.35*exp(-pow(lam-2210,2)/(2*22*22));
    }
    printf("Generating %zux%zu scene with %zu bands...\n",nrows,ncols,nbands);
    for(size_t r=0;r<nrows;r++){
        for(size_t c=0;c<ncols;c++){
            double a[4]={0,0,0,0};
            if(r<2&&c<4){a[0]=0.7;a[1]=0.2;a[3]=0.1;}
            else if(r<2){a[0]=0.1;a[1]=0.7;a[2]=0.1;a[3]=0.1;}
            else if(r<4&&c<3){a[2]=0.6;a[3]=0.2;a[0]=0.1;a[1]=0.1;}
            else if(r<4){a[3]=0.6;a[0]=0.2;a[1]=0.1;a[2]=0.1;}
            else{a[0]=a[1]=a[2]=a[3]=0.25;}
            size_t base=c*nrows*nbands+r*nbands;
            for(size_t b=0;b<nbands;b++){
                double val=0;
                for(size_t k=0;k<4;k++)val+=a[k]*em[k*nbands+b];
                val+=((double)rand()/RAND_MAX-0.5)*0.02;
                dc.data[base+b]=val;
            }
        }
    }
    printf("Running VCA endmember extraction...\n");
    hspec_vca_params_t vp;vp.n_endmembers=4;vp.max_iterations=5;vp.convergence_tol=1e-4;vp.use_snr_proj=1;
    hspec_unmixing_result_t vr=hspec_vca_extract(&dc,vp);
    printf("VCA extracted %zu endmembers\n",vr.n_endmembers);
    hspec_fcls_params_t fp;fp.max_iterations=100;fp.tolerance=1e-5;fp.regularization=0;fp.enforce_sum_to_one=1;
    hspec_unmixing_result_t ur=hspec_unmix_scene(&dc,vr.endmembers,vr.n_endmembers,fp);
    printf("RMSE: %.6f  SAD: %.4f deg  Sum-to-one err: %.6f\n",ur.rmse,ur.mean_sad*180/M_PI,ur.sum_to_one_error);
    hspec_pixel_t px;px.nbands=nbands;px.reflectance=malloc(nbands*sizeof(double));
    size_t cnt[4]={0,0,0,0};
    for(size_t r=0;r<nrows;r++)for(size_t c=0;c<ncols;c++){
        hspec_datacube_extract_pixel(&dc,r,c,&px);
        hspec_sam_result_t sr=hspec_sam_classify(&px,em,4,nbands);
        if(sr.best_class<4)cnt[sr.best_class]++;
        free(sr.angles_rad);
    }
    printf("\nClassification:\n");
    for(size_t k=0;k<4;k++)printf("  %-14s: %zu pixels\n",minerals[k],cnt[k]);
    free(px.reflectance);
    hspec_unmixing_result_free(&vr);hspec_unmixing_result_free(&ur);
    hspec_datacube_free(&dc);
    printf("\nMineral detection completed.\n");
    return 0;
}
