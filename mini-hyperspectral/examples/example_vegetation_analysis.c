/**
 * @file    example_vegetation_analysis.c
 * @brief   Vegetation health monitoring using hyperspectral indices
 * L6: Vegetation health monitoring | L7: Precision agriculture (NASA AVIRIS)
 */
#include "../include/hyperspectral_core.h"
#include "../include/hyperspectral_spectroscopy.h"
#include "../include/hyperspectral_radiometry.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(void) {
    printf("=== Vegetation Health Monitoring ===\n\n");
    size_t nrows=6,ncols=8,nbands=30;
    hspec_datacube_t dc=hspec_datacube_alloc(nrows,ncols,nbands);
    double wl[30],fw[30];
    hspec_create_uniform_bands(400.0,1000.0,nbands,wl,fw);
    memcpy(dc.wavelengths,wl,nbands*sizeof(double));
    memcpy(dc.fwhm,fw,nbands*sizeof(double));
    strcpy(dc.sensor_name,"AVIRIS-NG-Sim");
    dc.spatial_resolution=3.0;
    printf("Sensor: %s | %zux%zu @ %.1fm GSD\n",dc.sensor_name,nrows,ncols,dc.spatial_resolution);
    printf("Generating vegetation scene...\n");
    for(size_t r=0;r<nrows;r++){
        for(size_t c=0;c<ncols;c++){
            int vt;
            if(r<2&&c<4)vt=0;else if(r<2)vt=1;else if(r<4&&c<3)vt=2;else if(r<4)vt=3;else vt=4;
            size_t base=c*nrows*nbands+r*nbands;
            for(size_t b=0;b<nbands;b++){
                double lam=wl[b],refl=0;
                if(vt==0){if(lam<700)refl=0.03+(lam-400)*0.0005;else refl=0.03+0.45*(1-exp(-(lam-700)/80));}
                else if(vt==1){if(lam<700)refl=0.05+(lam-400)*0.0006;else refl=0.05+0.20*(1-exp(-(lam-700)/60));}
                else if(vt==2){if(lam<700)refl=0.04+(lam-400)*0.0005;else refl=0.04+0.30*(1-exp(-(lam-700)/70));refl-=0.15*exp(-pow(lam-970,2)/(2*25*25));}
                else if(vt==3)refl=0.1+(lam-400)/1000*0.3;
                else{if(lam<700)refl=0.08+(lam-400)*0.001;else refl=0.15+(lam-700)*0.0002;}
                refl+=((double)rand()/RAND_MAX-0.5)*0.005;
                if(refl<0)refl=0;
                if(refl>1)refl=1;
                dc.data[base+b]=refl;
            }
        }
    }
    printf("Computing vegetation indices...\n\n");
    const char *tn[]={"Healthy","Stressed","WaterStress","BareSoil","Senescent"};
    size_t rp[5][2]={{0,0},{0,6},{2,0},{2,4},{4,2}};
    hspec_pixel_t px;px.nbands=nbands;px.reflectance=malloc(nbands*sizeof(double));
    printf("%-14s %8s %8s %8s %8s\n","Type","NDVI","SAVI","EVI","PRI");
    for(int t=0;t<5;t++){
        hspec_datacube_extract_pixel(&dc,rp[t][0],rp[t][1],&px);
        hspec_spectral_indices_t idx;
        hspec_compute_spectral_indices(&px,wl,&idx);
        printf("%-14s %8.3f %8.3f %8.3f %8.3f\n",tn[t],idx.ndvi,idx.savi,idx.evi,idx.pri);
    }
    double ndvi_sum=0,ndvi_min=1e9,ndvi_max=-1e9;
    size_t veg=0;
    for(size_t r=0;r<nrows;r++)for(size_t c=0;c<ncols;c++){
        hspec_datacube_extract_pixel(&dc,r,c,&px);
        hspec_spectral_indices_t idx;
        hspec_compute_spectral_indices(&px,wl,&idx);
        ndvi_sum+=idx.ndvi;
        if(idx.ndvi<ndvi_min)ndvi_min=idx.ndvi;
        if(idx.ndvi>ndvi_max)ndvi_max=idx.ndvi;
        if(idx.ndvi>0.3)veg++;
    }
    double ndvi_mean=ndvi_sum/(nrows*ncols);
    printf("\nScene NDVI: mean=%.3f min=%.3f max=%.3f vegetated=%zu/%zu\n",ndvi_mean,ndvi_min,ndvi_max,veg,nrows*ncols);
    if(ndvi_mean>0.6)printf("Overall: EXCELLENT crop health\n");
    else if(ndvi_mean>0.4)printf("Overall: MODERATE crop health\n");
    else printf("Overall: POOR crop health\n");
    free(px.reflectance);hspec_datacube_free(&dc);
    printf("\nVegetation analysis completed.\n");
    return 0;
}
