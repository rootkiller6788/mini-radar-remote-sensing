#include "../include/cfar_detector.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int tr=0,tp=0;
#define T(n) do{tr++;printf("  %s... ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

int main(void){
    printf("=== cfar tests ===\n");

    T("cfar_config_init");
    cfar_config_t cfg;
    assert(cfar_config_init(CFAR_CA,16,4,1e-4,&cfg)==0);
    assert(cfg.num_reference_cells==16);
    assert(cfg.threshold_factor>0.0);
    P();

    T("ca_cfar_detect");
    double data[200];
    for(int i=0;i<200;i++) data[i]=1.0;
    data[100]=20.0;
    cfar_detection_t dets[200];
    size_t nd;
    assert(ca_cfar_detect(data,200,&cfg,dets,&nd)==0);
    assert(nd>50);
    int found=0;
    for(size_t i=0;i<nd;i++) if(dets[i].cell_index==100&&dets[i].detected) found=1;
    assert(found);
    P();

    T("go_cfar_detect");
    for(int i=0;i<50;i++) data[i]=5.0;
    assert(go_cfar_detect(data,200,&cfg,dets,&nd)==0);
    P();

    T("so_cfar_detect");
    assert(so_cfar_detect(data,200,&cfg,dets,&nd)==0);
    P();

    T("os_cfar_detect");
    cfg.type=CFAR_OS;
    cfg.os_rank=12;
    cfar_config_init(CFAR_OS,16,4,1e-4,&cfg);
    assert(os_cfar_detect(data,200,&cfg,dets,&nd)==0);
    P();

    T("cfar_2d_scan");
    double *rd = calloc(32*32,sizeof(double));
    rd[16+16*32]=50.0;
    cfar_detection_map_t map;
    cfar_detection_map_alloc(32,32,&map);
    cfg.type=CFAR_CA;
    cfar_config_init(CFAR_CA,8,2,1e-3,&cfg);
    assert(cfar_2d_scan(rd,32,32,&cfg,&map)==0);
    assert(map.num_detections>0);
    cfar_detection_map_free(&map);
    free(rd);
    P();

    printf("\n=== %d/%d passed ===\n",tp,tr);
    return (tp==tr)?0:1;
}
