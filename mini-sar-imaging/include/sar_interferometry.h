
#ifndef SAR_INTERFEROMETRY_H
#define SAR_INTERFEROMETRY_H
#include "sar_core.h"
#include "sar_geometry.h"

/* --- L1: Interferometric Definitions --- */
typedef struct {
    sar_image_t *master;
    sar_image_t *slave;
    double baseline_m;
    double baseline_parallel_m;
    double baseline_perp_m;
    double baseline_angle_rad;
    double height_ambiguity_m;
    double critical_baseline_m;
} sar_insar_pair_t;

/* --- L1: Coherence --- */
typedef struct {
    size_t nrows, ncols;
    double **coherence;
    double mean_coherence;
    double **interferometric_phase;
    double **phase_stddev;
} sar_coherence_map_t;

sar_insar_pair_t *sar_insar_pair_alloc(void);
void sar_insar_pair_free(sar_insar_pair_t *pair);
void sar_insar_set_baseline(sar_insar_pair_t *pair, double B, double alpha, double H, double R0, double lambda);

sar_coherence_map_t *sar_coherence_alloc(size_t nrows, size_t ncols);
void sar_coherence_free(sar_coherence_map_t *cmap);

/* L5: Coherence estimation via spatial averaging */
void sar_coherence_estimate(const sar_image_t *master, const sar_image_t *slave, sar_coherence_map_t *cmap, size_t win_rows, size_t win_cols);

/* L5: Interferometric phase computation */
void sar_interferogram_compute(const sar_image_t *master, const sar_image_t *slave, double **ifgram, size_t nrows, size_t ncols);

/* L5: Flat-earth phase removal */
void sar_flat_earth_removal(double **ifgram, size_t nrows, size_t ncols, double range_spacing, double lambda, double H, double look_angle);

/* L5: Phase unwrapping - Goldstein branch cut */
int sar_phase_unwrap_goldstein(double **wrapped_phase, double **unwrapped_phase, size_t nrows, size_t ncols);

/* L5: Phase unwrapping - quality guided */
int sar_phase_unwrap_quality(double **wrapped_phase, double **unwrapped_phase, size_t nrows, size_t ncols);

/* L5: Phase to height conversion */
void sar_phase_to_height(const double **unwrapped_phase, size_t nrows, size_t ncols, double B_perp, double R0, double lambda, double H, double theta, double **dem);

/* L7: DInSAR - differential interferometry for displacement */
int sar_dinsar_displacement(const sar_image_t *master, const sar_image_t *slave, double **displacement_m, size_t nrows, size_t ncols, double lambda, const double **topo_phase);

/* L7: Phase filtering - Goldstein filter */
void sar_goldstein_filter(double **ifgram, size_t nrows, size_t ncols, double alpha, int patch_size);

#endif
