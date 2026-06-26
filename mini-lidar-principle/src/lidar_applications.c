/**
 * @file    lidar_applications.c
 * @brief   LiDAR applications — DEM, forestry, building extraction, perception
 *
 * Knowledge covered:
 *   L6: DEM/DSM/CHM generation, building extraction, ground filtering
 *   L7: Forestry metrics (canopy height, cover, LAI, biomass),
 *       autonomous vehicle object detection,
 *       atmospheric LiDAR (aerosol backscatter, boundary layer height)
 *   L8: FMCW LiDAR beat frequency processing
 *
 * Reference:
 *   - Vosselman & Maas, *Airborne and Terrestrial Laser Scanning*, 2010.
 *   - Shan & Toth, *Topographic Laser Ranging and Scanning*, 2nd ed., 2018.
 *   - Lefsky et al., *BioScience* 52(1), pp.19-30, 2002.
 *   - Amann et al., *Optical Engineering* 40(1), pp.10-19, 2001.
 *   - Klett, J.D., *Applied Optics* 20(2), pp.211-220, 1981.
 */

#include "lidar_applications.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ─── Internal helpers ──────────────────────────────────────────────────── */


static int cmp_double(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) ? 1 : ((da < db) ? -1 : 0);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L6: DEM / DTM generation
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_dem_config_from_scan(const lidar_scan_t *scan,
                                double cell_size,
                                lidar_dem_config_t *config) {
    if (!scan || !config || cell_size <= 0.0 || scan->num_points == 0) return -1;

    double min_xyz[3], max_xyz[3];
    if (lidar_scan_bounding_box(scan, min_xyz, max_xyz) != 0) return -1;

    config->cell_size = cell_size;
    config->x_min = min_xyz[0];
    config->x_max = max_xyz[0];
    config->y_min = min_xyz[1];
    config->y_max = max_xyz[1];

    config->cols = (size_t)ceil((max_xyz[0] - min_xyz[0]) / cell_size) + 1;
    config->rows = (size_t)ceil((max_xyz[1] - min_xyz[1]) / cell_size) + 1;
    if (config->cols == 0) config->cols = 1;
    if (config->rows == 0) config->rows = 1;

    return 0;
}

int lidar_dem_generate(const lidar_scan_t *scan,
                        const lidar_dem_config_t *config,
                        lidar_dem_t *dem) {
    if (!scan || !config || !dem || scan->num_points == 0) return -1;

    size_t total_cells = config->cols * config->rows;
    if (total_cells == 0 || total_cells > 100000000) return -1;

    dem->data = (double*)malloc(total_cells * sizeof(double));
    if (!dem->data) return -1;

    dem->cols = config->cols;
    dem->rows = config->rows;
    dem->cell_size = config->cell_size;
    dem->x_min = config->x_min;
    dem->y_min = config->y_min;
    dem->nodata_value = NAN;

    /* Per-cell min Z accumulator (for DTM) */
    double *z_min = (double*)malloc(total_cells * sizeof(double));
    char *has_data = (char*)calloc(total_cells, 1);
    if (!z_min || !has_data) {
        free(z_min); free(has_data); free(dem->data);
        return -1;
    }
    for (size_t i = 0; i < total_cells; i++) z_min[i] = 1e100;

    for (size_t i = 0; i < scan->num_points; i++) {
        double x = scan->points[i].x;
        double y = scan->points[i].y;
        double z = scan->points[i].z;

        size_t col = (size_t)((x - config->x_min) / config->cell_size);
        size_t row = (size_t)((y - config->y_min) / config->cell_size);
        if (col < config->cols && row < config->rows) {
            size_t idx = col + row * config->cols;

            /* For DTM, use ground-classified (2) points if available,
               otherwise use all points */
            if (scan->points[i].class_label == 2 || !has_data[idx]) {
                if (z < z_min[idx]) {
                    z_min[idx] = z;
                    has_data[idx] = 1;
                }
            }
            /* Always update if no data yet */
            if (!has_data[idx]) {
                z_min[idx] = z;
                has_data[idx] = 1;
            }
        }
    }

    /* Fill DEM */
    for (size_t idx = 0; idx < total_cells; idx++) {
        dem->data[idx] = has_data[idx] ? z_min[idx] : dem->nodata_value;
    }

    free(z_min);
    free(has_data);
    return 0;
}

int lidar_dsm_generate(const lidar_scan_t *scan,
                        const lidar_dem_config_t *config,
                        lidar_dem_t *dsm) {
    if (!scan || !config || !dsm || scan->num_points == 0) return -1;

    size_t total_cells = config->cols * config->rows;
    if (total_cells == 0 || total_cells > 100000000) return -1;

    dsm->data = (double*)malloc(total_cells * sizeof(double));
    if (!dsm->data) return -1;

    dsm->cols = config->cols;
    dsm->rows = config->rows;
    dsm->cell_size = config->cell_size;
    dsm->x_min = config->x_min;
    dsm->y_min = config->y_min;
    dsm->nodata_value = NAN;

    double *z_max = (double*)malloc(total_cells * sizeof(double));
    char *has_data = (char*)calloc(total_cells, 1);
    if (!z_max || !has_data) {
        free(z_max); free(has_data); free(dsm->data);
        return -1;
    }
    for (size_t i = 0; i < total_cells; i++) z_max[i] = -1e100;

    for (size_t i = 0; i < scan->num_points; i++) {
        double x = scan->points[i].x;
        double y = scan->points[i].y;
        double z = scan->points[i].z;

        size_t col = (size_t)((x - config->x_min) / config->cell_size);
        size_t row = (size_t)((y - config->y_min) / config->cell_size);
        if (col < config->cols && row < config->rows) {
            size_t idx = col + row * config->cols;
            /* DSM uses first returns (return_num == 1) or all if no classification */
            if (scan->points[i].return_num == 1 || !has_data[idx]) {
                if (z > z_max[idx]) z_max[idx] = z;
                has_data[idx] = 1;
            }
            if (!has_data[idx]) {
                z_max[idx] = z;
                has_data[idx] = 1;
            }
        }
    }

    for (size_t idx = 0; idx < total_cells; idx++) {
        dsm->data[idx] = has_data[idx] ? z_max[idx] : dsm->nodata_value;
    }

    free(z_max);
    free(has_data);
    return 0;
}

int lidar_chm_compute(const lidar_dem_t *dsm,
                       const lidar_dem_t *dtm,
                       lidar_dem_t *chm) {
    if (!dsm || !dtm || !chm) return -1;
    if (dsm->cols != dtm->cols || dsm->rows != dtm->rows) return -1;

    size_t total = dsm->cols * dsm->rows;
    chm->data = (double*)malloc(total * sizeof(double));
    if (!chm->data) return -1;

    chm->cols = dsm->cols;
    chm->rows = dsm->rows;
    chm->cell_size = dsm->cell_size;
    chm->x_min = dsm->x_min;
    chm->y_min = dsm->y_min;
    chm->nodata_value = NAN;

    for (size_t i = 0; i < total; i++) {
        if (!isnan(dsm->data[i]) && !isnan(dtm->data[i])) {
            chm->data[i] = dsm->data[i] - dtm->data[i];
            if (chm->data[i] < 0.0) chm->data[i] = 0.0;
        } else {
            chm->data[i] = chm->nodata_value;
        }
    }

    return 0;
}

void lidar_dem_free(lidar_dem_t *dem) {
    if (!dem) return;
    free(dem->data);
    memset(dem, 0, sizeof(*dem));
}

int lidar_dem_hillshade(const lidar_dem_t *dem,
                         double sun_azimuth, double sun_zenith,
                         uint8_t *hillshade) {
    if (!dem || !dem->data || !hillshade) return -1;

    size_t cols = dem->cols, rows = dem->rows;
    double cell_size = dem->cell_size;

    for (size_t row = 0; row < rows; row++) {
        for (size_t col = 0; col < cols; col++) {
            size_t idx = col + row * cols;
            double z = dem->data[idx];
            if (isnan(z)) {
                hillshade[idx] = 0;
                continue;
            }

            /* Compute slope and aspect from 8 neighbors */
            double z_left  = (col > 0) ? dem->data[(col-1) + row * cols] : z;
            double z_right = (col + 1 < cols) ? dem->data[(col+1) + row * cols] : z;
            double z_up    = (row > 0) ? dem->data[col + (row-1) * cols] : z;
            double z_down  = (row + 1 < rows) ? dem->data[col + (row+1) * cols] : z;

            if (isnan(z_left)) z_left = z;
            if (isnan(z_right)) z_right = z;
            if (isnan(z_up)) z_up = z;
            if (isnan(z_down)) z_down = z;

            double dz_dx = (z_right - z_left) / (2.0 * cell_size);
            double dz_dy = (z_down - z_up) / (2.0 * cell_size);

            /* Slope */
            double slope = atan(sqrt(dz_dx * dz_dx + dz_dy * dz_dy));

            /* Aspect */
            double aspect = atan2(dz_dy, -dz_dx);
            if (aspect < 0.0) aspect += 2.0 * M_PI;

            /* Hillshade = 255 * (cos(zenith)*cos(slope) + sin(zenith)*sin(slope)*cos(azimuth-aspect)) */
            double hs = cos(sun_zenith) * cos(slope)
                       + sin(sun_zenith) * sin(slope) * cos(sun_azimuth - aspect);
            if (hs < 0.0) hs = 0.0;
            if (hs > 1.0) hs = 1.0;

            hillshade[idx] = (uint8_t)(hs * 255.0);
        }
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L6: Building extraction
 * ═══════════════════════════════════════════════════════════════════════════ */

size_t lidar_extract_buildings(const lidar_scan_t *scan,
                                 double min_height,
                                 double min_area,
                                 double cluster_dist,
                                 lidar_building_footprint_t *footprints,
                                 size_t max_footprints) {
    if (!scan || !footprints || max_footprints == 0 || scan->num_points == 0) return 0;

    /* Find minimum Z (ground reference) */
    double z_min = scan->points[0].z;
    for (size_t i = 1; i < scan->num_points; i++) {
        if (scan->points[i].z < z_min) z_min = scan->points[i].z;
    }

    /* Filter potential building points */
    size_t n_bldg = 0;
    size_t *bldg_indices = (size_t*)malloc(scan->num_points * sizeof(size_t));
    if (!bldg_indices) return 0;

    for (size_t i = 0; i < scan->num_points; i++) {
        double height = scan->points[i].z - z_min;
        if (height >= min_height) {
            bldg_indices[n_bldg++] = i;
        }
    }
    if (n_bldg < 3) { free(bldg_indices); return 0; }

    /* Simple connected-components clustering (flood-fill) */
    char *visited = (char*)calloc(n_bldg, 1);
    if (!visited) { free(bldg_indices); return 0; }

    size_t n_found = 0;
    double cluster_dist2 = cluster_dist * cluster_dist;

    for (size_t seed = 0; seed < n_bldg && n_found < max_footprints; seed++) {
        if (visited[seed]) continue;

        /* Start a new cluster (BFS with static array as queue) */
        size_t *queue = (size_t*)malloc(n_bldg * sizeof(size_t));
        if (!queue) break;
        size_t q_head = 0, q_tail = 0;
        queue[q_tail++] = seed;
        visited[seed] = 1;

        double x_sum = 0, y_sum = 0, z_sum = 0;
        double xmin = 1e100, xmax = -1e100, ymin = 1e100, ymax = -1e100;
        size_t cluster_size = 0;

        while (q_head < q_tail) {
            size_t ci = queue[q_head++];
            size_t pi = bldg_indices[ci];

            double x = scan->points[pi].x;
            double y = scan->points[pi].y;
            double z = scan->points[pi].z;

            x_sum += x; y_sum += y; z_sum += z;
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
            cluster_size++;

            /* Find neighbors */
            for (size_t nj = 0; nj < n_bldg; nj++) {
                if (visited[nj]) continue;
                size_t pj = bldg_indices[nj];
                double dx = scan->points[pj].x - x;
                double dy = scan->points[pj].y - y;
                if (dx * dx + dy * dy <= cluster_dist2) {
                    visited[nj] = 1;
                    queue[q_tail++] = nj;
                }
            }
        }
        free(queue);

        /* Validate cluster */
        double area = (xmax - xmin) * (ymax - ymin);
        if (area < min_area || cluster_size < 3) continue;

        /* Fill footprint */
        footprints[n_found].x_min = xmin;
        footprints[n_found].x_max = xmax;
        footprints[n_found].y_min = ymin;
        footprints[n_found].y_max = ymax;
        footprints[n_found].area = area;
        footprints[n_found].mean_elevation = z_sum / (double)cluster_size;
        footprints[n_found].height = footprints[n_found].mean_elevation - z_min;
        footprints[n_found].num_points = cluster_size;
        n_found++;
    }

    free(visited);
    free(bldg_indices);
    return n_found;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L7: Forestry metrics
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_forestry_metrics(const lidar_scan_t *scan,
                             double ground_elevation,
                             lidar_forestry_metrics_t *metrics) {
    if (!scan || !metrics || scan->num_points == 0) return -1;
    memset(metrics, 0, sizeof(*metrics));

    /* Compute height above ground for each point */
    size_t n_veg = 0, n_ground = 0, n_first = 0, n_first_above = 0;
    double max_h = -1e100;
    double sum_h = 0.0, sum_h2 = 0.0;

    double *heights = (double*)malloc(scan->num_points * sizeof(double));
    if (!heights) return -1;

    for (size_t i = 0; i < scan->num_points; i++) {
        double h = scan->points[i].z - ground_elevation;
        if (h < 0.0) h = 0.0;
        heights[i] = h;

        if (scan->points[i].class_label == 2) {
            n_ground++;
        } else {
            n_veg++;
            if (h > max_h) max_h = h;
            sum_h += h;
            sum_h2 += h * h;
        }

        if (scan->points[i].return_num == 1) {
            n_first++;
            if (h > 2.0) n_first_above++;
        }
    }

    metrics->canopy_height_max = (max_h > -1e99) ? max_h : 0.0;

    if (n_veg > 0) {
        metrics->canopy_height_mean = sum_h / (double)n_veg;
        double var = sum_h2 / (double)n_veg - metrics->canopy_height_mean * metrics->canopy_height_mean;
        if (var < 0.0) var = 0.0;
        metrics->canopy_height_std = sqrt(var);
    }

    if (n_first > 0) {
        metrics->canopy_cover = (double)n_first_above / (double)n_first;
    }

    /* Leaf Area Index via Beer-Lambert: LAI = -cos(theta) * ln(P_gap) */
    double gap_fraction = 1.0 - metrics->canopy_cover;
    if (gap_fraction < 0.01) gap_fraction = 0.01;
    if (gap_fraction > 0.99) gap_fraction = 0.99;
    double mean_scan_angle = 0.0; /* near-nadir assumption */
    metrics->leaf_area_index = -cos(mean_scan_angle) * log(gap_fraction);

    /* Percentile heights */
    if (n_veg > 0) {
        double *veg_heights = (double*)malloc(n_veg * sizeof(double));
        if (veg_heights) {
            size_t vi = 0;
            for (size_t i = 0; i < scan->num_points; i++) {
                if (scan->points[i].class_label != 2) {
                    veg_heights[vi++] = heights[i];
                }
            }
            qsort(veg_heights, n_veg, sizeof(double), cmp_double);

            metrics->percentile_h50 = veg_heights[(size_t)(0.50 * (double)(n_veg - 1))];
            metrics->percentile_h75 = veg_heights[(size_t)(0.75 * (double)(n_veg - 1))];
            metrics->percentile_h90 = veg_heights[(size_t)(0.90 * (double)(n_veg - 1))];
            metrics->percentile_h95 = veg_heights[(size_t)(0.95 * (double)(n_veg - 1))];
            metrics->percentile_h99 = veg_heights[(size_t)(0.99 * (double)(n_veg - 1))];

            free(veg_heights);
        }
    }

    /* Above-ground biomass (AGB) estimate — allometric equation
       AGB = a * H_mean^b * CC^c
       Typical a=2.5, b=1.5, c=0.8 for temperate mixed forest */
    double a = 2.5, b = 1.5, c = 0.8;
    double agb = a * pow(metrics->canopy_height_mean, b)
                 * pow(metrics->canopy_cover, c);
    metrics->biomass_estimate = agb; /* Mg/ha — assumes plot ~1 ha */

    metrics->total_returns = scan->num_points;
    metrics->ground_returns = n_ground;

    free(heights);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L7: Object detection for autonomous driving
 * ═══════════════════════════════════════════════════════════════════════════ */

size_t lidar_detect_objects_euclidean(const lidar_scan_t *scan,
                                        double cluster_dist,
                                        size_t min_cluster_size,
                                        size_t max_cluster_size,
                                        lidar_object_detection_t *detections,
                                        size_t max_detections) {
    if (!scan || !detections || max_detections == 0 || scan->num_points == 0) return 0;

    double cd2 = cluster_dist * cluster_dist;
    char *visited = (char*)calloc(scan->num_points, 1);
    if (!visited) return 0;

    size_t n_detected = 0;

    for (size_t seed = 0; seed < scan->num_points && n_detected < max_detections; seed++) {
        if (visited[seed]) continue;

        /* BFS cluster */
        size_t *queue = (size_t*)malloc(scan->num_points * sizeof(size_t));
        if (!queue) break;
        size_t q_head = 0, q_tail = 0;
        queue[q_tail++] = seed;
        visited[seed] = 1;

        double x_sum = 0, y_sum = 0, z_sum = 0;
        double xmin = 1e100, xmax = -1e100, ymin = 1e100, ymax = -1e100;
        double zmin = 1e100, zmax = -1e100;
        size_t csize = 0;

        while (q_head < q_tail) {
            size_t ci = queue[q_head++];
            double x = scan->points[ci].x;
            double y = scan->points[ci].y;
            double z = scan->points[ci].z;

            x_sum += x; y_sum += y; z_sum += z;
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
            csize++;

            for (size_t nj = 0; nj < scan->num_points; nj++) {
                if (visited[nj]) continue;
                double dx = scan->points[nj].x - x;
                double dy = scan->points[nj].y - y;
                double dz = scan->points[nj].z - z;
                if (dx * dx + dy * dy + dz * dz <= cd2) {
                    visited[nj] = 1;
                    queue[q_tail++] = nj;
                }
            }
        }
        free(queue);

        if (csize < min_cluster_size || csize > max_cluster_size) continue;

        lidar_object_detection_t *det = &detections[n_detected];
        det->x_center = x_sum / (double)csize;
        det->y_center = y_sum / (double)csize;
        det->width    = xmax - xmin;
        det->length   = ymax - ymin;
        det->height   = zmax - zmin;
        det->heading  = 0.0;  /* Cannot determine heading without shape fitting */
        det->velocity = 0.0;  /* Single-frame: no velocity estimate */
        det->class_id = 0;    /* Unclassified */
        det->confidence = (csize >= min_cluster_size * 2) ? 0.8 : 0.5;
        n_detected++;
    }

    free(visited);
    return n_detected;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L8: FMCW LiDAR
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_fmcw_beat_frequency(double B_chirp, double T_chirp,
                                double f0, double R, double v,
                                double *f_beat) {
    if (!f_beat || B_chirp <= 0.0 || T_chirp <= 0.0 || f0 <= 0.0) return -1;

    /* f_beat = 2*B*R/(c*T) + 2*v*f0/c */
    double f_range = 2.0 * B_chirp * R / (LIDAR_C * T_chirp);
    double f_doppler = 2.0 * v * f0 / LIDAR_C;
    *f_beat = f_range + f_doppler;

    return 0;
}

double lidar_fmcw_range_from_beat(double f_beat,
                                    double B_chirp, double T_chirp) {
    if (B_chirp <= 0.0 || T_chirp <= 0.0) return 0.0;
    /* R = f_beat * c * T_chirp / (2 * B_chirp) */
    return f_beat * LIDAR_C * T_chirp / (2.0 * B_chirp);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * L7: Atmospheric LiDAR — aerosol backscatter (Klett-Fernald method)
 * ═══════════════════════════════════════════════════════════════════════════ */

int lidar_aerosol_backscatter(const double *range_corrected,
                                const double *ranges,
                                size_t num_bins,
                                double lidar_ratio,
                                size_t ref_range_idx,
                                double ref_beta,
                                double *beta_aer) {
    if (!range_corrected || !ranges || !beta_aer || num_bins < 2) return -1;
    if (ref_range_idx >= num_bins) return -1;

    /* Klett backward inversion:
       beta_aer(R) = X(R)*exp(2*(S_aer - S_mol)*(tau_mol(R_ref) - tau_mol(R)))
                     / (X(R_ref)/beta(R_ref) + 2*S_aer*int_R^Rref X(r)*exp(...)) */

    /* For simplicity, implement the near-end solution with constant lidar ratio.
       Assume molecular contribution is negligible (valid at IR wavelengths). */

    double S_aer = lidar_ratio; /* aerosol lidar ratio [sr] */
    double beta_ref = ref_beta;
    double X_ref = range_corrected[ref_range_idx];

    /* Denominator integral: trapezoidal from current bin to reference bin.
       Working backwards from reference to near range. */

    for (size_t i = 0; i < num_bins; i++) {
        if (i == ref_range_idx) {
            beta_aer[i] = beta_ref;
            continue;
        }

        /* Integrate from i to ref_range_idx */
        double integral = 0.0;
        if (i < ref_range_idx) {
            for (size_t j = i; j < ref_range_idx; j++) {
                double dr = ranges[j + 1] - ranges[j];
                if (dr < 0.0) dr = -dr;
                integral += 0.5 * (range_corrected[j] + range_corrected[j + 1]) * dr;
            }
        } else {
            for (size_t j = ref_range_idx; j < i; j++) {
                double dr = ranges[j + 1] - ranges[j];
                if (dr < 0.0) dr = -dr;
                integral += 0.5 * (range_corrected[j] + range_corrected[j + 1]) * dr;
            }
            integral = -integral;
        }

        double denom = X_ref / beta_ref + 2.0 * S_aer * integral;
        if (fabs(denom) < 1e-30) {
            beta_aer[i] = 0.0;
        } else {
            beta_aer[i] = range_corrected[i] / denom;
        }
    }

    return 0;
}

double lidar_boundary_layer_height(const double *range_corrected,
                                     const double *altitudes,
                                     size_t num_bins,
                                     size_t smooth_window) {
    if (!range_corrected || !altitudes || num_bins < 2 * smooth_window + 3) return -1.0;

    /* Compute gradient of range-corrected signal */
    double *gradient = (double*)malloc(num_bins * sizeof(double));
    if (!gradient) return -1.0;

    /* Smoothed derivative */
    for (size_t i = smooth_window; i < num_bins - smooth_window; i++) {
        double dy = range_corrected[i + smooth_window] - range_corrected[i - smooth_window];
        double dx = altitudes[i + smooth_window] - altitudes[i - smooth_window];
        if (dx > 0.0) {
            gradient[i] = dy / dx;
        } else {
            gradient[i] = 0.0;
        }
    }

    /* Find altitude of minimum gradient (strongest negative slope = PBL top) */
    double min_grad = 1e100;
    double pbl_height = -1.0;
    for (size_t i = smooth_window; i < num_bins - smooth_window; i++) {
        if (gradient[i] < min_grad) {
            min_grad = gradient[i];
            pbl_height = altitudes[i];
        }
    }

    free(gradient);
    return pbl_height;
}