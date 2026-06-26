/**
 * @file    lidar_applications.h
 * @brief   LiDAR applications — DEM generation, forestry, building extraction
 *
 * Knowledge covered:
 *   L6: Digital Elevation Model (DEM) generation, canopy height model,
 *       building extraction from airborne LiDAR
 *   L7: Autonomous vehicle perception, terrain mapping, forestry metrics,
 *       urban 3D modeling, atmospheric sensing
 *   L8: Multi-spectral LiDAR, single-photon counting, 4D LiDAR concepts
 *
 * Reference:
 *   - Vosselman & Maas, *Airborne and Terrestrial Laser Scanning*,
 *     Whittles Publishing, 2010.
 *   - Shan & Toth, *Topographic Laser Ranging and Scanning: Principles
 *     and Processing*, 2nd ed., CRC Press, 2018.
 *   - Dong & Chen, "LiDAR Remote Sensing and Applications", CRC Press, 2018.
 */

#ifndef LIDAR_APPLICATIONS_H
#define LIDAR_APPLICATIONS_H

#include "lidar_core.h"
#include "lidar_geometry.h"
#include <stddef.h>

/* ─── L6: Digital Elevation Model (DEM) generation ──────────────────────── */

/**
 * @brief DEM / DTM configuration
 */
typedef struct {
    double cell_size;        /**< Grid cell size [m] */
    double x_min, x_max;    /**< Grid extent X [m] */
    double y_min, y_max;    /**< Grid extent Y [m] */
    size_t cols, rows;      /**< Grid dimensions */
} lidar_dem_config_t;

/**
 * @brief Digital Elevation Model (raster grid)
 *
 * Stores elevation values on a regular 2D grid.
 * Missing data (no points in cell) is marked with NAN.
 */
typedef struct {
    double *data;           /**< Elevation values, row-major: data[row*cols+col] */
    size_t  cols, rows;     /**< Grid dimensions */
    double  cell_size;      /**< Cell size [m] */
    double  x_min, y_min;   /**< Grid origin (lower-left corner) */
    double  nodata_value;   /**< Value for cells with no data */
} lidar_dem_t;

/**
 * @brief Compute DEM configuration from point cloud bounds
 *
 * @param scan     Input point cloud
 * @param cell_size Desired cell size [m]
 * @param config   Output DEM config
 * @return         0 on success
 */
int lidar_dem_config_from_scan(const lidar_scan_t *scan,
                                double cell_size,
                                lidar_dem_config_t *config);

/**
 * @brief Generate Digital Elevation Model (DEM) from point cloud
 *
 * For each grid cell, takes the minimum Z of ground-classified points
 * (class 2).  If no ground points exist in a cell, uses the lowest
 * point overall.
 *
 * This produces a Digital Terrain Model (DTM) when using ground points,
 * or a Digital Surface Model (DSM) when using all points (first returns).
 *
 * @param scan   Input point cloud
 * @param config DEM configuration
 * @param dem    Output DEM (must not be initialized — allocated internally)
 * @return       0 on success
 */
int lidar_dem_generate(const lidar_scan_t *scan,
                        const lidar_dem_config_t *config,
                        lidar_dem_t *dem);

/**
 * @brief Generate Digital Surface Model (DSM) from first returns
 *
 * Uses the maximum elevation in each cell from first-return points.
 * The DSM represents the top of the surface including buildings,
 * vegetation canopy, etc.
 *
 * @param scan   Input point cloud
 * @param config DEM configuration
 * @param dsm    Output DSM
 * @return       0 on success
 */
int lidar_dsm_generate(const lidar_scan_t *scan,
                        const lidar_dem_config_t *config,
                        lidar_dem_t *dsm);

/**
 * @brief Compute Canopy Height Model (CHM = DSM - DTM)
 *
 * CHM represents vegetation height above ground.
 *   CHM[i][j] = DSM[i][j] - DTM[i][j]
 *
 * @param dsm  Digital Surface Model
 * @param dtm  Digital Terrain Model
 * @param chm  Output Canopy Height Model (must be same dimensions)
 * @return     0 on success
 */
int lidar_chm_compute(const lidar_dem_t *dsm,
                       const lidar_dem_t *dtm,
                       lidar_dem_t *chm);

/**
 * @brief Free DEM memory
 */
void lidar_dem_free(lidar_dem_t *dem);

/**
 * @brief Compute hillshade (shaded relief) from DEM for visualization
 *
 * Illuminates the DEM from a given sun azimuth and elevation angle.
 * hillshade = 255 * (cos(zenith) * cos(slope) + sin(zenith) * sin(slope) * cos(azimuth - aspect))
 *
 * @param dem        Input DEM
 * @param sun_azimuth  Sun azimuth angle [rad] (0 = north, clockwise)
 * @param sun_zenith   Sun zenith angle [rad] (0 = overhead)
 * @param hillshade    Output hillshade array (same dimensions as DEM, 0-255)
 * @return             0 on success
 */
int lidar_dem_hillshade(const lidar_dem_t *dem,
                         double sun_azimuth, double sun_zenith,
                         uint8_t *hillshade);

/* ─── L6: Building extraction ───────────────────────────────────────────── */

/**
 * @brief Building footprint result
 */
typedef struct {
    double x_min, y_min;    /**< Bounding box [m] */
    double x_max, y_max;
    double height;          /**< Building height above ground [m] */
    double area;            /**< Footprint area [m²] */
    double mean_elevation;  /**< Mean roof elevation [m] */
    size_t num_points;      /**< Number of points in building */
} lidar_building_footprint_t;

/**
 * @brief Extract building footprints from classified point cloud
 *
 * Method:
 *   1. Filter building-class points (class 6) or non-ground points above
 *      a minimum height threshold
 *   2. Cluster points using Euclidean distance
 *   3. Compute bounding box, area, and height for each cluster
 *   4. Filter clusters by minimum area to remove noise
 *
 * @param scan             Input point cloud (should have ground classification)
 * @param min_height       Minimum height above ground for building [m]
 * @param min_area         Minimum building area [m²]
 * @param cluster_dist     Maximum distance for point clustering [m]
 * @param footprints       Output array of building footprints
 * @param max_footprints   Maximum number of footprints
 * @return                 Number of buildings found
 */
size_t lidar_extract_buildings(const lidar_scan_t *scan,
                                 double min_height,
                                 double min_area,
                                 double cluster_dist,
                                 lidar_building_footprint_t *footprints,
                                 size_t max_footprints);

/* ─── L7: Forestry metrics ──────────────────────────────────────────────── */

/**
 * @brief Forestry plot metrics from LiDAR point cloud
 */
typedef struct {
    double canopy_height_max;    /**< Maximum canopy height [m] */
    double canopy_height_mean;   /**< Mean canopy height [m] */
    double canopy_height_std;    /**< Canopy height standard deviation [m] */
    double canopy_cover;         /**< Canopy cover fraction [0-1] (returns above threshold / total) */
    double leaf_area_index;      /**< Estimated LAI (dimensionless) */
    double biomass_estimate;     /**< Above-ground biomass estimate [Mg/ha] */
    double percentile_h50;       /**< 50th percentile height [m] */
    double percentile_h75;       /**< 75th percentile height [m] */
    double percentile_h90;       /**< 90th percentile height [m] */
    double percentile_h95;       /**< 95th percentile height [m] */
    double percentile_h99;       /**< 99th percentile height [m] */
    size_t total_returns;        /**< Total LiDAR returns in plot */
    size_t ground_returns;       /**< Number of ground returns */
} lidar_forestry_metrics_t;

/**
 * @brief Compute forestry metrics from a plot-level point cloud
 *
 * Uses vegetation points (class 3, 4, 5) and ground points (class 2).
 * Height is relative to the DTM (z - z_ground).
 *
 * Canopy cover is the fraction of first returns above a height threshold
 * (typically 2 m) relative to total first returns.
 *
 * LAI is estimated using the Beer-Lambert law:
 *   LAI = -cos(θ) · ln( P_gap(θ) )
 * where P_gap is the gap fraction at scan angle θ.
 *
 * Biomass is estimated using an allometric relationship:
 *   AGB = a · (H_mean)ᵇ · (CC)ᶜ
 * with typical parameters for temperate forests.
 *
 * Reference: Lefsky et al., "LiDAR Remote Sensing for Ecosystem Studies",
 *            *BioScience* 52(1), pp.19-30, 2002.
 *
 * @param scan            Point cloud for a forest plot
 * @param ground_elevation Mean ground elevation [m] (from DTM)
 * @param metrics         Output forestry metrics
 * @return                0 on success
 */
int lidar_forestry_metrics(const lidar_scan_t *scan,
                             double ground_elevation,
                             lidar_forestry_metrics_t *metrics);

/* ─── L7: Autonomous vehicle perception ─────────────────────────────────── */

/**
 * @brief Object detection result for autonomous driving
 */
typedef struct {
    double x_center, y_center;  /**< Object center in sensor frame [m] */
    double width, length, height; /**< Bounding box dimensions [m] */
    double heading;             /**< Object heading angle [rad] */
    double velocity;            /**< Estimated radial velocity [m/s] */
    int    class_id;            /**< Object class (0=car, 1=pedestrian, 2=cyclist, ...) */
    double confidence;          /**< Detection confidence [0-1] */
} lidar_object_detection_t;

/**
 * @brief Simple Euclidean clustering for object detection
 *
 * Clusters points within a distance threshold, then fits oriented
 * bounding boxes to each cluster.  This is the first stage of
 * LiDAR-based perception pipelines for autonomous vehicles
 * (e.g., Waymo, Tesla, Cruise).
 *
 * @param scan             Input point cloud
 * @param cluster_dist     Maximum point-to-point distance for cluster [m]
 * @param min_cluster_size Minimum number of points per cluster
 * @param max_cluster_size Maximum points per cluster (reject large clusters)
 * @param detections       Output detection array
 * @param max_detections   Maximum detections
 * @return                 Number of detected objects
 */
size_t lidar_detect_objects_euclidean(const lidar_scan_t *scan,
                                        double cluster_dist,
                                        size_t min_cluster_size,
                                        size_t max_cluster_size,
                                        lidar_object_detection_t *detections,
                                        size_t max_detections);

/* ─── L8: Advanced LiDAR modalities ─────────────────────────────────────── */

/**
 * @brief Simulate FMCW LiDAR beat frequency
 *
 * In FMCW LiDAR, a chirped laser is transmitted.  The reflected signal
 * is mixed with the local oscillator, producing a beat frequency:
 *
 *   f_beat = (2 · B_chirp · R) / (c · T_chirp) + (2 · v · f_c) / c
 *
 * where:
 *   B_chirp = chirp bandwidth [Hz]
 *   T_chirp = chirp duration [s]
 *   R       = target range [m]
 *   v       = radial velocity [m/s] (positive = approaching)
 *   f_c     = carrier frequency [Hz]
 *
 * Range measurement:  R = f_beat · c · T_chirp / (2 · B_chirp)
 * Velocity measurement: v = (f_beat_up - f_beat_down) · c / (4 · f_c)
 *
 * Reference: Amann et al., "Laser Ranging: A Critical Review of Usual
 *            Techniques for Distance Measurement", *Optical Engineering*
 *            40(1), pp.10-19, 2001.
 *
 * @param B_chirp   Chirp bandwidth [Hz]
 * @param T_chirp   Chirp duration [s]
 * @param f0        Carrier frequency [Hz]
 * @param R         Target range [m]
 * @param v         Radial velocity (+ approaching) [m/s]
 * @param f_beat    Output: beat frequency [Hz]
 * @return          0 on success
 */
int lidar_fmcw_beat_frequency(double B_chirp, double T_chirp,
                                double f0, double R, double v,
                                double *f_beat);

/**
 * @brief Range extraction from FMCW beat frequency
 *
 * R = f_beat · c · T_chirp / (2 · B_chirp)
 *
 * @param f_beat    Beat frequency [Hz]
 * @param B_chirp   Chirp bandwidth [Hz]
 * @param T_chirp   Chirp duration [s]
 * @return          Range [m]
 */
double lidar_fmcw_range_from_beat(double f_beat,
                                    double B_chirp, double T_chirp);

/* ─── L7: Atmospheric LiDAR applications ────────────────────────────────── */

/**
 * @brief Compute aerosol backscatter coefficient from LiDAR signal
 *
 * Uses the LiDAR equation inversion (Fernald/Klett method):
 *
 *   β_aer(R) = X(R) · exp(2·(S_aer - S_mol)·∫ β_mol(r)·dr)
 *              / (C - 2·S_aer·∫ X(r)·exp(2·(S_aer-S_mol)·∫ β_mol·dr')·dr)
 *
 * where X(R) = P(R)·R² is the range-corrected signal,
 * S = α/β is the lidar ratio (extinction-to-backscatter ratio).
 *
 * Reference: Klett, J.D., "Stable Analytical Inversion Solution for
 *            Processing Lidar Returns", *Applied Optics* 20(2), pp.211-220, 1981.
 *            Fernald, F.G., "Analysis of Atmospheric Lidar Observations",
 *            *Applied Optics* 23(5), pp.652-653, 1984.
 *
 * @param range_corrected  Range-corrected signal P(R)·R²
 * @param ranges           Range values [m]
 * @param num_bins         Number of range bins
 * @param lidar_ratio      Aerosol lidar ratio S_aer [sr]
 * @param ref_range_idx    Reference range bin index for boundary condition
 * @param ref_beta         Backscatter coefficient at reference range [m⁻¹·sr⁻¹]
 * @param beta_aer         Output: aerosol backscatter coefficient [m⁻¹·sr⁻¹]
 * @return                 0 on success
 */
int lidar_aerosol_backscatter(const double *range_corrected,
                                const double *ranges,
                                size_t num_bins,
                                double lidar_ratio,
                                size_t ref_range_idx,
                                double ref_beta,
                                double *beta_aer);

/**
 * @brief Simple atmospheric boundary layer height detection
 *
 * Detects the top of the atmospheric boundary layer by finding
 * the altitude of maximum negative gradient in the range-corrected
 * LiDAR signal.  This is the standard gradient method for ceilometer
 * and aerosol LiDAR data.
 *
 * @param range_corrected  Range-corrected signal
 * @param altitudes        Altitude values [m] (AGL)
 * @param num_bins         Number of range bins
 * @param smooth_window    Smoothing window size
 * @return                 Boundary layer height [m], -1 if not detected
 */
double lidar_boundary_layer_height(const double *range_corrected,
                                     const double *altitudes,
                                     size_t num_bins,
                                     size_t smooth_window);

#endif /* LIDAR_APPLICATIONS_H */