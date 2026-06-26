/**
 * @file    sar_geometry.h
 * @brief   SAR Geometry -- L1 Definitions + L2 Core Concepts + L3 Math + L4 Laws
 *
 * @details Geometric foundations: range equations (hyperbolic/parabolic/quartic),
 *          Doppler centroid & rate, synthetic aperture time, Range Cell Migration,
 *          coordinate transformations (ECEF/geodetic, slant/ground),
 *          antenna beam pattern, effective velocity with Earth curvature,
 *          resolution and PRF constraint verification.
 *
 * Knowledge Mapping:
 *   L1: Slant range, ground range, Doppler centroid f_Dc, Doppler rate f_R,
 *       Range Cell Migration (RCM), synthetic aperture time T_sa
 *   L2: Hyperbolic range equation as SAR phase history, RCM correction,
 *       squint-Doppler relationship, antenna illumination pattern
 *   L3: Taylor expansion of range equation (parabolic to quartic),
 *       Chu relation (space-Doppler domain mapping)
 *   L4: SAR resolution (rho_r = c/(2B_r), rho_a = L_a/2),
 *       PRF Nyquist constraint (PRF >= 2v/L_a)
 *
 * Reference:
 *   Cumming & Wong, "Digital Processing of SAR Data" (2005)
 *   Curlander & McDonough, "Synthetic Aperture Radar" (1991)
 *   Raney, "Theory of Synthetic Aperture Radar" (1971)
 */

#ifndef SAR_GEOMETRY_H
#define SAR_GEOMETRY_H

#include "sar_core.h"

/** @brief Coordinate system types used in SAR processing. */
typedef enum {
    SAR_COORD_SLANT_RANGE = 0,  /**< (r, eta) natural radar coordinates */
    SAR_COORD_GROUND_RANGE = 1, /**< (g, x) projected to Earth surface */
    SAR_COORD_GEOGRAPHIC = 2,   /**< (lat, lon, h) WGS84 geodetic */
    SAR_COORD_MAP = 3           /**< (E, N, U) local map projection */
} sar_coord_type_t;

/** @brief Earth-Centered Earth-Fixed 3D coordinate [m]. */
typedef struct { double x_m, y_m, z_m; } sar_ecef_t;

/** @brief Geodetic coordinate on WGS84 ellipsoid. */
typedef struct { double lat_deg, lon_deg, alt_m; } sar_geodetic_t;

/* === L1: Range Equations === */

/** Exact hyperbolic: R(eta) = sqrt(R0^2 + v^2*(eta-eta_0)^2).
 *  Fundamental SAR range equation. At broadside dR/deta=0, R(eta_0)=R0. */
double sar_range_hyperbolic(double R0, double v, double eta, double eta_0);

/** Parabolic approx: R(eta) ~ R0 + v^2*(eta-eta_0)^2/(2*R0).
 *  Valid when v*|eta-eta_0| << R0 (small synthetic aperture). */
double sar_range_parabolic(double R0, double v, double eta, double eta_0);

/** Quartic (4th-order): R ~ R0 + v^2*d^2/(2R0) - v^4*d^4/(8R0^3).
 *  Corrects hyperbolic curvature for high-resolution / large-aperture SAR. */
double sar_range_quartic(double R0, double v, double eta, double eta_0);

/* === L1: Doppler Parameters === */

/** Doppler centroid: f_Dc = 2*v/lambda * sin(theta_sq).
 *  Zero for broadside, non-zero for squinted (causes range walk). */
double sar_doppler_centroid(double lambda, double v, double theta_sq);

/** Doppler rate: f_R = -2*v^2*cos^3(theta_sq)/(lambda*R0).
 *  Azimuth FM rate; varies inversely with R0 (range-dependent focusing). */
double sar_doppler_rate(double lambda, double v, double R0, double theta_sq);

/** Synthetic aperture time: T_sa = lambda*R0/(L_a*v*cos(theta_sq)).
 *  Duration a point target is illuminated by the azimuth antenna beam. */
double sar_synthetic_aperture_time(double lambda, double R0, double L_a,
                                   double v, double theta_sq);

/* === L1/L2: Range Cell Migration === */

/** RCM = R(eta)-R0 = linear range walk + quadratic range curvature. */
typedef struct {
    double R0, v, lambda, f_Dc, f_R;
    double rcm_linear_coeff, rcm_quadratic_coeff;
    double max_rcm_m, range_cell_size_m;
    size_t max_rcm_cells;
} sar_rcm_t;

/** Init RCM: C_lin=lambda*f_Dc/2, C_quad=-lambda*f_R/4. */
void sar_rcm_init(sar_rcm_t *rcm, double R0, double v, double lambda,
                  double theta_sq, double range_cell_size);

/** Evaluate RCM(eta) = C_lin*eta + C_quad*eta^2. */
double sar_rcm_at_eta(const sar_rcm_t *rcm, double eta);

/* === L1: Coordinate Transforms === */

/** Slant to ground: G = sqrt(R_slant^2 - H^2). */
double sar_slant_to_ground(double R_slant, double H);

/** Ground to slant: R_slant = sqrt(G^2 + H^2). */
double sar_ground_to_slant(double G, double H);

/** Incidence angle: theta_i = asin(H/R_slant). */
double sar_incidence_angle(double H, double R_slant);

/** ECEF to geodetic (WGS84) using Bowering iterative method. */
void sar_ecef_to_geodetic(const sar_ecef_t *ecef, sar_geodetic_t *geo);

/** Geodetic to ECEF (WGS84). */
void sar_geodetic_to_ecef(const sar_geodetic_t *geo, sar_ecef_t *ecef);

/* === L2: Antenna Beam Pattern === */

/** Two-way power pattern: G(theta)=sinc^4(L_a/lambda*sin(theta)). */
typedef struct {
    double antenna_length_m, wavelength_m;
    double beamwidth_rad, beamwidth_az_km;
} sar_antenna_t;

/** Init: theta_3dB = 0.886*lambda/L_a (uniform aperture). */
void sar_antenna_init(sar_antenna_t *ant, double L_a, double lambda, double R0);

/** Two-way gain at off-boresight angle. */
double sar_antenna_gain(const sar_antenna_t *ant, double theta);

/* === L3: Effective Velocity === */

/** v_eff = sqrt(v_plat * v_g), v_g = v_plat*Re/(Re+H).
 *  Accounts for Earth curvature in spaceborne SAR. */
double sar_effective_velocity(double v_plat, double H, double Re);

/* === L4: Resolution & PRF Verification === */

/** Verify rho_r=c/(2*B_r) and rho_a=L_a/2 within tolerance. */
int sar_resolution_check(const sar_params_t *sp, double tol);

/** Check PRF >= 2*v/L_a (azimuth Nyquist) and PRF <= c/(2*swath). */
int sar_prf_check(const sar_params_t *sp);

#endif /* SAR_GEOMETRY_H */
