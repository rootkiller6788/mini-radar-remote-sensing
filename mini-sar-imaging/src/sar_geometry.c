/**
 * @file    sar_geometry.c
 * @brief   SAR Geometry Implementation -- L1+L2+L3+L4
 * @details Range equations, Doppler parameters, RCM, coordinate transforms,
 *          antenna pattern, effective velocity, resolution/PRF verification.
 * Reference: Cumming & Wong (2005), Curlander & McDonough (1991)
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "sar_core.h"
#include "sar_geometry.h"

#define WGS84_A  6378137.0
#define WGS84_F  (1.0/298.257223563)
#define WGS84_E2 (2.0*WGS84_F - WGS84_F*WGS84_F)

/* L1: Exact hyperbolic range equation R(eta)=sqrt(R0^2+v^2*(eta-eta_0)^2) */
double sar_range_hyperbolic(double R0, double v, double eta, double eta_0)
{
    double d = eta - eta_0;
    return sqrt(R0*R0 + v*v*d*d);
}

/* L1: Parabolic approximation R(eta)~R0+v^2*(eta-eta_0)^2/(2*R0) */
double sar_range_parabolic(double R0, double v, double eta, double eta_0)
{
    double d = eta - eta_0;
    return R0 + (v*v*d*d)/(2.0*R0);
}

/* L1: 4th-order range approx for high-resolution SAR */
double sar_range_quartic(double R0, double v, double eta, double eta_0)
{
    double d = eta - eta_0;
    double d2 = d*d, d4 = d2*d2;
    return R0 + (v*v*d2)/(2.0*R0) - (v*v*v*v*d4)/(8.0*R0*R0*R0);
}

/* L1: Doppler centroid f_Dc = 2*v/lambda * sin(theta_sq) */
double sar_doppler_centroid(double lambda, double v, double theta_sq)
{
    return (2.0*v/lambda) * sin(theta_sq);
}

/* L1: Doppler rate f_R = -2*v^2*cos^3(theta_sq)/(lambda*R0) */
double sar_doppler_rate(double lambda, double v, double R0, double theta_sq)
{
    double cs = cos(theta_sq);
    return -2.0*v*v*cs*cs*cs / (lambda*R0);
}

/* L1: Synthetic aperture time T_sa = lambda*R0/(L_a*v*cos(theta_sq)) */
double sar_synthetic_aperture_time(double lambda, double R0, double L_a,
                                   double v, double theta_sq)
{
    double cs = cos(theta_sq);
    if (cs < 1e-6) cs = 1e-6;
    return lambda*R0/(L_a*v*cs);
}

/* L2: Initialize RCM parameters from geometry */
void sar_rcm_init(sar_rcm_t *rcm, double R0, double v, double lambda,
                  double theta_sq, double range_cell_size)
{
    if (!rcm) return;
    memset(rcm, 0, sizeof(*rcm));
    rcm->R0 = R0; rcm->v = v; rcm->lambda = lambda;
    rcm->f_Dc = sar_doppler_centroid(lambda, v, theta_sq);
    rcm->f_R  = sar_doppler_rate(lambda, v, R0, theta_sq);
    rcm->range_cell_size_m = range_cell_size;
    rcm->rcm_linear_coeff = lambda * rcm->f_Dc / 2.0;
    rcm->rcm_quadratic_coeff = -lambda * rcm->f_R / 4.0;
    double eta_max = 1.0;
    rcm->max_rcm_m = fabs(rcm->rcm_linear_coeff)*eta_max
                   + fabs(rcm->rcm_quadratic_coeff)*eta_max*eta_max;
    rcm->max_rcm_cells = (size_t)ceil(rcm->max_rcm_m / range_cell_size);
}

/* L2: RCM evaluation at azimuth time eta */
double sar_rcm_at_eta(const sar_rcm_t *rcm, double eta)
{
    if (!rcm) return 0.0;
    return rcm->rcm_linear_coeff*eta + rcm->rcm_quadratic_coeff*eta*eta;
}


/* L1: Slant to ground range G = sqrt(R_slant^2 - H^2) */
double sar_slant_to_ground(double R_slant, double H) {
    double r2 = R_slant*R_slant - H*H;
    return (r2 > 0.0) ? sqrt(r2) : 0.0;
}

/* L1: Ground to slant range */
double sar_ground_to_slant(double G, double H) {
    return sqrt(G*G + H*H);
}

/* L1: Incidence angle theta_i = asin(H/R_slant) */
double sar_incidence_angle(double H, double R_slant) {
    if (R_slant <= H || R_slant <= 0.0) return M_PI/2.0;
    double s = H/R_slant;
    if (s > 1.0) s = 1.0;
    return asin(s);
}

/* L1: ECEF to geodetic (WGS84) using Bowering iterative method
 * Convergence within 3-4 iterations for h < 1000 km */
void sar_ecef_to_geodetic(const sar_ecef_t *ecef, sar_geodetic_t *geo) {
    if (!ecef || !geo) return;
    double x=ecef->x_m, y=ecef->y_m, z=ecef->z_m;
    double p = sqrt(x*x + y*y);
    double lon = atan2(y, x);
    double lat = atan2(z, p*(1.0-WGS84_E2));
    double h = 0.0;
    for (int i=0; i<10; i++) {
        double sl = sin(lat), N = WGS84_A/sqrt(1.0-WGS84_E2*sl*sl);
        double hn = p/cos(lat) - N;
        double ln = atan2(z, p*(1.0-WGS84_E2*N/(N+h)));
        if (fabs(ln-lat)<1e-12 && fabs(hn-h)<1e-6) break;
        lat=ln; h=hn;
    }
    geo->lat_deg = lat*180.0/M_PI;
    geo->lon_deg = lon*180.0/M_PI;
    geo->alt_m = h;
}

/* L1: Geodetic to ECEF (WGS84)
 * N = a/sqrt(1-e^2*sin^2(lat)), x=(N+h)cos(lat)cos(lon),
 * y=(N+h)cos(lat)sin(lon), z=(N(1-e^2)+h)sin(lat) */
void sar_geodetic_to_ecef(const sar_geodetic_t *geo, sar_ecef_t *ecef) {
    if (!geo || !ecef) return;
    double lat=geo->lat_deg*M_PI/180.0, lon=geo->lon_deg*M_PI/180.0, h=geo->alt_m;
    double sl=sin(lat), cl=cos(lat);
    double N = WGS84_A/sqrt(1.0-WGS84_E2*sl*sl);
    ecef->x_m = (N+h)*cl*cos(lon);
    ecef->y_m = (N+h)*cl*sin(lon);
    ecef->z_m = (N*(1.0-WGS84_E2)+h)*sl;
}

/* L2: Antenna init - theta_3dB = 0.886*lambda/L_a */
void sar_antenna_init(sar_antenna_t *ant, double L_a, double lambda, double R0) {
    if (!ant) return;
    ant->antenna_length_m = L_a;
    ant->wavelength_m = lambda;
    ant->beamwidth_rad = (L_a>0.0 && lambda>0.0) ? 0.886*lambda/L_a : 0.0;
    ant->beamwidth_az_km = ant->beamwidth_rad*R0/1000.0;
}

/* L2: Two-way antenna gain G(θ) = sinc^4(pi*L_a/lambda*sin(theta))
 * Two-way because same antenna transmits and receives */
double sar_antenna_gain(const sar_antenna_t *ant, double theta) {
    if (!ant || ant->antenna_length_m<=0.0 || ant->wavelength_m<=0.0) return 0.0;
    double x = M_PI*ant->antenna_length_m*sin(theta)/ant->wavelength_m;
    if (fabs(x)<1e-12) return 1.0;
    double sv = sin(x)/x;
    double s2 = sv*sv;
    return s2*s2;
}

/* L3: Effective velocity with Earth curvature
 * v_g = v_plat*Re/(Re+H), v_eff = sqrt(v_plat*v_g) */
double sar_effective_velocity(double v_plat, double H, double Re) {
    if (Re <= 0.0) return v_plat;
    double vg = v_plat*Re/(Re+H);
    return sqrt(v_plat*vg);
}

/* L4: Resolution check -- verifies rho_r=c/(2*B_r) and rho_a=L_a/2 */
int sar_resolution_check(const sar_params_t *sp, double tol) {
    if (!sp || tol<0.0) return 0;
    double rr = SAR_C/(2.0*sp->bandwidth_hz);
    double ra = sp->antenna_length_m/2.0;
    double er = fabs(sp->range_resolution_m-rr)/rr;
    double ea = fabs(sp->azimuth_resolution_m-ra)/ra;
    return (er<=tol && ea<=tol) ? 1 : 0;
}

/* L4: PRF Nyquist check: PRF >= 2*v/L_a (azimuth) AND PRF <= c/(2*swath) (range) */
int sar_prf_check(const sar_params_t *sp) {
    if (!sp) return 0;
    double pmin = 2.0*sp->platform_velocity_ms/sp->antenna_length_m;
    double pmax = (sp->swath_width_m>0.0) ? SAR_C/(2.0*sp->swath_width_m) : 1e12;
    return (sp->prf_hz>=pmin && sp->prf_hz<=pmax) ? 1 : 0;
}


/* ==================================================================
 * Extended L2: Earth Curvature Models
 * ================================================================== */

/* Geocentric angle gamma between two points on Earth surface.
 * gamma = R/(Re+H) where R is slant range.
 * Used for spaceborne SAR geometry corrections. */
static double geocentric_angle(double R, double H, double Re) {
    return R / (Re + H);
}

/* Slant range correction for spherical Earth.
 * R_true = sqrt( (Re+H)^2 + Re^2 - 2*Re*(Re+H)*cos(gamma) )
 * Replaces flat-Earth approximation for spaceborne SAR. */
double sar_slant_range_spherical(double G, double H, double Re) {
    double gamma = G / Re;
    double R1 = Re + H;
    double cos_gamma = cos(gamma);
    return sqrt(R1*R1 + Re*Re - 2.0*Re*R1*cos_gamma);
}

/* Ground range from slant range on spherical Earth.
 * gamma = acos( ((Re+H)^2 + Re^2 - R^2) / (2*Re*(Re+H)) )
 * G = Re * gamma */
double sar_ground_range_spherical(double R_slant, double H, double Re) {
    double R1 = Re + H;
    double cos_gamma = (R1*R1 + Re*Re - R_slant*R_slant) / (2.0*Re*R1);
    if (cos_gamma > 1.0) cos_gamma = 1.0;
    if (cos_gamma < -1.0) cos_gamma = -1.0;
    return Re * acos(cos_gamma);
}

/* Local incidence angle on spherical Earth.
 * theta_local = asin( (Re+H)*sin(gamma) / R )
 * Accounts for Earth curvature (important for spaceborne InSAR). */
double sar_local_incidence_spherical(double R_slant, double H, double Re) {
    double R1 = Re + H;
    double gamma = sar_ground_range_spherical(R_slant, H, Re) / Re;
    double sin_theta = R1 * sin(gamma) / R_slant;
    if (sin_theta > 1.0) sin_theta = 1.0;
    if (sin_theta < -1.0) sin_theta = -1.0;
    return asin(sin_theta);
}

/* ==================================================================
 * L3: Doppler Parameter Validation
 * ================================================================== */

/* Compute azimuth FM rate from geometry without squint approximation.
 * For validation against the simplified formula.
 * Full expression: f_R = -2*v^2/(lambda*R0) * [1 + (v*eta)^2/R0^2]^(-3/2)
 * Evaluated at eta=0: simplifies to -2*v^2/(lambda*R0). */
double sar_doppler_rate_full(double lambda, double v, double R0, double eta) {
    double v2 = v * v;
    double factor = 1.0 + v2 * eta * eta / (R0 * R0);
    return -2.0 * v2 / (lambda * R0) * pow(factor, -1.5);
}

/* Compute the quadratic phase error from using parabolic approx.
 * QPE = pi * v^4 * T_sa^4 / (64 * lambda * R0^3)
 * Must be < pi/4 for RDA validity. Returns QPE in radians. */
double sar_quadratic_phase_error(double lambda, double v, double R0, double T_sa) {
    double num = M_PI * v*v*v*v * T_sa*T_sa*T_sa*T_sa;
    double den = 64.0 * lambda * R0*R0*R0;
    return num / den;
}

/* Check if RDA is valid for given SAR parameters.
 * Returns 1 if QPE < pi/4 (parabolic approx valid). */
int sar_rda_validity_check(double lambda, double v, double R0, double L_a, double theta_sq) {
    double T_sa = sar_synthetic_aperture_time(lambda, R0, L_a, v, theta_sq);
    double QPE = sar_quadratic_phase_error(lambda, v, R0, T_sa);
    return (QPE < M_PI/4.0) ? 1 : 0;
}

/* Compute range-Doppler coupling phase for squinted SAR.
 * Secondary range compression (SRC) compensates this coupling.
 * SRC phase: phi_SRC(f_tau, f_eta) = -pi * f_tau^2 * f_eta^2 * lambda / (K_r^2 * c) */
double sar_src_phase(double f_tau, double f_eta, double lambda, double K_r) {
    return -M_PI * f_tau*f_tau * f_eta*f_eta * lambda / (K_r*K_r * SAR_C);
}

/* ==================================================================
 * L5: Azimuth Compression Filter Design
 * ================================================================== */

/* Generate azimuth matched filter coefficients for a given range bin.
 * H_az[n] = exp(-j*pi*f_R*(n*d_eta)^2), n = -N/2...N/2-1
 * The filter length N = T_sa / d_eta (rounded to integer).
 * Output: h_I, h_Q arrays of length N. Caller must free. */
int sar_azimuth_matched_filter(double f_R, double PRF, double T_sa,
                               double **h_I_out, double **h_Q_out, size_t *N_out)
{
    double d_eta = 1.0 / PRF;
    size_t N = (size_t)ceil(T_sa / d_eta);
    if (N < 4) { *h_I_out = NULL; *h_Q_out = NULL; *N_out = 0; return -1; }
    if (N % 2 == 0) N++;

    double *h_I = (double *)malloc(N * sizeof(double));
    double *h_Q = (double *)malloc(N * sizeof(double));
    if (!h_I || !h_Q) { free(h_I); free(h_Q); return -2; }

    for (size_t n = 0; n < N; n++) {
        double eta = ((double)((int64_t)n - (int64_t)(N/2))) * d_eta;
        double phase = -M_PI * f_R * eta * eta;
        h_I[n] = cos(phase);
        h_Q[n] = sin(phase);
    }
    *h_I_out = h_I; *h_Q_out = h_Q; *N_out = N;
    return 0;
}
