/**
 * test_kalman_filter.c — Tests for Kalman filter variants
 *
 * Tests:
 *   1. Linear KF predict/update
 *   2. KF steady state (CV model)
 *   3. UKF sigma point generation
 *   4. Information filter update
 *   5. Adaptive R estimation
 *   6. Filter diagnostics (NIS)
 */

#include "kalman_filter.h"
#include "motion_models.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define EPS 1e-4

static int tests_passed = 0;
static int tests_total = 0;

#define TEST(name) do { tests_total++; printf("  %-50s", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

/* ============================================================================
 * Test 1: Linear KF predict/update (1D CV)
 * ============================================================================
 */

static void test_kf_predict_update(void)
{
    TEST("Linear KF predict + update (1D CV)");

    /* 1D CV model: state = [x, v], meas = [x] */
    kf_t *kf = kf_alloc(2, 1, 0);
    assert(kf != NULL);

    /* Model: F = [1, T; 0, 1], H = [1, 0], Q = q*[T^3/3, T^2/2; T^2/2, T], R = [sigma_r^2] */
    double T = 1.0;
    double q = 1.0;
    double F[4] = {1, T, 0, 1};
    double H[2] = {1, 0};
    double Q[4] = {T*T*T/3.0, T*T/2.0, T*T/2.0, T};
    for (int i = 0; i < 4; i++) Q[i] *= q;
    double R[1] = {0.01};

    kf_set_model(kf, F, H, Q, R, NULL, NULL);

    double x0[2] = {0.0, 1.0}; /* start at 0, moving at 1 m/s */
    double P0[4] = {1.0, 0, 0, 1.0};
    kf_set_state(kf, x0, P0);

    /* Predict */
    assert(kf_predict(kf) == 0);
    assert(fabs(kf->x_pred[0] - 1.0) < EPS); /* x = x + v*T = 1 */
    assert(fabs(kf->x_pred[1] - 1.0) < EPS); /* v unchanged */

    /* Update with measurement z = 1.2 */
    double z[1] = {1.2};
    assert(kf_update(kf, z) == 0);

    /* State should move toward measurement */
    assert(kf->x[0] > 0.0 && kf->x[0] < 2.0);
    assert(kf->P[0] < P0[0]); /* uncertainty reduced */

    kf_free(kf);
    PASS();
}

/* ============================================================================
 * Test 2: KF on CV 2D model
 * ============================================================================
 */

static void test_kf_cv2d(void)
{
    TEST("Linear KF with CV 2D model");

    motion_model_t *cv = motion_create_cv(DIM_2D, 1.0, 0.01);
    assert(cv != NULL);

    kf_t *kf = kf_alloc(cv->state_dim, 2, 0);
    assert(kf != NULL);

    /* H measures position only */
    double H[2 * 4]; /* 2 x 4 */
    memset(H, 0, sizeof(H));
    H[0] = 1.0; H[1 * 4 + 2] = 1.0; /* measure x and y */

    double R[4] = {1.0, 0, 0, 1.0};

    kf_set_model(kf, cv->F, H, cv->Q, R, NULL, NULL);

    double x0[4] = {0.0, 10.0, 0.0, 5.0}; /* start at origin, 10 m/s east, 5 m/s north */
    double P0[16];
    memset(P0, 0, sizeof(P0));
    for (int i = 0; i < 4; i++) P0[i * 4 + i] = 10.0;

    kf_set_state(kf, x0, P0);

    /* Run 10 steps */
    for (int k = 0; k < 10; k++) {
        assert(kf_predict(kf) == 0);

        /* Simulate measurement at true position */
        double true_x = 10.0 * (k + 1);
        double true_y = 5.0 * (k + 1);
        double z[2] = {true_x + 0.1 * ((double)rand()/RAND_MAX - 0.5),
                       true_y + 0.1 * ((double)rand()/RAND_MAX - 0.5)};
        assert(kf_update(kf, z) == 0);
    }

    /* After 10s, estimated velocity should be close to true (10, 5) */
    assert(fabs(kf->x[1] - 10.0) < 2.0);
    assert(fabs(kf->x[3] - 5.0) < 2.0);

    kf_free(kf);
    motion_free(cv);
    PASS();
}

/* ============================================================================
 * Test 3: UKF sigma point generation
 * ============================================================================
 */

static void test_ukf_sigma_points(void)
{
    TEST("UKF sigma point generation");

    kf_t *kf = kf_alloc(2, 1, 0);
    assert(kf != NULL);

    double x0[2] = {1.0, 2.0};
    double P0[4] = {4.0, 0.0, 0.0, 1.0};

    kf_set_state(kf, x0, P0);
    ukf_init_params(kf, 1.0, 2.0, 0.0);

    int ret = ukf_generate_sigma_points(kf);
    assert(ret == 0);
    assert(kf->sigma_points != NULL);
    assert(kf->sigma_weights != NULL);

    /* Check number of sigma points (2n+1 = 5) */
    assert(kf->num_sigma == 5);

    /* Check sigma point 0 is the mean */
    assert(fabs(kf->sigma_points[0] - 1.0) < EPS);
    assert(fabs(kf->sigma_points[1] - 2.0) < EPS);

    /* Weights should sum to 1 */
    double sum_weights = 0.0;
    for (int i = 0; i < 5; i++) sum_weights += kf->sigma_weights[i];
    assert(fabs(sum_weights - 1.0) < EPS);

    kf_free(kf);
    PASS();
}

/* ============================================================================
 * Test 4: Information filter update
 * ============================================================================
 */

static void test_info_filter(void)
{
    TEST("Information filter update");

    int n = 2, m = 1;
    double Y[4] = {1.0, 0.0, 0.0, 1.0}; /* initial information */
    double y[2] = {0.5, 0.5};

    double H[2] = {1.0, 0.0}; /* measure first state */
    double R_inv[1] = {10.0}; /* high precision */
    double z[1] = {2.0};

    int ret = info_filter_update(Y, y, H, R_inv, z, n, m);
    assert(ret == 0);

    /* Y should increase: Y += H' * R_inv * H = [[10, 0], [0, 0]] */
    assert(fabs(Y[0] - 11.0) < EPS);
    assert(fabs(Y[1] - 0.0) < EPS);

    /* y should increase: y += H' * R_inv * z = H' * 10 * 2 = [20, 0] */
    assert(fabs(y[0] - 20.5) < EPS);
    assert(fabs(y[1] - 0.5) < EPS);

    PASS();
}

/* ============================================================================
 * Test 5: Adaptive R estimation
 * ============================================================================
 */

static void test_adaptive_R(void)
{
    TEST("Adaptive R estimation");

    kf_t *kf = kf_alloc(2, 1, 0);
    assert(kf != NULL);

    /* Setup simple KF */
    double F[4] = {1, 1, 0, 1};
    double H[2] = {1, 0};
    double Q[4] = {0.33, 0.5, 0.5, 1.0};
    double R[1] = {1.0};
    kf_set_model(kf, F, H, Q, R, NULL, NULL);
    kf->forgetting_factor = 0.95;

    /* Manual setup of internal state for adaptation test */
    kf->R_est = (double *)calloc(1, sizeof(double));
    kf->R_est[0] = 1.0;

    /* Set innovation and P_pred */
    kf->nu[0] = 2.0;
    kf->P_pred[0] = 4.0; kf->P_pred[1] = 1.0;
    kf->P_pred[2] = 1.0; kf->P_pred[3] = 2.0;

    kf_adapt_R(kf);
    assert(kf->R_est[0] > 0.0);

    kf_free(kf);
    PASS();
}

/* ============================================================================
 * Test 6: Filter diagnostics (NIS)
 * ============================================================================
 */

static void test_kf_nis(void)
{
    TEST("KF NIS computation");

    kf_t *kf = kf_alloc(2, 1, 0);
    assert(kf != NULL);

    /* Setup innovation and S */
    kf->nu[0] = 1.0;
    kf->S[0] = 2.0; /* 1x1 S matrix */

    double nis = kf_nis(kf);
    /* NIS = nu' * S^{-1} * nu = 1 * 0.5 * 1 = 0.5 */
    assert(fabs(nis - 0.5) < EPS);

    kf_free(kf);
    PASS();
}

/* ============================================================================
 * Test 7: Filter consistency test
 * ============================================================================
 */

static void test_consistency_test(void)
{
    TEST("KF consistency chi-squared test");

    /* For dof=2, alpha=0.05, threshold = chi2inv(0.95, 2) ≈ 5.991
     * epsilon = 3.0 < 5.991 => consistent */
    int consistent = kf_consistency_test(3.0, 2, 0.05);
    assert(consistent == 1);

    /* epsilon = 10.0 > 5.991 => inconsistent */
    consistent = kf_consistency_test(10.0, 2, 0.05);
    assert(consistent == 0);

    PASS();
}

/* ============================================================================
 * Test 8: NEES computation
 * ============================================================================
 */

static void test_nees(void)
{
    TEST("KF NEES computation");

    kf_t *kf = kf_alloc(2, 1, 0);
    assert(kf != NULL);

    kf->x[0] = 1.0; kf->x[1] = 2.0;
    kf->P[0] = 4.0; kf->P[1] = 0.0;
    kf->P[2] = 0.0; kf->P[3] = 1.0;

    double x_true[2] = {2.0, 2.0};

    /* dx = [-1, 0], P^{-1} = [[0.25, 0], [0, 1]]
     * NEES = 1*0.25*1 + 0 + 0 + 0 = 0.25 */
    double nees = kf_nees(kf, x_true);
    assert(fabs(nees - 0.25) < EPS);

    kf_free(kf);
    PASS();
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void)
{
    printf("=== Test Suite: mini-target-tracking (Kalman Filter) ===\n\n");

    test_kf_predict_update();
    test_kf_cv2d();
    test_ukf_sigma_points();
    test_info_filter();
    test_adaptive_R();
    test_kf_nis();
    test_consistency_test();
    test_nees();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_total);
    return (tests_passed == tests_total) ? 0 : 1;
}
