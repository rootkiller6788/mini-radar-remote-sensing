/**
 * test_track_core.c — Tests for track management, gating, scoring, matrix ops
 *
 * Tests:
 *   1. Tracker initialization
 *   2. Track creation and lifecycle
 *   3. M/N confirmation logic
 *   4. Mahalanobis distance
 *   5. Chi-squared threshold computation
 *   6. Track scoring
 *   7. Matrix operations (multiply, inverse, determinant, solve)
 *   8. Vector operations
 *   9. Angle utilities
 */

#include "track_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define EPS 1e-6

static int tests_passed = 0;
static int tests_total = 0;

#define TEST(name) do { tests_total++; printf("  %-50s", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

/* ============================================================================
 * Test 1: Tracker initialization
 * ============================================================================
 */

static void test_tracker_init(void)
{
    TEST("Tracker initialization");
    tracker_t tracker;
    tracker_init(&tracker);

    assert(tracker.mgmt.M == 3);
    assert(tracker.mgmt.N == 5);
    assert(tracker.mgmt.coast_max == 5);
    assert(tracker.num_tracks == 0);
    assert(tracker.next_track_id == 1);
    assert(tracker.default_gate.shape == GATE_ELLIPSOIDAL);

    /* All slots should be FREE */
    for (int i = 0; i < TRACKER_MAX_TRACKS; i++) {
        assert(tracker.tracks[i].state == TRACK_STATE_FREE);
        assert(tracker.tracks[i].track_id == TRACK_INVALID_ID);
    }
    PASS();
}

/* ============================================================================
 * Test 2: Track creation and lifecycle
 * ============================================================================
 */

static void test_track_lifecycle(void)
{
    TEST("Track creation lifecycle");
    tracker_t tracker;
    tracker_init(&tracker);

    /* Create first track */
    track_t *t = tracker_create_track(&tracker, 4, 2, MEAS_CARTESIAN_2D);
    assert(t != NULL);
    assert(t->state == TRACK_STATE_TENTATIVE);
    assert(t->state_dim == 4);
    assert(t->meas_dim == 2);
    assert(t->track_id == 1);
    assert(tracker.num_tracks == 1);

    /* Create second track */
    track_t *t2 = tracker_create_track(&tracker, 6, 3, MEAS_CARTESIAN_3D);
    assert(t2 != NULL);
    assert(t2->track_id == 2);
    assert(tracker.num_tracks == 2);

    /* Delete first track */
    track_delete(t);
    assert(t->state == TRACK_STATE_DELETED);

    /* Cleanup: deleted track removed, only t2 remains */
    int removed = tracker_cleanup(&tracker);
    assert(removed == 1);
    assert(tracker.num_tracks == 1);
    assert(tracker.tracks[0].track_id == 2);

    PASS();
}

/* ============================================================================
 * Test 3: M/N confirmation logic
 * ============================================================================
 */

static void test_mn_logic(void)
{
    TEST("M/N confirmation logic");
    tracker_t tracker;
    tracker_init(&tracker);
    tracker.mgmt.M = 3;
    tracker.mgmt.N = 5;

    track_t *t = tracker_create_track(&tracker, 4, 2, MEAS_CARTESIAN_2D);
    assert(t != NULL);
    assert(t->state == TRACK_STATE_TENTATIVE);

    /* Simulate 3 detections in 5 scans */
    for (int scan = 0; scan < 5; scan++) {
        /* Every scan but the 2nd has a detection */
        if (scan != 2) {
            t->misses_since_last = 0; /* associated */
        } else {
            t->misses_since_last = 1; /* missed */
        }
        tracker_apply_mn_logic(&tracker);
    }

    /* After 5 scans with 4 detections ≥ M=3, track should be confirmed */
    assert(t->state == TRACK_STATE_CONFIRMED);
    tracker.total_tracks_confirmed++;

    PASS();
}

/* ============================================================================
 * Test 4: Mahalanobis distance
 * ============================================================================
 */

static void test_mahalanobis(void)
{
    TEST("Mahalanobis distance (2D example)");

    /* Simple 1D state + 1D measurement: H = [1], x = [5], P = [4], R = [1], z = [7]
     * Innovation = 7 - 5 = 2
     * S = 1*4*1 + 1 = 5
     * d² = 2 * (1/5) * 2 = 4/5 = 0.8
     */
    double z[1] = {7.0};
    double R[1] = {1.0};
    double x_pred[1] = {5.0};
    double P_pred[1] = {4.0};
    double H[1] = {1.0};

    double d2 = mahalanobis_distance_sq(z, R, x_pred, P_pred, H, 1, 1);
    assert(fabs(d2 - 0.8) < EPS);
    PASS();
}

/* ============================================================================
 * Test 5: Chi-squared threshold
 * ============================================================================
 */

static void test_chi2_threshold(void)
{
    TEST("Chi-squared threshold (Wilson-Hilferty)");

    /* At P_G = 0.95, n_z = 2: chi2inv(0.95, 2) ≈ 5.991 */
    double gamma = chi2_threshold_approx(0.95, 2);
    assert(fabs(gamma - 5.991) < 0.5); /* approximation tolerance */

    /* At P_G = 0.99, n_z = 2: chi2inv(0.99, 2) ≈ 9.21 */
    gamma = chi2_threshold_approx(0.99, 2);
    assert(fabs(gamma - 9.21) < 1.0);

    PASS();
}

/* ============================================================================
 * Test 6: Track scoring
 * ============================================================================
 */

static void test_track_scoring(void)
{
    TEST("Track scoring (LLR)");

    /* Build a simple innovation covariance S = [[2, 0], [0, 2]] */
    double S[4] = {2.0, 0.0, 0.0, 2.0};
    double d2 = 3.0;
    double P_D = 0.9;
    double beta_F = 1e-4;

    double inc = track_score_increment(P_D, beta_F, S, d2, 2);
    assert(isfinite(inc));

    double miss = track_score_miss_penalty(P_D);
    assert(miss < 0.0); /* penalty is negative */
    assert(fabs(miss - log(0.1)) < EPS); /* log(1-0.9) = log(0.1) */

    PASS();
}

/* ============================================================================
 * Test 7: Matrix operations
 * ============================================================================
 */

static void test_matrix_ops(void)
{
    TEST("Matrix operations");

    /* Matrix multiply: A = [1 2; 3 4], B = [5 6; 7 8]
     * C = A*B = [1*5+2*7, 1*6+2*8; 3*5+4*7, 3*6+4*8] = [19, 22; 43, 50]
     */
    double A[4] = {1, 2, 3, 4};
    double B[4] = {5, 6, 7, 8};
    double C[4];

    mat_mat_mul(A, B, C, 2, 2, 2);
    assert(fabs(C[0] - 19.0) < EPS);
    assert(fabs(C[1] - 22.0) < EPS);
    assert(fabs(C[2] - 43.0) < EPS);
    assert(fabs(C[3] - 50.0) < EPS);

    /* Matrix inverse: A = [2, 1; 1, 2] => A^{-1} = [2/3, -1/3; -1/3, 2/3] */
    double A2[4] = {2, 1, 1, 2};
    double A2_inv[4];
    double det = mat_inv_cholesky(A2, A2_inv, 2);
    assert(fabs(det - 3.0) < EPS);
    assert(fabs(A2_inv[0] - 0.666666) < EPS);
    assert(fabs(A2_inv[1] - (-0.333333)) < EPS);
    assert(fabs(A2_inv[2] - (-0.333333)) < EPS);
    assert(fabs(A2_inv[3] - 0.666666) < EPS);

    PASS();
}

/* ============================================================================
 * Test 8: Vector operations
 * ============================================================================
 */

static void test_vector_ops(void)
{
    TEST("Vector operations");

    double a[3] = {1.0, 2.0, 3.0};
    double b[3] = {4.0, 5.0, 6.0};
    double c[3];

    /* Dot product: 1*4 + 2*5 + 3*6 = 4+10+18 = 32 */
    double dot = vec_dot(a, b, 3);
    assert(fabs(dot - 32.0) < EPS);

    /* Norm: sqrt(1+4+9) = sqrt(14) */
    double norm = vec_norm2(a, 3);
    assert(fabs(norm - sqrt(14.0)) < EPS);

    /* Subtraction: c = a - b = [-3, -3, -3] */
    vec_sub(a, b, c, 3);
    assert(fabs(c[0] + 3.0) < EPS);
    assert(fabs(c[1] + 3.0) < EPS);
    assert(fabs(c[2] + 3.0) < EPS);

    PASS();
}

/* ============================================================================
 * Test 9: Angle utilities
 * ============================================================================
 */

static void test_angle_utils(void)
{
    TEST("Angle utilities");

    /* deg_to_rad */
    assert(fabs(deg_to_rad(180.0) - M_PI) < EPS);
    assert(fabs(deg_to_rad(90.0) - M_PI/2) < EPS);

    /* rad_to_deg */
    assert(fabs(rad_to_deg(M_PI) - 180.0) < EPS);

    /* wrap_angle */
    double a = wrap_angle(3.0 * M_PI); /* 3π should wrap to π */
    assert(fabs(a - M_PI) < EPS || fabs(a + M_PI) < EPS);

    a = wrap_angle_2pi(-0.5);
    assert(a >= 0.0 && a < 2.0 * M_PI);

    PASS();
}

/* ============================================================================
 * Test 10: Linear system solve via Cholesky
 * ============================================================================
 */

static void test_cholesky_solve(void)
{
    TEST("Cholesky linear solve");

    /* Solve [[4, 1], [1, 3]] * x = [5, 4]
     * Solution: x = [1, 1] (since 4*1 + 1*1 = 5, 1*1 + 3*1 = 4)
     */
    double A[4] = {4, 1, 1, 3};
    double b[2] = {5, 4};
    double x[2];

    int ret = mat_solve_cholesky(A, b, x, 2);
    assert(ret == 0);
    assert(fabs(x[0] - 1.0) < EPS);
    assert(fabs(x[1] - 1.0) < EPS);

    PASS();
}

/* ============================================================================
 * Test 11: Determinant calculation
 * ============================================================================
 */

static void test_determinant(void)
{
    TEST("Determinant via Cholesky");

    /* det([[3, 1], [1, 2]]) = 3*2 - 1*1 = 5 */
    double A[4] = {3, 1, 1, 2};
    double det = mat_det_cholesky(A, 2);
    assert(fabs(det - 5.0) < EPS);

    PASS();
}

/* ============================================================================
 * Test 12: Normal quantile (probit)
 * ============================================================================
 */

static void test_normal_quantile(void)
{
    TEST("Normal quantile approximation");

    /* Φ⁻¹(0.975) ≈ 1.96 */
    double z = normal_quantile(0.975);
    assert(fabs(z - 1.96) < 0.1);

    /* Φ⁻¹(0.5) = 0 */
    z = normal_quantile(0.5);
    assert(fabs(z) < 0.1);

    /* Φ⁻¹(0.025) ≈ -1.96 */
    z = normal_quantile(0.025);
    assert(fabs(z + 1.96) < 0.1);

    PASS();
}

/* ============================================================================
 * Test 13: Gate check
 * ============================================================================
 */

static void test_gate_check(void)
{
    TEST("Ellipsoidal gate check");

    assert(ellipsoidal_gate_check(4.0, 5.0) == 1);  /* inside */
    assert(ellipsoidal_gate_check(6.0, 5.0) == 0);  /* outside */
    assert(ellipsoidal_gate_check(-1.0, 5.0) == 0); /* invalid */

    PASS();
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void)
{
    printf("=== Test Suite: mini-target-tracking (Track Core) ===\n\n");

    test_tracker_init();
    test_track_lifecycle();
    test_mn_logic();
    test_mahalanobis();
    test_chi2_threshold();
    test_track_scoring();
    test_matrix_ops();
    test_vector_ops();
    test_angle_utils();
    test_cholesky_solve();
    test_determinant();
    test_normal_quantile();
    test_gate_check();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_total);
    return (tests_passed == tests_total) ? 0 : 1;
}
