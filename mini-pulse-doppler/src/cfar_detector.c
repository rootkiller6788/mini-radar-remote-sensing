#include "cfar_detector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int cfar_config_init(cfar_type_t type, uint32_t num_reference,
                     uint32_t num_guard, double pfa, cfar_config_t *config)
{
    if (config == NULL || num_reference == 0 || pfa <= 0.0 || pfa >= 1.0) return -1;
    config->type = type;
    config->num_reference_cells = num_reference;
    config->num_guard_cells = num_guard;
    config->pfa = pfa;
    config->os_rank = (double)(3 * num_reference / 4);
    config->vi_threshold = 4.76;
    config->mr_threshold = 1.806;

    double N = (double)num_reference;
    switch (type) {
    case CFAR_CA:
    case CFAR_GO:
    case CFAR_SO:
    case CFAR_CAGO:
    case CFAR_CASO:
        config->threshold_factor = N * (pow(pfa, -1.0/N) - 1.0);
        break;
    case CFAR_OS:
    case CFAR_OSCA:
        config->threshold_factor = pow(pfa, -1.0/config->os_rank) - 1.0;
        break;
    case CFAR_TM:
        config->threshold_factor = N * (pow(pfa, -1.0/(N-1.0)) - 1.0);
        break;
    case CFAR_VI:
        config->threshold_factor = N * (pow(pfa, -1.0/N) - 1.0);
        break;
    default: return -1;
    }
    return 0;
}

int ca_cfar_detect(const double *data, size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections, size_t *num_detections)
{
    if (data == NULL || config == NULL || detections == NULL || num_detections == NULL)
        return -1;

    uint32_t G = config->num_guard_cells;
    uint32_t R = config->num_reference_cells / 2;
    size_t count = 0;

    for (size_t i = R + G; i + R + G < length; i++) {
        double sum = 0.0;
        for (size_t j = i - R - G; j < i - G; j++)
            sum += data[j];
        for (size_t j = i + G + 1; j <= i + G + R; j++)
            sum += data[j];

        double noise_est = sum / (double)(2 * R);
        double threshold = config->threshold_factor * noise_est;

        detections[count].cell_index = (uint32_t)i;
        detections[count].cell_value = data[i];
        detections[count].estimated_noise_power = noise_est;
        detections[count].threshold_value = threshold;
        detections[count].detected = (data[i] > threshold) ? 1 : 0;
        count++;
    }
    *num_detections = count;
    return 0;
}

int go_cfar_detect(const double *data, size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections, size_t *num_detections)
{
    if (data == NULL || config == NULL || detections == NULL || num_detections == NULL)
        return -1;

    uint32_t G = config->num_guard_cells;
    uint32_t R = config->num_reference_cells / 2;
    size_t count = 0;

    for (size_t i = R + G; i + R + G < length; i++) {
        double sum_lead = 0.0, sum_lag = 0.0;
        for (size_t j = i - R - G; j < i - G; j++) sum_lag += data[j];
        for (size_t j = i + G + 1; j <= i + G + R; j++) sum_lead += data[j];

        double avg_lead = sum_lead / (double)R;
        double avg_lag = sum_lag / (double)R;
        double noise_est = (avg_lead > avg_lag) ? avg_lead : avg_lag;

        double threshold = config->threshold_factor * noise_est;
        detections[count].cell_index = (uint32_t)i;
        detections[count].cell_value = data[i];
        detections[count].estimated_noise_power = noise_est;
        detections[count].threshold_value = threshold;
        detections[count].detected = (data[i] > threshold) ? 1 : 0;
        count++;
    }
    *num_detections = count;
    return 0;
}

int so_cfar_detect(const double *data, size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections, size_t *num_detections)
{
    if (data == NULL || config == NULL || detections == NULL || num_detections == NULL)
        return -1;

    uint32_t G = config->num_guard_cells;
    uint32_t R = config->num_reference_cells / 2;
    size_t count = 0;

    for (size_t i = R + G; i + R + G < length; i++) {
        double sum_lead = 0.0, sum_lag = 0.0;
        for (size_t j = i - R - G; j < i - G; j++) sum_lag += data[j];
        for (size_t j = i + G + 1; j <= i + G + R; j++) sum_lead += data[j];

        double avg_lead = sum_lead / (double)R;
        double avg_lag = sum_lag / (double)R;
        double noise_est = (avg_lead < avg_lag) ? avg_lead : avg_lag;

        double threshold = config->threshold_factor * noise_est;
        detections[count].cell_index = (uint32_t)i;
        detections[count].cell_value = data[i];
        detections[count].estimated_noise_power = noise_est;
        detections[count].threshold_value = threshold;
        detections[count].detected = (data[i] > threshold) ? 1 : 0;
        count++;
    }
    *num_detections = count;
    return 0;
}

static int compare_double_asc(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

int os_cfar_detect(const double *data, size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections, size_t *num_detections)
{
    if (data == NULL || config == NULL || detections == NULL || num_detections == NULL)
        return -1;

    uint32_t G = config->num_guard_cells;
    uint32_t R = config->num_reference_cells / 2;
    size_t count = 0;
    double *window = (double *)malloc((2 * R) * sizeof(double));
    if (window == NULL) return -1;

    for (size_t i = R + G; i + R + G < length; i++) {
        uint32_t idx = 0;
        for (size_t j = i - R - G; j < i - G; j++) window[idx++] = data[j];
        for (size_t j = i + G + 1; j <= i + G + R; j++) window[idx++] = data[j];
        qsort(window, 2*R, sizeof(double), compare_double_asc);

        int k_rank = (int)config->os_rank;
        if (k_rank < 1) k_rank = 1;
        if (k_rank > (int)(2*R)) k_rank = (int)(2*R);
        double noise_est = window[k_rank - 1];

        double threshold = config->threshold_factor * noise_est;
        detections[count].cell_index = (uint32_t)i;
        detections[count].cell_value = data[i];
        detections[count].estimated_noise_power = noise_est;
        detections[count].threshold_value = threshold;
        detections[count].detected = (data[i] > threshold) ? 1 : 0;
        count++;
    }
    free(window);
    *num_detections = count;
    return 0;
}

int vi_cfar_detect(const double *data, size_t length,
                   const cfar_config_t *config,
                   cfar_detection_t *detections, size_t *num_detections)
{
    if (data == NULL || config == NULL || detections == NULL || num_detections == NULL)
        return -1;

    uint32_t G = config->num_guard_cells;
    uint32_t R = config->num_reference_cells / 2;
    size_t count = 0;

    for (size_t i = R + G; i + R + G < length; i++) {
        double sumA = 0.0, sumA2 = 0.0, sumB = 0.0, sumB2 = 0.0;
        for (size_t j = i - R - G; j < i - G; j++) {
            sumA += data[j]; sumA2 += data[j]*data[j];
        }
        for (size_t j = i + G + 1; j <= i + G + R; j++) {
            sumB += data[j]; sumB2 += data[j]*data[j];
        }

        double muA = sumA / (double)R, muB = sumB / (double)R;
        double varA = sumA2/(double)R - muA*muA;
        double varB = sumB2/(double)R - muB*muB;
        double VIA = 1.0 + varA/(muA*muA + 1e-30);
        double VIB = 1.0 + varB/(muB*muB + 1e-30);
        double MR = (muA + 1e-30)/(muB + 1e-30);

        int non_homog_A = (VIA > config->vi_threshold) ? 1 : 0;
        int non_homog_B = (VIB > config->vi_threshold) ? 1 : 0;
        int same_mean = (MR > 1.0/config->mr_threshold && MR < config->mr_threshold) ? 1 : 0;

        double noise_est;
        if (!non_homog_A && !non_homog_B && same_mean)
            noise_est = (sumA + sumB) / (double)(2*R);
        else if (!non_homog_A && !non_homog_B && !same_mean)
            noise_est = (muA > muB) ? muA : muB;
        else if (non_homog_A && !non_homog_B)
            noise_est = muB;
        else if (!non_homog_A && non_homog_B)
            noise_est = muA;
        else
            noise_est = (muA < muB) ? muA : muB;

        double threshold = config->threshold_factor * noise_est;
        detections[count].cell_index = (uint32_t)i;
        detections[count].cell_value = data[i];
        detections[count].estimated_noise_power = noise_est;
        detections[count].threshold_value = threshold;
        detections[count].detected = (data[i] > threshold) ? 1 : 0;
        count++;
    }
    *num_detections = count;
    return 0;
}

int cfar_2d_scan(const double *rd_map, uint32_t num_range, uint32_t num_doppler,
                 const cfar_config_t *config, cfar_detection_map_t *result)
{
    if (rd_map == NULL || config == NULL || result == NULL) return -1;
    if (result->detection_map == NULL || result->threshold_map == NULL) return -1;
    if (result->noise_power_map == NULL) return -1;

    result->num_range_bins = num_range;
    result->num_doppler_bins = num_doppler;
    result->num_detections = 0;

    uint32_t G = config->num_guard_cells;
    uint32_t R = config->num_reference_cells / 2;

    for (uint32_t r = 0; r < num_range; r++) {
        for (uint32_t d = 0; d < num_doppler; d++) {
            double val = rd_map[r + d * num_range];
            double sum = 0.0;
            uint32_t ref_count = 0;
            for (int dr = -(int)(R+G); dr <= (int)(R+G); dr++) {
                int nr = (int)r + dr;
                if (nr < 0 || nr >= (int)num_range) continue;
                if (abs(dr) <= (int)G) continue;
                for (int dd = -(int)(R+G); dd <= (int)(R+G); dd++) {
                    int nd = (int)d + dd;
                    if (nd < 0 || nd >= (int)num_doppler) continue;
                    if (abs(dd) <= (int)G) continue;
                    if (dr == 0 && dd == 0) continue;
                    sum += rd_map[nr + nd * num_range];
                    ref_count++;
                }
            }
            double noise_power = (ref_count > 0) ? sum / (double)ref_count : val;
            double threshold = config->threshold_factor * noise_power;

            uint32_t idx = r + d * num_range;
            result->noise_power_map[idx] = noise_power;
            result->threshold_map[idx] = threshold;
            if (val > threshold) {
                result->detection_map[idx] = 1;
                result->num_detections++;
            } else {
                result->detection_map[idx] = 0;
            }
        }
    }
    return 0;
}

int cfar_threshold_factor_compute(cfar_type_t type, double pfa,
                                  uint32_t num_reference, double os_rank,
                                  double *alpha)
{
    if (alpha == NULL || pfa <= 0.0 || pfa >= 1.0 || num_reference == 0) return -1;
    double N = (double)num_reference;
    switch (type) {
    case CFAR_CA: case CFAR_GO: case CFAR_SO:
        *alpha = N * (pow(pfa, -1.0/N) - 1.0); break;
    case CFAR_OS:
        *alpha = pow(pfa, -1.0/os_rank) - 1.0; break;
    default:
        *alpha = N * (pow(pfa, -1.0/N) - 1.0); break;
    }
    return 0;
}

int cfar_pd_compute(cfar_type_t type, double snr_db, double pfa,
                    uint32_t num_reference, double *pd)
{
    if (pd == NULL || pfa <= 0.0 || num_reference == 0) return -1;
    (void)type;
    double snr_linear = pow(10.0, snr_db / 10.0);
    double alpha;
    if (cfar_threshold_factor_compute(type, pfa, num_reference, 1.0, &alpha) != 0)
        return -1;
    *pd = 1.0 / (1.0 + (alpha / (1.0 + snr_linear)));
    if (*pd > 1.0) *pd = 1.0;
    if (*pd < 0.0) *pd = 0.0;
    return 0;
}

int cfar_detection_map_alloc(uint32_t num_range, uint32_t num_doppler,
                              cfar_detection_map_t *map)
{
    if (map == NULL || num_range == 0 || num_doppler == 0) return -1;
    size_t total = (size_t)num_range * num_doppler;
    map->detection_map = (uint8_t *)calloc(total, sizeof(uint8_t));
    map->threshold_map = (double *)calloc(total, sizeof(double));
    map->noise_power_map = (double *)calloc(total, sizeof(double));
    if (!map->detection_map || !map->threshold_map || !map->noise_power_map) {
        free(map->detection_map); free(map->threshold_map);
        free(map->noise_power_map); return -1;
    }
    map->num_range_bins = num_range;
    map->num_doppler_bins = num_doppler;
    map->num_detections = 0;
    return 0;
}

void cfar_detection_map_free(cfar_detection_map_t *map)
{
    if (map != NULL) {
        free(map->detection_map);
        free(map->threshold_map);
        free(map->noise_power_map);
        map->detection_map = NULL;
        map->threshold_map = NULL;
        map->noise_power_map = NULL;
    }
}
