#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef std::vector<double> DoubleVec;

double avs_to_norm(double val) { return (val - 16.0) / 219.0; }
double norm_to_avs(double norm_val) { return norm_val * 219.0 + 16.0; }

double sinc(double x) {
    if (std::abs(x) < 1e-12) return 1.0;
    double px = M_PI * x;
    return std::sin(px) / px;
}

double ud_kernel(double x, const DoubleVec& coeffs, double s) {
    if (std::abs(x) > s) return 0.0;
    double val = sinc(x);
    for (size_t k = 0; k < coeffs.size(); ++k) {
        double shift = double(k + 1);
        val += coeffs.at(k) * (sinc(x + shift) + sinc(x - shift));
    }
    return val;
}

double find_recommended_support(const DoubleVec& coeffs, double R2, double threshold = 0.002) {
    double max_scan = 16.0;
    double step = 0.01;
    double last_significant_x = 2.0;

    for (double x = 0.0; x <= max_scan; x += step) {
        double val = std::abs(ud_kernel(x, coeffs, max_scan));
        if (val >= threshold) {
            last_significant_x = x;
        }
    }

    double grid_step = 0.5 / R2;
    double recommended_s = std::ceil(last_significant_x / grid_step) * grid_step + 0.05;
    return std::max(2.0, recommended_s);
}

struct DiscreteKernel {
    int half_win;
    DoubleVec weights;
};

DiscreteKernel build_discrete_kernel_phased(double R, const DoubleVec& coeffs, double s, double phase = 0.5, int extra_pad = 10) {
    int hw = int(std::ceil(s * R)) + extra_pad;
    DiscreteKernel dk;
    dk.half_win = hw;
    dk.weights.assign(2 * hw + 1, 0.0);

    double sum = 0.0;
    for (int n = -hw; n <= hw; ++n) {
        double x = (double(n) - phase) / R;
        double w = ud_kernel(x, coeffs, s);
        dk.weights.at(n + hw) = w;
        sum += w;
    }

    if (std::abs(sum) > 1e-12) {
        for (size_t i = 0; i < dk.weights.size(); ++i) {
            dk.weights.at(i) /= sum;
        }
    }
    return dk;
}

DoubleVec compute_step_response(const DoubleVec& w) {
    DoubleVec step(w.size(), 0.0);
    double acc = 0.0;
    for (size_t i = 0; i < w.size(); ++i) {
        acc += w.at(i);
        step.at(i) = acc;
    }
    return step;
}

double get_negative_overshoot(const DoubleVec& step) {
    double min_val = 0.0;
    for (size_t i = 0; i < step.size(); ++i) {
        if (step.at(i) < min_val) min_val = step.at(i);
    }
    return min_val;
}

double get_rise_time_px(const DoubleVec& step, double R_total) {
    double pos_10 = -1.0, pos_90 = -1.0;
    for (size_t i = 1; i < step.size(); ++i) {
        double prev = step.at(i - 1);
        double curr = step.at(i);
        if (pos_10 < 0.0 && prev <= 0.1 && curr > 0.1) {
            double frac = (0.1 - prev) / (curr - prev);
            pos_10 = double(i - 1) + frac;
        }
        if (pos_90 < 0.0 && prev <= 0.9 && curr > 0.9) {
            double frac = (0.9 - prev) / (curr - prev);
            pos_90 = double(i - 1) + frac;
            break;
        }
    }
    if (pos_10 < 0.0 || pos_90 < 0.0) return 10.0;
    return (pos_90 - pos_10) / R_total;
}

double get_overshoot_width_px(const DoubleVec& step, double R_total) {
    double min_val = get_negative_overshoot(step);
    if (min_val > -0.005) return 0.0;
    double half_level = min_val * 0.5;
    double count = 0.0;
    for (size_t i = 1; i < step.size(); ++i) {
        if (step.at(i) <= half_level) count += 1.0;
    }
    return count / R_total;
}

double get_secondary_ringing(const DoubleVec& step, double R_total) {
    int center = int(step.size()) / 2;
    int guard = int(1.6 * R_total);
    double ring_err = 0.0;
    for (int i = 0; i < int(step.size()); ++i) {
        if (i < center - guard) {
            ring_err += step.at(i) * step.at(i);
        }
        else if (i > center + guard) {
            double d = step.at(i) - 1.0;
            ring_err += d * d;
        }
    }
    return ring_err / R_total;
}

DoubleVec compute_cascade_for_phase(
    int R1, const DiscreteKernel& stage1,
    double R2, const DoubleVec& coeffs2, double s2,
    double phase2, int total_half_win)
{
    DoubleVec cascade_w(2 * total_half_win + 1, 0.0);
    int k_max = int(std::ceil(s2 * R2)) + 2;

    double sum_w2 = 0.0;
    std::vector<int> tap_k;
    DoubleVec tap_w;

    for (int k = -k_max; k <= k_max; ++k) {
        double x = (double(k) - phase2) / R2;
        double w = ud_kernel(x, coeffs2, s2);
        if (std::abs(w) > 1e-14) {
            tap_k.push_back(k);
            tap_w.push_back(w);
            sum_w2 += w;
        }
    }

    if (std::abs(sum_w2) < 1e-12) return cascade_w;

    int hw1 = stage1.half_win;
    for (size_t t = 0; t < tap_k.size(); ++t) {
        int k = tap_k.at(t);
        double wk = tap_w.at(t) / sum_w2;
        for (int n = -hw1; n <= hw1; ++n) {
            int pos = R1 * k + n + total_half_win;
            if (pos >= 0 && pos < int(cascade_w.size())) {
                cascade_w.at(pos) += wk * stage1.weights.at(n + hw1);
            }
        }
    }
    return cascade_w;
}

DoubleVec synthesize_input_model(double input_ov_pct, double s_ref) {
    DoubleVec res(2, 0.0);
    if (input_ov_pct <= 0.1) {
        res.at(0) = avs_to_norm(112.0);
        res.at(1) = avs_to_norm(16.0);
        return res;
    }
    double target_ov = -std::abs(input_ov_pct) / 100.0;
    double b_avs = 78.0;
    double low = -65.0, high = 16.0;
    for (int iter = 0; iter < 45; ++iter) {
        double mid = 0.5 * (low + high);
        DoubleVec test_c(2, 0.0);
        test_c.at(0) = avs_to_norm(b_avs);
        test_c.at(1) = avs_to_norm(mid);
        DiscreteKernel dk = build_discrete_kernel_phased(16.0, test_c, s_ref, 0.5, 5);
        double ov = get_negative_overshoot(compute_step_response(dk.weights));
        if (ov > target_ov) high = mid;
        else low = mid;
    }
    res.at(0) = avs_to_norm(b_avs);
    res.at(1) = avs_to_norm(0.5 * (low + high));
    return res;
}

struct OptimizationContext {
    int R1;
    double R2;
    double s2_wide;
    DiscreteKernel stage1;
    int total_half_win;
    double target_overshoot;
    DoubleVec test_phases;
};

double evaluate_loss(const DoubleVec& coeffs2, const OptimizationContext& ctx) {
    double R_tot = double(ctx.R1) * ctx.R2;
    double total_loss = 0.0;

    for (size_t p = 0; p < ctx.test_phases.size(); ++p) {
        double ph = ctx.test_phases.at(p);
        DoubleVec cascade_w = compute_cascade_for_phase(ctx.R1, ctx.stage1, ctx.R2, coeffs2, ctx.s2_wide, ph, ctx.total_half_win);
        DoubleVec step = compute_step_response(cascade_w);

        double ov = get_negative_overshoot(step);
        double rise_time = get_rise_time_px(step, R_tot);
        double ov_width = get_overshoot_width_px(step, R_tot);
        double ringing = get_secondary_ringing(step, R_tot);

        double err_ov = (ov - ctx.target_overshoot);

        double phase_loss = 25000.0 * (err_ov * err_ov) +
            1.8 * (rise_time * rise_time) +
            0.9 * (ov_width * ov_width) +
            150.0 * ringing;

        if (ctx.target_overshoot > -0.002 && ov < -0.002) {
            phase_loss += 5000.0 * (ov * ov);
        }
        total_loss += phase_loss;
    }

    double tail_sum = 0.0;
    for (double x = 2.6; x <= 8.0; x += 0.2) { double w = ud_kernel(x, coeffs2, 8.0); tail_sum += w * w; }
    return (total_loss / double(ctx.test_phases.size())) + 80.0 * tail_sum;
}

struct SimplexVertex {
    DoubleVec x;
    double val;
    bool operator<(const SimplexVertex& other) const { return val < other.val; }
};

DoubleVec nelder_mead(const DoubleVec& start, const OptimizationContext& ctx, int max_iter = 2500) {
    size_t n = start.size();
    std::vector<SimplexVertex> simplex(n + 1);

    simplex.at(0).x = start;
    simplex.at(0).val = evaluate_loss(start, ctx);

    for (size_t i = 0; i < n; ++i) {
        DoubleVec pt = start;
        pt.at(i) += 0.06;
        simplex.at(i + 1).x = pt;
        simplex.at(i + 1).val = evaluate_loss(pt, ctx);
    }

    for (int iter = 0; iter < max_iter; ++iter) {
        std::sort(simplex.begin(), simplex.end());
        if (std::abs(simplex.at(n).val - simplex.at(0).val) < 1e-12) break;

        DoubleVec centroid(n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                centroid.at(j) += simplex.at(i).x.at(j);
            }
        }
        for (size_t j = 0; j < n; ++j) centroid.at(j) /= double(n);

        DoubleVec xr(n, 0.0);
        for (size_t j = 0; j < n; ++j) {
            xr.at(j) = centroid.at(j) + (centroid.at(j) - simplex.at(n).x.at(j));
        }
        double fr = evaluate_loss(xr, ctx);

        if (fr < simplex.at(0).val) {
            DoubleVec xe(n, 0.0);
            for (size_t j = 0; j < n; ++j) {
                xe.at(j) = centroid.at(j) + 2.0 * (xr.at(j) - centroid.at(j));
            }
            double fe = evaluate_loss(xe, ctx);
            if (fe < fr) {
                simplex.at(n).x = xe;
                simplex.at(n).val = fe;
            }
            else {
                simplex.at(n).x = xr;
                simplex.at(n).val = fr;
            }
        }
        else if (fr < simplex.at(n - 1).val) {
            simplex.at(n).x = xr;
            simplex.at(n).val = fr;
        }
        else {
            DoubleVec xc(n, 0.0);
            if (fr < simplex.at(n).val) {
                for (size_t j = 0; j < n; ++j) {
                    xc.at(j) = centroid.at(j) + 0.5 * (xr.at(j) - centroid.at(j));
                }
            }
            else {
                for (size_t j = 0; j < n; ++j) {
                    xc.at(j) = centroid.at(j) - 0.5 * (centroid.at(j) - simplex.at(n).x.at(j));
                }
            }
            double fc = evaluate_loss(xc, ctx);
            if (fc < std::min(fr, simplex.at(n).val)) {
                simplex.at(n).x = xc;
                simplex.at(n).val = fc;
            }
            else {
                for (size_t i = 1; i <= n; ++i) {
                    for (size_t j = 0; j < n; ++j) {
                        simplex.at(i).x.at(j) = simplex.at(0).x.at(j) + 0.5 * (simplex.at(i).x.at(j) - simplex.at(0).x.at(j));
                    }
                    simplex.at(i).val = evaluate_loss(simplex.at(i).x, ctx);
                }
            }
        }
    }
    std::sort(simplex.begin(), simplex.end());
    return simplex.at(0).x;
}
void print_ascii_step(const DoubleVec& step, double R_total) {
    std::cout << "\n=== Step response of the output frame ===\n";
    int center = int(step.size()) / 2;
    int span = int(2.5 * R_total);
    int step_px = std::max(1, int(std::round(R_total / 4.0)));
    for (int idx = center - span; idx <= center + span; idx += step_px) {
        if (idx < 0 || idx >= int(step.size())) continue;
        double px_pos = double(idx - center) / R_total;
        double v = step.at(idx);
        int col_zero = 8;
        int col_one = 36;
        int p = std::clamp(int(std::round(col_zero + v * (col_one - col_zero))), 0, 44);
        std::string line(45, ' ');
        line.at(col_zero) = '|';
        line.at(col_one) = '|';
        line.at(p) = '*';
        std::cout << std::setw(5) << std::setprecision(2) << px_pos << " px | " << line
            << " | " << std::setw(6) << (v * 100.0) << "%\n";
    }
}
int main(int argc, char* argv[]) {
    int mode = 1;
    double R2 = 2.0; // The downscale ratio at the current resize process
    double in_ov_pct = 12.0;
    double out_ov_pct = 12.0;
    if (argc >= 4) {
        R2 = std::stod(argv[1]);
        in_ov_pct = std::stod(argv[2]);
        out_ov_pct = std::stod(argv[3]);
    }
    else {
        std::cout << "Select operating mode:\n";
        std::cout << "  1. Preserve (Keep) (Keep overshoot amplitude, for example 12% -> 12%)\n";
        std::cout << "  2. Film -> Video (Flat film source (overshoot 0%) -> sharp video look/makeup (overshoot 12%)\n";
        std::cout << "  3. Video -> Film (Sharp video look (overshoot 12%) -> flat film look/makeup (0% no overshoots)\n";
        std::cout << "  4. Manual (Manual input of the source and target overshootin %)\n";
        std::cout << "Select mode [1-4]: ";
        std::cin >> mode;
        std::cout << "Enter downscale ratio R2 (integer or floating point  like 1.5 or 2.0): ";
        std::cin >> R2;
        if (mode == 1) {
            std::cout << "Estimation of the overshoot in the source in % (for example 12%): ";
            std::cin >> in_ov_pct;
            out_ov_pct = in_ov_pct;
        }
        else if (mode == 2) {
            in_ov_pct = 0.0;
            out_ov_pct = 12.0;
        }
        else if (mode == 3) {
            in_ov_pct = 12.0;
            out_ov_pct = 0.0;
        }
        else {
            std::cout << "Source overshoot  (%): ";
            std::cin >> in_ov_pct;
            std::cout << "Target overshoot (%): ";
            std::cin >> out_ov_pct;
        }
    }
    if (R2 < 1.0) R2 = 1.0;
    int R1 = 12; // Some 'very big' emulated first stage downscale ratio (like in the shooting optical camera or other source like supersampled 3D render)
    double s_ref = 8.0; // reference support size of the Avisynth resampler filter (UserDefined2Resize s-param)
    DoubleVec coeffs_in = synthesize_input_model(in_ov_pct, s_ref);
    OptimizationContext ctx;
    ctx.R1 = R1;
    ctx.R2 = R2;
    ctx.s2_wide = 8.0;
    ctx.stage1 = build_discrete_kernel_phased(double(R1), coeffs_in, s_ref, 0.5, 10);
    ctx.total_half_win = int(std::ceil(s_ref * R1 * R2)) + 30;
    ctx.target_overshoot = -std::abs(out_ov_pct) / 100.0;
    if (std::abs(R2 - std::round(R2)) < 1e-6) {
        int r_int = int(std::round(R2));
        ctx.test_phases.push_back((r_int % 2 == 0) ? 0.5 : 0.0);
    }
    else {
        double ph_arr[] = { 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875 }; // phases array for non-integer resizing ratio overshoot amplitude evaluation
        for (int i = 0; i < 8; ++i) ctx.test_phases.push_back(ph_arr[i]);
    }
    double b_starts[] = { 82.0, 90.0, 102.0 }; // starting b-params values
    double c_starts[] = { -22.0, -16.0, 4.0 }; // starting c-params values
    DoubleVec best_res(2, 0.0);
    best_res.at(0) = avs_to_norm(80.0);
    best_res.at(1) = avs_to_norm(-20.0);
    double best_loss = 1e18;
    for (int s_idx = 0; s_idx < 3; ++s_idx) {
        DoubleVec st(2, 0.0);
        st.at(0) = avs_to_norm(b_starts[s_idx]);
        st.at(1) = avs_to_norm(c_starts[s_idx]);
        DoubleVec cand = nelder_mead(st, ctx);
        double l = evaluate_loss(cand, ctx);
        if (l < best_loss) {
            best_loss = l;
            best_res = cand;
        }
    }
    double s_opt = find_recommended_support(best_res, R2, 0.002);
    double s_strict = find_recommended_support(best_res, R2, 0.0005);
    double ov_min = 1e9, ov_max = -1e9, ov_sum = 0.0;
    double rise_sum = 0.0, width_sum = 0.0;
    DoubleVec rep_step;
    double R_tot = double(R1) * R2;
    for (size_t i = 0; i < ctx.test_phases.size(); ++i) {
        DoubleVec cw = compute_cascade_for_phase(R1, ctx.stage1, R2, best_res, s_opt, ctx.test_phases.at(i), ctx.total_half_win);
        DoubleVec st = compute_step_response(cw);
        double ov_pct = -get_negative_overshoot(st) * 100.0;
        ov_min = std::min(ov_min, ov_pct);
        ov_max = std::max(ov_max, ov_pct);
        ov_sum += ov_pct;
        rise_sum += get_rise_time_px(st, R_tot);
        width_sum += get_overshoot_width_px(st, R_tot);
        if (i == ctx.test_phases.size() / 2) rep_step = st;
    }
    double n_ph = double(ctx.test_phases.size());
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n======================================================\n";
    std::cout << "Ready to use Avisynth+ filter call (recommended s param value):\n";
    std::cout << "UserDefined2Resize(Round(width/" << R2 << "), Round(height/" << R2
        << "), b=" << norm_to_avs(best_res.at(0)) << ", c=" << norm_to_avs(best_res.at(1))
        << ", s=" << s_opt << ")\n";
    std::cout << "======================================================\n";
    std::cout << "Coefficient R2         : " << R2 << " (" << (ctx.test_phases.size() > 1 ? "float (multi phase computing)" : "integer (single phase computing") << ")\n";
    std::cout << "Recommended support s= : " << s_opt << " (threshold 0.2% from peak)\n";
    std::cout << "More linear processing : s=" << s_strict << " (threshold 0.05% from peak)\n";
    std::cout << "Computed overshoot     : " << (ov_sum / n_ph) << "% (target: " << out_ov_pct << "%)";
    if (ctx.test_phases.size() > 1) {
        std::cout << " [overshoot at different plahses: " << ov_min << "% .. " << ov_max << "%]";
    }
    std::cout << "\n";
    std::cout << "Front duration (rise time/size)    : " << (rise_sum / n_ph) << " output samples (10% -> 90%)\n";
    std::cout << "Overshoot width (50%)   : " << (width_sum / n_ph) << " output samples\n";
    print_ascii_step(rep_step, R_tot);
    return 0;
}