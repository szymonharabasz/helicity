import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import helicity_model_1d
import helicity_model_3d


def covariance_fit_scipy(predictive_mean, predictive_lower, predictive_upper, best, hist_index, ax, bins, bounds,
                         dir_name):
    def gaus1d(x, a, mean_x, sigma_x):
        x = x - mean_x
        z = a * np.exp(-0.5 * (x / sigma_x) ** 2)
        return z

    def gaus1d_offset(x, a, mean_x, sigma_x, offset):
        x = x - mean_x
        z = offset + a * np.exp(-0.5 * (x / sigma_x) ** 2)
        return z

    def fit_1d(ax):
        curr_bin = bins[hist_index]
        title = (f"{curr_bin.m_min} < $M_{{ee}}$ < {curr_bin.m_max}, {curr_bin.z_min} < $cos(\\theta^{{CM}}_{{"
                 f"\gamma*}})$ < {curr_bin.z_max}")

        mean_x = best[0].item()

        scale1 = bounds[1][0].item()
        scale2 = bounds[1][0].item()
        proj_min_x = bounds[0][0].item()
        proj_max_x = bounds[1][0].item()

        proj = predictive_mean
        proj_lower = predictive_lower
        proj_upper = predictive_upper

        # xmin_ind = max(0,   bin_index(mean_x, proj_min_x, proj_max_x)-50)
        # xmax_ind = min(100, bin_index(mean_x, proj_min_x, proj_max_x)+50)
        xmin_ind = 0
        xmax_ind = 100
        xmin = proj_min_x + xmin_ind / 101. * (proj_max_x - proj_min_x)
        xmax = proj_min_x + (xmax_ind + 1) / 101. * (proj_max_x - proj_min_x)

        x = np.linspace(-2, 2, 101)

        proj1 = proj[xmin_ind:xmax_ind]
        print("range ", xmin_ind, xmax_ind, xmin, xmax)

        # initial_guess = (1.0, mean_x, 0.2*(proj_max_x-proj_min_x))
        initial_guess = (1.0, mean_x, 0.2 * (proj_max_x - proj_min_x), 0)
        eps = 0.001
        # param_bounds = ([0,mean_x-eps,0],[np.inf,mean_x+eps,2])
        # if mean_x >= 1.0:
        #     param_bounds = ([0,mean_x-eps,0],[np.inf,2,2])
        param_bounds = ([0, mean_x - eps, 0, -np.inf], [np.inf, mean_x + eps, 2, np.inf])
        if mean_x >= 1.0:
            param_bounds = ([0, mean_x - eps, 0, -np.inf], [np.inf, 2, 2, np.inf])
        # ax.plot(x, proj, label="Estimated values")
        ax.plot(x, proj)
        ax.fill_between(x, proj_lower, proj_upper, alpha=0.5)
        ax.plot(x, proj, label="Estimated values")
        ax.set_title(title)
        ax.set_xlabel(r"$\lambda_{\theta}$")
        ax.set_ylabel("ndf$/\chi^2$")
        ax.set_ylim(0, 1.5 * proj_upper.max())
        try:
            popt, pcov = opt.curve_fit(gaus1d_offset, x[xmin_ind:xmax_ind], proj1, p0=initial_guess,
                                       bounds=param_bounds, maxfev=2000)
            fit_result = gaus1d_offset(x, *(popt))

            ax.plot(x, fit_result, label="Gaussian fit")
            # ax.set_ylim([0,2*predictive_upper.max()])
            # plt.rc('axes', titlesize=8)
            # plt.rc('axes', labelsize=8)
            # plt.rc('xtick', labelsize=8)
            # plt.rc('ytick', labelsize=8)
            ax.legend()
        except RuntimeError as e:
            print(f"There was an exception {e}")
            popt, pcov = None, None
        return popt, pcov

    # params0, _ = fit_1d(ax[1][0])
    params0, _ = fit_1d(ax)
    plt.savefig(f"{dir_name}/chi2_best_{hist_index}.png", bbox_inches="tight")

    try:
        return params0[1], params0[2]
    except:
        return None


def fit_simple_model(histsData, get_hist_maker_mc, range_used, learn_norm, analyse_3d):
    parameters_all= []
    losses_all = [[]] * len(range_used)
    fit_simple = helicity_model_3d.fit_simple if analyse_3d else helicity_model_1d.fit_simple

    for HIST_INDEX in range_used:
        simple_model = helicity_model_3d.Helicity3d(learn_norm) if analyse_3d else helicity_model_1d.Helicity1d(
            learn_norm)
        hist_data_simple = histsData[0][HIST_INDEX]
        hist_mc_simple = get_hist_maker_mc("pp", HIST_INDEX).make_hists(0.0)[0][HIST_INDEX]
        losses = fit_simple(simple_model, hist_data_simple, hist_mc_simple, 10000, 0.01, learn_norm)
        parameters_all.append([param for param in simple_model.parameters()])
        losses_all[HIST_INDEX] = losses

        return parameters_all, losses_all
