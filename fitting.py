import math
from eventsreader import SpdmeEventsReader
from distributionbuilder import DistributionBuilder
from surrogatedistributionbuilder import SurrogateDistributionBuilder
from distributionbuilder_1d import DistributionBuilder_1d
from surrogatedistributionbuilder_1d import SurrogateDistributionBuilder_1d
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt


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
