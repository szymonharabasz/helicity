import logging

import torch

import helicity_model_1d
import helicity_model_3d

logging.basicConfig(level=logging.WARNING)

def errors_1d(x, data_y, model_y, sigmas, debug = False):

    logging.debug("model_y ", model_y.shape)

    mask1 = data_y > 0
    mask2 = model_y > 0
    mask = torch.logical_and(mask1, mask2)
    indices = torch.nonzero(mask)
    data_y = data_y[indices].flatten()
    sigmas = sigmas[indices].flatten()
    model_y = model_y[indices].flatten()
    x = x[indices].flatten()
    data_integral = data_y.sum()
    model_integral = model_y.sum()
    sigmas = sigmas / data_integral
    model_y = model_y / model_integral

    col1 = torch.ones(x.size()) * model_y
    col2 = x**2 * model_y
    logging.debug("x, cols ", x.shape, col1.shape, col2.shape)
    logging.debug("x**2 ", (x**2).shape)
    at = torch.stack([col1, col2])
    logging.debug("at ", at.shape)
    a = at.transpose(0,-1)
    logging.debug("a", a.shape)
    logging.debug("sigmas ", sigmas.shape)
    vm1 = torch.diag(torch.pow(sigmas, -2))
    logging.debug("vm1 ", vm1.shape)
    tmp = torch.matmul(vm1, a)
    logging.debug("tmp ", tmp)
    result = torch.matmul(at, tmp)
    return torch.inverse(result)


def errors_3d(cos_theta, phi, data_y, model_y, sigmas, debug = False):

    logging.debug("model_y ", model_y.shape)

    mask1 = data_y > 0
    mask2 = model_y > 0
    mask = torch.logical_and(mask1, mask2)
    indices = torch.nonzero(mask)
    data_y = data_y[indices].flatten()
    sigmas = sigmas[indices].flatten()
    model_y = model_y[indices].flatten()
    cos_theta = cos_theta[indices].flatten()
    phi = phi[indices].flatten()
    data_integral = data_y.sum()
    model_integral = model_y.sum()
    sigmas = sigmas / data_integral
    model_y = model_y / model_integral

    theta = torch.acos(cos_theta)

    col1 = torch.ones(cos_theta.size()) * model_y
    col2 = cos_theta**2 * model_y
    col3 = torch.sin(2*theta) * torch.cos(phi) * model_y
    col4 = torch.sin(theta) ** 2 * torch.sin(2*phi) * model_y
    logging.debug("x, cols ", cos_theta.shape, col1.shape, col2.shape)
    logging.debug("x**2 ", (cos_theta**2).shape)
    at = torch.stack([col1, col2, col3, col4])
    logging.debug("at ", at.shape)
    a = at.transpose(0,-1)
    logging.debug("a", a.shape)
    logging.debug("sigmas ", sigmas.shape)
    vm1 = torch.diag(torch.pow(sigmas, -2))
    logging.debug("vm1 ", vm1.shape)
    tmp = torch.matmul(vm1, a)
    logging.debug("tmp ", tmp)
    result = torch.matmul(at, tmp)
    return torch.inverse(result)


def errors_1d_hists(hist_data, hist_mc):
    t_data = helicity_model_1d.hist_to_tensor(hist_data)
    t_mc = helicity_model_1d.hist_to_tensor(hist_mc)
    return errors_1d(t_data[0], t_data[1], t_mc[1], t_data[2])

def errors_3d_hists(hist_data, hist_mc):
    t_data = helicity_model_3d.hist_to_tensor(hist_data)
    t_mc = helicity_model_3d.hist_to_tensor(hist_mc)
    return errors_3d(t_data[0], t_data[1], t_data[2], t_mc[2], t_data[3])
