import torch

from helicity_model_1d import hist_to_tensor


def errors_1d(x, data_y, model_y, sigmas):

    print("model_y ", model_y.shape)

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
    data_y = data_y / data_integral
    sigmas = sigmas / data_integral
    model_y = model_y / model_integral

    col1 = torch.ones(x.size()) * model_y
    col2 = x**2 * model_y
    print("x, cols ", x.shape, col1.shape, col2.shape)
    print("x**2 ", (x**2).shape)
    at = torch.stack([col1, col2])
    print("at ", at.shape)
    a = at.transpose(0,-1)
    print("a", a.shape)
    print("sigmas ", sigmas.shape)
    vm1 = torch.diag(torch.pow(sigmas, -2))
    print("vm1 ", vm1.shape)
    tmp = torch.matmul(vm1, a)
    print("tmp ", tmp)
    result = torch.matmul(at, tmp)
    return torch.inverse(result)


def errors_1d_hists(hist_data, hist_mc):
    t_data = hist_to_tensor(hist_data)
    t_mc = hist_to_tensor(hist_mc)
    return errors_1d(t_data[0], t_data[1], t_mc[1], t_data[2])

sigmas = torch.tensor([[5.4208e-05],
        [2.5426e-04],
        [5.2835e-04],
        [8.3908e-04],
        [1.0968e-03]])
sigmas = sigmas.flatten()
print(sigmas)
vm1 = torch.diag(sigmas)
print(vm1)
