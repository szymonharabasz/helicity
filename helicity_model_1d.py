import torch

from hist_utils import hist_to_tensor


class Helicity1d(torch.nn.Module):
    def __init__(self, learn_norm):
        """
        In the constructor we instantiate four parameters and assign them as
        member parameters.
        """
        super().__init__()
        self.learn_norm = learn_norm
        self.lambda_theta = torch.nn.Parameter(torch.ones(()))
        if self.learn_norm:
            self.norm = torch.nn.Parameter(torch.ones(()))

    def forward(self, x):
        """
        In the forward function we accept a Tensor of input data and we must return
        a Tensor of output data. We can use Modules defined in the constructor as
        well as arbitrary operators on Tensors.
        """
        cosTheta = x[0]
        weight = 1.0 + self.lambda_theta * cosTheta ** 2
        if self.learn_norm:
            return torch.vstack([
                cosTheta,
                self.norm * weight * x[1],
                self.norm * weight * x[2]])
        else:
            return torch.vstack([
                cosTheta, weight * x[1],
                weight * x[2]])

    def string(self):
        """
        Just like any class in Python, you can also define custom method on PyTorch modules
        """
        return f'y = (1 + {self.lambda_theta.item()} x^2'


def chi2_loss(pred_y, y):
    mask1 = y[1] > 0
    mask2 = pred_y[1] > 0
    mask = torch.logical_and(mask1, mask2)
    indices = torch.nonzero(mask)
    data = y[1][indices]
    data_err = y[2][indices]
    model = pred_y[1][indices]
    data_integral = data.sum()
    model_integral = model.sum()
    data = data / data_integral
    data_err = data_err / data_integral
    model = model / model_integral
    # print("y[1] ", y[1])
    # print("y[1][indices] ", y[1][indices])
    # return (((y[1][indices] - pred_y[1][indices])/y[2][indices]) ** 2).sum()
    return (((data - model) / data_err) ** 2).sum()


def chi2_loss_learn_norm(pred_y, y):
    mask1 = y[1] > 0
    mask2 = pred_y[1] > 0
    mask = torch.logical_and(mask1, mask2)
    indices = torch.nonzero(mask)
    data = y[1][indices]
    data_err = y[2][indices]
    model = pred_y[1][indices]
    return (((data - model) / data_err) ** 2).sum()


def fit_simple(model, hist_data, hist_mc, n_epochs, lr, learn_norm):
    # hist_mc = hists_mc[0][hist_index]
    # hist_data = histsData_np[0][hist_index]

    train_x = hist_to_tensor(hist_mc)
    train_y = hist_to_tensor(hist_data)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    losses = []

    for epoch in range(n_epochs):
        optimizer.zero_grad()

        output = model(train_x)
        if learn_norm:
            loss = chi2_loss_learn_norm(output, train_y)
        else:
            loss = chi2_loss(output, train_y)
        loss.backward()
        optimizer.step()

        losses.append(loss.item())

    return losses

    # print("Epoch: {}, loss: {}, param: {}".format(epoch, loss.item(), [param.item() for param in model.parameters()]))
