import torch
from hist_utils import calc_one_chi2, HistMaker1d, geom_avg1d, ratio_err, symmetrize
from eventsreader import Frame
from bins import Bins

class Helicity1d(torch.nn.Module):
    def __init__(self):
        """
        In the constructor we instantiate four parameters and assign them as
        member parameters.
        """
        super().__init__()
        self.lambda_theta = torch.nn.Parameter(torch.ones(()))
        self.norm = torch.nn.Parameter(torch.randn(()))

    def forward(self, x):
        """
        In the forward function we accept a Tensor of input data and we must return
        a Tensor of output data. We can use Modules defined in the constructor as
        well as arbitrary operators on Tensors.
        """
        return torch.vstack([
            x[0],
            self.norm * (1.0 + self.lambda_theta * x[0] ** 2) * x[1],
            self.norm * (1.0 + self.lambda_theta * x[0] ** 2) * x[2]])
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
    # print("y[1] ", y[1])
    # print("y[1][indices] ", y[1][indices])
    return (((y[1][indices] - pred_y[1][indices])/y[2][indices]) ** 2).sum()


def hist_to_tensor(hist):
    n = hist.GetNbinsX()
    columns = []
    for b in range(1, n+1):
        t = torch.tensor([
            hist.GetBinCenter(b),
            hist.GetBinContent(b),
            hist.GetBinError(b)
        ])
        columns.append(t)
    result = torch.stack(columns).transpose(0,1)
    return result

def fit_simple(model, hist_data, hist_mc, n_epochs, lr):
   # hist_mc = hists_mc[0][hist_index]
   # hist_data = histsData_np[0][hist_index]

    train_x = hist_to_tensor(hist_mc)
    train_y = hist_to_tensor(hist_data)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    for epoch in range(n_epochs):
        optimizer.zero_grad()

        output = model(train_x)
        loss = chi2_loss(output, train_y)

        loss.backward()
        optimizer.step()

       # print("Epoch: {}, loss: {}, param: {}".format(epoch, loss.item(), [param.item() for param in model.parameters()]))