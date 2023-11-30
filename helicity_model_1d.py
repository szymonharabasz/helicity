import torch


class Helicity1d(torch.nn.Module):
    def __int__(self):
        super().__int__()
        self.lambda_theta = torch.nn.Parameter(torch.randn())

    def forward(self, x):
        return torch.vstack(
            (1.0 + self.a * x[0] ** 2) * x[1],
            (1.0 + self.a * x[0] ** 2) * x[2])

    def string(self):
        return f'y = (1 + {self.a.item()} x^2'


def chi2_loss(pred_y, y):
    return (((y[0] - pred_y[0])/pred_y[1]) ** 2).sum()


def hist_to_tensor(hist):
    n = hist.getNbinsX()
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
