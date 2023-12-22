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
        self.lambda_theta = torch.nn.Parameter(torch.randn(()))
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
   # print("Loss:")
   # print(y[1])
   # print(pred_y[1])
   # print(pred_y[2])
    mask = pred_y > 0
    return (((y[1] - pred_y[1])/y[2]) ** 2).sum()


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

def fit_simple(model, hist_data, hist_mc):
   # hist_mc = hists_mc[0][hist_index]
   # hist_data = histsData_np[0][hist_index]

    train_x = hist_to_tensor(hist_mc)
    train_y = hist_to_tensor(hist_data)

    optimizer = torch.optim.Adam(model.parameters(), lr=0.05)
    N_EPOCHS = 200
    for epoch in range(N_EPOCHS):
        optimizer.zero_grad()

        output = model(train_x)
        loss = chi2_loss(output, train_y)

        loss.backward()
        optimizer.step()

        print("Epoch: {}, loss: {}, param: {}".format(epoch, loss.item(), [param.item() for param in model.parameters()]))

model = Helicity1d()
params = model.parameters()
for param in params:
    print("param: ", param)

ekin = 1580
frame = Frame.HX
bins = Bins.readFrom("ranges.yml")
HIST_INDEX = 4

filename_data_np = "063_088_ag158ag_3200A_accepted_np.dat"
filename_MC_rho_heli0 = "mar19_diele_inmedium_heli0cm_ag1580ag_jver22_np_newCuts.dat"

histMakerData_np = HistMaker1d(filename_data_np, "_data_np", bins, frame, ekin)
histMakerMC_rho_heli0 = HistMaker1d(filename_MC_rho_heli0, "_MC_rho_heli0", bins, frame, ekin)

histsData_np = histMakerData_np.make_hists  ()
simple_model = Helicity1d()
hist_data_simple = histsData_np[0][HIST_INDEX]
hist_mc_simple = histMakerMC_rho_heli0.make_hists(0.0)[0][HIST_INDEX]
fit_simple(simple_model, hist_data_simple, hist_mc_simple)