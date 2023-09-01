#!/usr/bin/env python
# coding: utf-8

# In[1]:


from bins import Bins
from utils import calcAllChi2, calcOneChi2, HistMaker, diffHist
from ROOT import TFile, TH1, TCanvas, TStyle, gStyle, Form
import math


# In[2]:


bins = Bins.readFrom("ranges.yml")
histMakerData = HistMaker("088_090_ag123ag_2500A_accepted_np.dat", "_data", bins)
#histsData = histMakerData.makeHists()
histsData = histMakerData.makeHists()
#histMakerMC = HistMaker("medium_isotropic_eff_ag1230ag_np_9deg.dat", "_MC", bins    )
#histMakerMC = HistMaker("mar19_diele_inmedium_heli0cm_np.dat", "_MC", bins    )
histMakerMC = HistMaker("mar19_diele_pi0_heli0cm_np.dat", "_MC", bins    )

allHistsMC = []


# In[3]:


import torch

torch.device("mps")
torch.manual_seed(0)

import gpytorch
import botorch

import matplotlib.pyplot as plt

plt.style.use("bmh")
plt.rcParams["figure.figsize"] = (8, 6)

from tqdm.notebook import tqdm

import warnings


# In[4]:

torch.backends.mps.is_available()
mps_device = torch.device("mps")

TH1.SetDefaultSumw2

N_PARAMS = 3


# In[5]:


#lb = -1
#ub = 1

#bounds = torch.tensor([[lb]*N_PARAMS, [ub]*N_PARAMS], dtype=torch.float)
bounds = torch.tensor([[0, 0,         0],
                       [1, 2*math.pi, 2]], dtype=torch.float)

grid_r   = torch.linspace(bounds[0][0], bounds[1][0], 101)
grid_phi = torch.linspace(bounds[0][1], bounds[1][1], 101)
grid_z   = torch.linspace(bounds[0][2], bounds[1][2], 101)

grid_x1, grid_x2, grid_x3 = torch.meshgrid(grid_r, grid_phi, grid_z, indexing="ij")

xs = torch.vstack([grid_x1.flatten(), grid_x2.flatten(), grid_x3.flatten()]).transpose(-2,-1)


# In[6]:


class GPModel(gpytorch.models.ExactGP, botorch.models.gpytorch.GPyTorchModel):
    _num_outputs = 1

    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=2.5, ard_num_dims=1)
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        mean_x = mean_x.to(mps_device)
        covar_x = covar_x.to(mps_device)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def fit_gp_model(train_x, train_y, num_train_iters=500):
    train_x = train_x.to(mps_device)
    train_y = train_y.to(mps_device)

    # declare the GP
    noise = 1e-4

    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = GPModel(train_x, train_y, likelihood)
    model.likelihood.noise = noise
    likelihood = likelihood.to(mps_device)
    model = model.to(mps_device)

    # train the hyperparameter (the constant)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

    model.train()
    likelihood.train()

    for i in range(num_train_iters):
        optimizer.zero_grad()

        output = model(train_x)
        loss = -mll(output, train_y)

        loss.backward()
        
        optimizer.step()

    model.eval()
    likelihood.eval()

    return model, likelihood


# In[7]:


num_queries = 175
num_repeats = 1
num_samples = 1


# In[8]:


def lambdas(xx):
   # return x[0], x[1], x[2]
    r, phi, z = xx[0], xx[1], xx[2]
    x = r*math.cos(phi)
    y = r*math.sin(phi)
    lambda_theta = 0.5 * (2*x + z)
    lambda_phi   = 0.25 * (-2 - 2*x + 3*z)
    lambda_theta_phi = y/math.sqrt(2.)
    return lambda_theta, lambda_phi, lambda_theta_phi

def all_lambdas(xx):
    def generator(xx):
        for x in xx:
            lambda_theta, lambda_phi, lambda_theta_phi = lambdas(x)
            yield torch.tensor([lambda_theta, lambda_phi, lambda_theta_phi])
           # yield torch.tensor([1./(chi2 / ndf)])
    return torch.stack([a for a in generator(xx)])
    
HIST_INDEX = 1

def objective(xx):
    def generator(xx):
        for x in xx:
            lambda_theta, lambda_phi, lambda_theta_phi = lambdas(x)

            histsMC = histMakerMC.makeHists(lambda_theta, lambda_phi, lambda_theta_phi)
            chi2, ndf = calcOneChi2(histsMC[0][HIST_INDEX], histsData[0][HIST_INDEX])
            allHistsMC.append(histsMC[0][HIST_INDEX])
            if not chi2 or not ndf:
                return torch.tensor([0])
            yield torch.tensor([1.0/(chi2 / ndf)])
           # yield torch.tensor([1./(chi2 / ndf)])
    return torch.stack([a for a in generator(xx)])


# In[9]:


# -2 is the default value when no feasible has been found
default_value = 0
feasible_incumbents = torch.ones((num_repeats, num_queries)) * default_value

best_fs = []

for trial in range(num_repeats):
    print("trial", trial)

    torch.manual_seed(trial)
    train_x = bounds[0] + (bounds[1] - bounds[0]) * torch.rand(num_samples, 3)
###                inequality_constraints=[(torch.tensor([0,2]),torch.tensor([-2.0,-1.0],dtype=torch.float),-2.0)],
    while -2.0*train_x[0][0] - 1.0*train_x[0][2] < -2:
        train_x = bounds[0] + (bounds[1] - bounds[0]) * torch.rand(num_samples, 3)

    train_utility = objective(train_x).to(mps_device)
    train_x = train_x.to(mps_device)

    for i in tqdm(range(num_queries)):
        feasible_incumbents[trial, i] = train_utility.max()
        print("train_x, train_utikity", train_x, train_utility)

        utility_model, utility_likelihood = fit_gp_model(
            train_x, train_utility.squeeze(-1)
        )
        print("utility model, utility likelihood = ", utility_model, utility_likelihood)
        best_f = train_utility.max()
        best_fs.append(best_f.item())
                
       # policy = botorch.acquisition.monte_carlo.qExpectedImprovement(
        policy = botorch.acquisition.analytic.LogExpectedImprovement(
       # policy = botorch.acquisition.analytic.ProbabilityOfImprovement(
            model=utility_model,
            best_f=train_utility.max(),
        ).to(mps_device)
        bounds = bounds.to(mps_device)

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            next_x, acq_val = botorch.optim.optimize_acqf(
                policy,
                bounds=bounds,
                q=1,
                num_restarts=40,
                raw_samples=100,
                inequality_constraints=[(torch.tensor([0,2],device=mps_device),torch.tensor([-2.0,-1.0],dtype=torch.float,device=mps_device),-2.0)],
            )

        next_utility = objective(next_x)

        train_x = torch.cat([train_x, next_x])
        train_utility = torch.cat([train_utility, next_utility])


# In[10]:


torch.save(feasible_incumbents, f"./incumbents.pth")
fout = TFile("out_"+str(HIST_INDEX)+".root","RECREATE")
fout.cd()
for hist in allHistsMC:
    print ("Writing hist: ", hist.GetName())
    hist.Write()
for j, hists in enumerate(histsData):
    for k, hist in enumerate(hists):
            hist.Write()
            print ("Writing hist: ", j, k, hist)
fout.Close()


# In[11]:


#t = torch.tensor( [ [-1.0000, -0.4678,  0.5433]])
#print(t)
#print(t.shape)


# In[12]:


with torch.no_grad():
    predictive_distribution = utility_likelihood(utility_model(xs))
    predictive_mean = predictive_distribution.mean
    predictive_lower, predictive_upper = predictive_distribution.confidence_region()
    acquisition_score = policy(xs.unsqueeze(1))


# In[109]:


from ipywidgets import interact, Layout, IntSlider
import numpy as np

axis_titles = [r"R of cone coordinates", r"$\Phi$ of cone coordinates", r"Z of cone coordinates"]

def oneplot(ax, tensor, index, cmap, title, transpose=None):
        
        global xs
        tensor_3d = torch.reshape(tensor, (101,101,101))
        xs_3d = torch.reshape(xs, (101,101,101,3))
        if transpose:
               tensor_3d = tensor_3d.transpose(transpose[0],transpose[1])
               xs_3d     = xs_3d    .transpose(transpose[0],transpose[1])
        if transpose is None:
                extent=[
                        xs_3d[0][0][0][2],
                        xs_3d[0][0][-1][2],
                        xs_3d[0][0][0][1],
                        xs_3d[0][-1][0][1],
                ]
        elif transpose == (0,1):
                extent=[
                        xs_3d[0][0][0][2],
                        xs_3d[0][0][-1][2],
                        xs_3d[0][0][0][0],
                        xs_3d[0][-1][0][0],
                ]
        elif transpose == (0,2):
                extent=[
                        xs_3d[0][0][0][0],
                        xs_3d[0][0][-1][0],
                        xs_3d[0][0][0][1],
                        xs_3d[0][-1][0][1],
                ]
        if transpose is None:
                x_title=axis_titles[2]
                y_title=axis_titles[1]
        elif transpose == (0, 1):
                x_title=axis_titles[2]
                y_title=axis_titles[0]
        elif transpose == (0,2):
                x_title=axis_titles[0]
                y_title=axis_titles[1]
        
        pos = ax.imshow(tensor_3d[index], cmap=cmap, interpolation="nearest", origin="lower", 
                vmin=0, vmax=tensor.max(), extent=extent)
        ax.set_aspect((extent[1]-extent[0])/(extent[3]-extent[2]))
        ax.set_title(title)
        ax.set_xlabel(x_title)
        ax.set_ylabel(y_title)
        plt.colorbar(pos,fraction=0.046, pad=0.04)

def ff(x, transpose):
        cmap = "gist_rainbow"
        fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True)
        oneplot(ax[0], predictive_mean, x, cmap, "value", transpose)
        oneplot(ax[1], predictive_upper-predictive_lower, x, cmap, "confidence int. width", transpose)
       # oneplot(ax[2], acquisition_score, x, cmap, "acquisition score", transpose)

def f(x):
        ff(x, None)

interact(f, x=IntSlider(50, 0, 100, 1, layout=Layout(width='500px')))


# In[110]:


def g(x):
        ff(x, (0,1))
           
interact(g, x=IntSlider(50, 0, 100, 1, layout=Layout(width='500px')))


# In[111]:


def h(x):
        ff(x, (0,2))
           
interact(h, x=IntSlider(50, 0, 100, 1, layout=Layout(width='500px')))


# In[112]:


lds = all_lambdas(train_x)
lambda_thetas     = lds.transpose(-2,-1)[0]
lambda_phis       = lds.transpose(-2,-1)[1]
lambda_theta_phis = lds.transpose(-2,-1)[2]
c = torch.stack((
    torch.arange(0,num_queries, dtype=int),
    lambda_thetas[0:num_queries],
    lambda_phis[0:num_queries],
    lambda_theta_phis[0:num_queries],
    train_utility.squeeze()[0:num_queries],
    (feasible_incumbents==feasible_incumbents.max())[0],
    feasible_incumbents[0]
),0).transpose(-2,-1)
torch.set_printoptions(precision=4,threshold=10_000, linewidth=120)
c_sorted = c[c[:, 4].sort()[1]]
print(c_sorted)


# In[113]:


gStyle.SetOptStat(0)

def plotComparison(histMC, histData):
    can = TCanvas("can","can",600,600)
    can.Divide(2,2)
    can.Draw()
    can.cd(1)
    histMC.Draw("COLZ")
    can.cd(2)
    histData.Draw("COLZ")
    pad = can.cd(3)
    pad.SetRightMargin(0.155)
    pad = can.cd(3)
    hdiff = diffHist(histMC, histData)
    hdiff.Draw("COLZ")
    can.Update()
    palette = hdiff.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.85)
    palette.SetX2NDC(0.9)
    can.Modified()
    can.Update()
    can.cd(4)
    nby = histMC.GetNbinsX()
    py = histMC.ProjectionY("py",1,nby)
    py.Draw("HIST")
    return can, hdiff, py

bestIndex = int(c_sorted[-1][0].item())
can1, hdiff1, py1 = plotComparison(allHistsMC[bestIndex], histsData[0][HIST_INDEX])


# In[114]:


worstIndex = int(c_sorted[0][0].item())
can2, hdiff2, py2 = plotComparison(allHistsMC[worstIndex], histsData[0][HIST_INDEX])


# In[115]:


base_hists = histMakerMC.builder.base_hists
can3 = TCanvas("c","C",900,1200)
can3.Divide(3,4)
can3.Draw()
for i, hist in enumerate(base_hists[0]):
    can3.cd(i+1)
    if not isinstance(hist, list):
        hist.Draw("COLZ")
    else:
        print(hist)


# In[116]:


can4 = TCanvas("c","C",900,1200)
can4.Divide(3,4)
can4.Draw()
for i, hist in enumerate(histsData[0]):
    can4.cd(i+1)
    if not isinstance(hist, list):
        hist.Draw("COLZ")
    else:
        print(hist)


# In[117]:


from ROOT import TH3F

pm = predictive_mean.reshape(101,101,101)
print(pm[100][100][100].item())
hpm = TH3F("hpm","predictive mean", 101, 0, 1, 101, 0, 2*math.pi, 101, 0, 2)
for binx in range(1,101):
    for biny in range(1,101):
        for binz in range(1,101):
            hpm.SetBinContent(binx,biny,binz,pm[binx-1][biny-1][binz-1].item())
hpm.SaveAs("pm.root")


# In[ ]:





# In[118]:


torch.log(-acquisition_score)
acquisition_score.max()


# In[119]:


tttt = torch.tensor([[1.0,5.0,1.0]])

def calcDiff2(histMC, histData, bx, by):
    contentMC = histMC.GetBinContent(bx+1,by+1)
    contentData = histData.GetBinContent(bx+1,by+1)
    errorData = histData.GetBinError(bx+1,by+1)
    if contentMC > 0 and contentData > 0:
        print(contentData, errorData)
        return pow((contentData - contentMC)/errorData, 2)
    else:
        return 0.0

def calcOneChi2_2(histMC, histData):
    chi2 = 0.
    ndf = 0
    for bx in range(histMC.GetXaxis().GetNbins()):
        for by in range(histMC.GetYaxis().GetNbins()):
            diff = calcDiff2(histMC, histData, bx, by)
            if diff > 0:
                ndf = ndf + 1
                chi2 = chi2 + diff
    return chi2, ndf

def objective2(xx):
    def generator(xx):
        for x in xx:
            lambda_theta, lambda_phi, lambda_theta_phi = lambdas(x)

            histsMC = histMakerMC.makeHists(lambda_theta, lambda_phi, lambda_theta_phi)
            chi2, ndf = calcOneChi2_2(histsMC[0][HIST_INDEX], histsData[0][HIST_INDEX])
            allHistsMC.append(histsMC[0][HIST_INDEX])
            if not chi2 or not ndf:
                return torch.tensor([0])
            yield torch.tensor([1.0/(chi2 / ndf)])
           # yield torch.tensor([1./(chi2 / ndf)])
    return torch.stack([a for a in generator(xx)])

objective2(tttt)


# In[120]:


best_fs


# In[ ]:





# In[ ]:




