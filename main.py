import imp
from bins import Bins
from utils import calcAllChi2, makeHists

bins = Bins.readFrom("ranges.yml")
histsMC = makeHists("medium_isotropic_eff_ag1230ag_nn_9deg.dat", "_MC", bins)
histsData = makeHists("apr12_diele_088_090_ag123ag_2500A_accepted_nn_2.dat", "_data", bins)

calcAllChi2(histsMC, histsData)