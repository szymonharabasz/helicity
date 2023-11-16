import math
from ROOT import TLorentzVector, TMath
from abc import ABC, abstractmethod
from enum import Enum

class Frame(Enum):
    HX = 1
    CS = 2
    GJ = 3

class Event:
    def __init__(self, weight):
        self.weight = weight
        self.four_momenta = []
        self.mass = 0.
        self.z = 0.
        self.theta = 0.
        self.phi = 0.

    def getFourMomentum(self, i):
        return self.four_momenta[i]

    def getPx(self, i):
        return self.getFourMomentum(i)[0]

    def getPy(self, i):
        return self.getFourMomentum(i)[1]

    def getPz(self, i):
        return self.getFourMomentum(i)[2]

    def getE(self, i):
        return self.getFourMomentum(i)[3]

class EventsReader(ABC):

    @abstractmethod
    def __init__(self, filename, frame, ekin):
        self.events = []
        self.filename = filename
        self.frame = frame

        self.mn = 9.39565612792968750e+02
        self.mp = 9.38272338867187528e+02
        self.mnucl = (79*self.mp + 118*self.mn)/197
        self.Ep = self.mnucl + 1230
        self.pp = TMath.Sqrt(self.Ep*self.Ep - self.mnucl*self.mnucl)

    @abstractmethod
    def readEvents(self):
        pass

    def getEvents(self):
        if not self.events:
            self.readEvents()
        return self.events

    def getEvent(self, i):
        return self.getEvents()[i]


class SpdmeEventsReader(EventsReader):

    def __init__(self, filenam):
        super().__init__(filenam)

    def readEvents(self):
        file = open(self.filename, "r")

        i = 0;
        for line in file:
            if line.startswith("==>"):
                continue
            numbers  = [float(token) for token in line.split()]
            if len(numbers) <= 0:
                continue
            if i % 4 == 0:
                if len(self.events) > 0:
                    self.setEventProperties(self.events[-1])

                try:
                    event_weight = numbers[1]
                except IndexError as e:
                    print("There was an error ", e)
                    print("Line read from file: ", line)
               # print(numbers)
                self.events.append(Event(event_weight))
               # print(self.events[-1].weight)
            else:
                self.events[-1].four_momenta.append(tuple(numbers))
            i = i + 1

    def getPhis(self, dilep, dilep_cm, ev4_dilep):
        beam = TLorentzVector(0., 0., self.pp, self.Ep)
        target = TLorentzVector(0., 0., 0., self.mp)
        beam.Boost(-dilep.BoostVector())
        target.Boost(-dilep.BoostVector())
        beam_dir = beam.Vect()
        target_dir = target.Vect()
        beam_dir.SetMag(1.)
        target_dir.SetMag(1.)

        zCS = beam_dir - target_dir
        yCS = beam_dir.Cross(target_dir)
        xCS = yCS.Cross(zCS)
        xGJ = yCS.Cross(beam_dir)
        xHeli = yCS.Cross(dilep_cm.Vect())
        zCS.SetMag(1.)
        yCS.SetMag(1.)
        xCS.SetMag(1.)
        xGJ.SetMag(1.)
        xHeli.SetMag(1.)
        ev4_dilep_y = 1000 * ev4_dilep.Vect().Dot(yCS);
        ev4_dilep_x_hl = 1000 * ev4_dilep.Vect().Dot(xHeli);
        ev4_dilep_x_cs = 1000 * ev4_dilep.Vect().Dot(xCS);
        ev4_dilep_x_gj = 1000 * ev4_dilep.Vect().Dot(xGJ);
        phiHeli = TMath.ATan2(ev4_dilep_y, ev4_dilep_x_hl) + TMath.Pi()
        phiCS = TMath.ATan2(ev4_dilep_y, ev4_dilep_x_cs) + TMath.Pi()
        phiGJ = TMath.ATan2(ev4_dilep_y, ev4_dilep_x_gj) + TMath.Pi()

        return phiHeli, phiCS, phiGJ

    def setEventProperties(self, event):
        event.mass = math.sqrt(
            math.pow(event.getE(0), 2) - 
            math.pow(event.getPx(0), 2) - 
            math.pow(event.getPy(0), 2) - 
            math.pow(event.getPz(0), 2)
        )

        gamma = event.getFourMomentum(0)
        lepton = event.getFourMomentum(2)
        neutron = event.getFourMomentum(1)
        vecGamma = TLorentzVector(gamma[0], gamma[1], gamma[2], gamma[3])
        vecLepton = TLorentzVector(lepton[0], lepton[1], lepton[2], lepton[3])
        vecNeutron = TLorentzVector(neutron[0], neutron[1], neutron[2], neutron[3])
       # print("masses: ", vecLepton.M(), vecNeutron.M())

        beam = TLorentzVector(0.,0.,self.pp,self.Ep)
        target = TLorentzVector(0.,0.,0.,self.mnucl)

        vecCM = target + beam
       # vecCM = vecGamma + vecNeutron
        boost = vecCM.BoostVector()
       # print("Boost vector: ", boost.X(), boost.Y(), boost.Z())
        vecGamma_cm = TLorentzVector(vecGamma)
        vecLepton_cm = TLorentzVector(vecLepton)
        vecGamma_cm.Boost(-vecCM.BoostVector())
        vecLepton_cm.Boost(-vecCM.BoostVector())
        vecLepton_gamma = TLorentzVector(vecLepton_cm)
        vecLepton_gamma.Boost(-vecGamma_cm.BoostVector())

       # print("lepton:       ", vecLepton.X(), vecLepton.Y(), vecLepton.Z(), vecLepton.T())
       # print("lepton_cm:    ", vecLepton_cm.X(), vecLepton_cm.Y(), vecLepton_cm.Z(), vecLepton_cm.T())
       # print("lepton_gamma: ", vecLepton_gamma.X(), vecLepton_gamma.Y(), vecLepton_gamma.Z(), vecLepton_gamma.T())

        beam = TLorentzVector(0., 0., self.pp, self.Ep)
        target = TLorentzVector(0., 0., 0., self.mp)
        beam.Boost(-vecGamma.BoostVector())
        target.Boost(-vecGamma.BoostVector())
        beam_dir = beam.Vect()
        target_dir = target.Vect()
        beam_dir.SetMag(1.)
        target_dir.SetMag(1.)

        z_cs = beam_dir - target_dir

        if self.frame == Frame.HX:
            event.theta = vecLepton_gamma.Angle(vecGamma_cm.Vect())
            event.phi, _, _ = self.getPhis(vecGamma, vecGamma_cm, vecLepton_gamma, 1230)
        elif self.frame == Frame.CS:
            event.theta = vecLepton_gamma.Angle(z_cs)
            _, event.phi, _ = self.getPhis(vecGamma, vecGamma_cm, vecLepton_gamma, 1230)
        else:
            event.theta = vecLepton_gamma.Angle(beam_dir)
            _, _, event.phi = self.getPhis(vecGamma, vecGamma_cm, vecLepton_gamma, 1230)

        event.z = TMath.Cos(vecGamma_cm.Theta())
       # print(event.theta, event.phi, event.z)


