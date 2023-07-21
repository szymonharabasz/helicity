import math
from ROOT import TLorentzVector, TMath
from abc import ABC, abstractmethod
from typing import List, Tuple

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
    def readEvents(self, filename):
        pass

    def getEvent(self, i):
        return self.events[i]

    events = []

class SpdmeEventsReader(EventsReader):

    def readEvents(self, filename):
        file = open(filename, "r")
        lines = file.readlines()

        i = 0;
        for line in lines:
            numbers  = [float(token) for token in line.split()]
            if i % 4 == 0:
                if len(self.events) > 0:
                    self.setEventProperties(self.events[-1])

                event_weight = numbers[1]
               # print(numbers)
                self.events.append(Event(event_weight))
               # print(self.events[-1].weight)
            else:
                self.events[-1].four_momenta.append(tuple(numbers))
            i = i + 1


    def getPhiHeli(self, dilep, dilep_cm, ev4_dilep, ekin):                                                                                                                                    

        mp = 9.38272338867187528e+02
        Ep = mp + ekin
        pp = TMath.Sqrt(Ep*Ep - mp*mp)
        beam = TLorentzVector(0.,0.,pp,Ep)
        target = TLorentzVector(0.,0.,0.,mp)
        beam.Boost(-dilep.BoostVector())
        target.Boost(-dilep.BoostVector())
        beam_dir = beam.Vect()
        target_dir = target.Vect()
        beam_dir.SetMag(1.)
        target_dir.SetMag(1.)

        yCS = beam_dir.Cross(target_dir)
        xHeli = yCS.Cross(dilep_cm.Vect())
        yCS.SetMag(1.)
        xHeli.SetMag(1.)
        ev4_dilep_y = 1000*ev4_dilep.Vect().Dot(yCS);
        ev4_dilep_x_hl = 1000*ev4_dilep.Vect().Dot(xHeli);
        phiHeli = TMath.ATan2(ev4_dilep_y, ev4_dilep_x_hl) + TMath.Pi()
        return phiHeli


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

        vecCM = vecGamma + vecNeutron
        vecGamma_cm = TLorentzVector(vecGamma)
        vecLepton_cm = TLorentzVector(vecLepton)
        vecGamma_cm.Boost(-vecCM.BoostVector())
        vecLepton_cm.Boost(-vecCM.BoostVector())
        vecLepton_gamma = TLorentzVector(vecLepton_cm)
        vecLepton_gamma.Boost(-vecGamma_cm.BoostVector())

        event.theta = vecLepton_gamma.Angle(vecGamma_cm.Vect())
        event.phi = self.getPhiHeli(vecGamma, vecGamma_cm, vecLepton_gamma, 1230)
        event.z = TMath.Cos(vecGamma_cm.Theta())
       # print("In this event: ", event.theta, event.phi, event.z)


