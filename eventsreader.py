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
        self.y = 0.
        self.theta = 0.
        self.phi = 0.

    def get_four_momentum(self, i):
        return self.four_momenta[i]

    def get_px(self, i):
        return self.get_four_momentum(i)[0]

    def get_py(self, i):
        return self.get_four_momentum(i)[1]

    def get_pz(self, i):
        return self.get_four_momentum(i)[2]

    def get_energy(self, i):
        return self.get_four_momentum(i)[3]


class EventsReader(ABC):

    @abstractmethod
    def __init__(self, filename, frame, ekin):
        self.events = []
        self.filename = filename
        self.frame = frame
        self.ekin = ekin

        self.mn = 9.39565612792968750e+02
        self.mp = 9.38272338867187528e+02
        self.mnucl = (79*self.mp + 118*self.mn)/197
        self.Ep = self.mnucl + self.ekin
        self.pp = TMath.Sqrt(self.Ep*self.Ep - self.mnucl*self.mnucl)

    @abstractmethod
    def read_events(self):
        pass

    def get_events(self):
        if not self.events:
            self.read_events()
        return self.events

    def get_event(self, i):
        return self.get_events()[i]


class SpdmeEventsReader(EventsReader):

    def __init__(self, filename, frame, ekin):
        super().__init__(filename, frame, ekin)

    def read_events(self):
        file = open(self.filename, "r")

        i = 0
        for line in file:
            if line.startswith("==>"):
                continue
            numbers = [float(token) for token in line.split()]
            if len(numbers) <= 0:
                continue
            if i % 4 == 0:
                if len(self.events) > 0:
                    self.set_event_properties(self.events[-1])

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

    def get_phis(self, dilep, dilep_cm, ev4_dilep):
        beam = TLorentzVector(0., 0., self.pp, self.Ep)
        target = TLorentzVector(0., 0., 0., self.mp)
        beam.Boost(-dilep.BoostVector())
        target.Boost(-dilep.BoostVector())
        beam_dir = beam.Vect()
        target_dir = target.Vect()
        beam_dir.SetMag(1.)
        target_dir.SetMag(1.)

        z_cs = beam_dir - target_dir
        y_cs = beam_dir.Cross(target_dir)
        x_cs = y_cs.Cross(z_cs)
        x_gj = y_cs.Cross(beam_dir)
        x_heli = y_cs.Cross(dilep_cm.Vect())
        z_cs.SetMag(1.)
        y_cs.SetMag(1.)
        x_cs.SetMag(1.)
        x_gj.SetMag(1.)
        x_heli.SetMag(1.)
        ev4_dilep_y = 1000 * ev4_dilep.Vect().Dot(y_cs)
        ev4_dilep_x_hl = 1000 * ev4_dilep.Vect().Dot(x_heli)
        ev4_dilep_x_cs = 1000 * ev4_dilep.Vect().Dot(x_cs)
        ev4_dilep_x_gj = 1000 * ev4_dilep.Vect().Dot(x_gj)
        phi_heli = TMath.ATan2(ev4_dilep_y, ev4_dilep_x_hl) + TMath.Pi()
        phi_cs = TMath.ATan2(ev4_dilep_y, ev4_dilep_x_cs) + TMath.Pi()
        phi_gj = TMath.ATan2(ev4_dilep_y, ev4_dilep_x_gj) + TMath.Pi()

        return phi_heli, phi_cs, phi_gj

    def set_event_properties(self, event):
        try:
            event.mass = math.sqrt(
                math.pow(event.get_energy(0), 2) -
                math.pow(event.get_px(0), 2) -
                math.pow(event.get_py(0), 2) -
                math.pow(event.get_pz(0), 2)
            )
        except ValueError:
            pass
           # print("ERROR! Wrong values: ", math.pow(event.get_energy(0), 2), math.pow(event.get_px(0), 2), math.pow(event.get_py(0), 2), math.pow(event.get_pz(0), 2))
           # print(event.get_energy(0), event.get_px(0), event.get_py(0), event.get_pz(0))

        gamma = event.get_four_momentum(0)
        lepton = event.get_four_momentum(2)
        # neutron = event.getFourMomentum(1)
        vec_gamma = TLorentzVector(gamma[0], gamma[1], gamma[2], gamma[3])
        vec_lepton = TLorentzVector(lepton[0], lepton[1], lepton[2], lepton[3])
        # vecNeutron = TLorentzVector(neutron[0], neutron[1], neutron[2], neutron[3])
        # print("masses: ", vec_lepton.M(), vecNeutron.M())

        beam = TLorentzVector(0., 0., self.pp, self.Ep)
        target = TLorentzVector(0., 0., 0., self.mnucl)

        vec_cm = target + beam
        # vec_cm = vec_gamma + vecNeutron
        # boost = vec_cm.BoostVector()
        # print("Boost vector: ", boost.X(), boost.Y(), boost.Z())
        vec_gamma_cm = TLorentzVector(vec_gamma)
        vec_lepton_cm = TLorentzVector(vec_lepton)
        vec_gamma_cm.Boost(-vec_cm.BoostVector())
        vec_lepton_cm.Boost(-vec_cm.BoostVector())
        vec_lepton_gamma = TLorentzVector(vec_lepton_cm)
        vec_lepton_gamma.Boost(-vec_gamma_cm.BoostVector())

        beam = TLorentzVector(0., 0., self.pp, self.Ep)
        target = TLorentzVector(0., 0., 0., self.mp)
        beam.Boost(-vec_gamma.BoostVector())
        target.Boost(-vec_gamma.BoostVector())
        beam_dir = beam.Vect()
        target_dir = target.Vect()
        beam_dir.SetMag(1.)
        target_dir.SetMag(1.)

        z_cs = beam_dir - target_dir

        if self.frame == Frame.HX:
            event.theta = vec_lepton_gamma.Angle(vec_gamma_cm.Vect())
            event.phi, _, _ = self.get_phis(vec_gamma, vec_gamma_cm, vec_lepton_gamma)
        elif self.frame == Frame.CS:
            event.theta = vec_lepton_gamma.Angle(z_cs)
            _, event.phi, _ = self.get_phis(vec_gamma, vec_gamma_cm, vec_lepton_gamma)
        else:
            event.theta = vec_lepton_gamma.Angle(beam_dir)
            _, _, event.phi = self.get_phis(vec_gamma, vec_gamma_cm, vec_lepton_gamma)

        event.z = TMath.Cos(vec_gamma_cm.Theta())
        event.y = vec_gamma_cm.Rapidity()
