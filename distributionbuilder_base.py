from abc import abstractmethod, ABC
from datetime import datetime

from ROOT import TH1F, TH2F

iter = 0


class DistributionBuilderBase(ABC):
    def __init__(self, histname_suffix, events, bins):
        self.histname_suffix = histname_suffix
        self.events = events
        self.bins = bins

    @abstractmethod
    def fill_angular(self, binIndex, event):
        pass

    @abstractmethod
    def add_angular(self):
        pass

    def binIndex(self, event, bins):
        index = 0
        for bin in bins:
            if event.mass >= bin.m_min and event.mass < bin.m_max and event.z >= bin.z_min and event.z < bin.z_max:
                return index
            index = index + 1
        return len(bins)

    def buildFromEvents(self):
        global iter

        now = datetime.now()
        current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print(f"[{current_time}]: iter ", iter)

        self.add_angular()

        hmassLowM = TH1F("hmassLowM" + self.histname_suffix, "hmassLowM" + self.histname_suffix, 100, 0, 1000)
        hzLowM = TH1F("hzLowM" + self.histname_suffix, "hzLowM" + self.histname_suffix, 100, -1, 1)
        hyLowM = TH1F("hyLowM" + self.histname_suffix, "hyLowM" + self.histname_suffix, 100, -2, 2)
        hmassHigM = TH1F("hmassHigM" + self.histname_suffix, "hmassHigM" + self.histname_suffix, 100, 0, 1000)
        hzHigM = TH1F("hzHigM" + self.histname_suffix, "hzHigM" + self.histname_suffix, 100, -1, 1)
        hyHigM = TH1F("hyHigM" + self.histname_suffix, "hyHigM" + self.histname_suffix, 100, -2, 2)
        counter = TH2F("counter" + self.histname_suffix, "counter" + self.histname_suffix, 4, -0.5, 3.5, 3, -0.5, 2.5)
        ievent = 0
        nevents = len(self.events)

        now = datetime.now()
        current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print(f"[{current_time}] Before processing events")
        print(f"Num events {nevents}")
        print(F"Num non-null events {len([evt for evt in self.events if evt is not None])}")
        for event in self.events:
            ievent = ievent + 1

            binIndex = self.binIndex(event, self.bins)

            if binIndex >= 0 and binIndex < len(self.angular_hists):
                self.fill_angular(binIndex, event)
                if binIndex < 3 and event.mass < self.bins[2].m_max:
                    hmassLowM.Fill(event.mass)
                    hzLowM.Fill(event.z)
                    hyLowM.Fill(event.y)
                if binIndex >= 3 and event.mass >= self.bins[2].m_max:
                    hmassHigM.Fill(event.mass)
                    hzHigM.Fill(event.z)
                    hyHigM.Fill(event.y)
                counter.Fill(binIndex // 3, binIndex % 3)

        now = datetime.now()
        current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print(f"[{current_time}] After processing events")
        for hist in [*self.angular_hists, hmassLowM, hzLowM, hyLowM]:
            if hist.Integral() > 0:
                pass
            # hist.Scale(1./hist.Integral())

        iter = iter + 1
        result = [self.angular_hists, [hmassLowM, hmassHigM], [hzLowM, hzHigM], [hyLowM, hyHigM], counter]
        print(f"#1 result length {len(result)}")
        return result

    def getHists(self):
        result = self.buildFromEvents()
        print(f"#2 result length {len(result)}")
        return result
