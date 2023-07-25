import math
import yaml

class Bins:
    def __init__(self, m_min, m_max, z_min, z_max):
        self.m_min = m_min
        self.m_max = m_max
        self.z_min = z_min
        self.z_max = z_max

    def __str__(self) -> str:
        return "m-min: %f, m-max: %f, z_min: %f, z-max: %f" % (self.m_min, self.m_max, self.z_min, self.z_max)
    
    def zToStr(self, value):
        result = "p" if value >= 0 else "m"
        if abs(value) < 1.0:
            result = result + "0"
        else:
            result = result + "1"
        result = result + str(int(10*abs(value)))
        return result

    def suffix(self):
        return "_" + str(int(self.m_min)) + "to" + str(int(self.m_max)) + "_" + self.zToStr(self.z_min) + self.zToStr(self.z_max)

    @staticmethod
    def readFrom(filename):
        with open(filename,"r") as stream:
            try:
                borders = yaml.safe_load(stream)
                m_borders = borders["m-borders"]
                z_borders = borders["z-borders"]

                bins = []
                for i in range(len(m_borders)-1):
                    for j in range(len(z_borders)-1):
                        bins.append(Bins(m_borders[i], m_borders[i+1], z_borders[j], z_borders[j+1]))

            except yaml.YAMLError as exc:
                print(exc)

        return bins