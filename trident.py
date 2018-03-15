from ROOT import TH1F, TH2F

import math


class Particle(object):

    def __init__(self, pdg_code, momentum, status):
        self.pdg = pdg_code
        self.momentum = momentum
        self.status = status
        self.mass = 0

        if abs(self.pdg) == 11:
            self.mass = 0.5109989461 / 1000. ## GeV
        elif abs(self.pdg) == 13:
            self.mass = 105.6583745 / 1000. ## GeV

    def mag_momentum(self):
        return math.sqrt(sum([i**2 for i in self.momentum]))

    def energy(self):
        return math.sqrt(self.mass**2 + self.mag_momentum()**2)


class EventMuMu(object):

    def __init__(self):
        self.muplus  = 0
        self.muminus = 0
        self.nuinc = 0
        self.nuout = 0

    def AddParticle(self, particle):

        if particle.pdg == 13:
            self.muminus = particle
            print "muminus"
        elif particle.pdg == -13:
            self.muplus = particle
            print "muplus"
        elif abs(particle.pdg) == 14:
            if particle.status == 1:
                self.nuout = particle
                print "nuout"
            else:
                self.nuinc = particle
                print "nuinc"

        return


class EventEE(object):

    def __init__(self):
        self.muplus  = 0
        self.muminus = 0
        self.nuinc = 0
        self.nuout = 0

    def AddParticle(self, particle):

        if particle.pdg == 11:
            self.muminus = particle
            print "muminus"
        elif particle.pdg == -11:
            self.muplus = particle
            print "muplus"
        elif abs(particle.pdg) == 14:
            if particle.status == 1:
                self.nuout = particle
                print "nuout"
            else:
                self.nuinc = particle
                print "nuinc"

        return



if __name__ == "__main__":

    #input_file = open("ArSMmumuDec11.dat", "r")
    input_file = open("ArSMeeDec11.dat", "r")

    ## Skip heading
    for line in input_file:
        if "<event>" in line:
            break

    ## Define histograms
    energy_muplus  = TH1F("energy_muplus", "", 40, 0., 20.)
    energy_muminus = TH1F("energy_muplus", "", 40, 0., 20.)
    energy_nuinc   = TH1F("energy_nuinc",  "", 40, 0., 20.)
    energy_nuout   = TH1F("energy_nuout",  "", 40, 0., 20.)
    corr_energy_mu = TH2F("h2_energies_mu", "", 20, 0., 10., 20, 0., 10.)

    ## Loop through all the events filling the corresponding histograms

    event = EventEE()

    for line in input_file:

        if "</event>" in line:
            ## End of event: Fill now the histograms.
            energy_muplus.Fill(event.muplus.energy())
            energy_muminus.Fill(event.muminus.energy())
            energy_nuinc.Fill(event.nuinc.energy())
            energy_nuout.Fill(event.nuout.energy())
            corr_energy_mu.Fill(event.muplus.energy(), event.muminus.energy())
            continue

        if "<event>" in line:
            ## New event: Reset
            event = EventEE()
            continue

        print line

        words = line.split()
        pdg_code = int(words[0])
        status   = int(words[1])
        momentum = (float(words[2]), float(words[3]), float(words[4])) ## GeV

        particle = Particle(pdg_code, momentum, status)
        event.AddParticle(particle)


    energy_muplus.Draw()
    energy_muminus.Draw("same")
    energy_nuinc.Draw("same")
    energy_nuout.Draw("same")

    raw_input("Press any key to continue...")

    corr_energy_mu.Draw()

    raw_input("Press any key to continue...")
