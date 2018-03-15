from ROOT import TH1F, TH2F

import math


class Particle(object):

    def __init__(self, pdg_code, momentum, energy, energy_f, start, end, tracklen):
        self.pdg = pdg_code
        self.momentum = momentum
        self.energydep = energy
        self.energy_f = energy_f
        self.start = start
        self.end = end
        self.tracklen = tracklen
        self.mass = 0
        self.status = 1

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
            #print "muminus"
        elif particle.pdg == -13:
            self.muplus = particle
            #print "muplus"
        elif abs(particle.pdg) == 14:
            if particle.status == 1:
                self.nuout = particle
                #print "nuout"
            else:
                self.nuinc = particle
                #print "nuinc"
                
        return

    def TrackSep(self):

        muplus_end = self.muplus.end
        muminus_end = self.muminus.end

        return math.sqrt((muplus_end[0]-muminus_end[0])**2 + (muplus_end[1]-muminus_end[1])**2 + (muplus_end[2]-muplus_end[2])**2)

    def TrackAngle(self):

        muplus_phi = math.atan(self.muplus.momentum[1]/self.muplus.momentum[0])
        muminus_phi = math.atan(self.muminus.momentum[1]/self.muminus.momentum[0])

        muplus_theta = math.atan((math.sqrt(self.muplus.momentum[0]**2+self.muplus.momentum[1]**2)/self.muplus.momentum[2]))
        muminus_theta = math.atan((math.sqrt(self.muminus.momentum[0]**2+self.muminus.momentum[1]**2)/self.muminus.momentum[2]))

        d_phi = abs(muplus_phi-muminus_phi)
        d_theta = abs(muplus_theta-muminus_theta)

        return d_phi


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

    input_file = open("ArSMmumuDec11.txt", "r")
    #input_file = open("ArSMeeDec11.txt", "r")

    ## Skip heading
    for line in input_file:
        if "<event>" in line:
            break

    ## Define histograms
    energy_muplus  = TH1F("energy_muplus", "", 40, 0., 20.)
    energy_muminus = TH1F("energy_muplus", "", 40, 0., 20.)
    #energy_nuinc   = TH1F("energy_nuinc",  "", 40, 0., 20.)
    energy_nuout   = TH1F("energy_nuout",  "", 40, 0., 20.)
    corr_energy_mu = TH2F("h2_energies_mu", "", 20, 0., 10., 20, 0., 10.)

    event_count = 0

    ## Loop through all the events filling the corresponding histograms

    event = EventMuMu()

    for line in input_file:

        if "</event>" in line:
            ## End of event: Fill now the histograms.
            energy_muplus.Fill(event.muplus.energy())
            energy_muminus.Fill(event.muminus.energy())
            #energy_nuinc.Fill(event.nuinc.energy())
            energy_nuout.Fill(event.nuout.energy())
            corr_energy_mu.Fill(event.muplus.energy(), event.muminus.energy())
            if event.muplus.energydep > 100 and event.muminus.energydep > 100 \
               and event.muplus.tracklen > 50 and event.muminus.tracklen > 50 \
               and event.TrackSep() > 2 \
               and event.TrackAngle() >= 0.04:
                event_count += 1
            continue

        if "<event>" in line:
            ## New event: Reset
            event = EventMuMu()
            continue

        #print line

        words = line.split()
        pdg_code = int(words[0])
        momentum = (float(words[1]), float(words[2]), float(words[3])) # GeV
        start    = (float(words[4]), float(words[5]), float(words[6])) # cm
        end      = (float(words[7]), float(words[8]), float(words[9])) # cm
        energy   = float(words[10])
        energy_f = float(words[11])
        tracklen = float(words[12])

        particle = Particle(pdg_code, momentum, energy, energy_f, start, end, tracklen)
        event.AddParticle(particle)


    energy_muplus.Draw()
    energy_muminus.Draw("same")
    #energy_nuinc.Draw("same")
    energy_nuout.Draw("same")

    print "Number of events which passed cuts: ", event_count

    raw_input("Press any key to continue...")

    corr_energy_mu.Draw()

    raw_input("Press any key to continue...")
