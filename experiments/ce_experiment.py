import os
import subprocess
import sys
from math import ceil
from numpy import arange

class Experiment:

    def __init__(self, path):
        self.path = path

        out_file = open("./" + path + "/results.csv", 'w')
        out_file.write(self.header_string())
        out_file.close()


        experiment_file = open("./"+path+"/experiment.txt", 'r')

        lines = experiment_file.readlines()

        self.settings = dict()

        self.m_key = 'm'

        for line in lines:
            vals = line.split(' ')

            if vals[0] == 'M':
                self.m_key = 'M'

            if vals[0] in ['n', 'r', 'm', 's', 'i']:
                if(len(vals) > 2):
                    self.settings[vals[0]] = range(int(vals[1]), int(vals[3]) + int(vals[2]), int(vals[2]))
                else:
                    self.settings[vals[0]] = [int(vals[1])]
            elif vals[0] in ['a', 'M', 't']:
                if(len(vals) > 2):
                    self.settings[vals[0]] = arange(float(vals[1]), float(vals[3]) + float(vals[2]), float(vals[2]))
                else:
                    self.settings[vals[0]] = [float(vals[1])]
            elif vals[0] == 'R':
                if(len(vals) > 2):
                    self.settings[vals[0]] = [pow(10, -1*x) for x in range(int(vals[1]), int(vals[3]) + int(vals[2]), int(vals[2]))]
                else:
                    self.settings[vals[0]] = [float(vals[1])]
            elif vals[0] == 'd':
                self.settings[vals[0]] = [s.strip() for s in vals[1:len(vals)]]

        experiment_file.close()
        print("read in file")
        print(self.settings)

    def begin(self):
        print("begin")
        
        for d in self.settings['d']:
            self.d = d
            for n in self.settings['n']:
                self.n = n
                for m in self.settings[self.m_key]:
                    self.m = m
                    if self.m_key == 'M':
                        self.m = ceil(m*n)
                    for r in self.settings['r']:
                        self.r = r
                        for R in self.settings['R']:
                            self.R = R
                            for a in self.settings['a']:
                                self.a = a;
                                self.run()

        

    def run(self):
        for file in os.listdir("./"+self.path+"/graphs/"):
            print(self.cl_string(file))
            subprocess.run(self.cl_string(file))

    def cl_string(self, file):
        return ["../cross_entropy.exe", 
                '-f', "./"+self.path+'/graphs/'+file,
                '-o', '-2', "./"+self.path+"/results.csv",
                '-n', str(self.n),
                '-m', str(self.m),
                '-r', str(self.r),
                '-R', str(self.R),
                '-a', str(self.a),
                '-d', str(self.d),
                '-s', str(self.settings['s'][0]),
                '-i', str(self.settings['i'][0]),
                '-t', str(self.settings['t'][0])]

    def header_string(self):
        return "graph, dom_type, vertices, edges, sample_sets, elite_sets, repeats_without_improvement, rho, alpha, best, time\n"

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("Please specify an experiment directory")
        exit(0)

    experiment = Experiment(sys.argv[1])
    experiment.begin()






