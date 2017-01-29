# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 18:31:35 2014

@author: Andreas
"""



class reader():
   
#   functie __init__ is soort van constructor, waar ook alle klassevariabelen
#   geinitialiseerd worden. xpositions en ypositions staan in die vorm omdat 
#   dat de makkelijkste vorm is voor de ax.plot particles om de data in te lezen.
    def __init__(self, f):
        
#       timeStep is de stap waarmee frames gegenereerd worden. gezien de variabele
#       tijdsstap is dit nodig.
        self.timeStep = 0.01
        self.f = f
        self.tMaster = 0.01
        
#       line is een array van waarden, wordt meteen zo ingelezen.
        self.line = self.f.readline().split('\t')

        self.t = self.line[0]
        self.xpositions = [float(self.line[1]), float(self.line[3]), float(self.line[5])]
        self.ypositions = [float(self.line[2]), float(self.line[4]), float(self.line[6])]
        self.energy = float(self.line[7])*1e7
        
        
    def update(self):
        self.line = self.f.readline().split('\t')
#       gezien het kan zijn dat update() aangeroepen wordt als de file al ten 
#       einde is, in een try-catch zetten, die alle waarden op 0 zet als er een
#       valueerror is.
        try:
            
            self.t = float(self.line[0])

#           controle dat timeStep gerespecteerd wordt
            while (self.t < self.tMaster):
                self.line = self.f.readline().split('\t')
                self.t = float(self.line[0])
    
#           er wordt nieuwe data doorgegeven voor een frame, tMaster moet dus 
#           aangepast worden met een nieuwe timeStep.
            self.tMaster += self.timeStep
            
#           updaten van de andere gegevens in reader
            self.xpositions = [float(self.line[1]), float(self.line[3]), float(self.line[5])]
            self.ypositions = [float(self.line[2]), float(self.line[4]), float(self.line[6])]
            self.energy = float(self.line[7])
            
        except ValueError:
            self.t = 0
            self.xpositions = [0, 0, 0]
            self.ypositions = [0, 0, 0]
            self.energy = 0
        
#   sluit file binnen deze klasse (wordt niet gebruikt in uiteindelijk programma)
    def kill(self):
        self.f.close()
        
#   getters        
    def getT(self):
        return self.t
        
    def getEnergy(self):
        return self.energy
        
    def getPositions(self, index):
        return self.xpositions[index], self.ypositions[index]
        
#code voor het testen van deze klasse

#f = open("solVarT0_01NoReg.dat", 'r')
#rd = reader(f)
#
#for i in range(10):
#    print rd.getT()
#    print rd.tMaster
#    print rd.getPositions()
#    rd.update()
#
#rd.kill()
