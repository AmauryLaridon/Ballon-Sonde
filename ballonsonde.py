"""
Started on Sun May  10 2020

@author: Amaury Laridon
"""
####################################################ENVIRONNEMENT#########################################################################
class Environment:
    """Classe décrivant l'environnement dans lequel va évolué le ballon"""
    def __init__(self, T_0=273, P_0=1, M_air=0.02896, rho_0=1.275, g=9.76):
        self.T_0   = T_0   #K
        self.P_0   = P_0   #bar
        self.M_air = M_air #kg/mol
        selg.g     = g     #m/s^2
########################################################BALLON############################################################################
class Ballon:
    """Classe décrivant le ballon atmosphérique"""
    def __init__(self, V_max, V_max_reached, gaz, M_gaz, P_gaz, M_str):
        self.V_max         = V_max
        self.V_max_reached = False
        self.gaz           = gaz
        self.M_gaz         = M_gaz
        self.P_gaz         = P_gaz
        self.M_str         = M_str

    def ballon_pdt(self, V_max, V_max_reached, gaz, M_gaz, P_gaz, M_str):
        """Permet de définir un ballon avec des valeurs par défaut"""
        self.V_max = 0.5
        self.gaz   = 'He'
        self.M_gaz = 0.06
        self.P_gaz = 1.2
        self.M_str = 0.3

######################################################  EDOSOLVER#########################################################################

import numpy as np
import csv
import matplotlib.pyplot as plt

class EDOSolver:
    """Fixe les conditions initiales et paramètre de la fonction à résoudre"""
    def __init__(self, function, t_0, t_f, Y_0, V_0, h):
        self.function = function
        self.t_0 = t_0
        self.t_f = t_f
        self.Y_0 = Y_0
        self.V_0 = V_0
        self.h   = h

        self.time_stored =[]
        self.pos_stored = []
        self.speed_stored = []
        self.method=""
        #self.tt=[]
        #self.yy=[]
        #self.position = [y[0] for y in self.Y_0]
        #self.vitesse = [y[1] for y in self.V_0]

    def solve(self, methode):
        """Implémente la méthode virtuelle afin de résoudre le problème"""
        raise NotImplemented

    def draw(self, file_name):
        #ECRITURE FICHIER
        with open(file_name,'w') as file:
            writer = csv.writer(file)
            writer.writerow(["temps","position", "vitesse"])
            for t, y, v in zip(self.time_stored, self.pos_stored, self.speed_stored):
                writer.writerow([t,y,v])
        #AFFICHAGE
        plt.subplot(1, 2, 1)
        plt.plot(self.time_stored,self.pos_stored, label='position')
        plt.title(self.method)
        plt.title("Graphique de la position \n en fonction du temps")
        plt.xlabel("Temps (s)")
        plt.ylabel("Position (m)")
        plt.legend()
        plt.grid()
        plt.subplot(1, 2, 2)
        plt.plot(self.time_stored, self.speed_stored, label='vitesse')
        plt.legend()
        plt.title("Graphique de la vitesse \n en fonction du temps")
        plt.xlabel("Temps (s)")
        plt.ylabel("Vitesse (m/s)")
        plt.grid()
        plt.show()

class Rk_4(EDOSolver):
    """Résolution d'une équation différentielle grâce à la méthode RK4"""
    def solve(self):
        #Conditions initiales
        self.Y = [self.Y_0]
        self.X = [self.t_0]
        self.V = [self.V_0]
        h = self.h
        t = self.t_0
        #Calcul
        while t < self.t_f :
        	Y_ = self.Y[-1]
        	V_ = self.V[-1]
        	V_inst, A_inst = self.function(Y_,V_)
        	K1, J1 = self.function(Y_,V_)
        	K2, J2 = self.function(Y_+(h*0.5*K1),V_+h*0.5*J1)
        	K3, J3 = self.function(Y_+(h*0.5*K2),V_+h*0.5*J2)
        	K4, J4 = self.function(Y_+h*K3,V_+h*J3)
        	y = Y_ + (h /6)*(K1 + 2*K2 + 2*K3 + K4)
        	v = V_ + (h/6)*(J1 + 2*J2 + 2*J3 + J4)
        	self.Y.append(y)
        	self.V.append(v)
        	t += h
        	self.X.append(t)
        #Variables intermédiaire pour l'affichage
        self.time_stored = self.X
        self.pos_stored = self.Y
        self.speed_stored = self.V
        self.method = "Rk4"
        #########################################
        return self.X, self.Y, self.V
