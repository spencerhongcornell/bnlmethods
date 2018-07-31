# This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#Author: Spencer Hong and Mehmet Topsakal

# WHAT DOES THIS CLASS DO?
# Class vasp_data reads the OSZICAR/OUTCAR/CONTCAR and provide relevant electronic
# information for the user.


import numpy as np
import os
import sys
import pickle

class vasp_data:

    
    def __init__(self, vasp_out = 'OUTCAR', vasp_geometry = 'CONTCAR', vasp_pos = 'POSCAR'):
        super(vasp_data, self).__init__() 
        
        # the items we need to collect
        self.elements = []
        self.energies = []
    
        
        
        #first load OUTCAR into memory     
        with open(vasp_out, mode='r') as f:
            lines = [line for line in f.readlines()]
            
            self.vasp_version = lines[0].split()[0]
            for i, line in enumerate(lines):

                if 'ENCUT' in line:
                    self.encut = float(lines[i].split()[2])
                if 'volume of cell' in line:
                    self.volume = float(lines[i].split()[4])
                
            for i, line in enumerate(lines):
                if 'free energy    TOTEN' in line:
                    self.energies.append(float(lines[i].split()[4])) ## in eV
            
            self.final_energy = self.energies[-1]
        
        self.natoms = 0

        #load CONTCAR into memory
        with open(vasp_geometry, mode = 'r') as f:
            lines = [line for line in f.readlines()]
            for j in lines[6].split():
                self.natoms += int(j)
            self.elementorder = []
            
            self.finalcrd = [[] for i in range(int(self.natoms))]
            
            for j in range(0, len(lines[6].split())):
                for m in range(0, int(lines[6].split()[j])):
                    self.elementorder.append(lines[5].split()[j])
                    
            for n in range(0, 3):
                for k in range(0, self.natoms):
                
                    self.finalcrd[k].append(float(lines[k + 8].split()[n]))
        

        self.finalcrd[0] = np.array(self.finalcrd[0])
        self.finalcrd[1] = np.array(self.finalcrd[1])
        self.finalcrd[2] = np.array(self.finalcrd[2])
        
        self.finalcrd = np.array(self.finalcrd).T
        
        ####################
        #
        #
        # coordinate syntax: 
        # data = vasp_data()
        # to get the final z coordinate of the 7th Iron (to find the order of elements, call data.elementorder):
        #     data.finalcrd[2][6]
        #
        ####################
        
        with open(vasp_pos, mode = 'r') as f:
            lines = [line for line in f.readlines()]
            for j in lines[6].split():
                self.natoms += int(j)
            self.elementorder = []
            
            self.initcrd = [[] for i in range(int(self.natoms))]

            for n in range(0, 3):
                for k in range(0, self.natoms):
                
                    self.initcrd[k].append(float(lines[k + 9].split()[n]))
        
        self.initcrd[0] = np.array(self.initcrd[0])
        self.initcrd[1] = np.array(self.initcrd[1])
        self.initcrd[2] = np.array(self.initcrd[2])
        
        self.initcrd = np.array(self.initcrd).T
              
        
        
        #finds the fractional to cartesian matrix to convert fractional coordinates to cartesian
        #Current version only works with standard cells
    def frac2cartmatrix(self, a, b, c, alpha, beta, gamma,
                                           angle_in_degrees=True):
        """
        Return the transformation matrix that converts fractional coordinates to
        cartesian coordinates.
        Parameters
        ----------
        a, b, c : float
            The lengths of the edges.
        alpha, gamma, beta : float
            The angles between the sides.
        angle_in_degrees : bool
            True if alpha, beta and gamma are expressed in degrees.
        Returns
        -------
        r : array_like
            The 3x3 rotation matrix. ``V_cart = np.dot(r, V_frac)``.
        """
        if angle_in_degrees:
            alpha = np.deg2rad(alpha)
            beta = np.deg2rad(beta)
            gamma = np.deg2rad(gamma)
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
        cosg = np.cos(gamma)
        sing = np.sin(gamma)
        volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)
        r = np.zeros((3, 3))
        r[0, 0] = a
        r[0, 1] = b * cosg
        r[0, 2] = c * cosb
        r[1, 1] = b * sing
        r[1, 2] = c * (cosa - cosb * cosg) / sing
        r[2, 2] = c * volume / sing
        return r

        
    def cart2fracmatrix(self, a, b, c, alpha, beta, gamma,
                                           angle_in_degrees=True):
        """
        Return the transformation matrix that converts cartesian coordinates to
        fractional coordinates.
        Parameters
        ----------
        a, b, c : float
            The lengths of the edges.
        alpha, gamma, beta : float
            The angles between the sides.
        angle_in_degrees : bool
            True if alpha, beta and gamma are expressed in degrees.
        Returns
        -------
        r : array_like
            The 3x3 rotation matrix. ``V_frac = np.dot(r, V_cart)``.
        """
        if angle_in_degrees:
            alpha = np.deg2rad(alpha)
            beta = np.deg2rad(beta)
            gamma = np.deg2rad(gamma)
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
        cosg = np.cos(gamma)
        sing = np.sin(gamma)
        volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)
        r = np.zeros((3, 3))
        r[0, 0] = 1.0 / a
        r[0, 1] = -cosg / (a * sing)
        r[0, 2] = (cosa * cosg - cosb) / (a * volume * sing)
        r[1, 1] = 1.0 / (b * sing)
        r[1, 2] = (cosb * cosg - cosa) / (b * volume * sing)
        r[2, 2] = sing / (c * volume)
        return r

        # find the distance given two elements
        # input must be of the following format: "elementtype element number" such as "Fe 23" which is the 23rd Iron (not 23rd atom)
    def distance(self, element1 = None, element2 = None):
            first = element1.split()[0]
            second = element2.split()[0]

            first_type = element1.split()[1]
            second_type = element2.split()[1]
            counter1 = 0
            counter2 = 0
            for i, element in enumerate(self.elementorder):
                if element == first:
                    counter1 = i  + int(first_type)
                    #print(element)
                    #print(counter1)
                    break

            for i, element in enumerate(self.elementorder):
                if element == second:
                    counter2 = i + int(second_type)
                    #print(element)
                    #print(counter2)
                    break

            x1 = self.finalcrd[0][counter1]
            x2 = self.finalcrd[0][counter2]
            y1 = self.finalcrd[1][counter1]
            y2 = self.finalcrd[1][counter2]
            z1 = self.finalcrd[2][counter1]
            z2 = self.finalcrd[2][counter2]

            coord_1 = [x1, y1, z1]
            coord_2 = [x2, y2, z2]

            trans_matrix = self.frac2cartmatrix( 10.3, 10.3, 10.3, 90, 90, 90)
            mat1 = np.matmul(trans_matrix, coord_1)
            mat2 = np.matmul(trans_matrix, coord_2)
            
            # print(coord_1)
            # print(coord_2)
            # print(mat1)
            # print(mat2)

            return np.sqrt((mat1[0] - mat2[0])**2.0 + (mat1[1] - mat2[1])**2.0 + (mat1[2] - mat2[2])**2.0)

    # converts the last coordinate captured in CONTCAR for a string format to be used as the next POSCAR
    def position_convert(self):
        crd_array = self.finalcrd
        output = ""
        for i in range(len(crd_array[0])):
            for j in range(0, 3):
                output = output + "  " + str(crd_array[j][i])
            output = output + "\n"
        return np.array(output).T