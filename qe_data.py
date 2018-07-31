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
# Class qe_data reads the SCF.in or relax.in file and provide relevant electronic
# information for the user. The function position_convert() provide with the last coordinate snapshot to be 
# re-run or manipulated for the next calculation. It is in a QE compatible format.

class qe_data:
  
    def __init__(self, qe_out=None, qe_in=None):
        super(qe_data, self).__init__() #there is no inherited parent class?
        
        
        
        # the items we need to collect
        self.qe_version = [] 
        self.iterations = [] 
        self.elements = []
        self.energies = []
        self.fermi = []
        
        #for coordinates, we create a dictionary with elements being keys
        # If we store  coordinates and lattice as array, we can do manupulations on them.
        self.crd = [1, 2 ,3]
        
        if qe_out is None:
            qe_out="scf.out"
            
            
        # first load qe_out into memory     
        with open(qe_out, mode='r') as f:
            lines = [line for line in f.readlines()]
            
            
            # read stuff that doesn't change at each iteration
            for i, line in enumerate(lines):
                
                if '     Program PWSCF' in line:
                    self.qe_version = lines[i].split()[2]
                    
                elif 'bravais-lattice index     =' in line:
                    self.ibrav  = lines[i].split()[3]
                    self.alat   = lines[i+1].split()[4]
                    self.volume = lines[i+2].split()[3]
                    self.natoms = lines[i+3].split()[4]
                    self.ntypes = lines[i+4].split()[5]                    
                    self.nelect = lines[i+5].split()[4]            
                    self.nstates= lines[i+6].split()[4]      
                    self.ecut   = lines[i+7].split()[3]  
                    self.ecutwfc= lines[i+8].split()[4]    

            self.crd[0] = [[] for i in range(int(self.natoms))]
            self.crd[1] = [[] for i in range(int(self.natoms))]
            self.crd[2] = [[] for i in range(int(self.natoms))]
            for i, line in enumerate(lines):
                if '     site n.' in line:
                    stop = i
            for j in range(1, int(self.natoms)+1):
                self.elements.append(lines[stop+j].split()[1])
                        
                    
#             # read stuff that does change at each iteration
            for i, line in enumerate(lines):                                  
                self.nit = 0
                if '!    total' in line:
                    self.energies.append(float(lines[i].split()[4])) #in Rydberg units
                    self.nit += 1
                    self.fermi.append(lines[i-2].split()[4]) #in electron-volts
                if 'ATOMIC_POSITIONS' in line:
                    for j in range(int(self.natoms)):
                        self.crd[0][j].append(float(lines[i+j+1].split()[1]))
                        self.crd[1][j].append(float(lines[i+j+1].split()[2]))
                        self.crd[2][j].append(float(lines[i+j+1].split()[3]))
            for i in range(len(self.crd[0])):
                self.crd[0][i] = np.array(self.crd[0][i])
                self.crd[1][i] = np.array(self.crd[1][i])
                self.crd[2][i] = np.array(self.crd[2][i])
                
            self.final_position= np.array([np.array([self.crd[0][i][-1] for i in range(int(self.natoms))]), 
                                  np.array([self.crd[1][i][-1] for i in range(int(self.natoms))]),
                                  np.array([self.crd[2][i][-1] for i in range(int(self.natoms))])])
            self.final_position = self.final_position.T
            
            self.final_energy = self.energies[-1]
            
  ##  crd has three lists, each for a coordinate system (x, y ,z)
  ##  each list has natom number of numpy arrays
  ##  each numpy array has the coordinates (in float) from the order of beginnning to end
  ##  final_position has a list of x, y, z arrays of the final positions
  
    #########################################
    # TESTED
    # Example: if you want to shift the final position by 1, it would be:
    #           new_pos = qe_data.final_position + 1
    #
    # TESTED
    # If you wanted to extract the z coordinate of 11th element at the 5th iteration:
    #           position = qe_data.crd[2][10][4]
    #########################################
    
    
    #Takes a qe_data class and converts the last final coordinates to an appropriate
    #string to be used for the next relax.in
    #example code: TESTED
    #data=qe_data(qe_out="relax.out")
    #output = position_convert(data)
    #print(output)
    def position_convert(self, data):
        crd_array = data.crd
        output = ""
        element_list = data.elements
        for i in range(len(crd_array[0])):
            output = output + element_list[i] 
            for j in range(0, 3):
                output = output + "  " + str(crd_array[j][i][-1])
            output = output + "\n"
        return output