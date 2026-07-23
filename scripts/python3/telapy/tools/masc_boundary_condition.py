"""
class of utilizing boundary conditions stated by case file
"""

import numpy as np

class BoundaryConditionCas:
    """
    main class
    """

    def __init__(self, chnl, bcID, nom):
        """
        Constructor
        @param chnl(obj) instance of class of study
        @param bcID(int) ID of the associated boundary condition
        @param nom(str) name of the associated boundary condition
        """

        self.bcID = bcID
        self.nom = nom
        self.size = 0

        self.type = chnl.masc.get('Model.Graph.Type', (self.bcID - 1))
        if self.type == 1:
            var1 = "Model.Graph.Discharge"
        elif self.type == 2:
            var1 = "Model.Graph.Level"
        elif self.type == 3:
            var1 = "Model.Graph.Discharge"
            var2 = "Model.Graph.Level"
        elif self.type == 7:
            var1 = "Model.Graph.InfLevel"
            var2 = "Model.Graph.SupLevel"

        _, self.size, _ = chnl.masc.get_var_size(var1, (self.bcID - 1))
        self.bc1_ini = np.ones([self.size], dtype=np.float64)
        if self.type == 3 or self.type == 7:
            self.bc2_ini = np.ones([self.size], dtype=np.float64)
        self.time_ini = np.ones([self.size], dtype=np.float64)
        for i in range(self.size):
            self.bc1_ini[i] = chnl.masc.get(var1, (self.bcID - 1), i)
            if self.type == 3 or self.type == 7:
                self.bc2_ini[i] = chnl.masc.get(var2, (self.bcID - 1), i)
            self.time_ini[i] = \
                chnl.masc.get('Model.Graph.Time', (self.bcID - 1), i)

    def interpolate_boundary_condtion(self, times):
        """
        Provides interpolated boundary conditions at given times
        @param times(1d array or float) given times
        @return bc1(, bc2)(1d array) interpolated boundary conditions
        """

        bc1 = np.interp(times, self.time_ini, self.bc1_ini)
        if self.type == 3 or self.type == 7:
            bc2 = np.interp(times, self.time_ini, self.bc2_ini)

        if self.type == 1 or self.type == 2:
            return bc1
        elif self.type == 3 or self.type == 7:
            return bc1, bc2