import lmfit
import numpy as np
import matplotlib.pyplot as plt


class SinusFitDegree:

    def __init__(self, t, y):
        """ Fits a sinus to the given data, where t is the angle in degree
        """
        self.rawT = t
        self.rawY = y
        
        # 1.) Filter out isnan 
        mask = np.isnan(y)
        tL = np.ma.array(t, mask=mask).compressed()        # list -> array 
        yL = np.ma.array(y, mask=mask).compressed()

        # 2.) Check lenth 
        if (len(tL) <= 4):
            print("Less than 5 elements")
            self.a   = 0
            self.phi = 0
            self.c   = 0
            
        # 3a.) set up estmated values:
        a    = np.sqrt(np.mean((yL - np.mean(yL))**2)) * np.sqrt(2)
        iMax = np.argmax(yL)
        phi  = 0 
        c    = np.mean(yL)

        # 3b.) Put the estimated values to the dictionary
        #     AND set plausible ranges
        params = lmfit.Parameters()
        params.add('a',   value=a, min = 0.5 *a , max = 2. * a)
        params.add('phi', value=phi, min=-np.pi, max=np.pi)
        params.add('c',   value=c)

        # 4.) Fit: use nelde
        # 4a.) calc residuum
        def res(params, t, data):
            """ Calc residuum """
            v = params.valuesdict()
            return self._function(v['a'],v['phi'],v['c'], t) - data

        # 4b.) Do optimization
        result = lmfit.minimize(res, params, args=(tL, yL), method='nelder')  
   
        # 5.) Extract fit
        self.a   = result.params['a'].value
        self.phi = result.params['phi'].value
        self.c   = result.params['c'].value


    def _function(self, a, phi, c, t):
        """ The model function """
        return  a * np.sin(t/180.0 * np.pi + phi) + c   


    def calc(self, t):
        """ Evaluate model function """
        return self._function(self.a, self.phi, self.c, t)

    
    def getRawData(self):
        """ Return the raw data """
        return np.column_stack(self.rawT, self.rawY)


    def getFitParameter(self):
        """ Return the fit parameter """
        return np.array([self.a, self.phi, self.c])

    
    def checkFit(self):
        # 1.) Plot raw data
        plt.plot(self.rawT, self.rawY, '-o', label="raw data")

        # 2.) Plot model
        x = np.linspace(self.rawT.min(), self.rawT.max(), 200)
        plt.plot(x, self.calc(x), '-', label = "model")
        plt.legend()
    
        

def unitTestSinusFitDegree():
    """ Unit test for the class SinusFit """
    x = np.linspace(0,360, 6)                                                   
    y = np.sin(x/180.0 * np.pi) + np.random.rand() * 0.1
    s = SinusFitDegree(x, y)
    plt.clf()
    s.checkFit() 
