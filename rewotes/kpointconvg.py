from __future__ import print_function

import numpy as np
import pandas as pd


class KPointConvgException(Exception):
    pass
        

class KPointConvg(object):
    def __init__(self, calculation, kmin=1, kmax=10):
        self.calculation = calculation
        self.kmin = kmin
        self.kmax = kmax
        self.data = pd.DataFrame(columns=['k1', 'k2', 'k3'])
        
    def find_convg(self, criteria='E_eV', tol=0.00001, debug=False, nproc=1):
        for i in range(self.kmin, self.kmax+1):
            print('running with k=[{} {} {}]...'.format(i, i, i))
            self.calculation.kpoints=[i, i, i]
            res = self.calculation.run(debug=debug, nproc=nproc, save_input='{}.{}.{}.in'.format(i, i, i), save_output='{}.{}.{}.out'.format(i, i, i))
            k = {'k1': i, 'k2': i, 'k3': i}
            res.update(k)
            self.data = self.data.append(res, ignore_index=True)
            if len(self.data) > 1:
                self.data['delta_{}'.format(criteria)] = self.data.diff()[criteria].abs()
                if np.any(self.data.loc[self.data['delta_{}'.format(criteria)] < tol]):
                    return list(self.data.loc[self.data['delta_{}'.format(criteria)] < tol, ['k1', 'k2', 'k3']].values[0])
        print('Did not converge')
        return None