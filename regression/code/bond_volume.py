from keras.models import model_from_json
import dill
from functional_fit_helper import FunctionalFitHelper
import os
import numpy as np

class BondVolume:
    
    #read in all of the model components
    def __init__(self, model_index, num_rbfs):
        
        #check that the model index is valid
        #only one model for now, but leaving room for more
        file_base = '../regression/code/bond_volume/model_{}'.format(model_index)
        assert os.path.exists(file_base)
        
        #the range of r allowed for prediction
        self.r_range = np.loadtxt('{}/r_range.txt'.format(file_base))
        self.A_range = np.loadtxt('{}/A_range.txt'.format(file_base))
        
        #the scalers for the input and r values to predict at
        with open('{}/r_scaler.dill'.format(file_base), 'rb') as f:
            self.r_scaler = dill.load(f)
        with open('{}/x_scaler.dill'.format(file_base), 'rb') as f:
            self.x_scaler = dill.load(f)
            
        #read in the regressor components
        with open('{}/model_cnts.json'.format(file_base), 'rb') as f:
            model_cnts_json = f.read()
        with open('{}/model_amps.json'.format(file_base), 'rb') as f:
            model_amps_json = f.read()
        with open('{}/model_efas.json'.format(file_base), 'rb') as f:
            model_efas_json = f.read()
        with open('{}/model_rbfs.json'.format(file_base), 'rb') as f:
            model_rbfs_json = f.read()

        #read in the regressor components
        with open('{}/model_cnts.json'.format(file_base), 'rb') as f:
            self.model_cnts = model_from_json(f.read())
        with open('{}/model_amps.json'.format(file_base), 'rb') as f:
            self.model_amps = model_from_json(f.read())
        with open('{}/model_efas.json'.format(file_base), 'rb') as f:
            self.model_efas = model_from_json(f.read())
        with open('{}/model_rbfs.json'.format(file_base), 'rb') as f:
            self.model_rbfs = model_from_json(f.read())
        self.model_cnts.load_weights('{}/model_cnts.h5'.format(file_base))
        self.model_amps.load_weights('{}/model_amps.h5'.format(file_base))
        self.model_efas.load_weights('{}/model_efas.h5'.format(file_base))
        self.model_rbfs.load_weights('{}/model_rbfs.h5'.format(file_base))
        
        #insert the components into final model
        self.functional_fit_helper = FunctionalFitHelper()
        self.functional_fit_helper.InitializeOldModel(self.model_cnts, self.model_amps,
                                                      self.model_efas, self.model_rbfs)
        self.functional_fit_helper.num_cnts = num_rbfs
        
        return None
    
    #set the discretization
    def InitializePredictionGrid(self, rs):
        
        #generate zero data for outside training
        num_l = sum((rs <= self.r_range[0]))
        self.f_l = np.zeros(num_l)
        num_r = sum((rs >= self.r_range[1]))
        self.f_r = np.zeros(num_r)
        
        #create a mask for focusing on r-range to actually predict
        self.mask_c = (rs >= self.r_range[0]) * (rs <= self.r_range[1])
        
        #only need to explicitly work with the central rs
        rs_c = rs[self.mask_c]
        
        #the rs must be normalized first
        rs_c = np.transpose(np.array([rs_c]))
        rs_c_norm = np.transpose(self.r_scaler.transform(rs_c))[0]
        self.model = self.functional_fit_helper.BuildModel(rs_c_norm)
        
        return None
    
    #predict... yay
    def Predict(self, A):
        
        #check that A is in training range
        assert A >= self.A_range[0] and A <= self.A_range[1]
        
        #central portion
        x = self.x_scaler.transform(np.array([[A]]))
        y = self.model.predict(x)[0]
        f_c = np.exp(y)-1.0
        
        #whole thing
        f = np.concatenate((self.f_l, f_c, self.f_r))
        
        return f