import dill
from keras.models import model_from_json
from functional_fit_helper import FunctionalFitHelper
import numpy as np

class NNRegressor:
        
    def __init__(self, num_rbfs):
        
        self.functional_helper = FunctionalFitHelper()
        self.functional_helper.num_cnts=num_rbfs
        
        file_base = '../regression/code'
        
        with open("{}/TrainedNNs/r_scaler.dill".format(file_base)) as r_scaler_file:
            self.r_scaler = dill.load(r_scaler_file)
        with open("{}/TrainedNNs/params_scaler.dill".format(file_base)) as params_scaler_file:
            self.params_scaler = dill.load(params_scaler_file)
            
        self.model_dict={}

        for mod in ['fe', 'cnts', 'amps', 'efas', 'rbfs']:
            modelfile='{}/TrainedNNs/model_'.format(file_base) + mod + '.json'
            json_file = open(modelfile, 'r')
            model_json = json_file.read()
            json_file.close()
            self.model_dict[mod] = model_from_json(model_json)
            weightsfile='{}/TrainedNNs/model_'.format(file_base) + mod + '.h5'
            self.model_dict[mod].load_weights(weightsfile)
            
        self.functional_helper.InitializeOldModel(self.model_dict['cnts'], self.model_dict['amps'], 
                                                  self.model_dict['efas'], self.model_dict['rbfs'])
                    
    def InitializeModel(self, r, rmax):
        if (max(r) - rmax) > 0.0001:
            print "I don't know anything about the RDF past r=", rmax, "; max(r)=", max(r)
            raise ValueError("invalid r range supplied")
        elif min(r) < 0.0:
            print "why are you trying to compute the RDF for negative r values, Ryan?? min(r)=", min(r)
            raise ValueError("invalid r range supplied")
        self.rnorm=np.transpose(self.r_scaler.fit_transform( np.transpose(np.array([r]))))[0]
        self.gr_model=self.functional_helper.BuildModel(self.rnorm)
        return None
    
    def predict_fe_gr(self, params, params_range):
        for i in range(0,5):
            if (params[i]-params_range[i][0] < -0.0001 or params[i]-params_range[i][1] > 0.0001):
                print "parameter", i, "out of trained range:", params[i], "not in", params_range[i]
                raise ValueError("parameter outside of trained range")
        eta_total = params[3]*(1.+params[2]**(-3.)*params[4]*params[0])
        if eta_total >= 0.5:
            print "model not trained total eta >= 0.5; eta=", eta_total
            raise ValueError("total system packing fraction is too high")
        params_norm = self.params_scaler.transform(np.array([params]))
        return self.model_dict['fe'].predict(params_norm)[0,0], self.gr_model.predict(params_norm)[0]

