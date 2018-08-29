from keras import backend as K
from keras.layers import Input, Dense, Lambda, concatenate
from keras.models import Model
import numpy as np

class FunctionalFitHelper:

    def InitializeNewModel(self, x_dim, num_cnts, hidden_dim, activation):
        
        #initialization details
        self.num_cnts = num_cnts
        self.x_dim = int(x_dim)
        input = Input(shape=(self.x_dim,)) 

        #create a single layer NN for the centers, amplitudes
        #and the exponential factors
        cnt_nns = []
        for i in range(self.num_cnts):
            cnt = Dense(hidden_dim, activation=activation)(input)
            cnt = Dense(1, activation='linear')(cnt)
            cnt_nns.append(cnt)

        amp_nns = []
        for i in range(self.num_cnts):
            amp = Dense(hidden_dim, activation=activation)(input)
            amp = Dense(1, activation='linear')(amp)
            amp_nns.append(amp)

        efa_nns = []
        for i in range(self.num_cnts):
            efa = Dense(hidden_dim, activation=activation)(input)
            efa = Dense(1, activation='linear')(efa)
            efa_nns.append(efa)

        #concatenate the amplitudes and standard deviations
        if num_cnts == 1:
            cnts = cnt_nns[0]
            amps = amp_nns[0]
            efas = efa_nns[0]
        else:
            cnts = concatenate(cnt_nns, axis=-1)
            amps = concatenate(amp_nns, axis=-1)
            efas = concatenate(efa_nns, axis=-1)

        #model for each NN for easy access in checkpoint
        self.model_cnts = Model(input, cnts, name='model_cnts')
        self.model_amps = Model(input, amps, name='model_amps')
        self.model_efas = Model(input, efas, name='model_efas')

        return None
    
    def InitializeOldModel(self, model):

        #get stored input dimenions
        self.x_dim = self.model_cnts.layers[0].input_shape[1]

        #get the various NNs
        self.model_cnts = model.get_layer('model_cnts')
        self.model_amps = model.get_layer('model_amps')
        self.model_efas = model.get_layer('model_efas')

        return None
    
    def BuildModel(self, rs):
        #extract relevant dimensions 
        input = Input(shape=(self.x_dim,))

        #input rs
        rs_ = K.variable(value=np.array([rs]))
        rs_ = K.repeat_elements(rs_, rep=self.num_cnts, axis=0)
        num_rs = len(rs)

        #centers, amplitude, and exponential factor NNs
        cnts = self.model_cnts(input)
        amps = self.model_amps(input)
        efas = self.model_efas(input)

        #centers
        cnts_ = Lambda(lambda x: K.expand_dims(x, axis=2), 
                       output_shape=(self.num_cnts,1,))(cnts)
        cnts_ = Lambda(lambda x: K.repeat_elements(x, rep=num_rs, axis=2), 
                       output_shape=(self.num_cnts,num_rs,))(cnts_)

        #amplitudes
        amps_ = Lambda(lambda x: K.expand_dims(x, axis=2), 
                       output_shape=(self.num_cnts,1,))(amps)
        amps_ = Lambda(lambda x: K.repeat_elements(x, rep=num_rs, axis=2), 
                       output_shape=(self.num_cnts,num_rs,))(amps_)

        #exponential factors
        efas_ = Lambda(lambda x: K.expand_dims(x, axis=2), 
                       output_shape=(self.num_cnts,1,))(efas)
        efas_ = Lambda(lambda x: K.repeat_elements(x, rep=num_rs, axis=2), 
                       output_shape=(self.num_cnts,num_rs,))(efas_)

        #custom function implementing the radial basis functions
        def rbfs(rs, input):
            cnts = input[0]
            amps = input[1]
            efas = input[2]
            return K.sum(amps*amps*K.exp(-efas*K.pow(rs-cnts, 2)), axis=1)

        #final transformation using rbfs wrapped in a lambda layer
        output = Lambda(lambda x: rbfs(rs_, x), name='rbfs', 
                        output_shape=(num_rs,))([cnts_, amps_, efas_])

        #create the model and return
        return Model(input, output)    