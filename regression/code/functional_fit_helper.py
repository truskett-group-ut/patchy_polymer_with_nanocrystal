from keras import backend as K
from keras.layers import Input, Dense, Lambda, concatenate
from keras.models import Model
import numpy as np

class FunctionalFitHelper:
    
    #build the various NN's used throughout and
    #store as separate models to employ later
    def InitializeNewModel(self, x_dim, num_cnts,
                           hidden_dim_rbf, activation_rbf,
                           hidden_dim_mod, activation_mod):
       
        #initialization details
        self.num_cnts = num_cnts
        self.x_dim = int(x_dim)
        input = Input(shape=(self.x_dim,))

        #create a single layer NN for the centers, amplitudes
        #and the exponential factors
        cnt_nns = []
        for i in range(self.num_cnts):
            cnt = Dense(hidden_dim_rbf, activation=activation_rbf, 
                        kernel_initializer='orthogonal')(input)
            cnt = Dense(1, activation='linear', 
                        kernel_initializer='orthogonal')(cnt)
            cnt_nns.append(cnt)

        amp_nns = []
        for i in range(self.num_cnts):
            amp = Dense(hidden_dim_rbf, activation=activation_rbf, 
                        kernel_initializer='orthogonal')(input)
            amp = Dense(1, activation='linear', 
                        kernel_initializer='orthogonal')(amp)
            amp_nns.append(amp)

        efa_nns = []
        for i in range(self.num_cnts):
            efa = Dense(hidden_dim_rbf, activation=activation_rbf, 
                        kernel_initializer='orthogonal')(input)
            efa = Dense(1, activation='linear', 
                        kernel_initializer='orthogonal')(efa)
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
        
        
        #use either a single NN to modify the output or a simple summation
        input = Input(shape=(self.num_cnts,))
        
        if hidden_dim_mod is None or activation_mod is None:
            rbfs = Lambda(lambda x: K.sum(x, axis=1, keepdims=True), name='rbfs_sum',
                        output_shape=(1,))(input)
        else:
            rbfs = Dense(hidden_dim_mod, activation=activation_mod)(input)
            rbfs = Dense(1, activation='linear')(rbfs)
       
        #model for the final NN, also for checkpointing
        self.model_rbfs = Model(input, rbfs, name='model_rbfs')

        return None
    
    #read in a saved model for re-instantiation
    def InitializeOldModel(self, 
                           model_cnts, 
                           model_amps, 
                           model_efas, 
                           model_rbfs):
        
        #get the various NNs
        self.model_cnts = model_cnts
        self.model_amps = model_amps
        self.model_efas = model_efas
        self.model_rbfs = model_rbfs

        #get stored input dimenions
        self.x_dim = self.model_cnts.layers[0].input_shape[1]

        return None
   
    #assemble the various NN's within the RBF functional form
    #as well as the final output layer which can be a NN too
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
            return K.abs(amps)*K.exp(-K.abs(efas)*K.pow(rs-cnts, 2))

        #rbfs wrapped in a lambda layer
        rbfs = Lambda(lambda x: rbfs(rs_, x), name='rbfs',
                        output_shape=(num_rs,))([cnts_, amps_, efas_])
       
        
        #custom function to extract the RBF signal data at a single r and flatten
        def ExtractRBFActsAndFlatten(i, input):
            return K.squeeze(input[:,:,i:i+1], axis=2)

        #extract the set of RBF activities for a given r in 
        #a flattened form (as row entries)
        rbfs_all = []
        for i in range(num_rs):
            rbfs_single = Lambda(lambda x: ExtractRBFActsAndFlatten(i, x),
                                         output_shape=(self.num_cnts,))(rbfs)
            rbfs_all.append(rbfs_single)

        #run each set through the final RBF NN
        rbfs_all_mod = []
        for rbfs_single in rbfs_all:
            rbfs_all_mod.append(self.model_rbfs(rbfs_single))

        #concatenate everything
        output = concatenate(rbfs_all_mod, axis=-1)

        #create the model and return
        return Model(input, output) 
