import numpy as np
import tensorflow as tf
from tensorflow import keras
from pickle import dump
import glob
import os 
cos, sin = tf.math.cos, tf.math.sin
from tensorflow.python.keras.engine import data_adapter
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

# configurate to allow memory dynamic allocation on gpu
config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

path = os.path.abspath(os.getcwd())
files = glob.glob(path+'/data_csv/'+'*.txt') # Get all .txt files in the folder

# Randomizing data
p = np.random.permutation(len(files))

# Creating train and tests sets. Train set is 80% of the whole set.
# The index of the closest to 80% (or something like that)
index280 = int(np.floor(.8*len(files)))
idx_cp_train_filepaths  = p[:index280]
idx_test_filepaths = p[index280:]
test_filepaths =  [files[i] for i in idx_test_filepaths]

# Creating train and validation sets with 80% of the whole train set.
index280train = int(np.floor(.8*len(idx_cp_train_filepaths)))
idx_train_filepaths = idx_cp_train_filepaths[:index280train]
idx_valid_filepaths = idx_cp_train_filepaths[index280train:]
train_filepaths = [files[i] for i in idx_train_filepaths]
valid_filepaths = [files[i] for i in idx_valid_filepaths]

## Creating a dataset 
def preprocess(line): # Function that preprocess the line
    defs = [0.,0.,0.,0.,0.,0.,0.,0.,tf.constant([], dtype=tf.float32)]
    fields = tf.io.decode_csv(line, record_defaults=defs)
    x = tf.stack(fields[:-5])
    y = tf.stack(fields[-1:])
    return x,y



def dataset_reader(filepaths, repeat=1, n_readers=8,
     shuffle_buffer_size=10000, batch_size=256):
    dataset = tf.data.Dataset.list_files(filepaths)
    dataset = dataset.interleave(lambda filepath: tf.data.TextLineDataset(filepath), 
        cycle_length=n_readers, num_parallel_calls=tf.data.AUTOTUNE)
    dataset = dataset.map(preprocess, num_parallel_calls=tf.data.AUTOTUNE)
    dataset = dataset.shuffle(shuffle_buffer_size).repeat(repeat)
    return dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

train_set = dataset_reader(train_filepaths)
valid_set = dataset_reader(valid_filepaths)
test_set = dataset_reader(test_filepaths)
# x_train = np.empty((1,4))
# y_train = np.empty((1,1))
# for train_filepath in train_filepaths:
#     if os.stat(train_filepath).st_size != 0:
#         data = np.genfromtxt(train_filepath, delimiter=',')
#         x_train = np.concatenate((x_train,data[:,0:4]),axis=0)
#         y_train = np.concatenate((y_train,data[:,-1:]),axis=0)

# x_test = np.empty((1,4))
# y_test = np.empty((1,1))
# for test_filepath in test_filepaths:
#     if os.stat(test_filepath).st_size != 0:
#         data = np.genfromtxt(test_filepath, delimiter=',')
#         x_test = np.concatenate((x_test,data[:,0:4]),axis=0)
#         y_test = np.concatenate((y_test,data[:,-1:]),axis=0)

# x_valid = np.empty((1,4))
# y_valid = np.empty((1,1))
# for valid_filepath in valid_filepaths:
#     if os.stat(valid_filepath).st_size != 0:
#         data = np.genfromtxt(valid_filepath, delimiter=',')
#         x_valid = np.concatenate((x_valid,data[:,0:4]),axis=0)
#         y_valid = np.concatenate((y_valid,data[:,-1:]),axis=0)


## Model of the system
def invert_pend(t,x,u,l,I,mb,mc,at,ar):
    x1, x2, x3, x4 = x[:,0], x[:,1], x[:,2], x[:,3]
    g = 9.8
    C1 = l*mb
    C2 = I+l**2*mb
    C3 = mb + mc
    x2p = (g*C1**2*cos(x3)*sin(x3)+C2*(u-at*x2)-ar*C1*cos(x3)*x4-C1*C2*sin(x3)*x4**2)/(C2*C3-C1**2*cos(x3)**2)
    x4p = (g*C1*C3*sin(x3)+C1*cos(x3)*(u-at*x2)-ar*C3*x4-C1**2*cos(x3)*sin(x3)*x4**2)/(C2*C3-C1**2*cos(x3)**2)
    dxdt = tf.stack([x2, x2p, x4, x4p])
    return dxdt

## Creating custom loss Function
def custom_loss(y_true, y_pred):
    alpha = 0.1
    msle = tf.keras.losses.MeanSquaredLogarithmicError()
    mse = tf.keras.losses.MeanSquaredError()
    return alpha*msle(y_true, y_pred)+(1-alpha)*mse(y_true, y_pred)

## Createing a custom model to allow custom loss function with gradient of the dnn

class CustomModel(keras.Model):
    @tf.function
    def train_step(self, data):
        # These are the only transformations `Model.fit` applies to user-input
        # data when a `tf.data.Dataset` is provided.
        data = data_adapter.expand_1d(data)
        x, y, sample_weight = data_adapter.unpack_x_y_sample_weight(data)

        with tf.GradientTape() as tape: # Allowing automatic differentiation
            with tf.GradientTape() as tape2: # Allowing automatic differentiation
                tape2.watch(x)
                y_out = self(x, training=True) # Predicts, aka Forward pass
                tape2.watch(y_out)
                y_pred = y_out[:,0] # The first output is the predicted control.
                tape2.watch(y_pred)
                V = y_out[:,1] # The second output is the lyapunov function
            ## Computing loss function
            # NOTE:
            # This overwrites the loss function configured in compile()!!!!!!!!
            # The first part of the loss function is the old loss function which
            # Tries to make the predicted control close to the optimal control with the 
            # convex combination of the mse and msle errors.
            # The second part of the loss function is the part that tries to make the forth
            # output of the dnn a lyapunov function for the system.
            ## Calculating Vdot
            # definindo constantes
            l, I, mb, mc, at, ar = 0.3, 2, 1, 3, .2, .2
            # Gradient of the V output wrt inputs
            dVdx = tape2.gradient(V,x)
            
            # Obtained by the chain rule.
            # Using tensordot to calculate the tensor product.
            # Since the formula is dVdx.dot(f), but there are batch_size number of lines in dVdx we need a for loop
            # to account for each element in batch_size.
            Vdot_list = []
            if x.shape[0] is None:
                i=0
                Vdot_list.append(tf.tensordot(dVdx[i:i+1,:],invert_pend(0,x[i:i+1,:],y_pred[i],l,I,mb,mc,at,ar),axes=1)[0])
            else:
                for i in range(x.shape[0]):
                    Vdot_list.append(tf.tensordot(dVdx[i:i+1,:],invert_pend(0,x[i:i+1,:],y_pred[i],l,I,mb,mc,at,ar),axes=1)[0])
            Vdot = tf.stack(Vdot_list)
            
            # Calculate the loss function
            mse = tf.keras.losses.MeanSquaredError() 
            relu = tf.keras.layers.ReLU()
            pen = 10e3
            loss = custom_loss(y,y_pred) + pen*relu(-V) + pen*relu(Vdot) # Adds the lyapunov penalty 

        # Run backwards pass.
        self.optimizer.minimize(loss, self.trainable_variables, tape=tape)
        self.compiled_metrics.update_state(y, y_pred, sample_weight)
        # Collect metrics to return
        return_metrics = {}
        for metric in self.metrics:
            result = metric.result()
            if isinstance(result, dict):
                return_metrics.update(result)
            else:
                return_metrics[metric.name] = result
        return return_metrics
## Creating the DNN with the functional api
layers = 14
inputs = keras.Input(shape=4)
mid = keras.layers.Dense(4, activation="linear")(inputs)
for i in range(layers):
    mid = keras.layers.Dense(10, activation="sigmoid")(mid)
    mid = keras.BatchNormalization()
outputs = keras.layers.Dense(2,  activation="linear")(mid)
model = CustomModel(inputs, outputs)



# # preload treined model 
# model = tf.keras.models.load_model('dnn_trained_segundo.h5',custom_objects={ 'custom_loss': custom_loss})
## Compiling the model
lr_scheduler = keras.callbacks.ReduceLROnPlateau(factor=0.9, patience=5, verbose=1, min_lr=0.000001, monitor='mse')
early_stop = tf.keras.callbacks.EarlyStopping(monitor='mse', patience=20, restore_best_weights=True)
model.compile( metrics=['mse'],optimizer=keras.optimizers.Nadam(learning_rate=0.001))
model.summary()

## Training the dnn
# from tensorflow.python.client import device_lib
# print(tf.config.list_physical_devices('GPU'))
with tf.device('/GPU:0'):
    history = model.fit(train_set, epochs=1000, validation_data=valid_set,callbacks=[lr_scheduler,early_stop])

# with tf.device('/GPU:0'):
#     history = model.fit(x=x_train, 
#         y=y_train, 
#         epochs=1000, 
#         validation_data=(x_valid,y_valid),
#         callbacks=[lr_scheduler,early_stop])


# ## Loading trained model
# model = tf.keras.models.load_model('dnn_trained_segundo.h5')

## Mean Square error on the test set
mse_test = model.evaluate(test_set)
# mse_test = model.evaluate(x=x_test,y=y_test)


## Saves the model 
model.save("dnn_trained_batchnorm.h5")

# ## Plot learning curves
# import pandas as pd
# import matplotlib.pyplot as plt
# pd.DataFrame(history.history).plot(figsize=(8, 5))
# plt.grid(True)
# plt.gca().set_ylim(0, 1) # set the vertical range to [0-1]
# plt.show()