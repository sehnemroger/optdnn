import numpy as np
import tensorflow as tf
from tensorflow import keras
from pickle import dump
import glob
import os 

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
     shuffle_buffer_size=10000, batch_size=512):
    dataset = tf.data.Dataset.list_files(filepaths)
    dataset = dataset.interleave(lambda filepath: tf.data.TextLineDataset(filepath), 
        cycle_length=n_readers, num_parallel_calls=tf.data.AUTOTUNE)
    dataset = dataset.map(preprocess, num_parallel_calls=tf.data.AUTOTUNE)
    dataset = dataset.shuffle(shuffle_buffer_size).repeat(repeat)
    return dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

train_set = dataset_reader(train_filepaths)
valid_set = dataset_reader(valid_filepaths)
test_set = dataset_reader(test_filepaths)


## Creating the DNN
layers = 14
inputs = keras.Input(shape=4)
mid = keras.layers.Dense(4, activation="linear")(inputs)
for i in range(layers):
    mid = keras.layers.Dense(10, activation="sigmoid")(mid)
    mid = keras.BatchNormalization()
outputs = keras.layers.Dense(2,  activation="linear")(mid)
model = tf.keras.Model(inputs, outputs)

## Creating custom loss Function
def custom_loss(y_true, y_pred):
    alpha = 1
    msle = tf.keras.losses.MeanSquaredLogarithmicError()
    mse = tf.keras.losses.MeanSquaredError()
    return alpha*msle(y_true, y_pred)+(1-alpha)*mse(y_true, y_pred)

# # preload treined model 
# model = tf.keras.models.load_model('dnn_trained_segundo.h5',custom_objects={ 'custom_loss': custom_loss})
## Compiling the model
lr_scheduler = keras.callbacks.ReduceLROnPlateau(factor=0.5, patience=2, verbose=1, min_lr=0.000001)
early_stop = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5, restore_best_weights=True)
model.compile(loss=custom_loss, metrics=['mse'],optimizer=keras.optimizers.Nadam(learning_rate=0.001))

## Training the dnn
# from tensorflow.python.client import device_lib
# print(tf.config.list_physical_devices('GPU'))
with tf.device('/GPU:0'):
    history = model.fit(train_set, epochs=1000, validation_data=valid_set,callbacks=[lr_scheduler,early_stop])


# ## Loading trained model
# model = tf.keras.models.load_model('dnn_trained_segundo.h5')

## Mean Square error on the test set
mse_test = model.evaluate(test_set)


## Saves the model 
model.save("dnn_trained_batchnorm.h5")

# ## Plot learning curves
# import pandas as pd
# import matplotlib.pyplot as plt
# pd.DataFrame(history.history).plot(figsize=(8, 5))
# plt.grid(True)
# plt.gca().set_ylim(0, 1) # set the vertical range to [0-1]
# plt.show()