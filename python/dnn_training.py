import numpy as np
import tensorflow as tf
from tensorflow import keras
from pickle import dump
import glob
import os
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
# RUN IN CPU
# files = glob.glob(path + "/python/Results/*/" + "*.dat")  # Get all .txt files in the folder
# RUN IN GPU
files = glob.glob(path + "/Results/*/" + "*.dat")

complete_list = []
for file in files:
    with open(file) as iter_file:
        for line in iter_file:
            local_list = []
            line_splited = line.split()
            for item in line_splited:
                local_list.append(float(item))
            complete_list.append(local_list)

data = np.array(complete_list)

p = np.random.permutation(data.shape[0])
index280 = int(np.floor(0.8 * data.shape[0]))

idx_cp_train_file = p[:index280]
idx_test_file = p[index280:]
test_file = [data[i] for i in idx_test_file]

# # Creating train and validation sets with 80% of the whole train set.
index280train = int(np.floor(0.8 * len(idx_cp_train_file)))
idx_train_file = idx_cp_train_file[:index280train]
idx_valid_file = idx_cp_train_file[index280train:]
train_file = [data[i] for i in idx_train_file]
valid_file = [data[i] for i in idx_valid_file]

train_set = np.array(train_file)[:, 1:]
valid_set = np.array(valid_file)[:, 1:]
test_set = np.array(test_file)[:, 1:]

x_train = train_set[:, :-1]
x_test = test_set[:, :-1]
x_valid = valid_set[:, :-1]

y_train = train_set[:, -1]
y_test = test_set[:, -1]
y_valid = valid_set[:, -1]

## Creating the DNN
layers = 10
inputs = keras.Input(shape=4)
mid = keras.layers.Dense(4, activation="linear")(inputs)
for i in range(layers):
    mid = keras.layers.Dense(50, activation='sigmoid', kernel_initializer='he_normal')(mid)
    mid = keras.layers.BatchNormalization()(mid)
outputs = keras.layers.Dense(1, activation="linear")(mid)
model = tf.keras.Model(inputs, outputs)
#tf.keras.layers.LeakyReLU(alpha=0.01)
model.summary()

# Creating custom loss Function
def custom_loss(y_true, y_pred):
    alpha = 0.1
    msle = tf.keras.losses.MeanSquaredLogarithmicError()
    mse = tf.keras.losses.MeanSquaredError()
    return alpha * msle(y_true, y_pred) + (1 - alpha) * mse(y_true, y_pred)


# # preload treined model
# model = tf.keras.models.load_model('dnn_trained_segundo.h5',custom_objects={ 'custom_loss': custom_loss})
# Compiling the model
lr_scheduler = keras.callbacks.ReduceLROnPlateau(
    monitor='loss', factor=0.95, patience=3, verbose=1, min_lr=0.000000001
)
early_stop = tf.keras.callbacks.EarlyStopping(
    monitor='loss', patience=50, restore_best_weights=True
)
model.compile(
    loss=custom_loss,
    metrics=["mse"],
    optimizer=keras.optimizers.Nadam(learning_rate=0.001),
)

# Training the dnn
# from tensorflow.python.client import device_lib
# print(tf.config.list_physical_devices('GPU'))
with tf.device("/GPU:0"):
    history = model.fit(
        x=x_train,
        y=y_train,
        epochs=10000,
        batch_size=256,
        validation_data=(x_valid, y_valid),
        callbacks=[lr_scheduler, early_stop],
    )

# ## Loading trained model
# model = tf.keras.models.load_model('dnn_trained_segundo.h5')

## Mean Square error on the test set
mse_test = model.evaluate(x=x_test, y=y_test)

## Saves the model
model.save("dnn_trained_batchnorm.h5")

# ## Plot learning curves

# import pandas as pd
# import matplotlib.pyplot as plt
# pd.DataFrame(history.history).plot(figsize=(8, 5))
# plt.grid(True)
# plt.gca().set_ylim(0, 1) # set the vertical range to [0-1]
# plt.show()
