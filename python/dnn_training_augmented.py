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
# files = glob.glob(path + "/python/Results/data_augmented_db/" + "*.txt")
files = glob.glob(path + "/Results/data_augmented_db/" + "*.txt")

# Randomizing data
p = np.random.permutation(len(files))

# Creating train and tests sets. Train set is 80% of the whole set.
# The index of the closest to 80% (or something like that)
index280 = int(np.floor(0.8 * len(files)))
idx_cp_train_filepaths = p[:index280]
idx_test_filepaths = p[index280:]
test_filepaths = [files[i] for i in idx_test_filepaths]

# Creating train and validation sets with 80% of the whole train set.
index280train = int(np.floor(0.8 * len(idx_cp_train_filepaths)))
idx_train_filepaths = idx_cp_train_filepaths[:index280train]
idx_valid_filepaths = idx_cp_train_filepaths[index280train:]
train_filepaths = [files[i] for i in idx_train_filepaths]
valid_filepaths = [files[i] for i in idx_valid_filepaths]

## Creating a dataset
def preprocess(line):  # Function that preprocess the line
    defs = [0.0, 0.0, 0.0, 0.0, 0.0, tf.constant([], dtype=tf.float32)]
    fields = tf.io.decode_csv(line, record_defaults=defs)
    x = tf.stack(fields[1:-1])
    y = tf.stack(fields[-1:])
    return x, y


def dataset_reader(
    filepaths, repeat=1, n_readers=8, shuffle_buffer_size=10000, batch_size=256
):
    dataset = tf.data.Dataset.list_files(filepaths)
    dataset = dataset.interleave(
        lambda filepath: tf.data.TextLineDataset(filepath),
        cycle_length=n_readers,
        num_parallel_calls=tf.data.AUTOTUNE,
    )
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

# Creating custom loss Function
def custom_loss(y_true, y_pred):
    alpha = 0.1
    msle = tf.keras.losses.MeanSquaredLogarithmicError()
    mse = tf.keras.losses.MeanSquaredError()
    return alpha * msle(y_true, y_pred) + (1 - alpha) * mse(y_true, y_pred)


# Creating the DNN
layers = 10
inputs = keras.Input(shape=4)
mid = keras.layers.Dense(4, activation="linear")(inputs)
for i in range(layers):
    mid = keras.layers.Dense(20, activation="sigmoid", kernel_initializer="he_normal")(
        mid
    )
    mid = keras.layers.BatchNormalization()(mid)
outputs = keras.layers.Dense(1, activation="linear")(mid)
model = tf.keras.Model(inputs, outputs)
# tf.keras.layers.LeakyReLU(alpha=0.01)


# # preload treined model
# model = tf.keras.models.load_model('dnn_trained_segundo.h5',custom_objects={ 'custom_loss': custom_loss})
## Compiling the model
lr_scheduler = keras.callbacks.ReduceLROnPlateau(
    monitor="loss", factor=0.95, patience=3, verbose=1, min_lr=0.000000001
)
early_stop = tf.keras.callbacks.EarlyStopping(
    monitor="loss", patience=50, restore_best_weights=True
)
model.compile(
    loss=custom_loss,
    metrics=["mse"],
    optimizer=keras.optimizers.Nadam(learning_rate=0.001),
)
model.summary()

## Training the dnn
# from tensorflow.python.client import device_lib
# print(tf.config.list_physical_devices('GPU'))
with tf.device("/GPU:0"):
    history = model.fit(
        train_set,
        epochs=10000,
        validation_data=valid_set,
        callbacks=[lr_scheduler, early_stop],
    )

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
model.save("dnn_trained_batchnorm_augmented.h5")

# ## Plot learning curves
# import pandas as pd
# import matplotlib.pyplot as plt
# pd.DataFrame(history.history).plot(figsize=(8, 5))
# plt.grid(True)
# plt.gca().set_ylim(0, 1) # set the vertical range to [0-1]
# plt.show()
