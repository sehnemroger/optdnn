import glob
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from sklearn.svm import l1_min_c
import tensorflow as tf
from tensorflow.python.keras.engine import data_adapter
cos, sin = tf.math.cos, tf.math.sin

def controle(t,x):
    u = model.predict(np.array([x]))[0,0]
    return u

def invert_pend(t,x,l,I,mb,mc,at,ar):
    x1,x2,x3,x4 = x
    u = controle(t,x)
    g = 9.8
    C1 = l*mb
    C2 = I+l**2*mb
    C3 = mb + mc
    x2p = (g*C1**2*cos(x3)*sin(x3)+C2*(u-at*x2)-ar*C1*cos(x3)*x4-C1*C2*sin(x3)*x4**2)/(C2*C3-C1**2*cos(x3)**2)
    x4p = (g*C1*C3*sin(x3)+C1*cos(x3)*(u-at*x2)-ar*C3*x4-C1**2*cos(x3)*sin(x3)*x4**2)/(C2*C3-C1**2*cos(x3)**2)
    dxdt = [x2, x2p, x4, x4p]
    return dxdt

class CustomModel(tf.keras.Model):
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
                V1 = y_out[:,1] # The second output is the square root of the target for the 
                # lyapunov function.
                V2 = y_out[:,2] # The third output is the negative square root of the target for the
                # lyapunov function derivative wrt time.
                V = y_out[:,3] # The fourth output is the lyapunov function.
            ## Computing loss function
            # NOTE:
            # This overwrites the loss function configured in compile()!!!!!!!!
            # The first part of the loss function is the old loss function which
            # Tries to make the predicted control close to the optimal control with the 
            # convex combination of the mse and msle errors.
            # The second part of the loss function is the part that tries to make the forth
            # output of the dnn a lyapunov function for the system.
            eps = 1 # Value for the target V and Vdot
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
            loss = custom_loss(y,y_pred) + mse(V1**2+eps,V) + mse(-V2**2-eps,Vdot)

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

def custom_loss(y_true, y_pred):
    alpha = 0.1
    msle = tf.keras.losses.MeanSquaredLogarithmicError()
    mse = tf.keras.losses.MeanSquaredError()
    return alpha*msle(y_true, y_pred)+(1-alpha)*mse(y_true, y_pred)

## Loading trained model
model = tf.keras.models.load_model('dnn_trained.h5',
    custom_objects={ 'custom_loss': custom_loss,
        'CustomModel': CustomModel})

def lyap(t,x):
    u = model.predict(np.array([x]))[0,3]
    return u

# Create cloud of data around origin 
size = 10
x1 = np.random.normal(loc=0,scale=1,size=size)
x2 = np.random.normal(loc=0,scale=1,size=size)
x3 = np.random.normal(loc=0,scale=1,size=size)
x4 = np.random.normal(loc=0,scale=1,size=size)
# An array that holds xs, V and Vdot
# xx1, xx2, xx3, xx4 = np.meshgrid(x1,x2,x3,x4)

V = np.zeros((size**4,3))
line = 0
for i in range(size):
    for j in range(size):
        for k in range(size):
            for l in range(size):
                V[line] = np.array([x1[i],x2[j],lyap(0,np.array([x1[i],x2[j],x3[k],x4[l]]))])
                line += 1

V = np.array(V)
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(projection='3d')
img = ax.scatter(V[:,0], V[:,1], V[:,2])
fig.colorbar(img)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

print(min(V[:,2]))
plt.show()