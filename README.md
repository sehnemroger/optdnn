# Optimal closed loop control through DNN

Optimal control problems (OCP) can only, in general, be solved in open loop, i.e. with the control and states
as a function of time.

This project aims to learn a function that maps from the states to the control "optimally" 
trought Deep Neural Network (DNN).

The OCP is solved using PSOPT and the DDN using TensorFlow.

# Dockerfile
There is a dockerfile in the docker folder that creates an enviroment to run PSOPT. 

## Build
To build the container on your local machine named as optdnn_env, just use: 

docker build . -t optdnn_env 

Note: don't forget the dot.

## Run 
To run something in the environment, for example, the example launch of PSOPT, you can use:

docker run -it --rm  optdnn_env:latest ./psopt/build/examples/launch/launch 
