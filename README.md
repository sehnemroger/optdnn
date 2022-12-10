# Optimal closed loop control through DNN

Optimal control problems (OCP) can only, in general, be solved in open loop, i.e. with the control and states
as a function of time.

This project aims to learn a function that maps from the states to the control "optimally" 
trought Deep Neural Network (DNN).

The OCP is solved using PSOPT and the DDN using TensorFlow.

# Docker
There is a dockerfile in the docker folder that creates an enviroment to run PSOPT. 

## Build
To build the container on your local machine named as optdnn_env, cd to /docker and use: 
```
docker build . -t optdnn_env 
```
Note: don't forget the dot.

## Run 
To run something in the environment, for example, the example launch of PSOPT, you can use:
```
docker run -it --rm  optdnn_env:latest ./psopt/build/examples/launch/launch 
```
Note: To compile your own example you should implement it in the file `/PSOPT/examples/user/user.cxx`. Then, cd to `/PSOPT/build/examples/user` and run 
```
make
```
This will create your compiled file in `/PSOPT/build/examples/user`.

### Run inside the docker container
To run inside the docker container we will have to bind mount two folders with the container, namely, `/build` to `~/build/examples/user` and `/src` to `~/examples/user`, this means that you will have to place your `.cxx` file inside `/src` and then, to run, just run in your terminal the following command:
```
docker run -it --rm -v $PWD/src:/usr/src/optdnn/psopt/examples/user/ -v $PWD/build:/usr/src/optdnn/psopt/build/examples/user optdnn_env:latest bash

```
In this command `-u$(id -u):$(id -g)` maps your user to the docker container user, `-it` runs interactively, `--rm` deletes the container after you leave, `-v $PWD/src:/usr/src/optdnn/psopt/examples/user/` bind mount `/src` to `/examples/user`, `-v $PWD/build:/usr/src/optdnn/psopt/build/examples/user` bind mount `/build` to `/build/examples/user`, `optdnn_env:latest` is the image and `bash` runs the bash.

With this your implemented `.cxx` file is already inside `/psopt/examples/user/` and you can just cd to `/psopt/build/examples/user` and run `make` and then the `user` executable.

### Run TensorFlow in Docker
To run tensorflow inside the docker container use the following command
```
docker run -it --rm -v $PWD:/tmp -w /tmp optdnn_env:latest python ./script.py
```
