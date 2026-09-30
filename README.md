# Double_Pendulum
This project is a 2D simulation of a `double pendulum`. 

## Build  and execution
In order to run the simulation you'll have to compile the project with `CMake` with the following instructions:

```bash
% cmake -S . -B build -G"Ninja Multi-Config"
% cmake --build build --config Debug
% cmake --build build --config Release
```
Now the project is compiled and we can run it with the following commands:
```bash
% ./build/Debug/progetto
```
or
```bash
% ./build/Release/progetto
```
Whend the execution will start, it'll appear a window with the double pendulum in its rest position. Is now possible to grab and drag with the mouse one of the two masses around the space. When the mouse will be released the evolution of the double pendulum will begin.

However, if you want to set some starting angle for the masses you can do it. The only thing you have to do is to write the values (in degrees) after the command of execution of the project.
Here is an example:
```bash
% ./build/Release/progetto 45 60
```
The first pendulum will now start with an angle `theta1 = 45°` and the second pendulum will start with an angle `theta2 = 60°`.

## How the evolution works
