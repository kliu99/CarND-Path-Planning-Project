# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

## Prediction

The model used to predict other vehicles on the highway is simply linear model assuming constant acceleration.

## Behavior Planning

Finite State Machines (FSM) with three states are considered.

- Lane Change Left (LCL)
- Keep Lane (KL)
- Lane Change Right (LCR)

The ego car prefers to maintain its current lane (KL) when it is not too close to the front vehicle (`s > 30m`). 

If the ego car is too close to the front vehicle, it will first consider if it is safe to change lane to the left. If LCL is not safe, it will check if LCR is safe. If both LCL and LCR are not safe, it will keep the current lane (KL) with reduced speed.

To avoid accelerate or decelerate to quick. The speed is gradually reduced by `0.224m/s` every cycle.

## Trajectory Generation

To allow smooth path generation, the `spline` library is used. The spline path is generated from three 30m evenly spaced points ahead of the vehicle starting point. 

For continuity, all the previously un-runned `previous_path` is reused. At each cycle, the trajectory generation will only generate `50 - previous_path.size()` number of points.

The newly generated points are queried from the spline with vehicle velocity and target lane as the parameters to allow smooth lane change without violation any acceleration and jerk constraints while maintaining desired velocity.

