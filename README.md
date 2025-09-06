## Params for Experiments
TODO: put all the params for experiments here

## PGO
Certifable 
- all datasets were set at 200 max iterations rest default 
- Kitti_00 and 05 needed 1e-15 for abs and rel error tolerance

Certifiable Odom
- running with 200 max iters
- With odom most datasets could be run at rel error at 1e-5 and abs at 1e-3
- Rim is only pgo dataset not working with odom ( not verifying )

Local
- ran with large max iterations (500)
- set with default error tolerances

## Landmark 
GT LMK &  Odom
- default for both datasets
Local
- default for both datasets 
Certifable 
- default for cityTrees
- Victoria park abs and rel tolerance was set to e-12

## Range 
Certifiable
- Goats 15 starting rank to 3
- Had toruble verifying on plaza1 with relative eta so set abs eta to 3e-3, had to have inital p = 3
- plaza 2 had absolute eta set to 2e-3, starting rank was 2

odom 
- Goats 15 starting rank to 3
- Had trouble verifying on plaza1 with relative eta so set abs eta to 3e-3, had to have inital p = 3
- plaza 2 had absolute eta set to 2e-3, starting rank was 2

local
- set max iters to 1000 (plaza 1 and 2 still has a sizeable cost change even at that high iteration count)



