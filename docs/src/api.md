## Summary

There are three main RAFF structures: 
1. *main functions:* called by user; 
1. *auxiliary functions:* used like internal auxiliary function but can be modify user;
1. *output type:* type defined to manipulate output information.

## Main functions
```@docs
lmlovo
raff
praff
setRAFFOutputLevel
setLMOutputLevel
```
## Auxiliary functions
```@docs
RAFF.eliminate_local_min!
RAFF.SortFun!
RAFF.update_best
RAFF.consume_tqueue
RAFF.check_and_close
RAFF.generateTestProblems
RAFF.get_unique_random_points
RAFF.generateNoisyData
```

## Output type
```@docs
RAFFOutput
```