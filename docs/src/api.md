## Summary

There are four main RAFF structures: 
1. *[Main functions](@ref):* directly called by user; 
1. *[Auxiliary functions](@ref):* used like internal auxiliary functions;
1. *[Random generation](@ref):* used to generate random sets of data, in order to test `RAFF`
1. *[Output type](@ref):* type defined to manipulate output information.

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
RAFF.setLMOutputLevel
RAFF.setRAFFOutputLevel
```

## Random generation
```@docs
RAFF.generateTestProblems
RAFF.get_unique_random_points
RAFF.generateNoisyData
```

## Output type
```@docs
RAFFOutput
```