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
set_raff_output_level
set_lm_output_level
```

## Auxiliary functions
```@docs
RAFF.eliminate_local_min!
RAFF.sort_fun!
RAFF.update_best
RAFF.consume_tqueue
RAFF.check_and_close
RAFF.check_ftrusted
RAFF.interval_rand!
```

## Random generation
```@docs
RAFF.generate_test_problems
RAFF.get_unique_random_points
RAFF.get_unique_random_points!
RAFF.generate_noisy_data!
RAFF.generate_noisy_data
RAFF.generate_clustered_noisy_data!
RAFF.generate_clustered_noisy_data
RAFF.model_list
```

## Output type
```@docs
RAFFOutput
```
