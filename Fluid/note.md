1. When using Template Class, don't forget to precompile with 

```
template class Nova::FluidQuantity<float, 2>;
template class Nova::FluidQuantity<float, 3>;
#ifdef COMPILE_WITH_DOUBLE_SUPPORT
template class Nova::FluidQuantity<double, 2>;
template class Nova::FluidQuantity<double, 3>;
#endif
```

2. 

```
void advect(T timestep, FluidQuantity *_v[d]);
```

usage
```
advect(T t, FluidQuantity _v);
```

3. OMG my offset2 index has bug!!!

The `y <- (os - x) / m` I made m as n.....

4. CalculateRHS order??? Nopes.... it was correct

5. 

```
        for (int x = 1; x <= m; x++)
        {
            _v[1]->at(T_INDEX{x, 1}) = 0.0;
            _v[1]->at(T_INDEX{x, n + 1}) = 0.0;
        }

        for (int y = 1; y <= n; y++)
        {
            _v[0]->at(T_INDEX{m + 1, y}) = 0.0;
            _v[0]->at(T_INDEX{1, y}) = 0.0;
        }
```






        for(int idx = 0; i < size; idx++)
        {
            T_INDEX index = offset2index(index);
            
            for(int axis = 0; axis < d; axis++)
            {
                /* code */
            }
            
        }