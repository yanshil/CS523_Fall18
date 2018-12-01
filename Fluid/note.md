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