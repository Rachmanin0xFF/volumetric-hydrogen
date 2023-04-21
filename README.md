# volumetric-hydrogen
A realtime volumetric visualization of the hydrogen-like atom's eigenfunctions.

I wrote this because I was tired of seeing atomic orbitals represented with hard isosurfaces. By default, I show wavefunction argument as hue and squared modulus as 'fog' density.
This code is a bit different from the majority of similar visualizations in that it efficiently calculates Laguerre polynomials and spherical harmonics on-the-fly using dynamic programming.

The visualization is a WIP; eventually I'll publish it to my website with some nice controls (maybe a writeup / educational video as well).

In the meantime, though, you can view a working copy [here](https://www.shadertoy.com/view/cdSSDw).

![example screenshot](./n5_l3_m2.png)
