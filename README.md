# Least Squares Conformal Maps for Automatic Texture Atlas Generation
Implementation of the research paper *[Least squares conformal maps for automatic texture atlas generation](https://doi.org/10.1145/566654.566590)*

## Requirements
Install [libigl](https://libigl.github.io/) on the project folder, one directory above or in the other paths listed in `cmake/FindLIBIGL.cmake`

Install [Eigen](http://eigen.tuxfamily.org/) v3.3
<!-- ### Mac:
```
brew install eigen
```
### Linux:
```
sudo apt-get install libeigen3-dev
``` -->

## Compiling
```
mkdir build
cd build
cmake ..
make
```

## Running

In `build/` folder:

### Default (Bunny Iso-(u,v) curves):

```
./lscm
```

### Custom mesh and mode:

```
./lscm <path to .obj mesh> <mode>
```

#### .obj meshes:
- Fish: `../data/blub_triangulated.obj`
- Bunny: `../data/LSCM_bunny.obj`
- Pumpkin: `../data/pumpkin_tall_10k.obj`
- Cow: `../data/spot_triangulated.obj`
- Teddy: `../data/teddy.obj`

#### modes:
- 0 : show initial features
- 1 : show feature curves
- 2 : show charts
- 3 : show iso-(u,v) curves
- 4 : show packing

#### examples:

```
./lscm ../data/blub_triangulated.obj 0
./lscm ../data/LSCM_bunny.obj 1
./lscm ../data/pumpkin_tall_10k.obj 2
./lscm ../data/spot_triangulated.obj 3
./lscm ../data/teddy.obj 4
```

