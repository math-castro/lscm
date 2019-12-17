# project-inf574
Least Squares Conformal Maps for Automatic Texture Atlas Generation Implementation

## Requirements
Install libigl on the project root folder or one directory above

Install Eigen v3.3
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

