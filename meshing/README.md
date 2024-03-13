# Requirements

Requires FreeSurfer.

> [!NOTE]  
> Script `extract-ventricles.sh` copied from mri2fem book (`chp4/extract-ventricles.sh`), updated location of input filename.

Requires `meshio`, `fenics` and `SVMTK`.
Can be installed using docker or conda.
Freesurfer is part of the [Dockerfile](./docker/Dockerfile), but the user has to supply a license file in the [docker-folder](./docker/) as `.license`.
A license can be obtained through [FreeSurfer License](https://surfer.nmr.mgh.harvard.edu/fswiki/License).

The docker image can be build from this folder with

```bash
docker build -t waterscales-brains -f ./docker/Dockerfile .
```

The docker image should be ran from the [parent-folder](./..) with

```bash
docker run -ti -v $(pwd):/root/shared -w /root/shared waterscales-brains
```

## Instructions

### Generate STL (intermediate surface triangulation) of the brain from MRI

To generate full brain mesh (without ventricles), with gray and white markers, from scratch, run:

```bash
python3 create_mesh.py --generate_stl --preprocess_surfaces --create 16
```

### Create volume mesh

To just create the mesh (after generating surfaces etc, run)

```bash
python3 create_mesh.py --create 8
```

> [!NOTE]  
> NB: The above command depends on FreeSurfer, meshio, SVMTK, but not on FEniCS, and will take perhaps 1-2 minutes to complete.

Adjust other parameters in [create_mesh](./create_mesh.py) functions to adjust smoothing, surface meshing parameters etc.

### Post-process mesh in FEniCS

Run the FEniCS processing of the mesh, here `abby_16`:

```bash
python3 create_mesh.py --postprocess 16
```

Check mesh with Dokken's topology check [script](https://gist.github.com/jorgensd/4ebce24192b75a060f342ccf6f4c6555).

```bash
python3 mesh_checker.py --infile abby_16 --topology="mesh/topology" --geometry="mesh/coordinates"
```

To clean everything and start a new, do

```bash
python3 create_mesh.py --clean
```

# Notes

```bash
abby_4: 37179 cells, 7848 vertices. # Too coarse for anything but prototyping
abby_6: 90367 cells, 17750 vertices.
abby_7: 117045 cells, 23024 vertices.
abby_8: 145195 cells, 28652 vertices.
```
