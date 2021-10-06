
def read_h5_mesh(hdf5file):
    # Read the mesh and mesh data from .h5 :

    import dolfin
    import numpy
    
    print("Importing mesh")
    mesh = dolfin.Mesh()
    hdf = dolfin.HDF5File(mesh.mpi_comm(), hdf5file, "r")
    hdf.read(mesh, "/mesh", False)

    d = mesh.topology().dim()
    subdomains = dolfin.MeshFunction("size_t", mesh ,d)
    hdf.read(subdomains, "/subdomains")

    boundaries = dolfin.MeshFunction("size_t", mesh, d-1)
    hdf.read(boundaries, "/boundaries")

    print(".. with %d cells and %d vertices." % (mesh.num_cells(), mesh.num_vertices()))
    print("Subdomain markers:", numpy.unique(subdomains.array()))
    print("Boundary markers", numpy.unique(boundaries.array()))
    
    hdf.close()

def preprocess_surface(f):

    import os.path
    import SVMTK as svmtk
    
    print("Loading %s" % f)
    surface = svmtk.Surface(f)
            
    print("Remeshing %s, this may take some time..." % f)
    L = 2.0 # mm
    m = 5   # Quantitative mesh size parameter 
    do_not_move_boundary_edges = False
    surface.isotropic_remeshing(L, m, do_not_move_boundary_edges)
    
    remeshed = "".join([f[:-3], "remesh.stl"])
    print("Saving remeshed surface as %s" % remeshed)
    surface.save(remeshed)
    
    print("Smoothing %s, this may take some time..." % f)
    n_its = 2
    surface.smooth_taubin(n_its)
    
    smooth = "".join([f[:-3], "smooth.stl"])
    print("Saving smooth surface as %s" % smooth)
    surface.save(smooth)

def convert_mesh_data(infile, outfile="abby.h5", tmpdir="tmp-xdmf"):
    import meshio
    import os.path

    # Read the mesh from the infile
    mesh = meshio.read(infile)
    points = mesh.points

    cells = {"tetra": mesh.cells_dict["tetra"]}
    facets = {"triangle": mesh.cells_dict["triangle"]}
    subdomains = {"subdomains": [mesh.cell_data_dict["medit:ref"]["tetra"]]}
    boundaries = {"boundaries": [mesh.cell_data_dict["medit:ref"]["triangle"]]}
    
    # Write the mesh, subdomains and boundaries to .xdmf format in the
    # directory 'tmpdir'
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)

    print("Writing mesh, subdomains and boundaries to %s/..." % tmpdir)
    meshio_mesh = meshio.Mesh(points, cells)
    meshio.write("%s/mesh.xdmf" % tmpdir, meshio_mesh)

    meshio_domains = meshio.Mesh(points, cells, cell_data=subdomains)
    meshio.write("%s/subdomains.xdmf" % tmpdir, meshio_domains)

    meshio_boundaries = meshio.Mesh(points, facets, cell_data=boundaries)
    meshio.write("%s/boundaries.xdmf" % tmpdir, meshio_boundaries)
    
    import dolfin
    print("Reading mesh, subdomains and boundaries back in from %s/..." % tmpdir)
    mesh = dolfin.Mesh()
    with dolfin.XDMFFile("%s/mesh.xdmf" % tmpdir) as infile:
        infile.read(mesh)

    d = mesh.topology().dim()
    subdomains = dolfin.MeshFunction("size_t", mesh, d)
    with dolfin.XDMFFile("%s/subdomains.xdmf" % tmpdir) as infile:
        infile.read(subdomains)

    boundaries = dolfin.MeshFunction("size_t", mesh, d-1)
    with dolfin.XDMFFile("%s/boundaries.xdmf" % tmpdir) as infile:
        infile.read(boundaries)

    # Write mesh and its data to a single .h5 file:
    ofile = dolfin.HDF5File(mesh.mpi_comm(), outfile, "w")
    ofile.write(mesh, "/mesh")
    ofile.write(subdomains, "/subdomains")
    ofile.write(boundaries, "/boundaries")
    ofile.close()
    
def create_brain_mesh(files, outbase, n=12):

    import SVMTK as svmtk

    # Load as surfaces
    surfaces = [svmtk.Surface(s) for s in files]
    surfaces[2].union(surfaces[3])    
    surfaces.pop(3)

    # Label the different regions
    tags = {"pial": 1, "white": 2, "ventricle": 3}
    smap = svmtk.SubdomainMap()
    smap.add("1000", tags["pial"])
    smap.add("0100", tags["pial"])
    smap.add("1010", tags["white"])
    smap.add("0110", tags["white"])
    smap.add("1110", tags["white"])
    smap.add("1011", tags["ventricle"])
    smap.add("0111", tags["ventricle"])
    smap.add("1111", tags["ventricle"])
    
    print("Creating mesh")
    domain = svmtk.Domain(surfaces, smap)
    domain.create_mesh(n)

    print("Removing ventricles")
    domain.remove_subdomain(tags["ventricle"])

    meshfile = "%s.mesh" % outbase
    print("Writing mesh to %s" % meshfile)
    domain.save(meshfile)

    if False:
        tmpdir = "tmp-xml"
        if not os.path.isdir(tmpdir):
            os.mkdir(tmpdir)
        xdmffile = os.path.join(tmpdir, "%s.xdmf" % outbase)
        xmlfile = os.path.join(tmpdir, "%s.xml" % outbase)
        print("Converting mesh also to %s via .xml" % xdmffile)
        subprocess.run(["meshio-convert", meshfile, xmlfile])
        subprocess.run(["meshio-convert", xmlfile, xdmffile])

def clean_all_outputs():
    import glob
    import shutil
    import os.path
    
    globs = []
    for ext in ["tmp*", "*.h5", "*.mesh", "*.xdmf", "*.xml", "*.stl", "stl-*"]:
        globs.extend(glob.glob(ext))

    for f in globs:
        print("Deleting %s" % f)
        if os.path.isfile(f):
            os.remove(f)
        else:
            shutil.rmtree(f)
        
if __name__ == "__main__":
    
    import argparse
    import subprocess
    import os.path
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean", action="store_true", help="Remove all generated files")
    parser.add_argument("--generate_stl", action="store_true", help="Copy and convert FreeSurfer files")
    parser.add_argument("--preprocess_surfaces", action="store_true", help="Remesh and smoothen STL surfaces")
    parser.add_argument("--create", default=0, help="Create mesh with given int resolution")
    parser.add_argument("--postprocess", default=0, help="Convert meshio data to FEniCS (HDF5) file")
    args = parser.parse_args()
    
    if args.clean:
        clean_all_outputs()
        
    # Where shall we put the stl files?
    stldir = "stl-surfaces"
    files = ["lh.pial", "rh.pial", "lh.white", "rh.white", "lh.ventricles"]

    if args.generate_stl:

        if not os.path.isdir(stldir):
            os.mkdir(stldir)

        # First run script to extract ventricles and create ./lh.ventricles.stl
        subprocess.run(["./extract-ventricles.sh"])
        outfile = "lh.ventricles.stl"
        print("Moving %s to %s/..." % (outfile, stldir))
        subprocess.run(["mv", outfile, os.path.join(stldir, outfile)])

        # Where are the FreeSurfer surface files located?
        base = "../freesurfer/abby/surf/"

        # Convert all other FreeSurfer surfaces to STL format and move to subdirectory
        for f in files[:-1]:
            outfile = f +".stl"
            subprocess.run(["mris_convert", base + f, "./" + outfile])
            print("Moving %s to %s/..." % (outfile, stldir))
            subprocess.run(["mv", outfile, os.path.join(stldir, outfile)])
    
    if args.preprocess_surfaces:
        print("Remeshing and smoothing all surfaces")
        stls = [os.path.join(stldir, "%s.stl" % f) for f in files]
        for f in stls:
            preprocess_surface(f)

    # Generate mesh
    if args.create:
        n = int(args.create)
        sfiles = [os.path.join(stldir, "".join([f, ".smooth.stl"])) for f in files[:-1]]
        sfiles.append(os.path.join(stldir, "".join([files[-1], ".stl"])))
        create_brain_mesh(sfiles, "abby_%d" % n, n=n)

    if args.postprocess:
        n = int(args.postprocess)
        convert_mesh_data("abby_%d.mesh" % n, "abby_%d.h5" % n)
        read_h5_mesh("abby_%d.h5" % n)
