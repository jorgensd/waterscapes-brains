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
    hdf.close()

    print(".. with %d cells and %d vertices." % (mesh.num_cells(), mesh.num_vertices()))
    print("Subdomain markers:", numpy.unique(subdomains.array()))
    print("Boundary markers", numpy.unique(boundaries.array()))

    return (mesh, subdomains, boundaries)
    
def preprocess_surface(f):

    import os.path
    import SVMTK as svmtk
    
    print("Loading %s" % f)
    surface = svmtk.Surface(f)

    print("Remeshing %s, this may take some time..." % f)
    L = 3.0 # mm This should be comparable with the mesh size of the mesh created
    m = 5   # This is an iteration number parameter
    do_not_move_boundary_edges = False
    surface.isotropic_remeshing(L, m, do_not_move_boundary_edges)
    
    remeshed = "".join([f[:-3], "remesh.stl"])
    print("Saving remeshed surface as %s" % remeshed)
    surface.save(remeshed)
    
    print("Filling holes in the surface %s" % f)
    surface.fill_holes()

    print("Keeping largest connected component")
    surface.keep_largest_connected_component()

    if "pial" in f:
        print("Adjusting boundary for %s" % f)
        surface.adjust_boundary(0.1)

    print("Separating close vertices")
    num_its = 5
    for it in range(num_its):

        print("Number of self-intersections:", surface.num_self_intersections() )
        if surface.num_self_intersections() == 0:
            break
        
        separation = 1.0 # Percentage/ratio 
        print("Separating close vertices with separation = %g" % separation)
        surface.separate_close_vertices(separation)

        print("Remeshing and smoothing again")
        surface.isotropic_remeshing(L, m, do_not_move_boundary_edges)
        surface.smooth_taubin(2)

        print("Number of self-intersections after:", surface.num_self_intersections() )
        # Default argument for separate_narrow_gaps is -0.33. 
        #print("Separating narrow gaps in the surface %s %d" % (f, it))
        #surface.separate_narrow_gaps(-0.5)

    print("Smoothing %s, this may take some time..." % f)
    n_its = 5
    surface.smooth_taubin(n_its)

    # Save final surface 
    final = "".join([f[:-3], "final.stl"])
    print("Saving final surface as %s" % final)
    surface.save(final)

    
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

    print(mesh.hmin())
    print(mesh.hmax())
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

def create_brain_mesh(files, outbase, n=12, wholebrain=False):

    import SVMTK as svmtk

    if wholebrain:
        print("Creating whole brain (pial + white - ventricles) mesh")
        # Load as surfaces
        surfaces = [svmtk.Surface(s) for s in files]
        surfaces[2].union(surfaces[3])    
        surfaces.pop(3)

        (lhpial, rhpial, white, ventricles) = surfaces
        # rhpial and lhpial are the surfaces that we try to separate
        # white argument: exclude stuff inside white
        # edge_movement is a relative factor
        print("Separating overlapping and close surfaces")
        svmtk.separate_overlapping_surfaces(rhpial, lhpial, white, edge_movement=-0.3,
                                            smoothing=0.3)
        svmtk.separate_close_surfaces(rhpial, lhpial, white, edge_movement=-0.3,
                                      smoothing=0.3)

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

        domain = svmtk.Domain(surfaces, smap)
        domain.create_mesh(n)

        print("Removing ventricles")
        domain.remove_subdomain(tags["ventricle"])
        
    else:
        print("Creating hemisphere (pial - ventricles) mesh")
        lhpial = svmtk.Surface(files[0])
        ventricles = svmtk.Surface(files[4])
        tags = {"pial": 1, "ventricle": 3}
        smap = svmtk.SubdomainMap()
        smap.add("10", tags["pial"])
        smap.add("11", tags["ventricle"])

        domain = svmtk.Domain([lhpial, ventricles], smap)
        domain.create_mesh(n)

        print("Removing ventricles")
        domain.remove_subdomain(tags["ventricle"])

    meshfile = "%s.mesh" % outbase
    print("Writing mesh to %s" % meshfile)
    domain.save(meshfile)

    if False:
        tmpdir = "current-xdmf-xml"
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
        print("Remeshing, smoothing and fixing all surfaces")
        stls = [os.path.join(stldir, "%s.stl" % f) for f in files]
        for f in stls:
            preprocess_surface(f)

    # Generate mesh
    if args.create:
        n = int(args.create)
        sfiles = [os.path.join(stldir, "".join([f, ".final.stl"])) for f in files[:-1]]
        sfiles.append(os.path.join(stldir, "".join([files[-1], ".stl"])))
        print("Creating brain mesh from %s" % sfiles)
        create_brain_mesh(sfiles, "abby_%d" % n, n=n)

    # Relies on FEniCS:
    if args.postprocess:
        n = int(args.postprocess)
        convert_mesh_data("abby_%d.mesh" % n, "abby_%d.h5" % n)
        read_h5_mesh("abby_%d.h5" % n)

    # Note to self distance between gaps need to be larger than
    # minimal cell size.
