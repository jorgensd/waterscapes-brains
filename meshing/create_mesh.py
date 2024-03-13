def read_h5_mesh(hdf5file):
    # Read the mesh and mesh data from .h5 :

    import dolfin
    import numpy
    
    print("Importing mesh from %s" % hdf5file)
    mesh = dolfin.Mesh()
    hdf = dolfin.HDF5File(mesh.mpi_comm(), hdf5file, "r")
    hdf.read(mesh, "/mesh", False)

    d = mesh.topology().dim()
    subdomains = dolfin.MeshFunction("size_t", mesh ,d)
    hdf.read(subdomains, "/subdomains")

    boundaries = dolfin.MeshFunction("size_t", mesh, d-1)
    hdf.read(boundaries, "/boundaries")
    hdf.close()

    volume = dolfin.assemble(1*dolfin.dx(domain=mesh))
    
    print(".. with %d cells and %d vertices." % (mesh.num_cells(), mesh.num_vertices()))
    print("Subdomain markers:", numpy.unique(subdomains.array()))
    print("Boundary markers:", numpy.unique(boundaries.array()))
    print("Domain volume (in dm^3 = L):", volume/(100**3))

    return (mesh, subdomains, boundaries)
    
def preprocess_surface(f, lhonly=False):
    import os.path
    import SVMTK as svmtk
    
    print("Loading %s" % f)
    surface = svmtk.Surface(f)

    print("Remeshing %s, this may take some time..." % f)
    L = 3.0 # (mm) This should be comparable with the mesh size of the
            # mesh to be created AND the resolution of the geometry
    m = 5   # This is an iteration number parameter, 5 is fine.
    do_not_move_boundary_edges = False
    surface.isotropic_remeshing(L, m, do_not_move_boundary_edges)
    
    remeshed = "".join([f[:-3], "remesh.stl"])
    print("Saving remeshed surface as %s" % remeshed)
    surface.save(remeshed)
    
    print("Filling holes in the surface %s" % f)
    surface.fill_holes()

    print("Keeping largest connected component")
    surface.keep_largest_connected_component()
    
    if "pial" in f and lhonly:
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
        #surface.separate_narrow_gaps(-0.33)

    print("Smoothing %s, this may take some time..." % f)
    n_its = 5
    surface.smooth_taubin(n_its)

    # Save final surface 
    final = "".join([f[:-3], "final.stl"])
    print("Saving final surface as %s" % final)
    surface.save(final)

def convert_mesh_data(infile, outfile, outdir="tmp-abby-meshes"):
    # Convert the .mesh file to FEniCS (legacy) .h5 format with
    # temporary .xdmf output in outdir

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
    # directory 'outdir'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    print("Writing mesh, subdomains and boundaries to %s/..." % outdir)
    stem = os.path.splitext(infile)[0]
    meshio_mesh = meshio.Mesh(points, cells)
    mesh_filename = "%s/%s_mesh.xdmf" % (outdir, stem)
    meshio.write(mesh_filename, meshio_mesh)

    meshio_domains = meshio.Mesh(points, cells, cell_data=subdomains)
    subdomains_filename = "%s/%s_subdomains.xdmf" % (outdir, stem)
    meshio.write(subdomains_filename, meshio_domains)

    meshio_boundaries = meshio.Mesh(points, facets, cell_data=boundaries)
    boundaries_filename = "%s/%s_boundaries.xdmf" % (outdir, stem)
    meshio.write(boundaries_filename, meshio_boundaries)
    
    import dolfin
    print("Reading mesh, subdomains and boundaries back in from %s/..." % outdir)
    mesh = dolfin.Mesh()
    with dolfin.XDMFFile(mesh_filename) as infile:
        infile.read(mesh)
    print("mesh.hmin = ", mesh.hmin())
    print("mesh.hmax = ", mesh.hmax())

    d = mesh.topology().dim()
    subdomains = dolfin.MeshFunction("size_t", mesh, d)
    with dolfin.XDMFFile(subdomains_filename) as infile:
        infile.read(subdomains)

    boundaries = dolfin.MeshFunction("size_t", mesh, d-1)
    with dolfin.XDMFFile(boundaries_filename) as infile:
        infile.read(boundaries)

    # Write mesh and its data to a single .h5 file:
    ofile = dolfin.HDF5File(mesh.mpi_comm(), outfile, "w")
    ofile.write(mesh, "/mesh")
    ofile.write(subdomains, "/subdomains")
    ofile.write(boundaries, "/boundaries")
    ofile.close()

def create_lh_mesh(pial, ventricles, outbase="tmp", n=12):

    import SVMTK as svmtk

    # Create hemisphere mesh given pial and ventricular surfaces
    print("Creating hemisphere mesh")
    pia = svmtk.Surface(pial)
    ventricle = svmtk.Surface(ventricles)

    tags = {"pial": 1, "ventricle": 3}
    smap = svmtk.SubdomainMap()
    smap.add("10", tags["pial"])
    smap.add("11", tags["ventricle"])
    
    domain = svmtk.Domain([pia, ventricle], smap)
    domain.create_mesh(n)
    
    print("Removing ventricles")
    domain.remove_subdomain(tags["ventricle"])
    
    meshfile = "%s.mesh" % outbase
    print("Saving mesh as %s" % meshfile)
    domain.save(meshfile)
    print("Hemisphere mesh generation complete.")
    
def create_brain_mesh(files, outbase, n=12):

    import SVMTK as svmtk

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
    svmtk.separate_overlapping_surfaces(rhpial, lhpial, white, edge_movement=-1.0,
                                        smoothing=0.3)
    svmtk.separate_close_surfaces(rhpial, lhpial, white, edge_movement=-1.0,
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
        
    meshfile = "%s.mesh" % outbase
    print("Saving mesh as %s" % meshfile)
    domain.save(meshfile)
    print("Brain mesh generation complete.")

def clean_all_outputs():
    # Danger zone: clean-up everything that can possibly have been
    # generated in this directory tree.
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
    parser.add_argument("--clean", action="store_true", help="Remove *all* generated files")
    parser.add_argument("--generate_stl", action="store_true", help="Copy and convert FreeSurfer files to STL")
    parser.add_argument("--preprocess_surfaces", action="store_true", help="Remesh and smoothen STL surfaces")
    parser.add_argument("--create", default=0, help="Create mesh with given int resolution")
    parser.add_argument("--create_lh", default=0, help="Create lh mesh with given int resolution")
    parser.add_argument("--postprocess", default=0, help="Convert meshio data to FEniCS (HDF5) file")
    parser.add_argument("--postprocess_lh", default=0, help="Convert lh meshio data to FEniCS (HDF5) file")
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

    # Preprocess surfaces. Note the L parameter in preprocess surface,
    # this adjusts the surface cell size. Adjust this size depending
    # on the target mesh size.
    if args.preprocess_surfaces:
        print("Remeshing, smoothing and fixing all surfaces")
        stls = [os.path.join(stldir, "%s.stl" % f) for f in files]
        for f in stls:
            preprocess_surface(f)

    # Generate brain mesh
    if args.create:
        n = int(args.create)
        sfiles = [os.path.join(stldir, "".join([f, ".final.stl"])) for f in files[:-1]]
        sfiles.append(os.path.join(stldir, "".join([files[-1], ".stl"])))
        print("Creating brain mesh from %s" % sfiles)
        create_brain_mesh(sfiles, "abby_%d" % n, n=n)

    # Generate left hemisphere mesh
    if args.create_lh:
        n = int(args.create_lh)
        sfiles = [os.path.join(stldir, "".join([f, ".final.stl"])) for f in files[:-1]]
        sfiles.append(os.path.join(stldir, "".join([files[-1], ".stl"])))
        print("Creating left hemisphere mesh from %s" % sfiles)
        create_lh_mesh(sfiles[0], sfiles[4], "abby_lh_%d" % n, n=n)
        
    # Relies on FEniCS:
    if args.postprocess:
        n = int(args.postprocess)
        # Convert brain mesh, subdomains, and boundary markers to FEniCS .h5 format
        convert_mesh_data("abby_%d.mesh" % n, "abby_%d.h5" % n)
        # Test that the data can be read back in.
        read_h5_mesh("abby_%d.h5" % n)

    if args.postprocess_lh:
        n = int(args.postprocess_lh)
        # Convert lh mesh, subdomains, and boundary markers to FEniCS .h5 format
        convert_mesh_data("abby_lh_%d.mesh" % n, "abby_lh_%d.h5" % n)
        # Test that the data can be read back in.
        read_h5_mesh("abby_lh_%d.h5" % n)
        
        
    # Step-by-step example: to create meshes 

    # Meshsize 32 below gives nice meshes, but you can make them
    # coarser as well, use e.g. 8, 16, or 24.
    
    # python3 create_mesh.py --generate_stl   
    # python3 create_mesh.py --preprocess_surface 
    # python3 create_mesh.py --create 32
    # python3 create_mesh.py --create_lh 32 

    # fenicsproject run
    # python3 create_mesh.py --postprocess 32
    # python3 create_mesh.py --postprocess_lh 32
    # exit

    # python3 mesh_checker.py --infile abby_32
    # python3 mesh_checker.py --infile abby_lh_32

    ## *Other notes to self and others*

    # Distance between gaps need to be larger than minimal cell size.

    # First, have a look at what L is set to in the code above. Adjust
    # at will (recommended range from 1.0 to 3.0 (mm): 1.0 is pretty
    # fine, 3.0 is pretty coarse). This has implications for the
    # largest meshsize you can consider.

    # Creating good left hemisphere meshes is pretty straight-forward
    # (increase the mesh size if you get issues with the mesh
    # checker), but creating good full brain meshes is not so easy due
    # to the merging of the multiple hemisphere surfaces. My best tip
    # thus far is to play with the processing, resolutions and mesh
    # sizes or consult with a brain meshing expert.
