# Copyright (C) 2021 Chris Richardson, Jørgen S. Dokken
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numba
import networkx
from collections import defaultdict
import h5py
import numpy as np
import datetime
import argparse
import meshio


def check_mesh(file_name: str, topology_path: str, geometry_path: str, outfile: str):
    """
    Check if cells in mesh in connected by edges or vertex only (not by facet)
    """
    data_file = h5py.File(f"{file_name}.h5", "r")

    topology = data_file[topology_path]

    # Create a list of all triangular facets, from tetrahedra

    ntet = topology.shape[0]

    print(f"{datetime.datetime.now()}: Number of tetrahedra: {ntet}")

    facet_array = np.empty((ntet * 4, 3), dtype=int)

    # Create triangles from tetrahedra
    print(f"{datetime.datetime.now()}: Facets 1")
    facet_array[:ntet, :] = topology[:, :3]
    print(f"{datetime.datetime.now()}: Facets 2")
    facet_array[ntet:2 * ntet, :] = topology[:, [0, 2, 3]]
    print(f"{datetime.datetime.now()}: Facets 3")
    facet_array[2 * ntet:3 * ntet, :] = topology[:, [0, 1, 3]]
    print(f"{datetime.datetime.now()}: Facets 4")
    facet_array[3 * ntet:, :] = topology[:, 1:]

    print(f"{datetime.datetime.now()}: Sort each row")
    facet_array.sort(axis=1)

    print(f"{datetime.datetime.now()}: Facet array: {facet_array.shape}")

    # Presort by first column, to speed up np.unique
    facet_array = facet_array[np.argsort(facet_array[:, 0])]
    print(f"{datetime.datetime.now()}: call np.unique (facets)")
    ftri, counts = np.unique(facet_array, axis=0, return_counts=True)
    bad_facets = np.where(counts > 2)[0]
    if len(bad_facets) > 0:
        print("The following facets exist in more than two cells")
        for i in bad_facets:
            print(f"Facet: {ftri[i]} in  {counts[i]} cells")
    del facet_array
    print(f"{datetime.datetime.now()}: Facets OK, taking exterior facets")

    # Take external facets only, and compute edges
    external = np.where(counts == 1)[0]
    ntri = len(external)
    print(f"{datetime.datetime.now()}: Number of exterior facets: {ntri}")
    ftri = ftri[external]

    print(f"{datetime.datetime.now()}: Edges 1")
    edge_array = np.empty((ntri * 3, 2), dtype=int)
    edge_array[:ntri, :] = ftri[:, :2]
    print(f"{datetime.datetime.now()}: Edges 2")
    edge_array[ntri:2 * ntri, :] = ftri[:, [0, 2]]
    print(f"{datetime.datetime.now()}: Edges 3")
    edge_array[2 * ntri:, :] = ftri[:, 1:]

    edge_array.sort(axis=1, kind='stable')
    print(datetime.datetime.now(), end=': ')
    print(f'Exterior edges array = {edge_array.shape}')

    # Presort by first column, to speed up np.unique
    edge_array = edge_array[np.argsort(edge_array[:, 0])]
    print(f"{datetime.datetime.now()}: call np.unique (edges)")
    edges, count_edge = np.unique(edge_array, axis=0, return_counts=True)
    bad_edges = np.where(count_edge != 2)[0]
    if len(bad_edges) > 0:
        print(f"The following {len(bad_edges)} edges exist in more than two surface facets:")
        for i in bad_edges:
            print(f"Edge: {edges[i]} in {count_edge[i]} surface facets")

    print(f"{datetime.datetime.now()}: Analysing surface nodes")

    # List of lists, collect edges on each vertex - should form a closed set
    d = defaultdict(list)
    for f in ftri:
        d[f[0]] += [[f[1], f[2]]]
        d[f[1]] += [[f[0], f[2]]]
        d[f[2]] += [[f[0], f[1]]]

    print(f'{datetime.datetime.now()}: Num surface nodes: {len(d)}')

    print(f"{datetime.datetime.now()}: Compute graph with networkx")
    for k in d:
        G = networkx.Graph(d[k])
        if networkx.number_connected_components(G) != 1:
            print(f"Bad vertex at {k}, {d[k]}")

    # Print bad vertex coordinates
    print(f"{datetime.datetime.now()}: Load mesh geometry")
    geometry = data_file[geometry_path]
    np.set_printoptions(precision=5)
    for edge in edges[bad_edges]:
        print(f"Bad edge with vertex coordinates:\n {geometry[edge]}")

    # Compile numba on small set
    extract_cell_index(
        np.array([[0, 1, 2, 7], [1, 2, 3, 4], [4, 1, 0, 2]], dtype=np.int32),
        np.array([[0, 1]], dtype=np.int32))
    print(f"{datetime.datetime.now()}: Extract cells of bad edges")
    bad_cells = extract_cell_index(np.asarray(topology, dtype=np.int32), edges[bad_edges])

    print(f"{datetime.datetime.now()}: Write new mesh to xdmf")
    c_vals = np.zeros(topology.shape[0])
    c_vals[bad_cells] = edges[bad_edges].shape[0] + 2
    mesh_out = meshio.Mesh(geometry, [("tetra", topology), ("line", edges[bad_edges])],
                           cell_data={"val": [c_vals,
                                              np.arange(1, edges[bad_edges].shape[0] + 1)]})
    meshio.write(f"{outfile}.xdmf", mesh_out)
    return edges[bad_facets], edges[bad_edges]


#@numba.njit
def extract_cell_index(topology, edges):
    """
    Given a mesh topology, and a set of edges defined by its vertices, find the cells
    """
    bad_cells = np.zeros(topology.shape[0], dtype=np.int32)
    i = 0
    for index, cell in enumerate(topology):
        for edge in edges:
            num_edges = 0
            for v0 in cell:
                for v1 in edge:
                    if v0 == v1:
                        num_edges += 1
            if num_edges == 2:
                bad_cells[i] = index
                i += 1
    return bad_cells[:i]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Mesh verification, testing if there are cells that are only connected through edges or vertices",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--infile", default="Benchmark1_order1", type=str, dest="fname",
                        help="Name of input file (without h5 extension)")
    parser.add_argument("--outfile", default="mesh_out", type=str, dest="outname",
                        help="Name of output file (without xdmf extension)")
    parser.add_argument("--topology", default="mesh/topology", type=str, dest="xpath_top",
                        help="Path to mesh topology in input file")
    parser.add_argument("--geometry", default="mesh/coordinates", type=str, dest="xpath_geom",
                        help="Path to mesh geometry in input file")

    args = parser.parse_args()
    check_mesh(args.fname, topology_path=args.xpath_top, geometry_path=args.xpath_geom, outfile=args.outname)
