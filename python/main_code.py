#z0tool
#Author: Christina Fratzeskou

from shapely.geometry import Polygon, Point, LinearRing
import time
import pywavefront
import numpy as np
import json, sys
from rtree import index

#Functions
def triangle_generator(v):
    n_v = len(v)
    for i in range(0, n_v, 9):
        yield [(v[i], v[i + 1]), (v[i + 3], v[i + 4]), (v[i + 6], v[i + 7])]

# Return bounding box for triangles
def tri2poly(v):
    polys = [Polygon(t) for t in triangle_generator(v)]
    return polys

# Return bounding boxes for all face centers
def pts_bbox(pts):
    pts_bounds = []
    for p in pts:
        # offset = 0.0001
        # pts_bounds.append((p[0]-offset, p[1]-offset, p[0]+offset, p[1]+offset))
        pts_bounds.append((p[0], p[1], p[0], p[1]))
    return pts_bounds

# Built rtree using stream for interleaved = True
def data_generator(z0, triangles):
    idx = 0
    for key in triangles.keys():
        for t in triangles[key]:
            yield idx, t.bounds, (t, z0[key])
            idx += 1

def build_rtree(z0, mesh):
    print("-----Start building the rtree-----")
    polygons = {}
    for key in z0.keys():
        vertices = mesh.materials[key].vertices
        polygons[key] = tri2poly(vertices)
    idx = index.Index(data_generator(z0, polygons))
    print("-----Finished building the rtree-----")
    return idx

# Returns point intersections with rtree of triangles bounding boxes
def find_intersections(pts, pts_bounds, rtree_idx):
    p_idx = 0
    for p_bbox in pts_bounds:
        hits = list(rtree_idx.intersection(p_bbox, objects=True)) # crosses or contains
        z0 = []
        distances = []
        if len(hits) == 0:
            nearest_hits = list(rtree_idx.nearest(p_bbox, 1, objects=True)) # nearest triangles to point
            for o in nearest_hits:
                tri_ext = LinearRing(o.object[0].exterior.coords)
                d = tri_ext.project(Point(p_bbox[0], p_bbox[1]))
                z0.append(o.object[1])
                distances.append(d)
            pts[p_idx] = z0[distances.index(min(distances))]
            p_idx += 1
            continue
        else:
           # Compute distances to hits and assign z0 of nearest triangle
            z0 = []
            distances = []
            for o in hits:
                tri_ext = LinearRing(o.object[0].exterior.coords)
                dist = tri_ext.project(Point(p_bbox[0], p_bbox[1]))
                z0.append(o.object[1])
                distances.append(dist)
            pts[p_idx] = z0[distances.index(min(distances))]
            p_idx += 1

    missed = []
    for p in pts:
        if p == -1:
            missed.append(p)
    print("Out of", len(pts), "face centers,", len(missed), "were not assigned a roughness length value.")
    return pts

def write_to_txt(list_z0, output):
    with open(output, "w") as f:
        print("-----Writing z0 values to", "z0.txt", "file-----")
        f.write("{0}\n".format(len(list_z0)))
        f.write("{0}\n".format("("))
        for i in list_z0:
            f.write("{0}\n".format(i))
        f.write("{0}\n".format(");"))

def write_to_xyz(Cf, z0):
    with open("xyz0.xyz", "w") as f:
        print("-----Writing face centers to", "z0.xyz", "file-----")
        f.write("{0} {1} {2}\n".format("x", "y", "z0"))
        for i in range(0, len(Cf)):
            f.write("{0} {1} {2}\n".format(Cf[i][0], Cf[i][1], z0[i]))

def main():
    # read parameters from the file 'setz0Dict.json'
    try:
        jparams = json.load(open('setZ0Dict.json'))
    except:
        print("ERROR: something is wrong with the setZ0Dict.json file.")
        sys.exit

    # load in memory the triangulated mesh
    if jparams["format"] == "obj":
        print("-----Loading obj file-----")
        scene = pywavefront.Wavefront(jparams["input_mesh"], strict=True, create_materials=True)  # ,collect_faces=True)
        print("-----Finished loading the obj file-----")

        # check if obj has the metadata for material or materials were created upon load
        k = "default0"
        if k not in scene.materials.keys():
            # load face centers
            face_centers = np.genfromtxt(jparams["input_pts"], delimiter=' ', skip_header=1, usecols=[0,1])
            n_faces = np.size(face_centers[:, :1])

            # Set face centers z0 to -1
            faces_z0 = [-1] * n_faces

            # Create and build rtree
            start_time = time.time()
            rtree_idx = build_rtree(jparams["landcover"], scene)
            time_rtree = time.time()
            print("Finished building rtree in: ", time_rtree - start_time,"s")

            pts_bounds = pts_bbox(face_centers)  # Find bounds of face centers
            finZ0_list = find_intersections(faces_z0, pts_bounds, rtree_idx)  # Find intersections face centers - mesh filtered triangles
            time_search_fin = time.time()
            print("Finished assigning z0 in: ", time_search_fin - time_rtree, "s")
            # write_to_xyz(face_centers, finZ0_list)
            if jparams["switch_txtWrite"]["switch"] == "on":
                write_to_txt(finZ0_list, jparams["switch_txtWrite"]["output_file"])
            else:
                print("The output is not written to file.")

        else:
            print("The input obj does not contain landcover information.")
            sys.exit

    else:
        print("The provided file format is not compatible.")
        sys.exit

if __name__ == '__main__':
    main()

