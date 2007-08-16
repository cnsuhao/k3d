#python

import k3d
import testing

document = k3d.new_document()

setup = testing.setup_mesh_reader_test("PLYMeshReader", "mesh.modifier.PGPRemesh.torus.ply")

modifier = setup.document.new_node("PGPRemesh")
modifier.use_smooth = True
modifier.smooth_4 = False
modifier.steps = 5
modifier.h = 1500
modifier.omega = 10
modifier.div = 2
document.set_dependency(modifier.get_property("input_mesh"), setup.reader.get_property("output_mesh"))

#print "source output: " + repr(source.output_mesh)
#print "triangles output: " + repr(triangles.output_mesh)
#print "modifier output: " + repr(modifier.output_mesh)

testing.mesh_comparison(document, modifier.get_property("output_mesh"), "mesh.modifier.PGPRemesh", 1)

