#python

import testing

setup = testing.setup_mesh_source_test("KnotCurve")
testing.mesh_comparison(setup.document, setup.source.get_property("output_mesh"), "mesh.source.KnotCurve", 300)
