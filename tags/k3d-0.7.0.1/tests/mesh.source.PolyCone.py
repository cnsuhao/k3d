#python

import testing

setup = testing.setup_mesh_source_image_test("PolyCone")
testing.mesh_comparison(setup.document, setup.source.get_property("output_mesh"), "mesh.source.PolyCone", 1)
