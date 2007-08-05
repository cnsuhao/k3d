#python

import k3d
import testing

setup = testing.setup_mesh_modifier_test("PolyTorus", "PGPRemesh")
testing.mesh_comparison(setup.document, setup.modifier.get_property("output_mesh"), "mesh.modifier.PGPRemesh", 1)

