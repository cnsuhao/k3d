#python

import k3d
import testing

setup = testing.setup_two_mesh_modifier_test("PolyTorus", "TriangulateFaces", "PGPRemesh")
testing.mesh_comparison(setup.document, setup.modifier2.get_property("output_mesh"), "mesh.modifier.PGPRemesh", 1)

