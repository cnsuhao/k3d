#python

doc = Document
doc.start_change_set()
try:
	material = doc.new_node("RenderManMaterial")
	material.name = "Patch Material"
	material.color = (1, 1, 1)

	frozen_mesh = doc.new_node("FrozenMesh")
	frozen_mesh.name = "Bilinear Patch"

	mesh = frozen_mesh.new_mesh()

	points = ((-5, 0, 5), (5, 0, 5), (0, -5, -5), (0, 5, -5))
	for position in points:
		mesh.new_point(position)

	patch = mesh.new_bilinear_patch()
	patch.material = material
	for i in range(4):
		patch.control_points[i] = mesh.points[i]

#	patch.uniform_data.set_color("Cs", (1, 0, 0))
#	mesh.points[0].vertex_data.set_color("Cs", (1, 0, 0))
#	mesh.points[1].vertex_data.set_color("Cs", (0, 1, 0))
#	mesh.points[2].vertex_data.set_color("Cs", (0, 0, 1))
#	mesh.points[3].vertex_data.set_color("Cs", (1, 1, 1))
	patch.varying_data[0].set_color("Cs", (1, 0, 0))
	patch.varying_data[1].set_color("Cs", (0, 1, 0))
	patch.varying_data[2].set_color("Cs", (0, 0, 1))
	patch.varying_data[3].set_color("Cs", (1, 1, 1))

	mesh_instance = doc.new_node("MeshInstance")
	mesh_instance.name = "Bilinear Patch Instance"
	doc.set_dependency(mesh_instance.get_property("input_mesh"), frozen_mesh.get_property("output_mesh"))

	doc.finish_change_set("Create Bilinear Patch")

except:
	doc.cancel_change_set()
	raise
