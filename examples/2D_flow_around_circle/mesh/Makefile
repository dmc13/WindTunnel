mesh:	
	gmsh -2 mesh.geo
	dolfin-convert -o xml mesh.msh mesh.xml
	python ../../../scripts/convert_to_new_xml.py
