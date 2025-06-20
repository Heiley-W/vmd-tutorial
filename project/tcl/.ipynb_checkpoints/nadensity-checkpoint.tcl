set na [atomselect top "name SOD"]

volmap density $na -allframes -combine avg -res 1 -o efnadens.dx

mol load dx efnadens.dx
