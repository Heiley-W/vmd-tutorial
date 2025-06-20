set cl [atomselect top "name CLA"]

volmap density $cl -allframes -combine avg -res 1 -o efcldens.dx

mol load dx efcldens.dx
