# Calculate the current for a trajectory. 
# Results are in "time(ns) current(nA)"

# Parameters
set startFrame 0
## dt in fs
set dt 50

# Output
set outFile curr.dat

# System size
set lx 75.7006384
set ly 75.7006384
set lz 115.125

# Load the system
mol load pdb traj-fit-nodt.pdb 
set sel [atomselect top "name SOD"]

# Load the trajectory
mol addfile traj-fit-nodt.xtc waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Open the output file
set out [open $outFile w]

# Summation loop
## puts $out "sum of q*v for $psf with trajectory $dcd"
## puts $out "t(ns) I(A)"

for {set i 0} {$i < 1} {incr i} {
     # Get the charge of each atom.
     set q [$sel get charge]

     # Get the position data for the first frame.
     molinfo top set frame $startFrame
     set z0 [$sel get z]
}

# Begin summation
set n 1
for {set f [expr $startFrame+1]} {$f < $nFrames && $n > 0} {incr f} {
     molinfo top set frame $f
     
     # Get the position data for the current frame.
     set z1 [$sel get z]
     
     # Find the displacements in the z-direction.
     set dz {}
     foreach r0 $z0 r1 $z1 {
          # Compensate for jumps across the periodic cell.
          set z [expr $r1-$r0]
          if {[expr $z > 0.5*$lz]} {set z [expr $z-$lz]}
          if {[expr $z <-0.5*$lz]} {set z [expr $z+$lz]}
          
          lappend dz $z
     }

     # Compute the average charge*velocity between the two frames.
     set qvsum [expr [vecdot $dz $q] / $dt]
          
     # We first scale by the system size to obtain the z-current in e/fs.
     set currentZ [expr $qvsum/$lz]

     # Now we convert to nanoamperes.
     set currentZ [expr $currentZ*1.60217733e5]

     # Get the time in nanoseconds for this frame.
     set t [expr ($f+0.5)*$dt*1.e-6]
               
     # Write the current.
     puts $out "$t $currentZ"
     puts -nonewline [format "FRAME %i: " $f]
     puts "$t $currentZ"
     
     # Store the postion data for the next computation.
     set z0 $z1
}
close $out
