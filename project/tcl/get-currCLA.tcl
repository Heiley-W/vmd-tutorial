proc mypermeation {frame} {
  global out current_out current upperEnd lowerEnd wat IdList LabelList EntryList CoorList num1 num2 dt Lz
  set oldLabelList $LabelList
  set oldEntryList $EntryList
  set oldCoorList $CoorList
  set LabelList {}
  set EntryList {}
  set CoorList {}
  set dt 5e-11
  set Lz 20
  set q -1.60217663e-19
  
# 2 = Crossed upper end 
# -2 = Crossed lower end 
# -1 = Crossed lower end previously 
# 1 = Crossed upper end previously 

  foreach z [$wat get z] oldLab $oldLabelList id $IdList oldEntry $oldEntryList oldCoor $oldCoorList {
    if {$z > $upperEnd} {
      set newLab 2
      set newEntry $frame
      set newCoor $z
      if {$oldLab == -1} {
          incr num1
          puts $out "wat $id permeated along +z direction at frame $frame in [expr $frame - $oldEntry] frames"
          set current [expr {$current + 1/($dt*$Lz)*$q*($newCoor-$oldCoor)}]
      }
    } elseif {$z < $lowerEnd} {
      set newLab -2
      set newEntry $frame
      set newCoor $z
      if {$oldLab == 1} {
          incr num2
          puts $out "wat $id permeated along -z direction at frame $frame in [expr $frame - $oldEntry] frames"
          set current [expr {$current + 1/($dt*$Lz)*$q*($newCoor-$oldCoor)}]
      }
    } elseif {abs($oldLab) > 1} {
      set newLab [expr $oldLab / 2]
      set newEntry $oldEntry
      set newCoor $z
      set current [expr {$current + 1/($dt*$Lz)*$q*($newCoor-$oldCoor)}]
    } else {
      set newLab $oldLab
      set newEntry $oldEntry
      set newCoor $z
    }
    lappend LabelList $newLab
    lappend EntryList $newEntry
    lappend CoorList $z
  }
  puts "Done with frame $frame"
  puts $current_out "$frame $current"
}
  
source ./bigdcd.tcl

set d 10;

set out [open permeationCLA-$d.dat w]
set out2 [open overallCLA-$d.dat w]
set current_out [open currentCLA-$d.dat w]

mol load pdb traj-fit-nodt.pdb
#animate delete beg 0 end 0 waitfor all  ;# get rid of the unaligned initial frame

set procen [atomselect top "protein and resid 1 to 32 and name CA"]
set cpro [measure center $procen]
set cproz [lindex $cpro 2]

set upperEnd [expr $cproz + $d]
set lowerEnd [expr $cproz - $d]

set wat [atomselect top "name CLA"]
set IdList [$wat get index]
set LabelList {}
set EntryList {}
set CoorList {}

foreach foo $IdList {
  lappend LabelList 0
  lappend EntryList 0
  lappend CoorList 0
}

set num1 0
set num2 0
set current 0

bigdcd mypermeation traj-fit-nodt-500ns.xtc
bigdcd_wait 
close $out

puts $out2 "$num1 CLA crossed along +z direction"
puts $out2 "$num2 CLA crossed along -z direction"
close $out2

close $current_out