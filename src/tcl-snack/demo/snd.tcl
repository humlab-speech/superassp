#!/usr/bin/env wize

package require snack
set dirname [file dirname [file normalize [info script]]]
set file [file join $dirname ex1]
set s [::snack::sound]
toplevel [set t .snd]
proc Play {t s ext} {
  global file
  wm title $t "SOUND DEMO: $ext"
  sounds op read $s $file.$ext
  sounds op play $s -start 6000
}

pack [canvas $t.c] -fill both -expand y
pack [frame $t.f] -fill x
pack [button $t.f.w -text WAV -command "Play $t $s wav"] -side left
pack [button $t.f.m -text MP3 -command "Play $t $s mp3"] -side left
pack [button $t.f.g -text OGG -command "Play $t $s ogg"] -side left
pack [button $t.f.e -text Source -command [list edit [info script]]] -side left
Play $t $s wav
$t.c create waveform "0 0" -anchor nw -sound $s


