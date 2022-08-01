extensions [vid array palette]
breed [ TCRs TCR ]
breed [ pLATs pLAT ]
pLATs-own [ phosCount linkCount leader followers direction speed stepsize pLATsThatTriedMe n_followers D_cluster]
patches-own [patch_LATden]

breed [ scalebars scalebar ]
breed [ banners banner ]
breed [ stopwatches stopwatch ]
breed [ circles circle ]

; Let's not make different breeds for different phosphorylation states, because that will be too many breeds.
globals [
  trialNidx
  time
  numXlinked
  top10leaders ; 여기서부터 이어서 짜보자!

  pLAT_phos_prob
  sm_pLAT_Phosp_Prob

  vidStartFlag
  N_uLAT2pLAT-inTCRpatch
  pLATsz_inPUnit
  TCRsz_inPUnit
  linkProb

;  fnameheader
  cluster_formed?
  timepoint_list
  Pcluster_list
  clusterFlagCnt_list

;  tlapse-dt
  next-tlapse-t
;  vid-dt
  next-vid-t
]

to adjust_world_size
  set-patch-size worldWidth_inPx / (max [pxcor] of patches - min [pxcor] of patches + 1)
end

to multiple-runs
  set timepoint_list (range 0 endtime timestep)           ; [0 0.01 0.02 ... endtime]
  set Pcluster_list n-values (length timepoint_list) [0]  ; [0 0 0 0 0 ... 0]
  set clusterFlagCnt_list n-values (length timepoint_list) [0]  ; [0 0 0 0 0 ... 0]

  let timepointN 0
;  let trialNidx 0
  repeat Ntrials [
    setup
    set timepointN 0
    while [time < endtime - timestep] [
      tick
      set time time + timestep
      ask stopwatches [ set label (word precision time 3 " s") ]

      Fn_uLATs2pLATs
      Fn_phos_pLATs
      Fn_formLinks
      Fn_move
      Fn_setcolors
      Fn_killAtEdges
      if calc-patchLATden? [ Fn_calc-patchLATdensity ]
      Fn_calc-Pcluster timepointN trialNidx
      display
      set-current-plot  "P_cluster-vs-time"
      set-plot-pen-mode 2
      plot (item timepointN Pcluster_list)
      set timepointN timepointN + 1
    ]
;    clear-plot
    set trialNidx trialNidx + 1
    ;    print trialNidx
  ]
  export-plot "P_cluster-vs-time" (word fnameheader " P_cluster.csv")
end






to Fn_calc-Pcluster [timepoint_n trial_n]
  if cluster_formed?
  [ set clusterFlagCnt_list replace-item timepoint_n clusterFlagCnt_list (item timepoint_n clusterFlagCnt_list + 1)  ]
  set Pcluster_list replace-item timepoint_n Pcluster_list (item timepoint_n clusterFlagCnt_list / (trial_n + 1))
;  print Pcluster_list
end

to-report num-date-time  ; current date in numerical format, yyyy-mm-dd
  let $dt substring date-and-time 16 27  ; "21-May-2013"
  let $dt2 substring date-and-time 0 8  ; 01:19
  report (word (substring $dt 7 11)           ; yyyy
    "-" (month-num substring $dt 3 6)  ; mm
    "-" (substring $dt 0 2)           ; dd
    "_" (substring $dt2 0 2) "_" (substring $dt2 3 5) "_" (substring $dt2 6 8))
end

to-report month-num [ #mon ]
  let $index 1 + position #mon
    ["Jan""Feb""Mar""Apr""May""Jun""Jul""Aug""Sep""Oct""Nov""Dec"]
  report substring (word (100 + $index)) 1 3  ; force 2-digit string
end

to setup
;  ca
  ct
;  set tlapse-dt 5
  set next-tlapse-t 0
;  set vid-dt 0.05
  set next-vid-t 0
  ask patches [ set patch_LATden denLAT ]
  set top10leaders no-turtles
  set-default-shape turtles "circle"
  set time 0
  set cluster_formed? false
  set fnameheader (word num-date-time " k_phos " k_phos " k_link-in-arm " k_link-in-arms-length " patchLen_inUm " patchLen_inUm " timestep " timestep  " denLAT " denLAT " D_LAT " D_LAT " armLen_inUm " armLen_inUm)

  set pLATsz_inPUnit pLATsz_inUm / patchLen_inUm
  set TCRsz_inPUnit TCRsz_inUm / patchLen_inUm
  set N_uLAT2pLAT-inTCRpatch k_phos * timestep * denLAT * patchLen_inUm ^ 2
  if N_uLAT2pLAT-inTCRpatch < 1 [ print (word "!!!!!! WARNING !!!!!! - N_uLAT2pLAT-inTCRpatch < 1: " N_uLAT2pLAT-inTCRpatch) ]
  set sm_pLAT_Phosp_Prob k_phos * timestep
  if sm_pLAT_Phosp_Prob > 1 [ print (word "!!!!!! WARNING !!!!!! - pLAT_phos_prob > 1: " sm_pLAT_Phosp_Prob)]
  set linkProb k_link-in-arms-length * timestep ;- When using the within-radius approach 2020/5/2 NKIM
  if linkProb >= 1 [user-message (word "!!!!!! WARNING !!!!!! - linkProb > 1 : " linkProb)] ;- when using the within-radius approach 2020.05.02 NKim

  create-scalebars 1 [
    set size 10
    setxy (max-pxcor - size) (min-pycor + size)
    set shape "line"
    set heading 90
    set color white
    set label word (size * patchLen_inUm) " um"
  ]

  create-stopwatches 1 [
    setxy (min-pxcor + 2 * [size] of one-of scalebars) (max-pycor - [size] of one-of scalebars)
    set label (word precision time 3 " s")
    set size 0
  ]

  create-TCRs 1 [
    setxy 0 0
    set color [0 200 255 150]
    set size TCRsz_inPUnit
    set heading 90
  ]
  reset-ticks
  set time 0

  if video? [ vidstart ]
  ;  Fn_saveimgvid
end

to slow-down
  let NtotalpLATs count pLATs with [self = leader]
  ask n-of (NtotalpLATs / 2) pLATs [ set D_cluster D_cluster / 10 ]
end


to test-Ds
  let Dlist [0.1 0.5 1 2 4 6 8]
  let index 0
  repeat length Dlist [
    setup
    let picked_D (item index Dlist)
    set D_LAT picked_D
    while [time < endtime] [
      tick
      set time time + timestep
      ask stopwatches [ set label (word precision time 3 " s") ]
      Fn_uLATs2pLATs
      Fn_phos_pLATs
      Fn_formLinks
      Fn_move
      Fn_setcolors
      Fn_saveimgvid
      Fn_killAtEdges
      if calc-patchLATden? [ Fn_calc-patchLATdensity ]
      if save-all-representations? [ Fn_saveAllRepresentation ]
      display
    ]
    set-current-plot  "xlink-vs-D"
    set-plot-pen-mode 2
    plotxy picked_D count links
    set index index + 1
  ]
end


to sel-dir
  set-current-directory user-directory
end

to go
  tick
  set time time + timestep
  ask stopwatches [ set label (word precision time 3 " s") ]

  ;if time = 5 [slow-down]

  Fn_uLATs2pLATs
  Fn_phos_pLATs
  Fn_formLinks
  Fn_move
  Fn_setcolors
  Fn_saveimgvid
  Fn_killAtEdges
  if calc-patchLATden? [ Fn_calc-patchLATdensity ]
  if save-all-representations? [ Fn_saveAllRepresentation ]
  display



  if time > endtime
  [ if vidStartFlag = true [ vidend ]
    if export-plots-and-interface? [ Fn_export-plots-and-interface ]
    stop ]


end

to Fn_export-plots-and-interface
  export-all-plots (word fnameheader ".csv")
  export-interface (word "Interface-" fnameheader ".png")
end

to Fn_saveAllRepresentation
  if time >= next-tlapse-t
  [ phosCountView
    export-view (word "CluSzView " fnameheader " t" (precision time 1)".png")
    ClusterSizeView
    export-view (word "PhosView " fnameheader " t" (precision time 1)".png")
    brightPatchView
    export-view (word "PatchView " fnameheader " t" (precision time 1)".png")
    set next-tlapse-t (next-tlapse-t + tlapse-dt)
  ]
end

to Fn_saveimgvid
  if tlapse? and time >= next-tlapse-t
  [ export-view (word fnameheader " t" (precision time 1)".png")
    set next-tlapse-t (next-tlapse-t + tlapse-dt)
  ]
  if vidStartFlag = true and time >= next-vid-t
  [ vid:record-view
    set next-vid-t (next-vid-t + vid-dt)
  ]
end

to record-xlink-vs-denLAT
  set-current-plot "xlink-vs-denLAT"
  set-plot-pen-mode 2
  plotxy denLAT count links
end

to record-xlink-vs-D
  set-current-plot "xlink-vs-D"
  set-plot-pen-mode 2
  plotxy D_LAT count links
end

to Fn_uLATs2pLATs
  let integer floor N_uLAT2pLAT-inTCRpatch
  let decimal remainder N_uLAT2pLAT-inTCRpatch 1
  create-pLATs integer
  [ set phosCount 1
    set linkCount 0
    set leader self
    set followers (turtle-set self)
    set color [255 0 0 100]
    set size pLATsz_inPUnit
    set pLATsThatTriedMe []
  ]
  if random-float 1 < decimal
  [ create-pLATs 1
    [ set phosCount 1
      set linkCount 0
      set leader self
      set followers (turtle-set self)
      set color [255 0 0 100]
      set size pLATsz_inPUnit
      set pLATsThatTriedMe []
      set D_cluster D_LAT
    ]
  ]
end

to Fn_phos_pLATs
  ask TCRs
  [ ask patch-here
    [ ask pLATs-here
      [ if phosCount < 3 and random-float 1 < sm_pLAT_Phosp_Prob
        [ set phosCount phosCount + 1 ] ] ]
  ]
end

to Fn_formLinks ;-  This was using 'within - radius' approach, which is slow. I'll make one with 'same patch' approach. [2020.05.02 NeilKim]
  ask pLATs [ set pLATsThatTriedMe [] ]

  foreach sort pLATs [ L1 -> ; sorting by who numbers
    ask L1 [let L1avail_pY phosCount - linkCount ;
      let L1pYidx 0
      while [L1pYidx < L1avail_pY] [
        if any? other pLATs in-radius (armLen_inUm / patchLen_inUm) ; Then pick a pLAT.
        [
          let linkable_pLATs other pLATs in-radius (armLen_inUm / patchLen_inUm) with [phosCount - linkCount > 0 and not link-neighbor? self and not member? self [pLATsThatTriedMe] of myself] ;and linkTestFlag = false ;  user-message (word "for who=" who ", linkable_pLATs: " linkable_pLATs)
          if any? linkable_pLATs
          [
            foreach sort linkable_pLATs [ L2 ->
              ask L2 [
                let L2avail_pY phosCount - linkCount
                let L2pYidx 0
                while [L2pYidx < L2avail_pY and L1pYidx < L1avail_pY and not link-neighbor? L1]
                [
                  set pLATsThatTriedMe (turtle-set L1 pLATsThatTriedMe)
                  if random-float 1 < linkProb
                  [
                    if [phosCount] of L1 - [linkCount] of L1 < 1 [ user-message "[phosCount] of L1 - [linkCount] of L1 < 1" print (word "L1pYidx: " L1pYidx "L1avail_pY: " L1avail_pY) ] ;  This is happening a lot now.
                    if [phosCount] of L2 - [linkCount] of L2 < 1 [ user-message "[phosCount] of L2 - [linkCount] of L2 < 1" ]
                    create-link-with L1 [set color white]
                    adjust_position L1 L2
                    set linkCount (linkCount + 1)
                    set L2avail_pY (L2avail_pY - 1)
                    ask L1
                    [ set linkCount (linkCount + 1)
                      set L1avail_pY (L1avail_pY - 1) ]
                    ask [leader] of L1
                    [
                      set followers (turtle-set followers ([followers] of [leader] of L2))
                      set n_followers count followers
                      ifelse count pLATs with [leader = self] < 10
                      [ set top10leaders pLATs with [leader = self] ]
                      [ set top10leaders max-n-of 10 pLATs with [leader = self] [n_followers] ]
                      ask followers
                      [ set leader [leader] of L1  ]

                      set D_cluster D_LAT / n_followers

                  ] set leader [leader] of L1

                  ]
                  set L2pYidx L2pYidx + 1  ]  ]  ]  ]  ]
        set L1pYidx L1pYidx + 1  ]  ]  ]
end

to adjust_position [L1 L2]
  let mover nobody
  let unmover nobody
  ifelse count [followers] of [leader] of L1 > count [followers] of [leader] of L1
  [ set unmover L1
    set mover L2  ]
  [ set unmover L2
    set mover L1 ]
  let d_x  [xcor] of mover - [xcor] of unmover
  let d_y  [ycor] of mover - [ycor] of unmover
  let dist_inPatches  sqrt(d_x ^ 2 + d_y ^ 2) ; in patchlength
  let armLen_inPatches armLen_inUm / patchLen_inUm
  ask mover [
    ifelse dist_inPatches = 0
    [ let randomdir random 360
      set xcor [xcor] of unmover + armLen_inPatches * cos(randomdir)
      set ycor [ycor] of unmover + armLen_inPatches * sin(randomdir)]
    [set xcor [xcor] of unmover + d_x * armLen_inPatches / dist_inPatches
    set ycor [ycor] of unmover + d_y * armLen_inPatches / dist_inPatches ]
  ]
  if abs(sqrt(([xcor] of mover - [xcor] of unmover) ^ 2 + ([ycor] of mover - [ycor] of unmover) ^ 2) / armLen_inPatches - 1) > 1e-5 [print "position not adjd"]
end

to phosCountView
  ask turtles [ st ]
  ask links [ show-link ]
  ask TCRs [ st ]
  ask patches [ set pcolor black ]
  ask pLATs [
    ifelse phosCount = 1 [ set color [255 0 0 100] set size pLATsz_inPUnit ][
      ifelse phosCount = 2 [ set color [0 255 0 150] set size pLATsz_inPUnit * 1 ][
        if phosCount = 3 [set color [255 255 0 255] set size pLATsz_inPUnit * 1 ]]]
  ]
end

to ClusterSizeView
  ask turtles [ st ]
  ask links [ show-link ]
  ask TCRs [ st ]
  ask patches [ set pcolor black ]
  let maxNfollowers 0
  ask pLATs with [leader = self]
  [ if count followers > maxNfollowers [ set maxNfollowers count followers ]  ]
  ask pLATs with [leader = self] [
    let nfr count followers
    ;let x1 map [[x] -> nfr * x / maxNfollowers  ] [255 0 0 150]
    let x1 map [[x] -> (nfr - 1) * x / maxNfollowers  ] [255 0 0 150]
    let x2 map [[x] -> (maxNfollowers - nfr + 1) * x / maxNfollowers ] [0 0 255 50]
    ask followers
    [ set color (map + x1 x2) ]
  ]
end

to Fn_calc-patchLATdensity
  set cluster_formed? false
  ask patches[
    set patch_LATden count pLATs-here / (patchLen_inUm ^ 2) + denLAT
    if patch_LATden > 2500 [ set cluster_formed? true  ]
  ]
end

to brightPatchView
  ;set cluster_formed? false
  ask patches[
    ifelse patch_LATden > 2500
    [ set pcolor [255 0 0]
      ;set cluster_formed? true
    ]
    [ set pcolor palette:scale-gradient palette:scheme-colors "Divergent" "RdYlBu" 7 patch_LATden 2500 0  ]
  ]
  ask pLATs [ ht ]
  ask links [ hide-link ]
;  ask TCRs [ ht ]
end

to Fn_move
  ask pLATs with [leader = self][
    rt random 360
    ;set stepsize sqrt (4 * D_LAT * timestep / count followers) / patchLen_inUm
    set stepsize sqrt (4 * D_cluster * timestep) / patchLen_inUm
    ask followers
    [ ; Followers 다 똑같이 움직이라는 것.
      set heading [heading] of myself
      set stepsize [stepsize] of myself
      fd stepsize
    ]
  ]
end

to Fn_setcolors
  ifelse color_option = "phCount"
  [ phosCountView
    fn_setvisibility  ]
  [ ifelse color_option = "clusterSz"
    [ ClusterSizeView
      fn_setvisibility ]
    [ brightPatchView ]   ]
end

to fn_setvisibility
  ifelse visibility_option = "all"
  [ ask pLATs [ st ]
    ask links [ show-link ]
  ]
  [ ifelse visibility_option = "xlinked"
    [ ask links [ show-link ]
      ask pLATs with [leader = self]
      [ ifelse count followers = 1
        [ ht ]
        [ st ]
    ]]
    ; big 10
    [
      ask links [ hide-link ]
      ask pLATs [ ht ]
      if any? top10leaders [
        ask top10leaders [
          ask followers [ st ]
          ask links [ hide-link ]  ]
      ]
    ]
  ]
end

to Fn_killAtEdges
  ask pLATs with [leader = self] [
    if any? followers with [[pxcor] of patch-here = min-pxcor or [pxcor] of patch-here = max-pxcor or [pycor] of patch-here = min-pycor or [pycor] of patch-here = max-pycor ]
    [ ask followers [ die ]   ]
  ]
end

to vidend
  let vidfname (word fnameheader ".mp4")
  vid:save-recording vidfname
  set vidStartFlag false
end

to vidstart
  vid:start-recorder
  set vidStartFlag true
  vid:record-view
end

; not used yet 2020/5/3
;to remove-links-between [ a b ]
;   ask a [
;    ask my-links with [ other-end = b ] [ die ]
;    set phosCount phosCount - 1
;    set linkCount linkCount - 1
;  ]
;  ask b [
;;    set phosCount phosCount - 1
;;    set linkCount linkCount - 1
;  ]
;end
@#$#@#$#@
GRAPHICS-WINDOW
528
23
1236
732
-1
-1
4.3478260869565215
1
15
1
1
1
0
0
0
1
-80
80
-80
80
0
0
1
ticks
30.0

BUTTON
25
20
88
53
NIL
setup
NIL
1
T
OBSERVER
NIL
S
NIL
NIL
1

BUTTON
25
57
88
90
go
go
T
1
T
OBSERVER
NIL
G
NIL
NIL
1

SLIDER
44
232
176
265
denLAT
denLAT
10
1500
300.0
10
1
NIL
HORIZONTAL

INPUTBOX
18
328
173
388
k_phos
33.34
1
0
Number

INPUTBOX
302
566
457
626
timestep
0.01
1
0
Number

INPUTBOX
18
398
173
458
k_link-in-arms-length
5.0
1
0
Number

INPUTBOX
300
498
455
558
patchLen_inUm
0.1
1
0
Number

BUTTON
16
276
85
309
NIL
vidend
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
208
15
521
108
fnameheader
2022-07-31_10_13_20 k_phos 33.34 k_link-in-arm 5 patchLen_inUm 0.1 timestep 0.01 denLAT 300 D_LAT 2 armLen_inUm 0.03
1
1
String

INPUTBOX
23
535
178
595
k_dephos
0.0
1
0
Number

INPUTBOX
303
637
455
697
D_LAT
2.0
1
0
Number

CHOOSER
180
541
272
586
dephos_type
dephos_type
"power" "inverse_proportional" "none" "equal"
3

TEXTBOX
118
335
190
353
/s
11
0.0
1

TEXTBOX
129
404
189
422
/s
11
0.0
1

INPUTBOX
211
127
299
187
pLATsz_inUm
0.03
1
0
Number

INPUTBOX
213
193
304
253
TCRsz_inUm
0.2
1
0
Number

PLOT
1290
64
1490
214
number-of-pLATs
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"pLAT" 1.0 0 -2674135 true "" "plot count pLATs with [phosCount = 1]"
"ppLAT" 1.0 0 -13840069 true "" "plot count pLATs with [phosCount = 2]"
"pppLAT" 1.0 0 -4079321 true "" "plot count pLATs with [phosCount = 3]"

PLOT
1297
247
1497
397
number-of-cross-links
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"links" 1.0 0 -16777216 true "" "plot count links"

INPUTBOX
386
421
470
481
armLen_inUm
0.03
1
0
Number

SLIDER
365
329
502
362
worldWidth_inPx
worldWidth_inPx
100
1200
700.0
50
1
NIL
HORIZONTAL

BUTTON
364
369
495
402
NIL
adjust_world_size
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
439
126
521
186
endtime
10.0
1
0
Number

CHOOSER
333
189
471
234
color_option
color_option
"phCount" "clusterSz" "patchSNR"
1

CHOOSER
334
245
472
290
visibility_option
visibility_option
"all" "xlinked" "big10"
0

MONITOR
94
470
290
515
NIL
max([patch_LATden] of patches)
5
1
11

MONITOR
204
406
338
451
NIL
linkProb
7
1
11

MONITOR
204
335
339
380
NIL
sm_pLAT_Phosp_Prob
7
1
11

BUTTON
1468
488
1589
521
NIL
record-xlink-vs-denLAT
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
1298
423
1458
543
xlink-vs-denLAT
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

PLOT
1298
558
1458
685
xlink-vs-D
NIL
NIL
0.0
2.0
0.0
1600.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

BUTTON
1269
724
1332
757
NIL
ca
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
22
97
147
130
NIL
record-xlink-vs-D
NIL
1
T
OBSERVER
NIL
1
NIL
NIL
1

SWITCH
113
149
203
182
video?
video?
1
1
-1000

MONITOR
40
678
125
723
NIL
time
7
1
11

MONITOR
102
272
185
317
NIL
count turtles
17
1
11

INPUTBOX
1525
61
1624
121
Ntrials
50.0
1
0
Number

PLOT
1529
131
1729
281
P_cluster-vs-time
NIL
NIL
0.0
3.0
0.0
1.1
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

BUTTON
1533
289
1635
322
NIL
multiple-runs
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
207
287
365
332
NIL
N_uLAT2pLAT-inTCRpatch
4
1
11

PLOT
1558
349
1718
469
max-patch-LAT-density
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot max([patch_LATden] of patches)"
"pen-1" 1.0 0 -2674135 true "" "plot 2500"

BUTTON
103
57
178
90
go once
go
NIL
1
T
OBSERVER
NIL
F
NIL
NIL
1

BUTTON
1621
699
1686
732
NIL
sel-dir
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
112
191
202
224
tlapse?
tlapse?
1
1
-1000

INPUTBOX
260
734
415
794
tlapse-dt
1.0
1
0
Number

INPUTBOX
80
737
235
797
vid-dt
0.025
1
0
Number

SWITCH
299
112
437
145
calc-patchLATden?
calc-patchLATden?
0
1
-1000

SWITCH
1486
764
1656
797
save-all-representations?
save-all-representations?
1
1
-1000

SWITCH
1486
731
1671
764
export-plots-and-interface?
export-plots-and-interface?
1
1
-1000

BUTTON
25
146
97
179
NIL
test-Ds
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1652
80
1710
125
NIL
trialNidx
1
1
11

BUTTON
1306
7
1396
40
NIL
slow-down
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
