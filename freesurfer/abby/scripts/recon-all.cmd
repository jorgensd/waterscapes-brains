

#---------------------------------
# New invocation of recon-all ma. 18. mars 16:07:35 +0100 2019 

 mri_convert /cluster/home/abby/abby-t1.mgz /work/jobs//abby/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor ma. 18. mars 16:08:45 +0100 2019

 cp /work/jobs//abby/mri/orig/001.mgz /work/jobs//abby/mri/rawavg.mgz 


 mri_convert /work/jobs//abby/mri/rawavg.mgz /work/jobs//abby/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /work/jobs//abby/mri/transforms/talairach.xfm /work/jobs//abby/mri/orig.mgz /work/jobs//abby/mri/orig.mgz 

#--------------------------------------------
#@# Talairach ma. 18. mars 16:09:44 +0100 2019

 mri_nu_correct.mni --no-rescale --i orig.mgz --o orig_nu.mgz --n 1 --proto-iters 1000 --distance 50 


 talairach_avi --i orig_nu.mgz --xfm transforms/talairach.auto.xfm 

talairach_avi log file is transforms/talairach_avi.log...

 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection ma. 18. mars 16:12:22 +0100 2019

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /cluster/software/VERSIONS/freesurfer-6.0.0/bin/extract_talairach_avi_QA.awk /work/jobs//abby/mri/transforms/talairach_avi.log 


 tal_QC_AZS /work/jobs//abby/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Nu Intensity Correction ma. 18. mars 16:12:27 +0100 2019

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 


 mri_add_xform_to_header -c /work/jobs//abby/mri/transforms/talairach.xfm nu.mgz nu.mgz 

#--------------------------------------------
#@# Intensity Normalization ma. 18. mars 16:15:32 +0100 2019

 mri_normalize -g 1 -mprage nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping ma. 18. mars 16:17:42 +0100 2019

 mri_em_register -rusage /work/jobs//abby/touch/rusage.mri_em_register.skull.dat -skull nu.mgz /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_withskull_2016-05-10.vc700.gca transforms/talairach_with_skull.lta 


 mri_watershed -rusage /work/jobs//abby/touch/rusage.mri_watershed.dat -T1 -brain_atlas /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_withskull_2016-05-10.vc700.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration ma. 18. mars 16:46:04 +0100 2019

 mri_em_register -rusage /work/jobs//abby/touch/rusage.mri_em_register.dat -uns 3 -mask brainmask.mgz nu.mgz /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize ma. 18. mars 17:05:08 +0100 2019

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg ma. 18. mars 17:06:47 +0100 2019

 mri_ca_register -rusage /work/jobs//abby/touch/rusage.mri_ca_register.dat -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_2016-05-10.vc700.gca transforms/talairach.m3z 

#--------------------------------------
#@# SubCort Seg ma. 18. mars 19:40:49 +0100 2019

 mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz transforms/talairach.m3z /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /work/jobs//abby/mri/transforms/cc_up.lta abby 

#--------------------------------------
#@# Merge ASeg ma. 18. mars 20:43:06 +0100 2019

 cp aseg.auto.mgz aseg.presurf.mgz 

#--------------------------------------------
#@# Intensity Normalization2 ma. 18. mars 20:43:06 +0100 2019

 mri_normalize -mprage -aseg aseg.presurf.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS ma. 18. mars 20:46:30 +0100 2019

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation ma. 18. mars 20:46:33 +0100 2019

 mri_segment -mprage brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.presurf.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill ma. 18. mars 20:48:33 +0100 2019

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh ma. 18. mars 20:49:18 +0100 2019

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Tessellate rh ma. 18. mars 20:49:42 +0100 2019

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh ma. 18. mars 20:49:49 +0100 2019

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Smooth1 rh ma. 18. mars 20:49:55 +0100 2019

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh ma. 18. mars 20:50:00 +0100 2019

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# Inflation1 rh ma. 18. mars 20:50:35 +0100 2019

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh ma. 18. mars 20:51:09 +0100 2019

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# QSphere rh ma. 18. mars 20:55:05 +0100 2019

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology Copy lh ma. 18. mars 20:59:19 +0100 2019

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 

#--------------------------------------------
#@# Fix Topology Copy rh ma. 18. mars 20:59:27 +0100 2019

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 

#@# Fix Topology lh ma. 18. mars 20:59:33 +0100 2019

 mris_fix_topology -rusage /work/jobs//abby/touch/rusage.mris_fix_topology.lh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 abby lh 

#@# Fix Topology rh ma. 18. mars 21:09:37 +0100 2019

 mris_fix_topology -rusage /work/jobs//abby/touch/rusage.mris_fix_topology.rh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 abby rh 


 mris_euler_number ../surf/lh.orig 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make White Surf lh ma. 18. mars 21:20:40 +0100 2019

 mris_make_surfaces -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs abby lh 

#--------------------------------------------
#@# Make White Surf rh ma. 18. mars 21:26:11 +0100 2019

 mris_make_surfaces -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs abby rh 

#--------------------------------------------
#@# Smooth2 lh ma. 18. mars 21:31:29 +0100 2019

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white.preaparc ../surf/lh.smoothwm 

#--------------------------------------------
#@# Smooth2 rh ma. 18. mars 21:31:35 +0100 2019

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white.preaparc ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh ma. 18. mars 21:31:40 +0100 2019

 mris_inflate -rusage /work/jobs//abby/touch/rusage.mris_inflate.lh.dat ../surf/lh.smoothwm ../surf/lh.inflated 

#--------------------------------------------
#@# Inflation2 rh ma. 18. mars 21:32:15 +0100 2019

 mris_inflate -rusage /work/jobs//abby/touch/rusage.mris_inflate.rh.dat ../surf/rh.smoothwm ../surf/rh.inflated 

#--------------------------------------------
#@# Curv .H and .K lh ma. 18. mars 21:32:51 +0100 2019

 mris_curvature -w lh.white.preaparc 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 lh.inflated 

#--------------------------------------------
#@# Curv .H and .K rh ma. 18. mars 21:34:15 +0100 2019

 mris_curvature -w rh.white.preaparc 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 rh.inflated 


#-----------------------------------------
#@# Curvature Stats lh ma. 18. mars 21:35:36 +0100 2019

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm abby lh curv sulc 


#-----------------------------------------
#@# Curvature Stats rh ma. 18. mars 21:35:41 +0100 2019

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm abby rh curv sulc 

#--------------------------------------------
#@# Sphere lh ma. 18. mars 21:35:45 +0100 2019

 mris_sphere -rusage /work/jobs//abby/touch/rusage.mris_sphere.lh.dat -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Sphere rh ma. 18. mars 22:27:37 +0100 2019

 mris_sphere -rusage /work/jobs//abby/touch/rusage.mris_sphere.rh.dat -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg lh ma. 18. mars 23:16:29 +0100 2019

 mris_register -curv -rusage /work/jobs//abby/touch/rusage.mris_register.lh.dat ../surf/lh.sphere /cluster/software/VERSIONS/freesurfer-6.0.0/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Surf Reg rh ti. 19. mars 00:08:04 +0100 2019

 mris_register -curv -rusage /work/jobs//abby/touch/rusage.mris_register.rh.dat ../surf/rh.sphere /cluster/software/VERSIONS/freesurfer-6.0.0/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh ti. 19. mars 00:56:50 +0100 2019

 mris_jacobian ../surf/lh.white.preaparc ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# Jacobian white rh ti. 19. mars 00:57:06 +0100 2019

 mris_jacobian ../surf/rh.white.preaparc ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh ti. 19. mars 00:57:18 +0100 2019

 mrisp_paint -a 5 /cluster/software/VERSIONS/freesurfer-6.0.0/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#--------------------------------------------
#@# AvgCurv rh ti. 19. mars 00:57:22 +0100 2019

 mrisp_paint -a 5 /cluster/software/VERSIONS/freesurfer-6.0.0/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh ti. 19. mars 00:57:31 +0100 2019

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 abby lh ../surf/lh.sphere.reg /cluster/software/VERSIONS/freesurfer-6.0.0/average/lh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.annot 

#-----------------------------------------
#@# Cortical Parc rh ti. 19. mars 00:58:00 +0100 2019

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 abby rh ../surf/rh.sphere.reg /cluster/software/VERSIONS/freesurfer-6.0.0/average/rh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf lh ti. 19. mars 00:58:36 +0100 2019

 mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs abby lh 

#--------------------------------------------
#@# Make Pial Surf rh ti. 19. mars 01:13:34 +0100 2019

 mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs abby rh 

#--------------------------------------------
#@# Surf Volume lh ti. 19. mars 01:28:08 +0100 2019
#--------------------------------------------
#@# Surf Volume rh ti. 19. mars 01:29:18 +0100 2019
#--------------------------------------------
#@# Cortical ribbon mask ti. 19. mars 01:29:49 +0100 2019

 mris_volmask --aseg_name aseg.presurf --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon abby 

#-----------------------------------------
#@# Parcellation Stats lh ti. 19. mars 01:46:54 +0100 2019

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab abby lh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.pial.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab abby lh pial 

#-----------------------------------------
#@# Parcellation Stats rh ti. 19. mars 01:49:59 +0100 2019

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab abby rh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.pial.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab abby rh pial 

#-----------------------------------------
#@# Cortical Parc 2 lh ti. 19. mars 01:52:43 +0100 2019

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 abby lh ../surf/lh.sphere.reg /cluster/software/VERSIONS/freesurfer-6.0.0/average/lh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Cortical Parc 2 rh ti. 19. mars 01:53:33 +0100 2019

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 abby rh ../surf/rh.sphere.reg /cluster/software/VERSIONS/freesurfer-6.0.0/average/rh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh ti. 19. mars 01:54:17 +0100 2019

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab abby lh white 

#-----------------------------------------
#@# Parcellation Stats 2 rh ti. 19. mars 01:56:16 +0100 2019

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab abby rh white 

#-----------------------------------------
#@# Cortical Parc 3 lh ti. 19. mars 01:57:31 +0100 2019

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 abby lh ../surf/lh.sphere.reg /cluster/software/VERSIONS/freesurfer-6.0.0/average/lh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Cortical Parc 3 rh ti. 19. mars 01:58:03 +0100 2019

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 abby rh ../surf/rh.sphere.reg /cluster/software/VERSIONS/freesurfer-6.0.0/average/rh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Parcellation Stats 3 lh ti. 19. mars 01:58:30 +0100 2019

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.DKTatlas.stats -b -a ../label/lh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab abby lh white 

#-----------------------------------------
#@# Parcellation Stats 3 rh ti. 19. mars 02:00:09 +0100 2019

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.DKTatlas.stats -b -a ../label/rh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab abby rh white 

#-----------------------------------------
#@# WM/GM Contrast lh ti. 19. mars 02:01:25 +0100 2019

 pctsurfcon --s abby --lh-only 

#-----------------------------------------
#@# WM/GM Contrast rh ti. 19. mars 02:02:58 +0100 2019

 pctsurfcon --s abby --rh-only 

#-----------------------------------------
#@# Relabel Hypointensities ti. 19. mars 02:03:33 +0100 2019

 mri_relabel_hypointensities aseg.presurf.mgz ../surf aseg.presurf.hypos.mgz 

#-----------------------------------------
#@# AParc-to-ASeg aparc ti. 19. mars 02:03:59 +0100 2019

 mri_aparc2aseg --s abby --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt 

#-----------------------------------------
#@# AParc-to-ASeg a2009s ti. 19. mars 02:09:20 +0100 2019

 mri_aparc2aseg --s abby --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt --a2009s 

#-----------------------------------------
#@# AParc-to-ASeg DKTatlas ti. 19. mars 02:15:02 +0100 2019

 mri_aparc2aseg --s abby --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /cluster/software/VERSIONS/freesurfer-6.0.0/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt --annot aparc.DKTatlas --o mri/aparc.DKTatlas+aseg.mgz 

#-----------------------------------------
#@# APas-to-ASeg ti. 19. mars 02:20:24 +0100 2019

 apas2aseg --i aparc+aseg.mgz --o aseg.mgz 

#--------------------------------------------
#@# ASeg Stats ti. 19. mars 02:20:39 +0100 2019

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --brainmask mri/brainmask.mgz --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --euler --ctab /cluster/software/VERSIONS/freesurfer-6.0.0/ASegStatsLUT.txt --subject abby 

#-----------------------------------------
#@# WMParc ti. 19. mars 02:25:10 +0100 2019

 mri_aparc2aseg --s abby --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject abby --surf-wm-vol --ctab /cluster/software/VERSIONS/freesurfer-6.0.0/WMParcStatsLUT.txt --etiv 

INFO: fsaverage subject does not exist in SUBJECTS_DIR
INFO: Creating symlink to fsaverage subject...

 cd /work/jobs/; ln -s /cluster/software/VERSIONS/freesurfer-6.0.0/subjects/fsaverage; cd - 

#--------------------------------------------
#@# BA_exvivo Labels lh ti. 19. mars 02:36:13 +0100 2019

 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA1_exvivo.label --trgsubject abby --trglabel ./lh.BA1_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA2_exvivo.label --trgsubject abby --trglabel ./lh.BA2_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA3a_exvivo.label --trgsubject abby --trglabel ./lh.BA3a_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA3b_exvivo.label --trgsubject abby --trglabel ./lh.BA3b_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA4a_exvivo.label --trgsubject abby --trglabel ./lh.BA4a_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA4p_exvivo.label --trgsubject abby --trglabel ./lh.BA4p_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA6_exvivo.label --trgsubject abby --trglabel ./lh.BA6_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA44_exvivo.label --trgsubject abby --trglabel ./lh.BA44_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA45_exvivo.label --trgsubject abby --trglabel ./lh.BA45_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.V1_exvivo.label --trgsubject abby --trglabel ./lh.V1_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.V2_exvivo.label --trgsubject abby --trglabel ./lh.V2_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.MT_exvivo.label --trgsubject abby --trglabel ./lh.MT_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.entorhinal_exvivo.label --trgsubject abby --trglabel ./lh.entorhinal_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.perirhinal_exvivo.label --trgsubject abby --trglabel ./lh.perirhinal_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA1_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA1_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA2_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA2_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA3a_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA3a_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA3b_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA3b_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA4a_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA4a_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA4p_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA4p_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA6_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA6_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA44_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA44_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.BA45_exvivo.thresh.label --trgsubject abby --trglabel ./lh.BA45_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.V1_exvivo.thresh.label --trgsubject abby --trglabel ./lh.V1_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.V2_exvivo.thresh.label --trgsubject abby --trglabel ./lh.V2_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.MT_exvivo.thresh.label --trgsubject abby --trglabel ./lh.MT_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.entorhinal_exvivo.thresh.label --trgsubject abby --trglabel ./lh.entorhinal_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/lh.perirhinal_exvivo.thresh.label --trgsubject abby --trglabel ./lh.perirhinal_exvivo.thresh.label --hemi lh --regmethod surface 


 mris_label2annot --s abby --hemi lh --ctab /cluster/software/VERSIONS/freesurfer-6.0.0/average/colortable_BA.txt --l lh.BA1_exvivo.label --l lh.BA2_exvivo.label --l lh.BA3a_exvivo.label --l lh.BA3b_exvivo.label --l lh.BA4a_exvivo.label --l lh.BA4p_exvivo.label --l lh.BA6_exvivo.label --l lh.BA44_exvivo.label --l lh.BA45_exvivo.label --l lh.V1_exvivo.label --l lh.V2_exvivo.label --l lh.MT_exvivo.label --l lh.entorhinal_exvivo.label --l lh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose 


 mris_label2annot --s abby --hemi lh --ctab /cluster/software/VERSIONS/freesurfer-6.0.0/average/colortable_BA.txt --l lh.BA1_exvivo.thresh.label --l lh.BA2_exvivo.thresh.label --l lh.BA3a_exvivo.thresh.label --l lh.BA3b_exvivo.thresh.label --l lh.BA4a_exvivo.thresh.label --l lh.BA4p_exvivo.thresh.label --l lh.BA6_exvivo.thresh.label --l lh.BA44_exvivo.thresh.label --l lh.BA45_exvivo.thresh.label --l lh.V1_exvivo.thresh.label --l lh.V2_exvivo.thresh.label --l lh.MT_exvivo.thresh.label --l lh.entorhinal_exvivo.thresh.label --l lh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose 


 mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.stats -b -a ./lh.BA_exvivo.annot -c ./BA_exvivo.ctab abby lh white 


 mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.thresh.stats -b -a ./lh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab abby lh white 

#--------------------------------------------
#@# BA_exvivo Labels rh ti. 19. mars 02:42:56 +0100 2019

 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA1_exvivo.label --trgsubject abby --trglabel ./rh.BA1_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA2_exvivo.label --trgsubject abby --trglabel ./rh.BA2_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA3a_exvivo.label --trgsubject abby --trglabel ./rh.BA3a_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA3b_exvivo.label --trgsubject abby --trglabel ./rh.BA3b_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA4a_exvivo.label --trgsubject abby --trglabel ./rh.BA4a_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA4p_exvivo.label --trgsubject abby --trglabel ./rh.BA4p_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA6_exvivo.label --trgsubject abby --trglabel ./rh.BA6_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA44_exvivo.label --trgsubject abby --trglabel ./rh.BA44_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA45_exvivo.label --trgsubject abby --trglabel ./rh.BA45_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.V1_exvivo.label --trgsubject abby --trglabel ./rh.V1_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.V2_exvivo.label --trgsubject abby --trglabel ./rh.V2_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.MT_exvivo.label --trgsubject abby --trglabel ./rh.MT_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.entorhinal_exvivo.label --trgsubject abby --trglabel ./rh.entorhinal_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.perirhinal_exvivo.label --trgsubject abby --trglabel ./rh.perirhinal_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA1_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA1_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA2_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA2_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA3a_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA3a_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA3b_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA3b_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA4a_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA4a_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA4p_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA4p_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA6_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA6_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA44_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA44_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.BA45_exvivo.thresh.label --trgsubject abby --trglabel ./rh.BA45_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.V1_exvivo.thresh.label --trgsubject abby --trglabel ./rh.V1_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.V2_exvivo.thresh.label --trgsubject abby --trglabel ./rh.V2_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.MT_exvivo.thresh.label --trgsubject abby --trglabel ./rh.MT_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.entorhinal_exvivo.thresh.label --trgsubject abby --trglabel ./rh.entorhinal_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /work/jobs//fsaverage/label/rh.perirhinal_exvivo.thresh.label --trgsubject abby --trglabel ./rh.perirhinal_exvivo.thresh.label --hemi rh --regmethod surface 


 mris_label2annot --s abby --hemi rh --ctab /cluster/software/VERSIONS/freesurfer-6.0.0/average/colortable_BA.txt --l rh.BA1_exvivo.label --l rh.BA2_exvivo.label --l rh.BA3a_exvivo.label --l rh.BA3b_exvivo.label --l rh.BA4a_exvivo.label --l rh.BA4p_exvivo.label --l rh.BA6_exvivo.label --l rh.BA44_exvivo.label --l rh.BA45_exvivo.label --l rh.V1_exvivo.label --l rh.V2_exvivo.label --l rh.MT_exvivo.label --l rh.entorhinal_exvivo.label --l rh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose 


 mris_label2annot --s abby --hemi rh --ctab /cluster/software/VERSIONS/freesurfer-6.0.0/average/colortable_BA.txt --l rh.BA1_exvivo.thresh.label --l rh.BA2_exvivo.thresh.label --l rh.BA3a_exvivo.thresh.label --l rh.BA3b_exvivo.thresh.label --l rh.BA4a_exvivo.thresh.label --l rh.BA4p_exvivo.thresh.label --l rh.BA6_exvivo.thresh.label --l rh.BA44_exvivo.thresh.label --l rh.BA45_exvivo.thresh.label --l rh.V1_exvivo.thresh.label --l rh.V2_exvivo.thresh.label --l rh.MT_exvivo.thresh.label --l rh.entorhinal_exvivo.thresh.label --l rh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose 


 mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.stats -b -a ./rh.BA_exvivo.annot -c ./BA_exvivo.ctab abby rh white 


 mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.thresh.stats -b -a ./rh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab abby rh white 

