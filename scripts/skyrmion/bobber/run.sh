#!/bin/bash -x
r_ratio_to_sample="16"
t_relax_state1_sec="3e-9"

magnum.af -g 1 -b opencl multilayer.cpp runs/string_to_ferromag/r4/case2/ 2 4 5e-9
magnum.af -g 2 -b opencl multilayer.cpp runs/string_to_ferromag/r8/case2/ 2 8 5e-9
magnum.af -g 3 -b opencl multilayer.cpp runs/string_to_ferromag/r16/case2/ 2 16 5e-9


#../../magnum.af -g 1 -b opencl multilayer.cpp runs/string_to_ferromag/r16/case1/ 1 16 3e-9
#../../magnum.af -g 2 -b opencl multilayer.cpp runs/string_to_ferromag/r16/case2/ 2 16 3e-9
#../../magnum.af -g 3 -b opencl multilayer.cpp runs/string_to_ferromag/r16/case3/ 3 16 3e-9

#../../magnum.af -S -g 1 -b opencl -o runs/string_to_ferromag/r16/case1/ multilayer.cpp 1 16 3e-9
#../../magnum.af -S -g 2 -b opencl -o runs/string_to_ferromag/r16/case2/ multilayer.cpp 2 16 3e-9
#../../magnum.af -S -g 3 -b opencl -o runs/string_to_ferromag/r16/case3/ multilayer.cpp 3 16 3e-9
#../../magnum.af -S -g 1 -b opencl multilayer.cpp runs/string_to_ferromag/r16/case1 1 "$r_ratio_to_sample" "$t_relax_state1_sec"
#../../magnum.af -S -g 2 -b opencl multilayer.cpp runs/string_to_ferromag/r16/case2 2 "$r_ratio_to_sample" "$t_relax_state1_sec"
#../../magnum.af -S -g 3 -b opencl multilayer.cpp runs/string_to_ferromag/r16/case3 3 "$r_ratio_to_sample" "$t_relax_state1_sec"
#../../magnum.af -S -g 1 -b opencl multilayer_case1.cpp runs/r16/case1
#../../magnum.af -S -g 2 -b opencl multilayer_case2.cpp runs/r16/case2
#../../magnum.af -S -g 3 -b opencl multilayer_case3.cpp runs/r16/case3

#string_case2_to_case3
magnum.af -f -g 1 -b opencl multilayer.cpp runs/string_case2_to_case3/ 2 16 5e-9
magnum.af -f -g 2 -b opencl multilayer.cpp runs/test_king_barrier_case2_image_51_53/ 2 16 5e-9

magnum.af -f -g 2 -b opencl multilayer.cpp runs/run1_case2_to_case3_image_55_59/ 2 16 5e-9

#magnum.af -f -g 3 -b opencl multilayer_case2_to_case3.cpp runs/v2_case3_to_case2/r16/nxy100 -1 16 5e-9
#magnum.af -f -g 2 -b opencl multilayer_case2_to_case3.cpp runs/v2_case3_to_case2/r16/nxy80 -1 16 5e-9
magnum.af -f -g 1 -b opencl multilayer_case2_to_case3.cpp runs/v2_case3_to_case2/r16/nxy120 -1 16 5e-9

magnum.af -f -g 2 -b opencl multilayer_case2_to_case3.cpp runs/v3_case3_to_case2/r16/nxy100 -1 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_case2_to_case3.cpp runs/v3_case3_to_case2/r16/nxy80 -1 16 5e-9 80
magnum.af -f -g 1 -b opencl multilayer_case2_to_case3.cpp runs/v3_case3_to_case2/r16/nxy120 -1 16 5e-9 120

magnum.af -f -g 0 -b opencl multilayer_case2_to_case3.cpp runs/v3_case3_to_case2/r16/nxy90 -1 16 5e-9 90
magnum.af -f -g 1 -b opencl multilayer_case2_to_case3.cpp runs/v3_case3_to_case2/r16/nxy110 -1 16 5e-9 110
magnum.af -f -g 2 -b opencl multilayer_case2_to_case3.cpp runs/v3_case3_to_case2/r16/nxy96 -1 16 5e-9 96
magnum.af -f -g 3 -b opencl multilayer_case2_to_case3.cpp runs/v3_case3_to_case2/r16/nxy104 -1 16 5e-9 104

magnum.af -g 0 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy90  2 16 5e-9 90
magnum.af -g 1 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy100 2 16 5e-9 100
magnum.af -g 2 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy110 2 16 5e-9 110
magnum.af -g 3 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy120 2 16 5e-9 120

magnum.af -g 0 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy70  2 16 5e-9 70
magnum.af -g 1 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy80  2 16 5e-9 80
magnum.af -g 2 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy130 2 16 5e-9 130
magnum.af -g 3 -b opencl multilayer.cpp runs/v3_case2_to_ferro/r16/nxy140 2 16 5e-9 140

magnum.af -g 0 -b opencl multilayer.cpp runs/v4_to_ferro_r16_nxy100/case1/ 1 16 5e-9 100
magnum.af -g 1 -b opencl multilayer.cpp runs/v4_to_ferro_r16_nxy100/case2/ 2 16 5e-9 100
magnum.af -g 2 -b opencl multilayer.cpp runs/v4_to_ferro_r16_nxy100/case3/ 3 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer.cpp runs/v5_to_ferro_r16_nxy100/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer.cpp runs/v5_to_ferro_r16_nxy100/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer.cpp runs/v5_to_ferro_r16_nxy100/case3/ 3 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer.cpp runs/v6_to_ferro_r16_nxy100/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer.cpp runs/v6_to_ferro_r16_nxy100/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer.cpp runs/v6_to_ferro_r16_nxy100/case3/ 3 16 5e-9 100

magnum.af -f -g 0 -b opencl multilayer.cpp runs/todel_memtest/case2/ 2 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer_resume.cpp runs/v6.1_to_ferro_r16_nxy100/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer_resume.cpp runs/v6.1_to_ferro_r16_nxy100/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_resume.cpp runs/v6.1_to_ferro_r16_nxy100/case3/ 3 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer_resume.cpp runs/v7_to_ferro_r16_nxy100_dt1em14/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer_resume.cpp runs/v7_to_ferro_r16_nxy100_dt1em14/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_resume.cpp runs/v7_to_ferro_r16_nxy100_dt1em14/case3/ 3 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer.cpp runs/v8_to_ferro_r16_nxy100_dt1em13_images120/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer.cpp runs/v8_to_ferro_r16_nxy100_dt1em13_images120/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer.cpp runs/v8_to_ferro_r16_nxy100_dt1em13_images120/case3/ 3 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer.cpp runs/v9_dmi_il_pos_to_ferro_r16_nxy100_dt1em13_images60/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer.cpp runs/v9_dmi_il_pos_to_ferro_r16_nxy100_dt1em13_images60/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer.cpp runs/v9_dmi_il_pos_to_ferro_r16_nxy100_dt1em13_images60/case3/ 3 16 5e-9 100
#magnum.af -f -g 4 -b opencl multilayer.cpp runs/v9_dmi_il_pos_to_ferro_r16_nxy100_dt1em13_images60_global_neg_dmi/case1/ 1 16 5e-9 100

magnum.af -f -g 0 -b opencl hys_multilayer.cpp runs/hys_v1/case1/ 1 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer.cpp runs/v9_dmi_il_pos_to_ferro_r16_nxy100_dt1em13_images60/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer.cpp runs/v9_dmi_il_pos_to_ferro_r16_nxy100_dt1em13_images60/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer.cpp runs/v9_dmi_il_pos_to_ferro_r16_nxy100_dt1em13_images60/case3/ 3 16 5e-9 100

magnum.af -f -g 0 -b opencl hys_multilayer.cpp runs/hys_v1/case1/ 1 16 5e-9 100
magnum.af -f -g 0 -b opencl hys_multilayer.cpp runs/hys_v1/case2/ 2 16 5e-9 100
magnum.af -f -g 0 -b opencl hys_multilayer.cpp runs/hys_v1/case3/ 3 16 5e-9 100

#sanity#magnum.af -f -g 0 -b opencl hys_multilayer.cpp runs/hys_v1/neg_ILdmi_case1/ 1 16 5e-9 100
magnum.af -f -g 0 -b opencl hys_multilayer.cpp runs/hys_v1/neg_ILdmi_case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl hys_multilayer.cpp runs/hys_v1/neg_ILdmi_case3/ 3 16 5e-9 100

magnum.af -g 1 -b opencl hys_multilayer.cpp runs/hys_v1/global_neg_dmi_case1/ 1 16 5e-9 100
magnum.af -g 2 -b opencl hys_multilayer.cpp runs/hys_v1/global_neg_dmi_case2/ 2 16 5e-9 100
magnum.af -g 3 -b opencl hys_multilayer.cpp runs/hys_v1/global_neg_dmi_case3/ 3 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run1/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run1/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run1/case3/ 3 16 5e-9 100

magnum.af -f -g 2 -b opencl multilayer_demag_above_top_layer_test_stability_of_switched_bottom_top_layers.cpp runs/demag_above_top_layer_test_stability_of_switched_bottom_top_layers/run1/case2/ 2 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run2/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run2/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run2/case3/ 3 16 5e-9 100

magnum.af -f -g 1 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run3/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run3/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run3/case3/ 3 16 5e-9 100

magnum.af -f -g 0 -b opencl multilayer_demag_above_top_layer_test_stability_of_switched_bottom_top_layers.cpp runs/demag_above_top_layer_test_stability_of_switched_bottom_top_layers/run2/case2/ 2 16 5e-9 100

#magnum.af -f -g 1 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run4_all_dmi_switched/case1/ 1 16 5e-9 100
magnum.af -f -g 2 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run4_all_dmi_switched/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run4_all_dmi_switched/case3/ 3 16 5e-9 100

#TODO
magnum.af -f -g 0 -b opencl multilayer.cpp runs/string/run1/case2/ 2 16 5e-9 100
magnum.af -f -g 1 -b opencl multilayer.cpp runs/string/run1/case3/ 3 16 5e-9 100

magnum.af -f -g 2 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run4/dmivec_switched/case2/ 2 16 5e-9 100
magnum.af -f -g 3 -b opencl multilayer_demag_above_top_layer.cpp runs/demag_above_top_layer/run4/dmivec_switched/case3/ 3 16 5e-9 100

magnum.af -f -g 0 -b opencl hys_multilayer.cpp runs/hys_v2/case2/ 2 16 5e-9 100
magnum.af -f -g 1 -b opencl hys_multilayer.cpp runs/hys_v2/case3/ 3 16 5e-9 100
