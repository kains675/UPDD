# Tail-Mechanism Probe Report

Date: 2026-07-10

Read-only analysis over existing Track B `.out`, manifest, and JSON artifacts.

## Cohort Diagnostics

| cohort | highlight | top \|ddG\| seed | \|ddG\| | top p99 seed | p99 | top frac>150 seed | frac>150 | min crossing seed | min crossing |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| w23a_gateA | s101 | s101 | 16.005 | s101 | 190.829 | s7 | 0.157 | s101 | 150 |
| w4a_c1_carved | s127 | s127 | 12.563 | s127 | 191.301 | s101 | 0.138 | s23 | 14 |
| w4a_c1_uncarved | s127,s199 | s199 | 4.408 | s199 | 191.024 | s199 | 0.139 | s163 | 10 |

## w23a_gateA

| seed | flag | ddG | bound dgb | free dgb | max p99 cell | p99 | max frac cell | frac>150 | min crossing cell | min crossing | RT sum | gate fail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| s7 |  | -7.806 | -7.922 | -0.116 | free/dplus:state0 | 190.065 | free/dplus:state0 | 0.157 | bound/dminus:state14 | 254 | 5 | False |
| s19 |  | -8.736 | -9.219 | -0.483 | free/dminus:state14 | 190.367 | free/dminus:state14 | 0.157 | free/dplus:state0 | 228 | 7 | False |
| s23 |  | -8.881 | -7.672 | 1.209 | free/dplus:state0 | 188.883 | free/dplus:state0 | 0.140 | bound/dplus:state0 | 284 | 4 | False |
| s42 |  | -7.611 | -9.589 | -1.978 | bound/dminus:state14 | 189.150 | free/dplus:state0 | 0.131 | free/dplus:state0 | 242 | 8 | False |
| s101 | * | -16.005 | -13.667 | 2.338 | bound/dplus:state0 | 190.829 | bound/dminus:state14 | 0.140 | bound/dminus:state14 | 150 | 4 | False |

Correlation diagnostics:

| metric | value |
| --- | --- |
| Pearson \|ddG\| vs max p99 | 0.652 |
| Pearson \|ddG\| vs max frac>150 | -0.186 |
| Pearson \|ddG\| vs min crossing | -0.880 |

Highlight minus non-highlight median:

| seed | abs ddG | p99 | frac>150 | min crossing |
| --- | --- | --- | --- | --- |
| s101 | 7.734 | 1.221 | -0.008 | -98.000 |

## w4a_c1_carved

| seed | flag | ddG | bound dgb | free dgb | max p99 cell | p99 | max frac cell | frac>150 | min crossing cell | min crossing | RT sum | gate fail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| s7 |  | -5.278 | -4.010 | 1.269 | bound/dminus:state11 | 189.947 | bound/dminus:state11 | 0.129 | free/dminus:state11 | 26 | 0 | True |
| s23 |  | -3.271 | -1.736 | 1.535 | bound/dplus:state0 | 187.655 | free/dplus:state0 | 0.116 | bound/dplus:state0 | 14 | 1 | True |
| s101 |  | -3.798 | -1.821 | 1.977 | free/dplus:state0 | 189.678 | free/dplus:state0 | 0.138 | bound/dminus:state11 | 18 | 1 | True |
| s127 | * | -12.563 | -6.667 | 5.895 | free/dplus:state0 | 191.301 | bound/dminus:state11 | 0.122 | free/dminus:state11 | 42 | 0 | True |
| s163 |  | 1.863 | 0.468 | -1.395 | free/dminus:state11 | 188.619 | bound/dminus:state11 | 0.121 | free/dplus:state0 | 18 | 3 | False |
| s199 |  | -3.289 | -0.770 | 2.520 | free/dminus:state11 | 187.973 | free/dplus:state0 | 0.115 | bound/dminus:state11 | 62 | 1 | False |

Correlation diagnostics:

| metric | value |
| --- | --- |
| Pearson \|ddG\| vs max p99 | 0.834 |
| Pearson \|ddG\| vs max frac>150 | 0.052 |
| Pearson \|ddG\| vs min crossing | 0.320 |

Highlight minus non-highlight median:

| seed | abs ddG | p99 | frac>150 | min crossing |
| --- | --- | --- | --- | --- |
| s127 | 9.273 | 2.681 | 0.001 | 24.000 |

## w4a_c1_uncarved

| seed | flag | ddG | bound dgb | free dgb | max p99 cell | p99 | max frac cell | frac>150 | min crossing cell | min crossing | RT sum | gate fail |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| s7 |  | -0.912 | 1.183 | 2.094 | bound/dplus:state0 | 188.491 | bound/dplus:state0 | 0.130 | free/dminus:state11 | 28 | 1 | True |
| s23 |  | -2.141 | 0.169 | 2.310 | bound/dplus:state0 | 188.330 | free/dplus:state0 | 0.117 | free/dminus:state11 | 58 | 3 | False |
| s101 |  | 0.679 | 0.150 | -0.530 | bound/dplus:state0 | 187.677 | bound/dminus:state11 | 0.107 | bound/dplus:state0 | 42 | 1 | False |
| s127 | * | 3.655 | -2.004 | -5.659 | bound/dplus:state0 | 189.510 | bound/dplus:state0 | 0.123 | bound/dminus:state11 | 30 | 2 | True |
| s163 |  | 1.214 | -0.845 | -2.059 | free/dminus:state11 | 189.216 | bound/dminus:state11 | 0.136 | bound/dplus:state0 | 10 | 0 | True |
| s199 | * | -4.408 | -2.869 | 1.539 | bound/dplus:state0 | 191.024 | free/dplus:state0 | 0.139 | free/dminus:state11 | 24 | 0 | True |

Correlation diagnostics:

| metric | value |
| --- | --- |
| Pearson \|ddG\| vs max p99 | 0.853 |
| Pearson \|ddG\| vs max frac>150 | 0.424 |
| Pearson \|ddG\| vs min crossing | -0.083 |

Highlight minus non-highlight median:

| seed | abs ddG | p99 | frac>150 | min crossing |
| --- | --- | --- | --- | --- |
| s127 | 2.593 | 1.099 | -0.000 | -5.000 |
| s199 | 3.345 | 2.613 | 0.015 | -11.000 |

## Top State-Localized pertE Tails

| cohort | seed | leg | dir | state | p99 | frac>150 | mean | n |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| w23a_gateA | s42 | bound | dminus | 14 | 198.003 | 0.917 | 171.172 | 800 |
| w23a_gateA | s23 | free | dplus | 0 | 197.589 | 0.973 | 173.736 | 800 |
| w23a_gateA | s101 | free | dplus | 0 | 197.544 | 0.964 | 173.665 | 800 |
| w4a_c1_uncarved | s163 | bound | dminus | 11 | 197.452 | 0.935 | 171.546 | 400 |
| w23a_gateA | s101 | bound | dplus | 0 | 197.385 | 0.930 | 177.145 | 800 |
| w4a_c1_uncarved | s7 | free | dminus | 11 | 197.324 | 0.848 | 167.436 | 400 |
| w4a_c1_carved | s127 | bound | dminus | 11 | 197.124 | 0.960 | 175.040 | 400 |
| w23a_gateA | s19 | free | dminus | 14 | 197.047 | 0.949 | 174.763 | 800 |
| w4a_c1_carved | s127 | free | dplus | 0 | 197.017 | 0.900 | 174.781 | 400 |
| w4a_c1_uncarved | s199 | free | dminus | 11 | 196.932 | 0.932 | 170.514 | 400 |
| w23a_gateA | s101 | bound | dminus | 14 | 196.884 | 0.929 | 171.562 | 800 |
| w4a_c1_uncarved | s127 | bound | dminus | 11 | 196.877 | 0.953 | 173.589 | 400 |
| w23a_gateA | s7 | bound | dminus | 14 | 196.830 | 0.904 | 169.690 | 800 |
| w4a_c1_carved | s23 | free | dminus | 11 | 196.801 | 0.927 | 170.934 | 400 |
| w4a_c1_uncarved | s23 | bound | dplus | 0 | 196.781 | 0.953 | 173.158 | 400 |
| w23a_gateA | s7 | free | dminus | 14 | 196.770 | 0.955 | 174.890 | 800 |
| w4a_c1_carved | s127 | free | dminus | 11 | 196.709 | 0.915 | 170.687 | 400 |
| w4a_c1_carved | s101 | free | dplus | 0 | 196.678 | 0.993 | 174.707 | 400 |
| w4a_c1_uncarved | s127 | bound | dplus | 0 | 196.650 | 0.985 | 174.918 | 400 |
| w23a_gateA | s23 | bound | dminus | 14 | 196.640 | 0.919 | 170.743 | 800 |
| w4a_c1_carved | s163 | bound | dminus | 11 | 196.634 | 0.912 | 171.135 | 400 |
| w4a_c1_uncarved | s163 | bound | dplus | 0 | 196.469 | 0.985 | 173.546 | 400 |
| w23a_gateA | s19 | bound | dminus | 14 | 196.452 | 0.944 | 173.536 | 800 |
| w4a_c1_uncarved | s199 | bound | dplus | 0 | 196.424 | 0.958 | 175.763 | 400 |
| w4a_c1_carved | s163 | free | dminus | 11 | 196.382 | 0.988 | 175.128 | 400 |
| w23a_gateA | s101 | free | dminus | 14 | 196.355 | 0.912 | 169.514 | 800 |
| w4a_c1_carved | s7 | free | dminus | 11 | 196.320 | 0.910 | 170.911 | 400 |
| w23a_gateA | s19 | free | dplus | 0 | 196.275 | 0.953 | 172.452 | 800 |
| w23a_gateA | s42 | free | dminus | 14 | 196.270 | 0.940 | 170.959 | 800 |
| w4a_c1_uncarved | s23 | free | dminus | 11 | 196.214 | 0.897 | 169.409 | 400 |
