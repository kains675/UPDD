# Apex Bridge P0v2 Reference Summary

Generated: 2026-07-17 08:16:33 +0900

Status: **REJECT_BRIDGE**

P0v2 used OpenMM Reference with zero MD/minimization steps and no GPU. It did not estimate a free energy.

Inventory digest: `8a5c2eabd6125bf70a8dbf6b334d29ded9d00091c6aab597e0c4af1c22cd3615`

| cell | status | endpoint gap (kJ/mol) | max force gap (kJ/mol/nm) | elapsed (s) |
|---|---|---:|---:|---:|
| w4a_uncarved_s101_bound | P0V2_PASS | 1599.915304657537 | 84923.013507088210 | 255.162 |
| w4a_uncarved_s127_bound | P0V2_PASS | 2042.957528679632 | 132518.308204555040 | 255.843 |
| w4a_uncarved_s163_bound | REJECT_BRIDGE | 1331583646.915235996246 | 276209950501.097167968750 | 268.418 |

Pending: `w4a_uncarved_s101_free`, `w4a_uncarved_s127_free`, `w4a_uncarved_s163_free`

Rejected: `w4a_uncarved_s163_bound`
