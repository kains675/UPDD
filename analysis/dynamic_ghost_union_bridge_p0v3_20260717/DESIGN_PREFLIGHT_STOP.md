# P0v3 Design Preflight Stop

Date: 2026-07-17 KST

Status: **SUPERSEDED BEFORE OFFICIAL SOURCE OR P0v3 OUTPUT**

The frozen P0v3 design correctly preserved the s101-free target counts during
an in-memory development smoke:

```text
final particles = 29992
water            = 9846
Na               = 27
Cl               = 27
total solvent    = 9900
placeholders     = 9 added, 9 removed
target drift     = 0.0 nm
```

However, `Modeller.addSolvent(numAdded=9900)` produced a cubic edge of
`6.8498 nm`, compared with `6.7890 nm` in the padding-built parent. A gate in
the frozen design compared total solvent molecules per total box volume and
required a relative difference within 1%. The observed total-box ratio changed
by 2.64%.

That metric is not a valid static bulk-density gate. Its denominator contains
the solute and the deliberately larger union excluded volume, so the same
solvent count can yield a different total-box `N/V` even when OpenMM places the
bulk water from the same pre-equilibrated template. Equilibrium density also
cannot be certified without MD/NPT, which this CPU-only stage forbids.

No official union source XML/PDB, inventory, KEEPER artifact, Reference cell,
MD step, minimization, GPU run, or free-energy estimate was produced under this
freeze. The preregistration and protocol remain unchanged as rejected design
evidence.

A separate revision must retain exact water/ion/particle counts, placeholder
removal, and both-image geometry gates, but report box volume without treating
static total-box `N/V` as equilibrium density. Density belongs to a future
explicit NPT equilibration gate after the revised CPU mechanics audit passes.
