# Guided mode selection at a recorded TS guess

```bash
python examples/ts_mode_selection/replay.py
```

This self-contained example uses the public `mappa` imports to load recorded R/P endpoints, a saved endpoint mapping and bond events, and one holdout TS guess with its Hessian modes. It performs the two endpoint-to-guess core searches and selects and scores the mode with `analyze_transition_state`; no external dataset or electronic-structure program is needed.

Expected selection: mode index **0** (one-based mode 1), **1129.9292i cm⁻¹**, overlap **0.869866**, progress **0.535756**, score **0.466036**. The same O/H/N target assignment is supported by both endpoint searches. The input includes the other recorded modes, not only the selected displacement.

[Watch the film and inspect the evidence](../../manuscript/animations/ts_mode_selection/README.md). The stable comparison mode is illustrative and is not eligible for imaginary-mode selection. The geometry remains a TS guess, not an optimized TS.
