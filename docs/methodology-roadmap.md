# Methodology roadmap

This is a living map of methods that may improve interpretation. Inclusion does not mean a method is
appropriate for every data source or is already implemented.

## Reference-population inference

- Overlap and positivity diagnostics before matching or weighting
- Exact, propensity, entropy-balancing, and calibration-weighting approaches
- Doubly robust estimation when outcome and selection models are defensible
- Negative controls, quantitative bias analysis, and partial identification for unmeasured confounding
- Transportability methods when the target cancer population differs from the available reference

## Expression and composition

- Donor-aware pseudobulk models for single-cell comparisons
- Hierarchical models for repeated tissues or samples per donor
- Bulk deconvolution and composition-aware sensitivity analysis
- Latent-factor and unwanted-variation methods, with diagnostics for removed biological signal
- Cross-platform calibration separated from biological population adjustment

## State, trajectory, and lineage

Snapshot-derived state geometry, inferred developmental trajectory, somatic phylogeny, and
experimental lineage tracing are different evidence classes. The software should label them
explicitly rather than use “lineage” as a generic synonym for similarity.

[Bonsai](https://www.nature.com/articles/s41587-026-03220-2) reconstructs a tree representation that
preserves structure in noisy high-dimensional single-cell measurements. It can support exploratory
state and differentiation hypotheses, but an expression-derived Bonsai tree is not direct lineage
observation. Integration of lineage-tracing mutations or barcodes is described as an extension.

Future lineage-aware work should distinguish:

- Evolving or static experimental barcodes
- Naturally occurring somatic or mitochondrial mutations
- Multi-region tumor phylogenies
- Longitudinal sampling and metastatic migration histories
- Transcriptomic trajectories inferred without lineage observations

## Evaluation standard

Every proposed method should be reviewed for its estimand, required design, independence unit,
missing-data assumptions, uncertainty model, failure modes, computational maturity, and compatibility
with public datasets. New methods should enter production only with simulations or benchmark data
that demonstrate both expected behavior and known failure cases.
