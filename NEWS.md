# SpaMTP 0.99.1

* Standardised all exported function names to lower camel case for the
  Bioconductor API. This is an intentional breaking change on the submission
  branch; the published-work development branch retains its historical API.
* Removed direct S4 slot access from package and test code. Internal container
  access now uses public SeuratObject, SummarizedExperiment, S4Vectors,
  Cardinal, and MSnbase accessors. Annotation and alignment provenance is
  stored with `SeuratObject::Misc()` while legacy `Tool()` entries remain
  readable.
* Added an internal container-access layer for Seurat and
  SummarizedExperiment-derived objects, including SingleCellExperiment and
  SpatialExperiment, as the foundation for broader Bioconductor container
  interoperability.
* Removed the unsupported manual mutation of Seurat modality-weight internals
  from `multiOmicIntegration()`; weighted-nearest-neighbour construction now
  delegates entirely to `Seurat::FindMultiModalNeighbors()`.

# SpaMTP 0.99.0

* Split versioned resources from the analysis code: `SpaMTPdb` now supplies
  the pruned RaMP annotation and pathway snapshot, and `SpaMTPData` provides
  named access to large experiment/vignette objects. `loadSpaMTPDatabase()`
  and `spaMTPDatabaseInfo()` form the public database interface. The software
  package retains only the small adduct and reaction-style constant tables.
* Updated the package citation to the peer-reviewed *Nature Methods* paper,
  designated Tianyao Lu as the current package maintainer, and added direct
  links to the developmental documentation and source branch. Andrew Causer
  remains credited as an original author and former maintainer.
* Added `applySpatialAlignment()` to apply existing SMINT coordinate outputs or
  homogeneous transforms, estimate landmark-based affine alignment in R, and
  run the SMINT-compatible STalign LDDMM workflow through an optional Python
  backend. Alignment provenance and nearest-target diagnostics are retained in
  the returned SpaMTP object without storing large Python velocity tensors.
* Pathway enrichment, pathway-assay construction, Fisher analysis of m/z
  inputs, and interactive pathway networks now consume the scored `Ramp_IDs`
  produced by the indexed annotation engine. `annotateSM()` resolves the
  versioned RaMP 3.0.7 chemical-property table through `SpaMTPdb`, records
  annotation provenance through `SeuratObject::Misc()`, with a compatibility
  fallback for legacy serialized objects.
  Legacy annotations require an explicit `annotation_source` fallback.
* Annotation candidates can now be stored once with `annotateSM(min_score =
  0)` and filtered later with `annotation_score_threshold`. Interactive
  pathway networks add `metabolite_detection = "annotated"` to display
  score-filtered pathway metabolites without requiring DE significance, while
  retaining separate leading-edge evidence in node styling and tooltips.
* Moved the pruned RaMP snapshot from 2.5.4 to versioned SpaMTPdb 3.0.7
  resources, with explicit upstream/source metadata and reproducible resource
  staging scripts.
* Added a chemically validated adduct rule table and reusable sorted m/z index.
* Added ppm pruning, proton/charge bounds, mass-defect scoring, isotope checks,
  and contextual adduct-family scoring while retaining the legacy annotation
  result columns.
* Added dependency-free SMILES graph decomposition, functional-group and atom
  site reporting, mode-specific protonation/deprotonation and alkali-binding
  priors, and `predictAdductsFromSMILES()`. Aromatic functional groups are
  perceived consistently from both lowercase aromatic and alternating-bond
  six-membered Kekule SMILES. Full RaMP runs join an independent precomputed
  `SpaMTPdb::smiles_features` resource; small custom databases are inferred
  automatically at runtime.
* `chem_props` can now be used directly as a RaMP-backed annotation database.
* Rebuilt `pathwayNetworkPlots()` around cached topology lookup, precomputed
  edges, selective sparse-matrix extraction, and JSON serialization. The new
  responsive viewer adds zoom/pan, focused/full networks, label controls,
  shared reaction markers, spatial inspection, and SVG export.
* Made legacy pathway annotation robust to unequal `Isomers_IDs` and
  `IsomerNames` counts, and added a cached Mouse Brain/DHB network demo.
* Pathway network metabolite labels now remain visible in the default label
  mode; dragged nodes stay pinned until released, and edge legends persist
  when the focused view changes.
* Added interactive repulsion, balanced-force, radial, and gene/metabolite
  pathway layouts, with a more widely spaced repulsion layout as the default.
* Harmonised all versioned pathway-graph node IDs with RaMP-DB 3.0.7, added a
  reproducible cross-version graph updater/auditor, and added stable source-ID
  labels when RaMP itself has no common name.

# SpaMTP 1.1.0

***Updated SpaMTP Release (Oct 2025)***

Additional Features:

* Updated `loadSM()` function for compatability with **Cardinal V3.8**
* Implementation of GraphPCA in R.
* Additional functions for handling large datasets - `annotateBigData` and `selectROIs`.
* Function for refining m/z annotations based on correlated pathway activity or MS/MS profiles.
* Package wide update for compatability with **Cardinal V3.8** and **Seurat V5.3**.


# SpaMTP 1.0.0

* Initial SpaMTP Release (March 2025)
