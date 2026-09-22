# SpaMTP 0.99.9

* Clarified native data access across function help pages: modalities use
  altExp(), expression matrices use assay(), and retained slot arguments
  select matrix names rather than S4 storage slots. Corrected return types,
  reduction locations and normalization descriptions. Added a runnable
  experimentAccess guide. Assay selectors now require exact character names;
  modality merging rejects duplicate primary aliases and prefers native
  scaled values over the legacy scale.data matrix.

* All vignette R chunks now execute by default, including full published
  human-brain, FMP10/Visium and DHB/Visium analyses. Resource versions are pinned;
  downloads and long builds are expected when verified caches are absent.
  Original Space Ranger scale factors supply the registered FMP10 pair's spot
  geometry. Synthetic file-I/O and replicate-design examples remain labelled.
  Installed mouse-brain recipes can resolve public resources when no local
  directory is supplied, while explicit directories retain offline defaults.

* Annotation statistics now select the requested MSI modality's candidate store
  for both single-feature and batch ranking, without borrowing another modality's
  latest root-level annotations. Historical and assay-local stores remain usable.
* Regional marker scores, spatial-block sensitivity and biological-replicate
  contrasts preserve matrix dimensions for single-feature inputs.

* Replaced the installed mouse-brain demo and preparation recipes with native
  SpatialExperiment workflows. Installed examples no longer load a legacy
  checkout, read S4 slots directly or require Seurat. Compound networks reuse
  scored annotations, scran-based regional effects and a shared pathway index.
  Optional optical images and spot geometry are supplied explicitly; FMP10
  mapping refuses a missing radius rather than guessing it from pixel spacing.
  Existing output paths are preserved and script sourcing has no side effects.

* Fixed pathway interaction corruption in the graphite-derived RaMP 3.0.7
  graphs. A pinned, checksum-guarded correction restores per-edge source labels,
  directions, and parallel interactions for KEGG, Reactome, WikiPathways and
  SMPDB. Interaction styles now use label semantics across databases, and
  undirected edges have no arrowhead. The original published resource bytes
  retain their checksums; SpaMTP applies corrections when loading them.
* Added the Pathway Database Integration vignette, a reproducible correction
  builder, source checksums, and regression tests for the Cell Cycle report.


* Reorganized the workflow around measurement quality, correspondence,
  representation comparisons, regional markers, replicate inference,
  conditional cross-omic association and shared-identity pathway analysis.
  Reports connect native heatmaps, H&E feature maps, effect plots, 3D views
  and pathway networks with their actual analytical inputs.
* Native PCA accepts explicit features and records the fitted feature set.
  Graph PCA uses a sparse Laplacian and all retained observations; its
  embeddings now feed spatial joint integration and optional clustering.
  Comparisons use matched feature selection/scaling, multiple requested K,
  external-reference ARI, approximate silhouettes and block-omission
  clustering sensitivity conditional on the fitted embeddings.
* findAllDEMs() supports descriptive scran marker effects/AUC and explicit
  biological-replicate limma contrasts. Spatial-block marker sensitivity
  carries no invented P values. The historical technical-pool mode warns
  that its pools are not biological replicates. Native heatmaps preserve
  the actual analysis matrix and rank numeric FDR correctly.
* findCorrelatedFeatures() supports generic feature IDs, explicit target
  sets and raw versus covariate-residual correlations. Exact RNA/protein
  features no longer require m/z metadata for native plots.
* findRegionalPathways() accepts complete unfiltered identity-level ranks.
  Competitive enrichment, scores, GESECA and network views reuse the same
  membership index. Compound identity assays support sparse inputs and
  audit ambiguity exclusions and exact aggregation weights. Annotation
  stores are preserved separately for multiple MS modalities.
* plotRegionalPathways() displays BH-adjusted significance and orders
  pathways using 1 minus Jaccard similarity, including single-pathway input.
  The region explorer adds directional AUC and marker-stability filters.
* Reports audit pathways with identical measured members across regions without
  shrinking the BH family. Native networks embed a pinned, licensed D3
  distribution for offline use, retain missing effects as unavailable, and
  label descriptive effects in the declared workflow units.
* Added portable FMP10/Visium and DHB/Visium recipes using pinned
  SpaMTPData/SpaMTPdb resources, plus a rewritten analytical tutorial with
  three synthetic modalities. Mouse genes are not mapped to human HGNC
  by changing case. Equal-weight joint embeddings remain distinct from
  the historical Seurat WNN method.

# SpaMTP 0.99.8

* Added analyzeSpaMTPRegions() and a portable offline region/DE explorer.
  Modality, region and replicate contrast selectors link feature tables,
  effect/volcano plots, spatial maps, regional means and distributions.
  Search, effect/direction/FDR filters, point inspection, zoom/pan and filtered
  CSV export operate on an explicitly labelled bounded preview; full results
  are exported separately. Blank region labels are excluded and audited.
* Reports now show native plotSpatialFeature(), mzViolinPlot(),
  findSpatiallyVariableMetabolites() and runSpatialGraphPCA() outputs, plus
  regional pathway score heatmaps when a compatible index is configured.
  Native spatial analyses record selected features and positions per sample.
  Descriptive region effects never receive pixel-based DE P values.
* Method descriptions use smaller serif text and collapsible detail panels.
  Capability cards identify computed and unavailable native modules.
  Existing workflow results can be extended without changing their original
  contrasts; core, extension and renderer versions remain distinguishable.
* Expanded the generic workflow tutorial with paired mouse brain, independent
  single-omic inputs, arbitrary region metadata, preview controls and native
  analysis settings. No study-specific feature names or species conversion
  are embedded in the workflow or report renderer.
* Single-modality workflows now apply requested clustering to their PCA and
  expose those exploratory labels to the region browser. Previously, clusters
  were only computed after joint integration and a single-omic request was
  silently ignored. Cluster input and settings are recorded.

# SpaMTP 0.99.7

* Added runSpaMTPWorkflow() for explicit paired or independent spatial
  modalities: pinned companion data acquisition, QC, registration/mapping,
  normalization, per-modality PCA, joint representation, exploratory
  associations, biological-replicate contrasts and optional annotation/pathways.
  Unmatched positions are audited and excluded from joint analysis.
* Added renderSpaMTPReport() with embedded figures, full CSV tables, native
  result RDS, registration weights, resource provenance and file checksums.
  Rendering needs neither Pandoc nor network access and preserves existing files.
* Added an offline end-to-end vignette demonstrating two/three modalities,
  known-transform registration and portable reporting, with explicit boundaries
  for raw preprocessing, biological inference and platform-specific settings.

# SpaMTP 0.99.6

* Added buildPathwayIndex() to share reconciled gene identities and pathway
  memberships across enrichment, RaMP expression assays, pathway scores,
  GESECA, named plots and network coverage. Pathway IDs distinguish different
  sets with the same display name. Reuse validates membership fingerprints.
* createPathwayAssay() now maps human genes with the audited HGNC index and
  preserves the selected expression layer. Identical duplicate gene rows are
  counted once; different rows require explicit mean or sum aggregation.
  createPathwayObject() accepts original gene identifiers or mapped assays.
* Pathway outputs retain database sizes before/after reconciliation, measured
  coverage, actual used members and excluded conflicting RaMP records.
  Zero-coverage and size-excluded pathways remain visible in coverage audits.
  Scoring and enrichment use the same measured membership sets. Regional and
  GESECA size filters use members actually supplied to their test engines.

# SpaMTP 0.99.5

* Added `buildGeneMappingIndex()`, `mapGeneIdentifiers()` and native
  `annotateGeneIdentifiers()` using SpaMTPdb's versioned HGNC reference.
  Approved/previous/alias symbols and stable IDs have explicit resolution
  status; ambiguous aliases and conflicting RaMP records are never expanded
  into multiple genes. Gene-specific IDs can resolve shared protein IDs.
* Fisher and regional pathway analyses merge RaMP records belonging to the
  same HGNC gene and unite their pathway memberships before counting. Network
  displays reuse the same identity mapping. Conflicting differential values
  from multiple features of one gene require resolution before analysis.
  Gene mapping and reference provenance are retained in outputs and rowData.
* Gene indices now record a normalized source-table fingerprint and reject
  incompatible reuse even when RaMP record IDs are unchanged. Older index
  schemas must be rebuilt. Shared crossreferences retain unresolved and
  conflicting raw candidates instead of implying a unique match. Conflict
  membership is indexed once and repeated queries are resolved once while
  preserving the original input order and per-row audit.
* Official database workflows use HGNC mapping automatically for human genes;
  custom fixtures remain offline, and `gene_mapping = "ramp"` supports explicit
  historical reproduction or curated non-human resources. HGNC mapping does
  not perform cross-species orthology conversion.

* Fixed the Fisher pathway contingency table: foreground non-members are
  K - overlap, not max(0, K - pathway_size). Added an explicit measured
  `universe`, deduplicated RaMP memberships and validated foreground inclusion.
  Background size is fixed before pathway-size filtering. Every eligible
  pathway, including zero-overlap sets, now belongs to the BH testing family.
  `pval_cutoff` now filters FDR as documented (previously it filtered raw p).
  Results include foreground/background counts, internal pathway IDs and an
  `enrichment` audit attribute; missing display metadata no longer drops tests.
  These corrections intentionally change previous p-values and FDRs.

# SpaMTP 0.99.4

* Updated the companion workflow for SpaMTPData >= 0.99.4 and its published
  native resource release 1.1.0. Historical archive examples explicitly select
  1.0.0; native downloads require no Seurat conversion.
* Preserve non-syntactic feature and pixel metadata names during Seurat
  conversion, including labelled metabolite intensities in the human brain
  archive. Added forward/reverse conversion regression tests.
* Support the companion data package's one-time native resource preparation:
  resulting SpatialExperiment RDS files are used directly without Seurat.

# SpaMTP 0.99.3

* Synchronized companion-package calls with the camelCase SpaMTPdb >= 0.99.2
  and SpaMTPData >= 0.99.1 APIs. Added an evaluated three-package workflow
  using a native SpatialExperiment example without downloads or Seurat.
* Passed local-file checksum verification through loadSpaMTPDatabase(), and
  isolated in-session cache entries by configured resource directory and
  verification setting. Explicit custom database bundles remain supported.
* Updated the full annotation vignette to use the current registry's bytes
  field rather than the obsolete serialized_bytes field.

# SpaMTP 0.99.2

* Migrated the complete indexed annotation tutorial to native containers and
  added a runnable public-accessor guide for expression, metadata, coordinates,
  images and paired modalities. Static regression checks cover actual R code
  in both functions and vignettes, including unevaluated examples, without
  mistaking SMILES chirality or `slot=` argument names for slot operations.
* Seurat spatial conversion now prioritizes a single image/FOV accessor over
  historical metadata coordinates and records the selected coordinate source.
  Multiple images require explicit selection; finite-radius centroids are read
  as pixel centres rather than expanded polygon vertices. Explicitly supplied
  coordinates remain supported, and metadata is the fallback for image-free
  objects only. Input objects and their archival metadata are preserved.
* Native CSV and METASPACE imports store coordinates only in `spatialCoords()`,
  removing redundant metadata copies that could become stale after alignment.
* Made `SpatialExperiment` the default container for aligned/binned MSI and
  retained Cardinal `MSImagingArrays` / `MSImagingExperiment` for raw, continuous, or file-backed
  spectra. Added S4 conversion generics, registered Cardinal/SPE coercions,
  and explicit optional Seurat conversion functions.
* Added raw-array S4 preprocessing methods, execute deferred Cardinal
  processing before binning/conversion, and cover the mass range of all raw
  spectra rather than inferring it only from the first pixel.
* Mapped MSI feature annotations, pixel metadata, coordinates, optical images,
  and paired transcriptomes to `rowData()`, `colData()`, `spatialCoords()`,
  `imgData()`, and `altExp()`, respectively. Seurat was removed from dependency
  declarations; SeuratObject is optional and used only by explicit converters.
* Added S4 methods for binning, normalization, multi-omic integration, spatial
  plotting, images, and paired transcriptomes. Native workflows reuse
  Cardinal, SpatialExperiment, SingleCellExperiment, scater, edgeR, limma,
  DropletUtils and fgsea functionality.
* Removed unused direct Imports of EBImage, sp, shinyjs, RColorBrewer, zeallot
  and matter after migrating the image/alignment workflows. Raw MSI continues
  to use matter through Cardinal; optical images use SpatialExperiment.
* Pathway and merged-modality matrices are stored as alternative experiments,
  allowing their feature dimensions to differ from the primary MSI assay.
* Migrated affine alignment, image attachment, pixel/spot mapping, ROI selection,
  spatial correlations, Moran's I and pathway-score plots to native containers.
  Mapping and spatial graphs separate sample identities; graphs also use
  `colPairs()` so standard subsetting updates their indices.
* Added regression coverage for reordered pixels, multi-sample coordinates,
  altExp normalization, intensity-conserving matrix bins, image scale factors,
  constant features and RaMP feature IDs. Fixed limma treatment-table extraction
  and preserved explicit sample labels during technical pooling with edgeR.
* Native integration defaults to main MSI plus the paired transcriptome; derived
  pathway and merged assays are not integrated automatically. Equal-weight
  concatenated PCA is explicitly distinguished from historical Seurat WNN.
* Migrated MSI feature/annotation plots, optical overlays, 3D views, mass
  spectra, interactive mass windows, density export and pathway-network
  container extraction to SpatialExperiment/altExp and imgData.
  Plotting rejects unsupported Seurat-only arguments instead of ignoring them.
  Mass windows now sum all in-range features and overlapping windows count
  each peak once; m/z axes no longer require specially formatted feature IDs.
  Fixed row-major raster colour ordering for transparent optical overlays and
  3D image points; S4 container classes are now imported explicitly for coercion.
* Replaced pathway-based annotation ranking's temporary Seurat/Cardinal
  conversion with direct paired Pearson correlations. Reused existing pathway
  scoring and current RaMP candidate resolution; respected requested assays,
  removed duplicate candidate IDs, and defined constant/missing-data handling.
  The retained weighted z-score tail probabilities are explicitly documented
  as heuristic ranks, not calibrated identification p-values.
* Fixed annotation candidate sorting when optional Score/Error columns are
  absent. Retired four unused internal plotting helpers and their help pages.
* Removed non-conversion Seurat branches and implicit Seurat importer output.
  Added `seuratToSingleCellExperiment()` for non-spatial input, explicit named
  coordinate conversion, strict Seurat layer matching and identity transfer.
* Added `scaleSMData()` for feature-wise centring/scaling, and corrected PCA
  and expression plotting to use the requested alternative experiment.
  QC plot inputs and `runDE()` now use `data`, not Seurat-named parameters;
  `subsetSPM()` no longer accepts retired object/slot-upgrade arguments.
* Native annotation, curated FMP10 annotation and feature subsetting preserve
  arbitrary feature IDs and paired modalities. Pathway entry points use
  native assays instead of Seurat-style double-bracket indexing.
* Added a fresh-process native workflow test that requires neither Seurat
  namespace, plus adapter-boundary and conversion regression tests.
  See `BIOCONDUCTOR_MIGRATION.md` for scientific limits and remaining validation.

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
