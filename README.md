# Access to germinal center IL-4 microniches drives tissue-divergence of IgE memory responses

## Article information

**Title:** Access to germinal center IL-4 microniches drives tissue-divergence of IgE memory responses

**Authors:** Sergio Villazala-Merino 1\*, Lucas Bertoia 1, Romain Fenouil 1, Mohammed Zghaebi 2, Myriam Moussa 1, Samuel Origlio 1, Claude Gregoire 1, Laura Almada 3, Mara Esposito 4, Stefano Colombo 5, Adriana Gruppi 3, Judith E Allen 5, Andrew MacDonald 5,6, Nicolas Fazilleau 7, Pierre Bruhns 8, Julia Eckl-Dorna 2, Carolyn G King 4, & Mauro Gaya 1,9\*

1. Centre d'Immunologie de Marseille-Luminy (CIML), Aix Marseille Université, INSERM, CNRS; Marseille, France.
2. Vienna Airway Lab, Department of Otorhinolaryngology, General Hospital and Medical University of Vienna, Vienna, Austria.
3. Centro de Investigaciones en Bioquímica Clínica e Inmunología (CIBICI-CONICET), Córdoba, Argentina.
4. Immune Cell Biology Laboratory, Department of Biomedicine, University of Basel, University Hospital Basel; Basel CH-4031, Switzerland.
5. Lydia Becker Institute of Immunology and Inflammation, School of Biological Sciences, Faculty of Biology, Medicine and Health, Manchester Academic Health Science Centre, University of Manchester; Manchester M13 9PT, United Kingdom.
6. Institute of Immunology and Infection Research, School of Biological Sciences, University of Edinburgh; Edinburgh, United Kingdom.
7. Infinity, Université de Toulouse, INSERM, CNRS; Toulouse, France.
8. Institut Pasteur, Université Paris Cité, INSERM UMR1222, Antibodies in Therapy and Pathology, 75015 Paris, France.
9. Lead contact

\* Correspondence: villazala@ciml.univ-mrs.fr (S.V.M.), gaya@ciml.univ-mrs.fr (M.G.)

**Summary:** Immunoglobulin E (IgE) drives allergy, yet how memory B cells (MBCs) reactivate to produce IgE, and how tissue localization shapes recall responses, remains unclear. Using mouse models of airborne exposure to house dust mites and Alternaria, we found that allergen sensitization generates lymphoid- and lung-resident MBCs. Upon allergen re-exposure, these populations followed distinct differentiation trajectories: lymph node MBCs engaged a germinal center (GC)-dependent pathway that generated IgG1+ and IgE+ plasma cells (PCs), whereas lung MBCs followed a GC-independent route producing mainly IgG1+ PCs. GC re-entry granted MBCs access to an IL-4-rich microniche formed by Tfh cells, which was essential for IgE production. Disrupting GC re-entry, IL-4 signaling, or Tfh-derived IL-4 during recall markedly reduced allergen-specific IgE. These findings reveal a spatially and cytokine-restricted mechanism that confines IgE memory to lymphoid organs, positioning GC IL-4 microniches as anatomical safeguards against IgE production at barrier sites frequently exposed to environmental antigens.

---

## Goal of this repository

This repository contains the source code (analysis scripts and Dockerfiles) used to produce the single-cell analyses reported in the article.

The three components of the study are distributed as follows:

| What | Where |
|---|---|
| Source code and container recipes | **this repository** |
| Raw sequencing data (FASTQ) | **ENA** [PRJEB105421](https://www.ebi.ac.uk/ena/browser/view/PRJEB105421) |
| Reference files, analysis results, HTML reports and container images | **Zenodo** [10.5281/zenodo.21504424](https://doi.org/10.5281/zenodo.21504424) |

---

## Dataset

A single 10x Genomics run (`230321_VH00228_159_AACKKJ7M5`, Illumina NextSeq 2000) of AID-YFP+ B cells from `AidCre(ERT2)-Rosa(eYFP)` mice exposed to house dust mite (HDM).

Cells were FACS-sorted from mediastinal lymph node (LN) and lung (LG) of 6 mice, before rechallenge (week 11, "day 0") and after rechallenge (week 12, "day 7"), and multiplexed with TotalSeq-C anti-mouse hashtag antibodies (anti-CD45 + anti-MHC class I):

| Hashtag | Sample | Hashtag | Sample |
|---|---|---|---|
| C0301, C0302, C0303 | mice 1, 2, 3 — day 0 | C0307 | mediastinal lymph node |
| C0304, C0305, C0306 | mice 4, 5, 6 — day 7 | C0308 | lung |

Sorted plasma cells (YFP+ IgD- CD45iv- CD39+ CD98+) and sorted non-plasma cells (CD39- CD98-) were loaded into **two separate 10x capture wells**, with a separate gene expression library for each. Hashtag (HTO) and BCR libraries were prepared for both wells.

### Points to be aware of when reusing the raw data

* The two HTO libraries were sequenced under a single sample identifier, and so were the two BCR libraries. Each therefore comes as **one FASTQ file pair containing both capture wells**. The i7 index recorded in each read header distinguishes them:

  | Library | Non-plasma-cell well | Plasma-cell well |
  |---|---|---|
  | HTO | `AGGCAGAA` | `TCCTGAGC` |
  | BCR | `GGACTCCT` | `TAGGCATG` |

* Consequently, `cellranger multi` was run **once over both wells**, as declared in `01_Reference/01_CellRanger/config_multi.csv`.

---

## Repository content

```
230321_VH00228_159_AACKKJ7M5/
├── 02_Container/     Dockerfiles of the execution environments
└── 03_Script/        analysis steps, numbered in execution order
    ├── globalParams.R    project paths, shared by every step
    ├── sampleParams.R    sample-level variables (mouse and tissue colours)
    └── <NN>_<name>/      one folder per analysis step
```

Each analysis step follows the same layout:

| File | Role |
|---|---|
| `analysisParams.R` | parameters of the step, **including the path to its input** |
| `Report.Rmd` | structure of the HTML report; the code lives in the numbered `.R` files it sources |
| `NN_*.R` | the analysis code itself (chunks read with `knitr::read_chunk`) |
| `launch_reports_compilation.R` | loads the parameter files and renders the report |

Some steps are plain scripts instead (`subsetSeurat.R`, `groupClusters.R`, `runVelocyto.sh`, `execute_cell_ranger_multi.sh`) or Jupyter notebooks (`scVeloAnalysis.ipynb`).

---

## Analysis workflow

Steps are numbered in execution order. Suffixes denote the subset analysed: `a` = lymph node, `b` = lung, `c` = both tissues.

| Step | What it does | Reads from |
|---|---|---|
| `01_CellRanger_FeatureBarcoding` | `cellranger multi`: alignment (mm10-2020-A), hashtag counting, BCR assembly (GRCm38 VDJ 7.0.0) | raw FASTQ |
| `02_Demux` | Seurat object, double hashtag demultiplexing (mouse, then tissue), QC and first clustering | `01` |
| `04_separateCells_LymphNode_Lung` | splits the object by tissue hashtag (LN / LG) | `02` |
| `06a`, `06b`, `06c` | main clustering, excluding T cells and excluding variable Ig genes and Ig heavy chains from the variable features | `04` (a, b), `02` (c) |
| `08_velocyto` | spliced / unspliced counts (`.loom`) | Cell Ranger BAM |
| `09c`, `09d` | merge the Seurat object with the loom counts, export `.h5ad` | `06a`, `06b`, `08` |
| `10a`, `10b` | differential expression, IgE+ cells versus other isotypes | `06a`, `06b` |
| `13a`, `13b`, `13c` | grouping and naming of clusters | `06a`, `06b`, `06c` |
| `14a`, `14b`, `14c` | analyses with the final named clusters | `06*`, `13*` |
| `18c`, `18d` | RNA velocity (scVelo), projected on the global UMAP | `09c`, `09d` |
| `20a` | projection on the global UMAP, cluster proportions (sccomp), IgE signatures | `14c`, `01` |
| `21b` | clonotype overlap between clusters and tissues (scRepertoire) | `14c`, `01` |
| `23` | **high-resolution figure export for the article** | `14c` |
| `24` | clonotype overlap, updated cluster names and colours | `14c`, `01` |
| `25` | IgE figures, updated cluster names and colours | `14c` |

Analysis steps that were purely exploratory are not included.

---

## Execution environments

Environments were built with Docker. The Dockerfiles are in `02_Container/`; the corresponding **binary images are on Zenodo**.

Rebuilding an image from its Dockerfile does not guarantee an identical environment (package versions have since moved), so downloading the images is the reliable route.

| Image (Zenodo) | Contents | Steps |
|---|---|---|
| `01_r411_tidyverse_seurat4.tar.gz` | R 4.1.1, Seurat 4.1.0, scRepertoire 1.4.0 | `02`, `04`, `06a/b/c`, `10a/b`, `13a/b/c`, `14a/b/c`, `23` |
| `02_r363_seurat_scvelo021_jupyterlab.tar.gz` | R 3.6.3, Python 3.7.3, scVelo 0.2.1, velocyto 0.17.17, JupyterLab | `08`, `09c/d`, `18c/d` |
| `03_r42_sccomp_shazam.tar.gz` | R 4.2.0, Seurat 5.0.0, sccomp 1.7.6, shazam 1.2.0, alakazam 1.3.0 | `20a`, `25` |
| `04_r44_scRepertoire.tar.gz` | R 4.4.3, Seurat 5.0.0, scRepertoire 2.5.3 | `21b`, `24` |

Alignment (`01`) used the 10x Genomics **Cell Ranger 7.0.1** singularity image, see `03_Script/01_CellRanger_FeatureBarcoding/execute_cell_ranger_multi.sh`.

---

## Reproducing the analyses

**1. Get the code and the data.**

```bash
git clone https://github.com/CIML-bioinformatic/MGlab_BmemGC_HDMrechallenge.git
```

Download from Zenodo the results archive (`230321_VH00228_159_AACKKJ7M5_noRDATA.tar.gz`) and the container images you need, then unpack the archive somewhere outside the repository:

```bash
tar xzf 230321_VH00228_159_AACKKJ7M5_noRDATA.tar.gz -C /path/to/data
docker load -i 01_r411_tidyverse_seurat4.tar.gz
```

The archive provides `01_Reference/` (gene lists, hashtag reference, cell selections) and `05_Output/` (the results of every step). Code and data are kept in separate trees, mirroring each other.

**2. Set the paths.**

Edit `03_Script/globalParams.R` so that `PATH_PROJECT` points to the folder where the data were unpacked. Every step derives its input and output paths from that file, so this is the only place to change.

**3. Run a step.**

Each step is run individually, inside the container listed above. For an R step, run its `launch_reports_compilation.R`: it loads `globalParams.R` and the step's `analysisParams.R`, then renders the report next to the results.

```bash
docker run --rm -v /path/to/data:/path/to/data -v $PWD:$PWD \
           -e PASSWORD=Pass -e USER=$(whoami) \
           -e USERID=$(id -u) -e GROUPID=$(id -g) \
           rfenouil/r411_tidyverse_seurat4 /init s6-setuidgid $(whoami) \
           Rscript "$PWD/230321_VH00228_159_AACKKJ7M5/03_Script/02_Demux/launch_reports_compilation.R"
```

Adapt the image name to the step, and mount both the data folder and the repository. The velocity steps (`09c`, `09d`, `18c`, `18d`) are Jupyter notebooks and are run in the JupyterLab image instead.

Because each step reads the output of earlier ones, a step can be re-run on its own as long as the results of its inputs are present — which is the case if the Zenodo archive has been unpacked.

---

## Notes and caveats

* **`08_velocyto` needs the Cell Ranger BAM file**, which is not included in the Zenodo archive because of its size. It can be regenerated by re-running step `01` from the ENA FASTQ files. The resulting `.loom` file *is* provided, so the downstream velocity steps can be run without it.
* **`18c` and `18d` read their input from the `09c` and `09d` output folders.** Their `analysisParams.R` sets `ANALYSIS_STEP_NAME` to the `09*` step on purpose: the `.h5ad` file lives there, and only the final HTML reports are written to the `18*` folders.
* Reports embed interactive plots and can be large (up to 1.5 GB for step `23`); opening them requires a browser with enough memory.
