# Functional-gut-dysbiosis-is-associated-with-overall-survival-after-solid-organ-transplantation
Codes for the project "Functional gut dysbiosis is associated with overall survival after solid organ transplantation" based on TransplantLines biobank and cohort study

                    1,231 metagenomes
              1,008 SOTRs + 223 healthy controls
                           │
                           ▼
              ┌─────────────────────────┐
              │  A. READ PREPROCESSING  │
              │ KneadData + Bowtie2     │
              │ + FastQC                │
              └────────────┬────────────┘
                           │
                    High-quality reads
                           │
             ┌─────────────┴─────────────┐
             ▼                           ▼
    ┌──────────────────┐       ┌──────────────────────┐
    │ B. MAG           │       │ C. GENE-centric      │
    └────────┬─────────┘       └──────────┬───────────┘
             │                            │
       Assembly                        Assembly
             │                            │
       MAG binning                   BAKTA CDS prediction
             │                            │
       MAG refinement                Gene filtering
             │                            │
       MAG abundance                 MMseqs2 catalogue
             │                            │
       CheckM quality                KofamScan
             │                            │
       GTDB-Tk taxonomy              Bowtie2 + CoverM
             │                            │
       GraPhlAn                      KO abundance
             │                            │
       anvi'o metabolism                  │
             │                            │
             └──────────────┬─────────────┘
                            ▼
                  Functional profiles
                            │
             ┌──────────────┼──────────────┐
             ▼              ▼              ▼
       D. Dysbiosis    E. Medication    F. Mortality
           score          regimen          analysis
             │              │              │
             └──────────────┼──────────────┘
                            ▼
                    Results & Figures 
