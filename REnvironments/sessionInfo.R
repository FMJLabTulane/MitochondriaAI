> sessionInfo()
R version 4.5.1 (2025-06-13 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] torch_0.16.0                      gprofiler2_0.2.3                  ragg_1.5.0                       
 [4] Pando_1.1.1                       DESeq2_1.48.1                     BiocParallel_1.42.2              
 [7] cicero_1.3.9                      Gviz_1.52.0                       chromVAR_1.30.1                  
[10] motifmatchr_1.30.0                BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.76.0                  
[13] rtracklayer_1.68.0                BiocIO_1.18.0                     Biostrings_2.76.0                
[16] XVector_0.48.0                    TFBSTools_1.46.0                  JASPAR2020_0.99.10               
[19] qs_0.27.3                         R.utils_2.13.0                    R.oo_1.27.1                      
[22] R.methodsS3_1.8.2                 devtools_2.4.5                    usethis_3.2.1                    
[25] ggVennDiagram_1.5.6               ggvenn_0.1.16                     DropletUtils_1.28.1              
[28] Nebulosa_1.18.0                   circlize_0.4.16                   ComplexHeatmap_2.24.1            
[31] viridis_0.6.5                     viridisLite_0.4.2                 EnrichmentBrowser_2.38.1         
[34] graph_1.86.0                      escape_2.4.0                      dittoSeq_1.20.0                  
[37] DOSE_4.2.0                        clusterProfiler_4.16.0            MeSHDbi_1.44.0                   
[40] AnnotationHub_3.16.1              BiocFileCache_2.16.2              dbplyr_2.5.1                     
[43] org.Hs.eg.db_3.21.0               GOSemSim_2.34.0                   glmGamPoi_1.20.0                 
[46] EnhancedVolcano_1.26.0            DoubletFinder_2.0.6               future_1.67.0                    
[49] patchwork_1.3.2                   clustree_0.5.1                    ggraph_2.2.2                     
[52] plotly_4.11.0                     EnsDb.Hsapiens.v86_2.99.0         ensembldb_2.32.0                 
[55] AnnotationFilter_1.32.0           GenomicFeatures_1.60.0            AnnotationDbi_1.70.0             
[58] scDblFinder_1.22.0                Signac_1.15.0                     harmony_1.2.3                    
[61] monocle3_1.4.26                   SingleCellExperiment_1.30.1       SummarizedExperiment_1.38.1      
[64] GenomicRanges_1.60.0              GenomeInfoDb_1.44.3               IRanges_2.42.0                   
[67] S4Vectors_0.46.0                  MatrixGenerics_1.20.0             matrixStats_1.5.0                
[70] Biobase_2.68.0                    BiocGenerics_0.54.0               generics_0.1.4                   
[73] Seurat_5.3.0                      SeuratObject_5.2.0                sp_2.2-0                         
[76] reticulate_1.43.0                 data.table_1.17.8                 lubridate_1.9.4                  
[79] forcats_1.0.1                     purrr_1.1.0                       readr_2.1.5                      
[82] tidyr_1.3.1                       tibble_3.3.0                      tidyverse_2.0.0                  
[85] dplyr_1.1.4                       ggridges_0.5.7                    Matrix_1.7-4                     
[88] cowplot_1.2.0                     Rcpp_1.1.0                        SoupX_1.6.2                      
[91] hdf5r_1.3.12                      stringr_1.5.2                     leiden_0.4.3.1                   
[94] ggrepel_0.9.6                     ggplot2_4.0.0                    

loaded via a namespace (and not attached):
  [1] igraph_2.1.4                Formula_1.2-5               ica_1.0-3                   scater_1.36.0              
  [5] maps_3.4.3                  tidyselect_1.2.1            bit_4.6.0                   doParallel_1.0.17          
  [9] clue_0.3-66                 lattice_0.22-7              rjson_0.2.23                urlchecker_1.0.1           
 [13] blob_1.2.4                  S4Arrays_1.8.1              parallel_4.5.1              dichromat_2.0-0.1          
 [17] seqLogo_1.74.0              png_0.1-8                   cli_3.6.5                   ggplotify_0.1.3            
 [21] ProtGenerics_1.40.0         goftest_1.2-3               textshaping_1.0.3           bluster_1.18.0             
 [25] BiocNeighbors_2.2.0         uwot_0.2.3                  curl_7.0.0                  evaluate_1.0.5             
 [29] mime_0.13                   tidytree_0.4.6              stringi_1.8.7               desc_1.4.3                 
 [33] backports_1.5.0             XML_3.99-0.19               httpuv_1.6.16               magrittr_2.0.4             
 [37] rappdirs_0.3.3              splines_4.5.1               RcppRoll_0.3.1              mclust_6.1.1               
 [41] RApiSerialize_0.1.4         jpeg_0.1-11                 DT_0.34.0                   sctransform_0.4.2          
 [45] ggbeeswarm_0.7.2            sessioninfo_1.2.3           DBI_1.2.3                   HDF5Array_1.36.0           
 [49] withr_3.0.2                 systemfonts_1.3.0           reformulas_0.4.1            rprojroot_2.1.1            
 [53] enrichplot_1.28.4           xgboost_1.7.11.1            lmtest_0.9-40               GSEABase_1.70.1            
 [57] tidygraph_1.3.1             BiocManager_1.30.26         htmlwidgets_1.6.4           fs_1.6.6                   
 [61] biomaRt_2.64.0              SparseArray_1.8.1           h5mread_1.0.1               annotate_1.86.1            
 [65] VariantAnnotation_1.54.1    zoo_1.8-14                  knitr_1.50                  UCSC.utils_1.4.0           
 [69] AUCell_1.30.1               TFMPvalue_0.0.9             timechange_0.3.0            foreach_1.5.2              
 [73] ggpointdensity_0.2.0        caTools_1.18.3              ggtree_3.16.3               rhdf5_2.52.1               
 [77] pwalign_1.4.0               RSpectra_0.16-2             irlba_2.3.5.1               ggrastr_1.0.2              
 [81] ellipsis_0.3.2              fastDummies_1.7.5           gridGraphics_0.5-1          lazyeval_0.2.2             
 [85] yaml_2.3.10                 survival_3.8-3              SpatialExperiment_1.18.1    scattermore_1.2            
 [89] BiocVersion_3.21.1          crayon_1.5.3                RcppAnnoy_0.0.22            mapproj_1.2.12             
 [93] RColorBrewer_1.1-3          progressr_0.16.0            tweenr_2.0.3                later_1.4.4                
 [97] Rgraphviz_2.52.0            base64enc_0.1-3             profvis_0.4.0               codetools_0.2-20           
[101] GlobalOptions_0.1.2         KEGGREST_1.48.1             Rtsne_0.17                  shape_1.4.6.1              
[105] limma_3.64.3                Rsamtools_2.24.1            filelock_1.0.3              foreign_0.8-90             
[109] pkgconfig_2.0.3             KEGGgraph_1.68.0            xml2_1.4.0                  spatstat.univar_3.1-4      
[113] GenomicAlignments_1.44.0    aplot_0.2.9                 biovizBase_1.56.0           spatstat.sparse_3.1-0      
[117] ape_5.8-1                   xtable_1.8-4                interp_1.1-6                msigdb_1.16.0              
[121] plyr_1.8.9                  httr_1.4.7                  rbibutils_2.3               tools_4.5.1                
[125] globals_0.18.0              pkgbuild_1.4.8              checkmate_2.3.3             htmlTable_2.4.3            
[129] beeswarm_0.4.0              nlme_3.1-168                lme4_1.1-37                 digest_0.6.37              
[133] farver_2.1.2                tzdb_0.5.0                  reshape2_1.4.4              ks_1.15.1                  
[137] yulab.utils_0.2.1           rpart_4.1.24                DirichletMultinomial_1.50.0 glue_1.8.0                 
[141] cachem_1.1.0                polyclip_1.10-7             Hmisc_5.2-3                 mvtnorm_1.3-3              
[145] parallelly_1.45.1           pkgload_1.4.1               statmod_1.5.0               here_1.0.2                 
[149] RcppHNSW_0.6.0              ScaledMatrix_1.16.0         minqa_1.2.8                 pbapply_1.7-4              
[153] httr2_1.2.1                 fields_17.1                 spam_2.11-1                 gson_0.1.0                 
[157] coro_1.1.0                  dqrng_0.4.1                 gtools_3.9.5                graphlayouts_1.2.2         
[161] gridExtra_2.3               shiny_1.11.1                GSVA_2.2.0                  GenomeInfoDbData_1.2.14    
[165] pals_1.10                   rhdf5filters_1.20.0         RCurl_1.98-1.17             memoise_2.0.1              
[169] rmarkdown_2.30              pheatmap_1.0.13             scales_1.4.0                RANN_2.6.2                 
[173] stringfish_0.17.0           spatstat.data_3.1-8         rstudioapi_0.17.1           cluster_2.1.8.1            
[177] spatstat.utils_3.2-0        hms_1.1.3                   fitdistrplus_1.2-4          colorspace_2.1-2           
[181] rlang_1.1.6                 DelayedMatrixStats_1.30.0   sparseMatrixStats_1.20.0    dotCall64_1.2              
[185] ggforce_0.5.0               scuttle_1.18.0              ggtangle_0.0.7              xfun_0.53                  
[189] remotes_2.5.0               iterators_1.0.14            abind_1.4-8                 treeio_1.32.0              
[193] Rhdf5lib_1.30.0             ps_1.9.1                    bitops_1.0-9                Rdpack_2.6.4               
[197] promises_1.3.3              RSQLite_2.4.3               qvalue_2.40.0               fgsea_1.34.2               
[201] DelayedArray_0.34.1         GO.db_3.21.0                compiler_4.5.1              prettyunits_1.2.0          
[205] boot_1.3-31                 distributional_0.5.0        beachmat_2.24.0             listenv_0.9.1              
[209] edgeR_4.6.3                 BiocSingular_1.24.0         tensor_1.5.1                progress_1.2.3             
[213] MASS_7.3-65                 UCell_2.12.0                spatstat.random_3.4-2       R6_2.6.1                   
[217] fastmap_1.2.0               fastmatch_1.1-6             vipor_0.4.7                 ROCR_1.0-11                
[221] nnet_7.3-20                 rsvd_1.0.5                  ggdist_3.3.3                gtable_0.3.6               
[225] KernSmooth_2.23-26          latticeExtra_0.6-31         miniUI_0.1.2                deldir_2.0-4               
[229] htmltools_0.5.8.1           RcppParallel_5.1.11-1       bit64_4.6.0-1               spatstat.explore_3.5-3     
[233] lifecycle_1.0.4             S7_0.2.0                    processx_3.8.6              callr_3.7.6                
[237] nloptr_2.2.1                restfulr_0.0.16             vctrs_0.6.5                 VGAM_1.1-13                
[241] spatstat.geom_3.6-0         scran_1.36.0                ggfun_0.2.0                 future.apply_1.20.0        
[245] pracma_2.4.4                pillar_1.11.1               magick_2.9.0                metapod_1.16.0             
[249] locfit_1.5-9.12             jsonlite_2.0.0              GetoptLong_1.0.5  
