#Shell script to run MAGMA pathway analysis on all kegg legacy pathways using DIAMANTE sex stratified T2D GWAS data. (change naming conventions to male/female depending on which sex you are using) The GWAS data, gene locations, and pathways are all publically available.
  

#First annotate the genome using snp locations from the gwas and MAGMA provided NCBI37.3 gene locations
  magma \
  --annotate \
  --snp-loc /path/to/diamante/gwas/female/snploc.txt \
  --gene-loc /path/to/NCBI37.3.gene.loc \
  --out NCBI_t2d_female
  
#Next run gene level enrichment analysis
 for chr in {1..22}; do
    magma \
      --bfile /path/to/european/1000G/reference/chromosome \
      --pval /path/to/Diamante/GWAS use=rsid,Pvalue ncol=Neff \
      --gene-annot NCBI_t2d_female.genes.annot \
      --out ncbi_t2d_female_gene_chr${chr}
done
#Make sure to merge individual chromosome results into one file before next step

#Use gene level results to run pathway level results for kegg legacy pathways

magma \
  --gene-results /path/to/chromosome/merged/raw/gene/results \
  --set-annot /path/to/kegg_legacy.v2025.gmt \
  --out pathway_level_allkegg_female
