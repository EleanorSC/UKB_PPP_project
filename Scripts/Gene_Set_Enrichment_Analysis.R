## ---------------------------
##
## Script Purpose: GSEA
## forked from https://benjacobs123456.github.io/ukb_proteomics/#GSEA_(MRI)
## -----------------------------

# read in DE proteins 
de = read_csv("~/ukb_proteomics/outputs/limma_results_full_normalised_t2_lesion_vol.csv")

# sort in descending order 
de_sorted = de %>%
  mutate(rank = sign(beta) * -log10(pval)) %>%
  arrange(desc(rank))

# create gene list 
protein_list = de_sorted$rank
names(protein_list) = de_sorted$protein

# fgsea 
library(msigdbr)
library(fgsea)

# get gene sets
kegg_gene_sets = msigdbr(species = "Homo sapiens", category = "C2") %>% filter(gs_subcat=="CP:KEGG")
kegg_list = split(x = kegg_gene_sets$gene_symbol, f = kegg_gene_sets$gs_name)

reactome_gene_sets = msigdbr(species = "Homo sapiens", category = "C2") %>% filter(gs_subcat=="CP:REACTOME")
reactome_list = split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)

hallmark_gene_sets = msigdbr(species = "Homo sapiens", category = "H") 
hallmark_list = split(x = hallmark_gene_sets$gene_symbol, f = hallmark_gene_sets$gs_name)

go_gene_sets = msigdbr(species = "Homo sapiens", category = "C5") %>% filter(gs_subcat=="GO:BP")
go_list = split(x = go_gene_sets$gene_symbol, f = go_gene_sets$gs_name)

# kegg
gsea_res_kegg = fgseaMultilevel(pathways = kegg_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)

# reactome 
gsea_res_reactome = fgseaMultilevel(pathways = reactome_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)

# hallmark
gsea_res_hallmark = fgseaMultilevel(pathways = hallmark_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)

# GO
gsea_res_go = fgseaMultilevel(pathways = go_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)

gsea_res_kegg %>%
  filter(NES > 0) %>%
  arrange(padj)

# plot 
plot_dat = gsea_res_kegg %>%
  mutate(pathway = str_remove_all(pathway,"KEGG_")) %>%
  mutate(pathway = str_replace_all(pathway,"_"," ")) %>%
  arrange((NES)) %>%
  mutate(sig = ifelse(padj < 0.05,"*"," "))
plot_dat$pathway=factor(plot_dat$pathway,levels=plot_dat$pathway,ordered=T)

p=ggplot(plot_dat,
       aes(NES,pathway,fill=NES,label=sig))+
  geom_col(color="black")+
  geom_text(size=5)+
  scale_fill_gradient2(low="purple",high="orange2",midpoint=0)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Normalised enrichment score (NES)",y="Pathway (KEGG)",fill="NES")
png("/data/home/hmy117/ukb_proteomics/outputs/mri_lesions_gsea.png",res=900,units="in",height=5,width=7)
p
dev.off()

# get leading edge proteins 
plots = list()
edge_proteins = plot_dat[plot_dat$pathway=="COMPLEMENT AND COAGULATION CASCADES",]$leadingEdge
for(i in c(1:length(edge_proteins[[1]]))){
  
  this_prot = edge_proteins[[1]][i]
  plots[[i]] = ggplot(filtered_pheno %>%
           filter(MS_status=="prevalent"),aes(.data[[this_prot]],volume_of_brain_greywhite_matter_normalised_for_head_size_f25009_2_0))+
    geom_point()+
    geom_smooth(method="lm",linetype="dashed",alpha=0.1,se = F)+
    theme_bw()+
    labs(y="Brain volume")
  plots <<- plots
  
}

png("~/ukb_proteomics/outputs/lesion_prot_plots.png",res=900,units="in",width=16,height=16)
print(gridExtra::grid.arrange(grobs = plots))
dev.off()

edge_proteins1 = edge_proteins

# brain vol 
de = read_csv("~/ukb_proteomics/outputs/limma_results_full_volume_of_brain_greywhite_matter_normalised_for_head_size_f25009_2_0.csv")
# sort in descending order 
de_sorted = de %>%
  mutate(rank = sign(beta) * -log10(pval)) %>%
  arrange(desc(rank))

# create gene list 
protein_list = de_sorted$rank
names(protein_list) = de_sorted$protein

gsea_res_kegg = fgseaMultilevel(pathways = kegg_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)

# reactome 
gsea_res_reactome = fgseaMultilevel(pathways = reactome_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)

# hallmark
gsea_res_hallmark = fgseaMultilevel(pathways = hallmark_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)

# GO
gsea_res_go = fgseaMultilevel(pathways = go_list,
      stats = protein_list,
      minSize = 30, 
      maxSize = 1000, 
      eps = 0, 
      nPermSimple = 10000)


gsea_res_kegg %>%
  filter(NES > 0) %>%
  arrange(padj)

# plot 
plot_dat = gsea_res_kegg %>%
  mutate(pathway = str_remove_all(pathway,"KEGG_")) %>%
  mutate(pathway = str_replace_all(pathway,"_"," ")) %>%
  arrange((NES)) %>%
  mutate(sig = ifelse(padj < 0.05,"*"," "))
plot_dat$pathway=factor(plot_dat$pathway,levels=plot_dat$pathway,ordered=T)

p=ggplot(plot_dat,
       aes(NES,pathway,fill=NES,label=sig))+
  geom_col(color="black")+
  geom_text(size=5)+
  scale_fill_gradient2(low="purple",high="orange2",midpoint=0)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Normalised enrichment score (NES)",y="Pathway (KEGG)",fill="NES")
png("/data/home/hmy117/ukb_proteomics/outputs/mri_bvol_gsea.png",res=900,units="in",height=5,width=7)
p
dev.off()


# get leading edge proteins 
plots = list()
edge_proteins = plot_dat[plot_dat$pathway=="COMPLEMENT AND COAGULATION CASCADES",]$leadingEdge
for(i in c(1:length(edge_proteins[[1]]))){
  
  this_prot = edge_proteins[[1]][i]
  plots[[i]] = ggplot(filtered_pheno %>%
           filter(MS_status=="prevalent"),aes(.data[[this_prot]],volume_of_brain_greywhite_matter_normalised_for_head_size_f25009_2_0))+
    geom_point()+
    geom_smooth(method="lm",linetype="dashed",alpha=0.1,se = F)+
    theme_bw()+
    labs(y="Brain volume")
  plots <<- plots
  
}

png("~/ukb_proteomics/outputs/bvol_prot_plots.png",res=900,units="in",width=16,height=16)
print(gridExtra::grid.arrange(grobs = plots))
dev.off()




# find all complement / clotting genes 
combo_res = combo_res %>% 
  mutate(complement_or_clotting = ifelse(
    protein %in% edge_proteins[[1]] & protein %in% edge_proteins1[[1]],
    "*",""))

write_csv(combo_res,"~/ukb_proteomics/outputs/mri_results_limma_res.csv")


p=ggplot(filtered_pheno %>% filter(MS_status=="prevalent"),
       aes(F11,normalised_t2_lesion_vol))+
    geom_point()+
  theme_bw()+
  labs(x="Plasma Factor XI",y="Normalised T2 lesion volume")+
  geom_smooth(method="lm",linetype="dashed",alpha=0.1,color="red",se=F)
png("~/ukb_proteomics/outputs/f11_t2lesion_vol.png",res=900,units="in",width=3,height=3)
p
dev.off()
