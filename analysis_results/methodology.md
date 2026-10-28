# Methodology: seed and seedling microbiome meta-analysis

Published 16S rRNA amplicon datasets (eight studies, 16 host species, 10 countries) were compiled into one phyloseq object (McMurdie & Holmes, 2013), following Simonin *et al.* (2022), with developmental stage (seeds versus seedlings) as the primary contrast. Non-bacterial, chloroplast, mitochondrial, and cyanobacterial ASVs, soil samples, samples with fewer than 500 reads, and ASVs present in fewer than three samples or with fewer than 10 reads were removed, leaving 1,886 samples and 2,823 ASVs; no contaminant removal was possible, as no negative controls were deposited.

Alpha diversity (observed richness, Shannon diversity, and Pielou's evenness) was calculated on data rarefied to 500 reads (1,000 as a sensitivity check) and compared with Kruskal–Wallis tests, with effect size (H − k + 1)/(n − k), across all samples for each grouping variable (stage, plant family, and host genus), and between stages within each study, plant family, and host genus that sampled both stages (at least three samples per stage for families and genera). The stage effect within plant families and host genera was also expressed as Hedges' g and tested across them by random-effects meta-analysis, with Cochran's Q for heterogeneity among strata. Beta diversity was analysed as Bray–Curtis dissimilarities (Bray & Curtis, 1957) of counts normalised by cumulative sum scaling (Paulson *et al.*, 2013), at the ASV and genus levels, visualised by principal coordinates analysis, and tested by PERMANOVA (Anderson, 2001; 999 permutations), for each variable alone and in a marginal model of stage, host genus, and study; seed and seedling samples were also ordinated separately, at both levels, and host genus, plant family, and study were tested within each stage by PERMANOVA. Composition was summarised at class, family, and genus levels as the mean of per-sample relative abundances, and co-occurrence among the 50 most abundant ASVs of each developmental stage as networks of Spearman correlations, built separately for seeds and seedlings (|ρ| ≥ 0.6, adjusted p < 0.05 within each stage). The mean genus-level profiles of each host plant genus at each stage (at least three samples) were clustered by average linkage on Bray–Curtis dissimilarities, separately for seeds and seedlings.

Associations between taxa and stage were examined in two steps. ASVs were screened by DESeq2 (Love *et al.*, 2014; `poscounts` size factors; adjusted p < 0.05, absolute log2 fold change > 1) and by indicator analysis (IndVal.g; De Cáceres & Legendre, 2009; 999 permutations, uncorrected p < 0.05), and the genera holding an ASV significant in both were retained (171 ASVs in 23 genera). Genus counts were then transformed to centred log-ratios (Aitchison, 1986; pseudocount 0.5), and the difference between seedlings and seeds was expressed as Hedges' g (Hedges, 1981) within each plant family and each host genus holding at least three samples per stage (six families and eight genera), where the genus occurred in at least 10% of the samples. Effects were pooled by random-effects meta-analysis (restricted maximum likelihood; Viechtbauer, 2010), with the Knapp–Hartung adjustment (Knapp & Hartung, 2003) for plant families and a random effect of study for host genera, and each pooled estimate was refitted with every stratum omitted in turn.

Host phylogenetic signal in alpha diversity was tested on a time-calibrated phylogeny of the host species (TimeTree 5; Kumar *et al.*, 2022), using the mean of each index per host species and stage (at least three samples), with Blomberg's *K* (Blomberg *et al.*, 2003; 999 permutations) and Pagel's λ (Pagel, 1999), and a bootstrap 95% confidence interval for λ (B = 500). Phylosymbiosis (Brooks *et al.*, 2016) was tested, within each stage, by relating the mean Bray–Curtis dissimilarity between host species to their divergence time, with Mantel tests (Mantel, 1967) and multiple regression on distance matrices (Lichstein, 2007; 999 permutations), with and without the proportion of between-study sample pairs as a covariate. The contributions of host (as genus or plant family) and stage were quantified by marginal PERMANOVA and variance partitioning (Peres-Neto *et al.*, 2006).

Genera were classified as seed specialists, seedling specialists, or stage generalists with the multinomial species classification method (CLAM; Chazdon *et al.*, 2011; specialisation threshold 2/3, α = 0.005), using the five studies that sampled both stages.

Host generalists and specialists were identified among ASVs and genera from the host species presence/absence matrix, in two steps. First, a breadth threshold was derived from a null model in which each host species received a random draw of taxa equal to its observed richness (1,000 iterations); the threshold was the smallest number of host species, above the expected breadth, at which more taxa were observed than the null mean (six host species for ASVs and seven for genera). Second, taxa meeting the threshold were tested for phylogenetic clustering of their hosts, by comparing the mean pairwise phylogenetic distance (Webb *et al.*, 2002) among their hosts with that of every possible set of the same number of host species. Taxa meeting the threshold without significant clustering were classed as generalists, those occurring on a single host species as specialists, and the remainder as intermediate; this gave 241 generalist and 1,120 specialist ASVs, and 53 generalist and 141 specialist genera.

Maximum-likelihood trees were built for the generalist ASVs (MUSCLE v5.3, Edgar, 2022; trimAl v1.5, Capella-Gutiérrez *et al.*, 2009; and IQ-TREE v3.0.1, Wong *et al.*, 2026, with ModelFinder, Kalyaanamoorthy *et al.*, 2017, and 1,000 ultrafast bootstraps, Hoang *et al.*, 2018) and for the specialist ASVs and both classes together (MAFFT v7.526, Katoh & Standley, 2013; trimAl; and FastTree v2.2.0, Price *et al.*, 2010, under GTR with gamma-distributed rates); alignments were trimmed with the `automated1` setting of trimAl. Trees were midpoint-rooted with phangorn (Schliep, 2011) and drawn with ggtree (Yu, 2022), with rings for taxonomy, predicted oxygen use, host plant family, and host breadth added with ggtreeExtra (Xu *et al.*, 2021). The two classes were compared by Fisher's exact tests at each taxonomic rank, and host breadth was related to non-zero relative abundance, within the five phyla richest in ASVs, by linear regression of log10 abundance and one-sided Spearman correlation.

Traits were predicted per ASV with BugBase (v0.1.0; Ward *et al.*, 2017), FAPROTAX (v1.2.12; Louca *et al.*, 2016), and metaTraits (GTDB genus-level summary; Podlesny *et al.*, 2026); for BugBase, ASVs were matched to Greengenes 13.8 97% OTUs (McDonald *et al.*, 2012) with vsearch at ≥ 97% identity (Rognes *et al.*, 2016). Binary traits with at least 100 annotated ASVs were tested against log2(n_hosts + 1) with binomial generalised linear models, and numeric traits by Spearman correlation with n_hosts. KEGG orthologue (Kanehisa & Goto, 2000) content was predicted with PICRUSt2 (v2.6.2; Douglas *et al.*, 2020; maximum NSTI 0.15); the 7,883 orthologues varying among ASVs were tested in the same logistic framework, classed as generalist- or specialist-enriched at FDR < 0.05 and an absolute log2 odds ratio of at least 0.5, and annotated to KEGG level-2 categories and to the plant-growth-promotion functions of the PGPg_finder ontology (Pellegrinetti *et al.*, 2024). Three groups (predicted anaerobic, nitrogen-fixing, and small-genome bacteria, the last being the smallest 1% of genome sizes) were described by ASV number, relative abundance, taxonomic and generalist–specialist composition, and host loading; their member taxa detected on at least three host genera were tested among host plant orders with Kruskal–Wallis tests (five anaerobic and four nitrogen-fixing taxa; too few small-genome taxa qualified). Traits, orthologue categories, and special groups were also compared between stages: the share of each sample's annotated relative abundance held by ASVs carrying a trait, the share of its predicted orthologue copies in a category (ASV counts divided by predicted 16S copy number), and the relative abundance of each group were expressed as Hedges' g (seedlings minus seeds) within each plant family holding at least three samples per stage, and pooled by random-effects meta-analysis with the Knapp–Hartung adjustment.

Analyses were performed in R (v4.3.2; R Core Team, 2023) with phyloseq, vegan (Oksanen *et al.*, 2025), metagenomeSeq, DESeq2, indicspecies, igraph (Csardi & Nepusz, 2006), phytools (Revell, 2024), ape (Paradis & Schliep, 2019), phangorn, metafor, ggplot2 (Wickham, 2016), ggtree, and ggtreeExtra, and in Python (v3.14), with settings read from `analysis_config.json`. Significance was set at α = 0.05, multiple comparisons were corrected by the Benjamini–Hochberg method (Benjamini & Hochberg, 1995), and headline statistics were checked by an independent reimplementation.

## References

Aitchison J. 1986. *The statistical analysis of compositional data*. London, UK: Chapman & Hall.

Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. *Austral Ecology* 26: 32–46.

Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B* 57: 289–300.

Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. *Evolution* 57: 717–745.

Bray JR, Curtis JT. 1957. An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs* 27: 325–349.

Brooks AW, Kohl KD, Brucker RM, van Opstal EJ, Bordenstein SR. 2016. Phylosymbiosis: relationships and functional effects of microbial communities across host evolutionary history. *PLoS Biology* 14: e2000225.

Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics* 25: 1972–1973.

Chazdon RL, Chao A, Colwell RK, Lin S-Y, Norden N, Letcher SG, Clark DB, Finegan B, Arroyo JP. 2011. A novel statistical method for classifying habitat generalists and specialists. *Ecology* 92: 1332–1343.

Csardi G, Nepusz T. 2006. The igraph software package for complex network research. *InterJournal, Complex Systems* 1695.

De Cáceres M, Legendre P. 2009. Associations between species and groups of sites: indices and statistical inference. *Ecology* 90: 3566–3574.

Douglas GM, Maffei VJ, Zaneveld JR, Yurgel SN, Brown JR, Taylor CM, Huttenhower C, Langille MGI. 2020. PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology* 38: 685–688.

Edgar RC. 2022. Muscle5: high-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. *Nature Communications* 13: 6968.

Hedges LV. 1981. Distribution theory for Glass's estimator of effect size and related estimators. *Journal of Educational Statistics* 6: 107–128.

Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. 2018. UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35: 518–522.

Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14: 587–589.

Kanehisa M, Goto S. 2000. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 28: 27–30.

Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30: 772–780.

Knapp G, Hartung J. 2003. Improved tests for a random effects meta-regression with a single covariate. *Statistics in Medicine* 22: 2693–2710.

Kumar S, Suleski M, Craig JM, Kasprowicz AE, Sanderford M, Li M, Stecher G, Hedges SB. 2022. TimeTree 5: an expanded resource for species divergence times. *Molecular Biology and Evolution* 39: msac174.

Lichstein JW. 2007. Multiple regression on distance matrices: a multivariate spatial analysis tool. *Plant Ecology* 188: 117–131.

Louca S, Parfrey LW, Doebeli M. 2016. Decoupling function and taxonomy in the global ocean microbiome. *Science* 353: 1272–1277.

Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15: 550.

Mantel N. 1967. The detection of disease clustering and a generalized regression approach. *Cancer Research* 27: 209–220.

McDonald D, Price MN, Goodrich J, Nawrocki EP, DeSantis TZ, Probst A, Andersen GL, Knight R, Hugenholtz P. 2012. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. *The ISME Journal* 6: 610–618.

McMurdie PJ, Holmes S. 2013. phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE* 8: e61217.

Oksanen J, Simpson GL, Blanchet FG, Kindt R, Legendre P, Minchin PR, O'Hara RB, Solymos P, Stevens MHH, Szoecs E *et al.* 2025. *vegan: community ecology package*. R package version 2.7-1.

Pagel M. 1999. Inferring the historical patterns of biological evolution. *Nature* 401: 877–884.

Paradis E, Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35: 526–528.

Paulson JN, Stine OC, Bravo HC, Pop M. 2013. Differential abundance analysis for microbial marker-gene surveys. *Nature Methods* 10: 1200–1202.

Pellegrinetti TA, Monteiro G, Lemos LN, Santos RAC, Barros A, Mendes L. 2024. PGPg_finder: a comprehensive and user-friendly pipeline for identifying plant growth-promoting genes in genomic and metagenomic data. *Rhizosphere* 30.

Peres-Neto PR, Legendre P, Dray S, Borcard D. 2006. Variation partitioning of species data matrices: estimation and comparison of fractions. *Ecology* 87: 2614–2625.

Podlesny D, Kim CY, Robbani SM, Schudoma C, Fullam A, Reimer LC, Koblitz J, Schober I, Iyappan A, Van Rossum T *et al.* 2026. metaTraits: a large-scale integration of microbial phenotypic trait information. *Nucleic Acids Research* 54: D835–D841.

Price MN, Dehal PS, Arkin AP. 2010. FastTree 2 – approximately maximum-likelihood trees for large alignments. *PLoS ONE* 5: e9490.

R Core Team. 2023. *R: a language and environment for statistical computing*. Vienna, Austria: R Foundation for Statistical Computing.

Revell LJ. 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). *PeerJ* 12: e16505.

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. 2016. VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4: e2584.

Schliep KP. 2011. phangorn: phylogenetic analysis in R. *Bioinformatics* 27: 592–593.

Simonin M, Briand M, Chesneau G, Rochefort A, Marais C, Sarniguet A, Barret M. 2022. Seed microbiota revealed by a large-scale meta-analysis including 50 plant species. *New Phytologist* 234: 1448–1463.

Viechtbauer W. 2010. Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software* 36: 1–48.

Ward T, Larson J, Meulemans J, Hillmann B, Lynch J, Sidiropoulos D, Spear JR, Caporaso G, Blekhman R, Knight R *et al.* 2017. BugBase predicts organism-level microbiome phenotypes. *bioRxiv*: 133462.

Webb CO, Ackerly DD, McPeek MA, Donoghue MJ. 2002. Phylogenies and community ecology. *Annual Review of Ecology and Systematics* 33: 475–505.

Wickham H. 2016. *ggplot2: elegant graphics for data analysis*. New York, NY, USA: Springer.

Wong TKF, Ly-Trong N, Ren H, Demotte P, Baños H, Roger AJ, Susko E, Bielow C, De Maio N, Goldman N *et al.* 2026. IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Molecular Biology and Evolution*. doi: 10.1093/molbev/msag117.

Xu S, Dai Z, Guo P, Fu X, Liu S, Zhou L, Tang W, Feng T, Chen M, Zhan L *et al.* 2021. ggtreeExtra: compact visualization of richly annotated phylogenetic data. *Molecular Biology and Evolution* 38: 4039–4042.

Yu G. 2022. *Data integration, manipulation and visualization of phylogenetic trees*. Boca Raton, FL, USA: Chapman & Hall/CRC.

