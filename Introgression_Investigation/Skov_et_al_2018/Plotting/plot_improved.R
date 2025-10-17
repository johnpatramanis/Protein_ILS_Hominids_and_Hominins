library(dplyr)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

map_to_range <- function(input_value) {
  
  input_value <- pmin(1000, input_value)
  output_value = round(input_value * 0.018, 0) + 7
  
  return(output_value)
}

#### Default to MMP20

GENE_NAME = "MMP20"
GENE_CHROM = "chr11"
GENE_START = 102576832
GENE_END = 102625332

#### If other Genes should be explored

if (args){
  if (size(args)==4){
    GENE_NAME = args[1]
    GENE_CHROM = args[2]
    GENE_START = as.integer(args[3]) 
    GENE_END = as.integer(args[4]) 
  }
}


# ---------------------------------------------
# You can change these:
REGION_TO_INCULDE = c(
  'MIDDLE_EAST',
  'AMERICA',
  'EAST_ASIA',
  'EUROPE',
  'CENTRAL_SOUTH_ASIA',
  'Genes')


WINDOW_SIZE =250000 
PLOT_EDGES = 250000
DOWNSAMPLE_TO = 1 # downsample data to 50%
SEGMENT_TYPES_TO_INCLUDE = c(
  'Neanderthal', 
  'Denisova', 
  'none',
  'Both')

# ---------------------------------------------

# Load genes
genes = read.table('GenCODE_uniq_hg38.txt', col.names = c('chrom', 'start', 'end', 'genelength', 'genename')) %>%
  filter(start < GENE_END + PLOT_EDGES, end > GENE_START - PLOT_EDGES, chrom == GENE_CHROM) %>%
  mutate(plot_rownumber = 1:n()) %>%
  mutate(region = 'Genes', population = 'Genes') %>%
  mutate(region = factor(region, REGION_TO_INCULDE)) 




# Load segments
data = read.table(paste0(GENE_NAME, '.segments'), header = T)  %>%
  filter(start < GENE_END + WINDOW_SIZE, 
         end > GENE_START - WINDOW_SIZE, 
         Archaic %in% SEGMENT_TYPES_TO_INCLUDE, 
         region %in% REGION_TO_INCULDE) %>%
  mutate(region = factor(region, REGION_TO_INCULDE)) 


order = data %>%
  group_by(individual_name, region, population) %>%
  summarize(min_start = min(start), max_end = max(end)) %>%
  ungroup() %>%
  group_by(region, population) %>%
  sample_n(round(n() * DOWNSAMPLE_TO,0)) %>%
  arrange(min_start, desc(max_end - min_start)) %>%
  
  mutate(plot_rownumber = 1:n()) %>%
  ungroup() %>%
  select(individual_name, plot_rownumber, min_start, max_end) 


n_ind = length(unique(order$individual_name)) 
PLOT_HEIGHT = map_to_range(n_ind)

cat('n ind = ',n_ind, 'plot height=', PLOT_HEIGHT, '\n')



# load snps
snp_data = read.table(paste0(GENE_NAME, '.snps'), header = T) %>%
  inner_join(order, by = c('individual_name'))  %>%
  filter(
    Archaic %in% SEGMENT_TYPES_TO_INCLUDE, 
    region %in% REGION_TO_INCULDE) %>%
  mutate(region = factor(region, REGION_TO_INCULDE)) 


PLOT_LEFT_END = GENE_START - PLOT_EDGES
PLOT_RIGHT_END = GENE_END + PLOT_EDGES


data %>%
  inner_join(order, by = c('individual_name')) %>%
  ggplot() +
  
  # add genes
  geom_segment(data = genes, aes( x = start, xend = end, y = -plot_rownumber-10, yend = -plot_rownumber-10), color = "#CC79A7", size = 1) +
  geom_text(data = genes, aes( x = end, y = -plot_rownumber-8, label = genename), hjust = 0, size = 4) +
  
  # segments 
  geom_segment(aes(x = start, xend = end, y = plot_rownumber, yend = plot_rownumber), color = 'lightgrey', size = 1) +
  
  # snps
  geom_point(data = snp_data, aes(x = pos, y = plot_rownumber, color = snptype, size = snptype), size = 0.5) + 
  
  
  # make plot pretty
  theme_bw() +
  theme(plot.title = element_text(size=42,face="bold",margin=margin(25,0,40,0)),
        axis.ticks.length=unit(0.5,"cm"),
        axis.text = element_text(size=25),
        axis.title = element_text(size=32,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.title.x = element_text(hjust=0.75,margin = margin(t = 25, r = 0, b = 0, l = 0)),
        legend.text=element_text(size=20), legend.title=element_text(size=25),
        legend.key.size = unit(3,"line")) +
  
  scale_color_manual(values = c("ND10" = "#56B4E9", 
                                "ND01" = "#009E73", 
                                "linkedDAV" = "black", 
                                "ND11" = "#E69F00", 
                                "unlinked" = "darkgrey")) +
  guides(color = guide_legend(override.aes = list(size = 7.5))) +
  facet_grid(region~., scale = 'free', space = 'free') +
  coord_cartesian(xlim = c(PLOT_LEFT_END, PLOT_RIGHT_END)) +
  scale_y_continuous('Introgressed haplotypes') +
  xlab('Genomic coordinates on Chromosome 11') +
  ggtitle(paste0('Introgressed Regions around ', GENE_NAME)) +
  theme(strip.text.y = element_text(angle = 0,size=20,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing=unit(0, "lines"),
        axis.text.x = element_text(hjust=-0.01,margin = margin(t = 10, r = 0, b = 2.5, l = 0)))


OUT_PLOT_FILE = paste0(GENE_CHROM, '_', GENE_START,'_', GENE_END, '_', GENE_NAME, '.pdf')
OUT_PLOT_FILE_2 = paste0(GENE_CHROM, '_', GENE_START,'_', GENE_END, '_', GENE_NAME, '.svg')
ggsave(OUT_PLOT_FILE, width = 22, height = PLOT_HEIGHT)
ggsave(OUT_PLOT_FILE_2, width = 22, height = PLOT_HEIGHT)

