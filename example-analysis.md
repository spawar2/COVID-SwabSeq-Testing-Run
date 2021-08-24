# Shrikant Pawar, 05/14/2020, SwabSeq COVID-19 RNA detection analysis
# Libraries to use
library(ggbeeswarm) 
library(MASS) 
library(speedglm)
library(furrr)
library(readxl) 
library(magrittr)
library(tidyverse)
select = dplyr::select

# Plotting style variable
theme_pub <- function(base_size = 11, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +# %+replace%
    theme(
      # grid lines
      panel.grid.major.x = element_line(colour="#ECECEC", size=0.5, linetype=1),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(colour="#ECECEC", size=0.5, linetype=1),
      panel.background   = element_blank(),
      
      # axis options
      axis.ticks.y   = element_blank(),
      axis.title.x   = element_text(size=rel(2), vjust=0.25),
      axis.title.y   = element_text(size=rel(2), vjust=0.35),
      axis.text      = element_text(color="black", size=rel(1)),
      
      # legend options
      legend.title    = element_text(size=rel(1.5)),
      legend.key      = element_rect(fill="white"),
      legend.key.size = unit(1, "cm"),
      legend.text     = element_text(size=rel(1.5)),
      
      # facet options
      strip.text = element_text(size=rel(2)),
      strip.background = element_blank(),
      
      # title options
      plot.title = element_text(size=rel(2.25), vjust=0.25, hjust=0.5)
    )
}
theme_set(theme_pub(base_size=8))

# loading barcode counts, metadata and link for barcode to amplicons
options(future.fork.enable = TRUE)
plan(multicore)
set.seed(42)

guess_max <- 100000
run_id = 'example'

#counts <- read_csv(paste0('/Users/yalegenomecenter/Desktop/SwabSeq/pipeline/', run_id, '/starcode.csv'))
counts <- read.csv(file="/Users/yalegenomecenter/Desktop/SwabSeq/pipeline/run-2/starcode.csv", header=TRUE, sep=",")
well.total <- counts %>%
  distinct(Sample_ID, Centroid, Count)  %>%
  count(Sample_ID, wt=Count, name = 'Well_Total') 

#cond <- read_csv(paste0('/Users/yalegenomecenter/Desktop/SwabSeq-master/pipeline/', run_id, '/conditions.csv'), guess_max=guess_max) 
cond <- read.csv(file="/Users/yalegenomecenter/Desktop/SwabSeq/pipeline/run-2/conditions.csv", header=TRUE, sep=",")

#bc.map <- read_csv('/Users/yalegenomecenter/Desktop/SwabSeq-master/data/barcode-map.csv') 
bc.map <- read.csv(file="/Users/yalegenomecenter/Desktop/SwabSeq/pipeline/run-2/bc-map.csv", header=TRUE, sep=",")


# Data visualization 
cond %>%
  distinct(bc_set) %>%
  inner_join(bc.map) %>%
  arrange(target)

counts %>%
  filter(Sample_ID == 'Plate1-A01') %>%
  distinct(Sample_ID, Centroid, Count) %>%
  mutate(bc_set = 'N1_S2_RPP30') %>%
  left_join(bc.map %>% rename(Centroid = sequence)) %>%
  head(n=10)

counts %>%
  distinct(Sample_ID, Centroid, Count)  %>%
  count(Sample_ID, wt=Count, name = 'Well_Total') %>%
  inner_join(cond) %>%
  separate(Sample_ID, into = c('Sample_Plate', 'Well'), sep = '-', remove=F) %>%
  mutate(
    Row = factor(str_sub(Well, 1, 1), levels = rev(LETTERS[1:16])),
    Col = str_sub(Well, 2)
  ) %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Well_Total))) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~paste(Sample_Plate, nCoV_amplicon, sep = ' - ')) +
  scale_fill_viridis_c(option = 'plasma')
  
well.total %>%
  separate(Sample_ID, into = c('Sample_Plate', 'Well'), sep = '-', remove=F) %>%
  mutate(
    Row = factor(str_sub(Well, 1, 1), levels = rev(LETTERS[1:16])),
    Col = str_sub(Well, 2)
  ) %>%
  inner_join(cond) %>%
  ggplot(aes(x=Col, y=Row, fill=lysate)) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~Sample_Plate)
  
# Adding explicit zeros to barcodes that drop out  

explicit.zeros <- function(df, bc.map) {
  # take only assays and targets from the current run
  # assumes df has been joined with condition sheet
  bc.map %>%
    filter(
      bc_set %in% unique(df$bc_set),
    ) %>%
    left_join(df, by = c('sequence', 'bc_set')) %>%
    replace_na(list(Count = 0))
}

df <- counts %>%
  select(-Centroid) %>%
  rename(sequence=barcode) %>% 
  inner_join(select(cond, Sample_ID, bc_set), by = 'Sample_ID') %>% 
  group_by(Sample_ID) %>%
  group_nest() %>%
  mutate(foo = future_map(data, ~explicit.zeros(.x, bc.map))) %>%
  select(-data) %>%
  unnest(foo) %>%
  inner_join(cond) %>%
  mutate(
    Row = factor(str_sub(Sample_Well, 1, 1), levels = rev(LETTERS)),
    Col = str_sub(Sample_Well, 2),
    expected_amplicon = if_else(nCoV_amplicon == 'N1', "N1 Expected", "S2 Expected")
  ) %>%
  select(-nCoV_amplicon)

# Calculating two different spike-in across the two different plates

df %>%
  filter(str_detect(amplicon, "spike")) %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Count+1))) +
  geom_raster() +
  coord_equal() +
  facet_grid(expected_amplicon ~ amplicon) +
  scale_fill_viridis_c(option = 'plasma')
  
# Expression Relative to Spike-in’s

df %>%
  ggplot(aes(x=Col, y=Row, fill=log10(RNA_copies+0.1))) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~expected_amplicon) +
  scale_fill_viridis_c()
  
# Filter out any of barcodes that aren’t expected for that condition (e.g. remove N1 reads from the S2 plate)

df.wide <- df %>%
  select(Sample_ID, Plate_ID, Row, Col, bc_set, lysate, expected_amplicon, RNA_origin, RNA_copies, amplicon, Count) %>%
  filter(amplicon == 'RPP30' | str_detect(expected_amplicon, str_sub(amplicon, end=2)))  %>%
  mutate(amplicon = case_when(amplicon == 'RPP30' ~ 'RPP30',
                              str_detect(amplicon, 'spike') ~ 'Spike',
                              TRUE ~ 'RNA')
  ) %>%
  spread(amplicon, Count)

# Null Distribution

df.wide %>%
  filter(RNA_copies == 0) %>%
  ggplot(aes(x=Col, y=Row, fill=lysate)) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~expected_amplicon)

nulls <- df.wide %>%
  filter(RNA_copies == 0) %>%
  select(-RNA_origin) %>%
  nest(null.df = c(-expected_amplicon, -lysate))

df.wide.nulls <- df.wide %>%
  filter(RNA_copies != 0) %>%
  nest(data = c(-expected_amplicon, -lysate, -RNA_origin)) %>%
  inner_join(nulls) %>%
  mutate(combo = map2(data, null.df, bind_rows)) %>%
  select(-data, -null.df) %>%
  unnest(combo)
  
# Detection Plots
# HEK293 Lysate

df.wide.nulls %>%
  filter(lysate == 'NP') %>%
  inner_join(well.total) %>%
  mutate(RNA_copies = if_else(RNA_copies == 0, 0.1, RNA_copies)) %>%
  ggplot(aes(x=RNA_copies, y=(RNA+1)/(Spike+1), group=RNA_copies)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(alpha=0.4, aes(color=log10(Well_Total))) +
  scale_x_log10(breaks = c(10^(-1:4)), labels = c(0,10^(0:4))) +
  scale_y_log10() +
  scale_color_viridis_c(option = 'plasma', direction = -1) +
  annotation_logticks() +
  facet_grid(expected_amplicon ~ RNA_origin)
  
# Simple Classifier

test.df <- df.wide.nulls %>%
  inner_join(well.total) %>%
  filter(
    lysate == 'HEK293',
    expected_amplicon == 'S2 Expected',
    RNA_origin == 'ATCC_RNA',
    Well_Total >= 1000
  )

test.df %>%
  inner_join(well.total) %>%
  mutate(RNA_copies = if_else(RNA_copies == 0, 0.1, RNA_copies)) %>%
  ggplot(aes(x=RNA_copies, y=(RNA+1)/(Spike+1), group=RNA_copies)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(alpha=0.4, aes(color=log10(Well_Total))) +
  scale_x_log10(breaks = c(10^(-1:4)), labels = c(0,10^(0:4))) +
  scale_y_log10() +
  scale_color_viridis_c(option = 'plasma', direction = -1) +
  annotation_logticks()