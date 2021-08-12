# Boxplot script

source("src/load_packages.R")

df.interactions <-
  read_excel("files/caltech.PD mice.040820.xlsx", 
             sheet = "interaction.regcoef") %>% 
  janitor::clean_names() %>%
  dplyr::filter(imputation == "imputed") %>% 
  filter(variable == "MicrobiotaSPF:GenotypeASO") 


# COlon Betaine
df.meta <-
  read_excel("files/Metadata_Metabolomics 2019 _checked 2020_updated.xlsx", 
             sheet = "Sheet2") %>% 
  select(Animal, group) %>% 
  mutate(Animal = as.character(Animal))

df.raw <-
  read_excel("files/5370_Cal Tech PD Mouse Model_Gut brain axis study_09112020.xlsx", 
             sheet = "q500 data")
plt <- 
  df.meta %>%
  left_join(df.raw, by = "Animal") %>% 
  filter(`Additional Information` == "colon") %>% 
  ggplot(aes(x=group, y=Betaine)) +
  geom_boxplot() +
  geom_point()
  
df.meta %>%
  left_join(df.raw, by = "Animal") %>% 
  filter(`Additional Information` == "duodenum content") %>%
  ggplot(aes(x=group, y=`TG(17:1_38:5)`)) +
  geom_boxplot() +
  geom_point()
