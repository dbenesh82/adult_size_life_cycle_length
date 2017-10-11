## ---- message=FALSE, warning=FALSE---------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

options(stringsAsFactors = FALSE) # never have species lists as factors; always as character vectors

# set theme
theme.o <- theme_update(axis.text = element_text(colour="black", size = 15),
                        axis.title = element_text(colour="black", size = 18, face = "bold", lineheight=0.25),
                        axis.ticks = element_line(colour="black"),
                        panel.border = element_rect(colour = "black",fill=NA),
                        panel.grid.minor=element_blank(),
                        panel.grid.major=element_line(color="gray",linetype = "dotted"),
                        panel.background= element_rect(fill = NA))

# import the life cycle database data tables
dataH <- read.csv(file="data/CLC_database_hosts.csv", header = TRUE, sep=",")
dataL <- read.csv(file="data/CLC_database_lifehistory.csv", header = TRUE, sep=",")

## ---- message=FALSE, warning=FALSE---------------------------------------
# filter to adult stages
dataL <- filter(dataL, Stage == 'adult', (Sex == 'f' | is.na(Sex)) ) # remove adult males

## ---- message=FALSE, warning=FALSE---------------------------------------
dataL <- mutate(dataL, biovolume = 
                  if_else(Shape %in% c("cylinder", "thread-like", "whip"), 
                          pi * (Width/2)^2 * Length, # calculate volume as a cylinder
                          if_else(Shape %in% c("coiled", "sphere", "ellipsoid"),
                                  4/3 * pi * Length/2 * Width/4, # calculate volume as a ellipsoid
                                  Length * Width # calculate volume as area for remaining ribbon, leaf shapes
                                  )),
                biovolume = biovolume * 1.1) # covert to biomass with assumed 1.1. g/cm3 tissue density 

# species averages
dataL.sp <- group_by(dataL, Parasite.species)%>%
  summarize(Biovolume = mean(biovolume, na.rm=T))

## ---- message=FALSE, warning=FALSE---------------------------------------
maxLCL <- group_by(dataH, Parasite.species)%>%summarize(maxLCL = max(Host.no))
minLCL <- filter(dataH, Facultative == "no")%>%
  group_by(Parasite.species)%>%summarise(minLCL = length(unique(Host.no)))
dataL.sp <- left_join(dataL.sp, maxLCL)
dataL.sp <- left_join(dataL.sp, minLCL)
rm(minLCL, maxLCL)

## ---- message=FALSE, warning=FALSE---------------------------------------
#let's try a scatter plot first
ggplot(data = dataL.sp,
       aes(x = maxLCL, y = log10(Biovolume))) +
  geom_point(alpha = 0.3, position = position_jitter(height = 0, width = 0.1)) +
  geom_smooth(se = F, color = 'darkgrey') +
  labs(y="Log(Adult biomass)\n", x="\nLife cycle length") 

## ---- message=FALSE, warning=FALSE---------------------------------------
mdl1 <- (lm(log10(Biovolume) ~ maxLCL, data = dataL.sp))
summary(mdl1)

## ---- message=FALSE, warning=FALSE---------------------------------------
# make life cycle length a factor and pool the few parasites with life cycles longer than three hosts
dataL.sp <- mutate(dataL.sp, maxLCL.fac = if_else(maxLCL > 3, "4", as.character(maxLCL)))%>%
  mutate(maxLCL.fac = factor(maxLCL.fac, labels = c("1", "2", "3", ">3")))

## ---- message=FALSE, warning=FALSE---------------------------------------

# boxplot with data points as overlay
outfig <- ggplot(data = dataL.sp,
                 aes(x = maxLCL.fac, y = log10(Biovolume))) + 
  geom_boxplot(outlier.color = "white", width = 0.9) +
  geom_jitter(width = 0.2, color = "red", alpha = 0.2) +
  labs(y="Log(Adult biomass)\n", x="\nLife cycle length") + 
  theme(panel.grid.major.x = element_blank())
outfig

# save svg and png of fig for a word doc
ggsave(filename = "figs/adultsize_vs_lcl.png", width = 5, height = 4.5, units = "in")
ggsave(filename = "figs/adultsize_vs_lcl.svg", width = 5, height = 4.5, units = "in")

## ---- message=FALSE, warning=FALSE---------------------------------------
mdl2 <- lm(log10(Biovolume) ~ maxLCL.fac, data = dataL.sp)
anova(mdl1, mdl2)

## ---- message=FALSE, warning=FALSE---------------------------------------
# host mass dataset
host.size <- read.csv(file="data/collated_host_mass_data.csv", header = TRUE, sep=",")

# host size either a length or a mass; restrict to just masses
host.size <- filter(host.size, !is.na(body.mass))%>%
  select(binomial, body.mass)%>%
  group_by(binomial)%>%
  summarize(body.mass = mean(body.mass)) # some species have multiple entries, so avg them

# filter host list to just adult stages
dataH <- filter(dataH, Def.int == 'def', Typical.host == 'typical') # just typical definitive hosts

# add host mass data to host species
dataH <- left_join(dataH, host.size, by = c("Host.species" = "binomial"))

# average host mass for each parasite species
dataH <- group_by(dataH, Parasite.species)%>%
  summarise(host.mass = mean(body.mass, na.rm=T))

# add host mass to life cycle data
dataL.sp <- left_join(dataL.sp, dataH)

## ---- message=FALSE, warning=FALSE---------------------------------------
outfig <- ggplot(dataL.sp,
                 aes(x = maxLCL.fac, y = log10(host.mass))) + 
  geom_boxplot(outlier.color = "white", width = 0.9) +
  geom_jitter(width = 0.2, color = "red", alpha = 0.2) +
  labs(y="Log(Definitive host mass)\n", x="\nLife cycle length") + 
  theme(panel.grid.major.x = element_blank())
outfig

# save svg and png of fig for a word doc
ggsave(filename = "figs/defhostmass_vs_lcl.png", width = 5, height = 4.5, units = "in")
ggsave(filename = "figs/defhostmass_vs_lcl.svg", width = 5, height = 4.5, units = "in")

## ---- message=FALSE, warning=FALSE---------------------------------------
ggplot(dataL.sp,
       aes(x = log10(host.mass), y = log10(Biovolume))) + 
  geom_point(alpha = 0.5) +
  geom_smooth(se = F, color = 'darkgrey') +
  labs(x = "Log(Definitive host mass)", y = "Log(Adult worm size)")

## ---- message=FALSE, warning=FALSE---------------------------------------
# trophic level data from my tropic vacuum study
host.tl <- read.csv(file = "data/TV_dryad.csv", header = TRUE, sep = ",")

# reduce to needed variables
host.tl <- select(host.tl, Parasite.species = Species, Min_TL, Max_TL, Avg_TL)

#add host trophic levels to life history data
dataL.sp <- left_join(dataL.sp, host.tl)

## ---- message=FALSE, warning=FALSE---------------------------------------
ggplot(dataL.sp, aes(x = Avg_TL, y = log10(host.mass))) + 
  geom_point(alpha = 0.4, position = position_jitter(width = 0.03, height = 0)) +
  geom_smooth(se = F, color = 'darkgray') +
  labs(x = "Definitive host trophic level", y = "Log(Definitive host mass)")

## ---- message=FALSE, warning=FALSE---------------------------------------
#correlation between def host TL and worm size
ggplot(dataL.sp, aes(x = Avg_TL, y = log10(Biovolume))) + 
  geom_point(alpha = 0.4, position = position_jitter(width = 0.03, height = 0)) +
  geom_smooth(se = F, color = 'darkgray') +
  labs(x = "Definitive host trophic level", y = "Log(Adult worm biovolume)")

## ---- message=FALSE, warning=FALSE---------------------------------------
summary(lm(log10(Biovolume) ~ Avg_TL, data = dataL.sp)) # sig but r2 of 2%

## ---- message=FALSE, warning=FALSE---------------------------------------
library(MuMIn) # load library for mult-model inference

data.nem <- na.omit(dataL.sp) # restrict to complete data, i.e. just nematodes

global.mod <- lm(log10(Biovolume) ~ log10(host.mass) + Avg_TL + maxLCL +
                 I(Avg_TL^2) + I(log10(host.mass)^2) + maxLCL.fac,
                 data = data.nem, na.action = "na.fail") # global model

## ---- message=FALSE, warning=FALSE---------------------------------------
# dredge the model set
model.set <- dredge(global.mod, subset = 
                      dc(log10(host.mass), I(log10(host.mass)^2)) && # only include quadratic if linear included
                      dc(Avg_TL, I(Avg_TL^2)) && # only include quadratic if linear included
                      xor(maxLCL, maxLCL.fac)) # do not include continuous and categorical lcl var in same model

subset(model.set, cumsum(weight) <= .95) # 95% ci set

