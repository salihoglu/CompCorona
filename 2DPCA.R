library("factoextra")
library("FactoMineR")
library("dplyr")
library("plotly")

clean_data <- read.csv("all_type_for_pca.csv")


# Select relevant columns
clean <- clean_data %>%
  select(-c(gene_name)) %>%
  replace(is.na(.), 0)

# Set row names
rownames(clean) <- clean_data$gene_name

# PCA
pca.data <- PCA(clean[, -5], scale.unit = TRUE, ncp = 2, graph = FALSE)

# Plotly PCA Biplot
# color options "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
biplot <- fviz_pca_biplot(pca.data, label = "var", 
                          # Fill individuals by groups
                          geom.ind = "point",
                          pointshape = 16,
                          pointsize = 2.2,
                          fill.ind = clean$Type,
                          col.var = "darkblue",
                          palette = "ucscgb",
                          addEllipses = TRUE, ellipse.level = 0.95,
                          legend.title = list(fill = "Infections"))

# Convert ggplot to plotly
biplot_plotly <- ggplotly(biplot)

# Show the interactive plot
biplot_plotly








fviz_pca_biplot(pca.data, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                addEllipses = TRUE,
                fill.ind = clean$Type,
                col.ind = "black",
                # Color variable by groups
                col.var = factor(c("Sars-CoV","Mers-CoV","Sars-CoV-2","Pbmc-Sars-CoV-2")),
                
                legend.title = list(fill = "Type", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  
  ggpubr::fill_palette("rickandmorty") +    # Indiviual fill color
  ggpubr::color_palette("rickandmorty")      # Variable colors
