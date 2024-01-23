plot_richness(Workshop, measures=c("Observed", "Shannon"), x="Site") +
  geom_boxplot() +
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  stat_compare_means(method = "anova")
ggsave("Alpha.png",height =6, width =12)

  
