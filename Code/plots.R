# together with 12_06_2018 file
## plots of PFOA scatter, top cpg female
p1
p2
grid.newpage()
grid.draw(g)

## test smooth, colour, title font size. 
p1 <-  ggplot(data = pfoa_M_Beta_female, mapping = aes(x = PFOA_ng_ml, y = M)) + 
  geom_point(colour = "gray15", size = 2, shape = 20) +
  geom_smooth(method = "loess", se = FALSE, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "tomato") +
  xlab("PFOA Concentration ng/ml") +
  ylab("M Values of cg19425295 (Female)") +
  theme_bw() +
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.x = element_text(face="bold", color="black", 
                               size=10, angle=0),
    axis.text.y = element_text(face="bold", color="black", 
                               size=10, angle=0)
  )

p1

p2 <-  ggplot(data = pfoa_M_Beta_female, mapping = aes(x = PFOA_ng_ml, y = Beta)) + 
  geom_point() +
  xlab("PFOA Concentration ng/ml") +
  geom_smooth(method = "loess", se = FALSE, colour = "steelblue") +
  ylab("% Methylation of cg19425295 (Female)") +
  theme_bw()
p1
