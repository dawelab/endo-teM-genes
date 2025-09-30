#endosperm tem genes
endo_tem_var<-mat_pref_var%>%
  filter(b73_gene%in%endosperm_teM_genes$V1)%>%
  group_by(b73_gene)%>%
  summarise(mean=mean(maternal_preference,na.rm=T))%>%
  mutate(group="endo teM genes")

#FMEGs
fmegs_var<-mat_pref_var%>%
  filter(b73_gene%in%MEG_FMEGs$x)%>%
  group_by(b73_gene)%>%
  summarise(mean=mean(maternal_preference,na.rm=T))%>%
  mutate(group="FMEG")

#core genes
coregene_var<-mat_pref_var%>%
  filter(b73_gene%in%core_gene$x)%>%
  group_by(b73_gene)%>%
  summarise(mean=mean(maternal_preference,na.rm=T))%>%
  mutate(group="Core Gene")


combined_var<-rbind(bmegs_var,fmegs_var,coregene_var)%>%
  mutate(group=factor(group,levels=c("BMEG","FMEG","Core Gene")))


#Fig 4A
ggplot(combined_var, aes(x = group, y = mean,color=group)) +
  geom_violin() +
  #geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  theme_minimal() +
  theme(text = element_text(size = 18))+
  labs(title = "Maternal expression  preference  ",
       x = "",
       y = "Maternal Expression Preference") +
  geom_hline(yintercept = 0.86, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.67, linetype = "dashed", color = "black")+
  geom_jitter(data = combined_var %>% filter(group %in% c("BMEG", "FMEG")),
              width = 0.2, alpha = 0.6, size = 1.5)



#Fig 4C
