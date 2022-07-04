## preapare HI titer data
## filter from influenza subtype and collection year

HI_titer_raw<-read.csv("uga_titers_2010_onward_H3.csv")
HI_titer_raw %>% 
  filter(ag_ha_type=="H3")%>%
  mutate(ag_collection_date2=as_date(ag_collection_date,format="%m/%d/%y"))%>%
  filter(ag_collection_date2>="2010-9-1" & ag_collection_date2 <= "2020-5-1")-> uga_titer_2020_H3


uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  group_by(ag_strain_name)%>%
  summarise(sr_count=n())->ag_srcount
  



uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  spread(key = "sr_strain_name", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->test_matrix2


uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  group_by(sr_strain_name)%>%
  summarise(count=n())->sr_strain_count





uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  group_by(ag_strain_name)%>%
  summarise(count=n())%>%
  ggplot(aes(x=count))+
  geom_histogram()+
  scale_x_continuous(limits = c(0,50,2))

## get strain name 

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,.keep_all= TRUE)%>%
  mutate(ag_strain_taxa=paste0(ag_epi_isolate_id,"|",ag_strain_name))->uga_titer_strain_2020_H3

write.table(uga_titer_2020_H3$ag_strain_taxa,"HI_titer_strain.txt",sep="\t",row.names=FALSE)

write.csv(test_matrix2,"2010_2020_raw_titer.csv",row.names = FALSE)

################################################################
## filter the sequence with high accuracy

uga_titer_strain_2020_H3%>%
  merge(ag_srcount,by=c("ag_strain_name","ag_strain_name"))%>%
  mutate(ag_collection_year=year(ag_collection_date2))%>%
  arrange(desc(titer_value,sr_count))%>%
  group_by(ag_collection_year)%>%
  slice(1:30)->uga_titer_strain_filteryear

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  filter(ag_strain_name %in% uga_titer_strain_filteryear$ag_strain_name)%>%
  
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  spread(key = "sr_strain_name", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->tes_matrix3

## get sequence to test


write.csv(tes_matrix3,"2010_2020_raw_titer_filter.csv",row.names = FALSE)
write.table(uga_titer_strain_filteryear$ag_strain_taxa,"HI_titer_yearfilter_strain.txt",sep="\t",row.names=FALSE)

