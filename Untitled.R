## preapare HI titer data
## filter from influenza subtype and collection year

HI_titer_raw<-read.csv("uga_titers_2010_onward_H3.csv")
HI_titer_raw %>% 
  filter(ag_ha_type=="H3")%>%
  mutate(ag_collection_date2=as_date(ag_collection_date,format="%m/%d/%y"))%>%
  filter(ag_collection_date2>="2010-9-1" & ag_collection_date2 <= "2020-5-1")-> uga_titer_2020_H3

## further filter by sr assay 

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  group_by(ag_strain_name)%>%
  summarise(count=n())%>%
  ggplot(aes(x=count))+
  geom_histogram()+
  scale_x_continuous( breaks = seq(0, 50, 2),
                      limits = c(0,50))
  


uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  group_by(ag_strain_name)%>%
  summarise(sr_count=n())->ag_srcount


%>%
  merge(uga_titer_2020_H3, by=c("ag_strain_name","ag_strain_name"))%>%
  distinct(ag_strain_name,.keep_all= TRUE)->test2

%>%
  mutate(ag_collection_year=year(ag_collection_date2))%>%
  arrange(sr_count)%>%
  group_by(ag_collection_year)%>%
  slice(1:50)->filter_seq



duplicated(uga_titer_2020_H3, by=c("ag_strain_name", "sr_strain_name"))->test
uga_titer_2020_H3[test,]

subset(uga_titer_2020_H3,duplicated(by=c("ag_strain_name","sr_strain_name")))-> duplicate_records

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  #select(ag_strain_name,sr_strain_name,titer_value)%>%
  #distinct(ag_strain_name,sr_strain_name)%>%
  group_by(ag_strain_name,sr_strain_name,ag_entry_id,sr_entry_id)%>%
  mutate(dupe = n()>1)%>%
  filter(dupe=="TRUE")->test2


uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  mutate(ag_strain_name2=paste0(ag_strain_name,"_",ag_entry_id))%>%
  mutate(sr_strain_name2=paste0(sr_strain_name,"_",sr_entry_id))%>%
  select(ag_strain_name2,sr_strain_name2,titer_value)%>%
  spread(key = "sr_strain_name2", value = "titer_value")->test_matrix



uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  spread(key = "sr_strain_name", value = "titer_value")->test_matrix2
