## preapare HI titer data
## filter from influenza subtype and collection year

HI_titer_raw<-read.csv("uga_titers_2010_onward_H3.csv")
HI_titer_raw %>% 
  filter(ag_ha_type=="H3")%>%
  mutate(ag_collection_date2=as_date(ag_collection_date,format="%m/%d/%y"))%>%
  filter(ag_collection_date2>="2010-9-1" & ag_collection_date2 <= "2020-5-1")-> uga_titer_2020_H3
  

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
  spread(key = "sr_strain_name", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->test_matrix2


uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  group_by(ag_strain_name)%>%
  summarise(count=n())-> HI_strain_count





uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  group_by(sr_strain_name)%>%
  summarise(count=n())-> sr_strain_count

## get strain name 

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,.keep_all= TRUE)%>%
  mutate(ag_strain_taxa=paste0(ag_epi_isolate_id,"|",ag_strain_name))->uga_titer_2020_H3

write.table(uga_titer_2020_H3$ag_strain_taxa,"HI_titer_strain.txt",sep="\t",row.names=FALSE)

write.csv(test_matrix2,"2010_2020_raw_titer.csv",row.names = FALSE)

################################################################
## filter the sequence with high accuracy