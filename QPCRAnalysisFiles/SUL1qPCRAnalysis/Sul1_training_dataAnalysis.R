#install.packages("tidyverse")#
#library(tidyverse)#
qpcr_data <- read.csv("Sul1_training_dataExtracted.csv", skip = 7, stringsAsFactors = FALSE)
head(qpcr_data)
##Separating "well" column into row and column##
tidy_data2 <- separate(qpcr_data, Well, into = c("row", "column"), 
                       sep = 1, convert = TRUE)
view(tidy_data2)
##Representing the data in form of the well plate##
ggplot(tidy_data2, aes(x = column, y = row, label = Sample.Name)) +
  geom_tile(colour = "black") +
  geom_text()
ggplot(tidy_data2, aes(x = column, y = row, fill = Task, label = Sample.Name)) +
  geom_tile(colour = "black") +
#geom_text(aes(size =x)) +#
  geom_text() +
  scale_y_discrete(limits = c("H", "G", "F", "E", "D", "C", "B", "A")) +
  scale_x_continuous(position = "bottom", guide = guide_axis(position = "top"), breaks = 1:12)
view(tidy_data2)
##Calculating average CT Values## 
view(tidy_data2)
tidy_data2 <- filter(tidy_data2, Task != "NTC", Task != "")
view(tidy_data2)
head(tidy_data2$Ct)
class(tidy_data2$Ct)
summarised_data <- tidy_data2 %>%
  mutate(Ct = as.numeric(Ct)) %>%
  group_by(Sample.Name) %>%
  summarise(mean_Ct = mean(Ct))
head(summarised_data)
##Calculating average gene concentration Values##
class(tidy_data2$Quantity)
summarised_data <- tidy_data2 %>%
  group_by(Sample.Name) %>%
  summarise(mean_Quantity = mean(Quantity,na.rm = TRUE))
##Calculating average gene concentration Values##
tidy_data2$z <- NA
view(tidy_data2)
tidy_data2$z <- ((tidy_data2$Quantity.Mean/10) * (10/3) * 100* (10/25))
view(tidy_data2)
colnames(tidy_data2)[colnames(tidy_data2) == "z"] <- "Absolute Abundance (GC/mL)"
view(tidy_data2)
