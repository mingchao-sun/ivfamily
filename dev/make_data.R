# Builds the example data files from the AER package (GPL). Run from the module root: Rscript dev/make_data.R
suppressMessages(library(AER))
data("CigarettesSW")
cg <- subset(CigarettesSW, year == "1995")
out <- data.frame(state = as.character(cg$state),
                  lq = log(cg$packs),
                  lp = log(cg$price / cg$cpi),
                  li = log(cg$income / cg$population / cg$cpi),
                  tdiff = (cg$taxs - cg$tax) / cg$cpi,
                  tax_r = cg$tax / cg$cpi)
write.csv(out, "data/cigarettes_sw.csv", row.names = FALSE)
data("CollegeDistance")
cd <- CollegeDistance[, c("wage", "education", "distance", "tuition", "gender", "ethnicity", "unemp")]
cd$gender <- as.character(cd$gender); cd$ethnicity <- as.character(cd$ethnicity)
write.csv(cd, "data/college_distance.csv", row.names = FALSE)
cat("cigarettes_sw.csv:", nrow(out), "rows; college_distance.csv:", nrow(cd), "rows\n")

# Angrist and Evans (1998): family size and labour supply. Binary treatment morekids (more than two children), binary
# instrument samesex (the first two children have the same sex), outcome weeks (weeks worked last year).
data("Fertility")
ft <- Fertility
out <- data.frame(weeks = ft$work, morekids = as.integer(ft$morekids == "yes"),
                  samesex = as.integer(ft$gender1 == ft$gender2), age = ft$age, boy1 = as.integer(ft$gender1 == "male"),
                  afam = as.integer(ft$afam == "yes"), hispanic = as.integer(ft$hispanic == "yes"), other = as.integer(ft$other == "yes"))
write.csv(out, "data/fertility.csv", row.names = FALSE)
cat("fertility.csv:", nrow(out), "rows\n")
