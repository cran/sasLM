require(sasLM)

f1 = yield ~ block + N*P*K
GLM(f1, npk)
REG(f1, npk)
ANOVA(f1, npk)
EMS(f1, npk)

f2 = weight ~ feed 
GLM(f2, chickwts)
REG(f2, chickwts)
ANOVA(f2, chickwts)
EMS(f2, chickwts)

m = mtcars
m[1, "mpg"] = NA
m[2, "disp"] = NA
m[3, "hp"] = NA
m[4, "drat"] = NA
m[5, "qsec"] = NA
m[6, "wt"] = NA

Cor.test(m)
Pcor.test(m, c("mpg", "hp", "disp", "qsec"), c("drat", "wt"))

tsum(lh)
t(tsum(CO2))
t(tsum(uptake ~ Treatment, CO2))
tsum(uptake ~ Type + Treatment, CO2)
tsum(uptake ~ conc + Type + Treatment, CO2)
