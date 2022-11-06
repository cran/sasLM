library(sasLM)

f1 = yield ~ block + N*P*K
GLM(f1, npk)
REG(f1, npk)
EMS(f1, npk)
lr(f1, npk)
lr0(f1, npk)
aov1(f1, npk)
aov2(f1, npk)
aov3(f1, npk)

f1b = yield ~ block + N*P*K - 1
GLM(f1b, npk[-1, ])
REG(f1b, npk[-1, ])
EMS(f1b, npk[-1, ])
lr(f1b, npk[-1, ])
lr0(f1b, npk[-1, ])
aov1(f1b, npk[-1, ])
aov2(f1b, npk[-1, ])
aov3(f1b, npk[-1, ])

f2 = weight ~ feed
GLM(f2, chickwts)
REG(f2, chickwts)
EMS(f2, chickwts)
lr(f2, chickwts)
lr0(f2, chickwts)
aov1(f2, chickwts)
aov2(f2, chickwts)
aov3(f2, chickwts)

f3 = uptake ~ conc - 1
GLM(f3, CO2)
REG(f3, CO2)
EMS(f3, CO2)
lr(f3, CO2)
lr0(f3, CO2)
aov1(f3, CO2)
aov2(f3, CO2)
aov3(f3, CO2)

Coll(mpg ~ disp + hp + drat + wt + qsec, mtcars)
Coll(mpg ~ disp + hp + drat + wt + qsec - 1, mtcars)
