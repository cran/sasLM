require(sasLM)

f1 = yield ~ block + N*P*K
GLM(f1, npk)
REG(f1, npk)
ANOVA(f1, npk)
T3MS(f1, npk)

f2 = weight ~ feed 
GLM(f2, chickwts)
REG(f2, chickwts)
ANOVA(f2, chickwts)
T3MS(f2, chickwts)
