# QTR

## INSTALLATION

Assuming you have `git`, the code can be made available locally as usual with
```bash
git clone https://github.com/jinghuazhao/QTR
```

## NOTES

### P values from `lmer()`

This is illustrated with the documentation example,
```{r}
require(lme4)
f <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
s <- summary(f)
s
class(with(s,coefficients))
t <- with(s,coefficients)[,3]
p <- 2*(1-pnorm(abs(t)))
p
```
