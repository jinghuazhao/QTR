# QTR

## INSTALLATION

The code can be made available locally as usual with
```bash
git clone https://github.com/jinghuazhao/QTR
```

## NOTES

### P values from Linear mixed models -- documentation example
```{r}
require(lme4)
f <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
s <- summary(f)
s
t <- with(s,coefficients)[,3]
class(with(s,coefficients))
p <- 2*(1-pnorm(abs(t)))
p
```
