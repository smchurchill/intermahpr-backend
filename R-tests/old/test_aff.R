
z <- RInterMAHP::rr_zhao
z$extrapolation = TRUE

tz <- data.frame(t(z))
flist <- lapply(tz,rr_chooser)

pts <- seq(1,250,0.5)

f <- flist[[1]][[1]](pts)
f

plot(pts,flist[[52]][[1]](pts))



compute_interval_aaf <- function(
    R1, R2,         # Real numbers, ratios
    P, R, B,        # Functions, Reals to Reals
    bb, lb, ub      # Real numbers, integral bounds
) {

    integrand1 <- function(x) {P(x) * (R(x) - 1)}
    integrand2 <- function(x) {P(x) * (B(x) - 1)}
    integral <- -1
    if((lb < bb)&&(bb < ub)){
        integral =
            (R1      * integrate(integrand1,lb,bb)$value) +
            (R2      * integrate(integrand2,lb,bb)$value) +
            integrate(integrand2,bb,ub)$value
    }
    else {
        integral =
            (R1      * integrate(integrand1,lb,ub)$value) +
            (R2      * integrate(integrand2,lb,ub)$value)
    }
    return(integral)
}



row30 = z[30,]
row30
betas30 = data.matrix(row30[8:23])
betas30

ihd_f_morb = ihd_rr(betas30,1)

smallpts = seq(0.03,60.03,0.1)

plot(smallpts, ihd_f_morb(smallpts))

numbers = 1:250
cut = ((0<numbers)&(numbers<=150))
print(cut*numbers)

row1 = RInterMAHP::rr_zhao[1,]
row1
row42 = z[42,]
betas = data.matrix(row1[8:23])
betas

simple <- simple_rr(betas)

s <- simple(1:250)
s
plot(s)
z = RInterMAHP::rr_zhao


nc = integrate(function(x) dgamma(x,shape=0.729266,scale=48.24401),0.03,250)$value
p_cd = 0.833608
p_bd = 0.443621
binge = 68.3
print(nc)

gamma_m_bc_2015_young <- function(x) {
    return(dgamma(x,shape=0.729266,scale=48.24401)*p_cd/nc)
}

gamma_m_bc_2015_young(13)

gamma_m_bc_2015_young(200)


p_bat = integrate(gamma_m_bc_2015_young,binge,250)$value

R1t = (p_cd - p_bd) / (p_cd - p_bat)
R2t = (p_bd - p_bat) / (p_cd - p_bat)

R1t
R2t

lb = 0.03
ub = 250

rr_male_cirrhosis_morb <- simple_rr(data.matrix(row42[8:23]))

rr_male_cirrhosis_morb(13.5)
rr_male_cirrhosis_morb(40)
rr_male_cirrhosis_morb(80)
rr_male_cirrhosis_morb(200)

v <- rr_male_cirrhosis_morb(c(10,40,80,200))
v

integrand1 <- function(x) {gamma_m_bc_2015_young(x) * (rr_male_cirrhosis_morb(x) - 1)}
integrand2 <- function(x) {gamma_m_bc_2015_young(x) * (rr_male_cirrhosis_morb(x) - 1)}


g1 <- integrate(integrand1,lb,binge)$value
g1
g2 <- integrate(integrand2,lb,binge)$value
g2
g3 <- integrate(integrand2,binge,ub)$value
g3

i1 <- R1t*g1
i2 <- R2t*g2

i1
i2

s <- i1 + i2 + g3
s

aaf_numerator <- compute_interval_aaf(
    R1t, R2t,
    gamma_m_bc_2015_young, rr_male_cirrhosis_morb, rr_male_cirrhosis_morb,
    binge, 20.2, 40.4
)

aaf_denominator <- compute_interval_aaf(
    R1t, R2t,
    gamma_m_bc_2015_young, rr_male_cirrhosis_morb, rr_male_cirrhosis_morb,
    binge, lb, ub
)

aaf_numerator

p_fd = 0.062952
rr_fd = 3.26
aaf_fd_numerator = p_fd * (rr_fd-1)
aaf_fd_numerator

aaf = (aaf_numerator) / (1 + aaf_fd_numerator + aaf_denominator)

aaf

