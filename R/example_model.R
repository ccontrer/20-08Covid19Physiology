library(brms)
dat <- read.csv("dat.csv")

###############################
###### define the model #######
###############################


b_mod1 <- bf(recall ~ betaMu + (alphaMu - betaMu) * exp(-exp(gammaMu) * trial),
             betaMu ~ 1 + (1 | c | subject),
             alphaMu ~ 1 + (1 | c | subject),
             gammaMu ~ 1 + (1 | c | subject),
             nl = TRUE) +
  nlf(sigma ~ betaSc + (alphaSc - betaSc) * exp(-exp(gammaSc) * trial),
      alphaSc ~ 1 + (1 | c | subject),
      betaSc ~ 1 + (1 | c | subject),
      gammaSc ~ 1 + (1 | c | subject))

#################################
###### define priors ############
#################################
prior1 <- c(set_prior("normal(25, 2)", nlpar = "betaMu"),
            set_prior("normal(8, 2)",  nlpar = "alphaMu"),
            set_prior("normal(-1, 2)",  nlpar = "gammaMu"),
            set_prior("normal(0, 1)", nlpar = "betaSc"),
            set_prior("normal(1, 1)",  nlpar = "alphaSc"),
            set_prior("normal(0, 1)",  nlpar = "gammaSc"),
            set_prior("student_t(3, 0, 5)", class = "sd", nlpar = "alphaMu"),
            set_prior("student_t(3, 0, 5)", class = "sd", nlpar = "betaMu"),
            set_prior("student_t(3, 0, 5)", class = "sd", nlpar = "gammaMu"),
            set_prior("student_t(3, 0, 1)", class = "sd", nlpar = "alphaSc"),
            set_prior("student_t(3, 0, 1)", class = "sd", nlpar = "betaSc"),
            set_prior("student_t(3, 0, 1)", class = "sd", nlpar = "gammaSc"))



# fit the model started 19:40, ended 23:10
b_fit2 <- brm(b_mod1, prior = prior1, data = dat,
              control = list(adapt_delta = .9999, max_treedepth = 12),
              chains = 1, iter = 100,
              inits = 0, cores = 4)

summary(b_fit1)

# we encourge the exploration of a variety of non-linear functions.
