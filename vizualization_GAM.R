library(gratia)
library(dplyr)
library(ggplot2)

## 1) Smooth estimates for s(Altitude_scaled) ----
sm <- smooth_estimates(
  mod_gam1,
  select = "s(Altitude_scaled)",  # <— use `select`, not `smooth`
  n = 200,
  unconditional = TRUE
)

## 2) Get inverse link and intercept ----
ilink <- family(mod_gam1)$linkinv
b0    <- coef(mod_gam1)[["(Intercept)"]]

## 3) Build both LINK- and RESPONSE-scale summaries ----
sm <- sm %>%
  mutate(
    eta   = b0 + estimate,                 # add intercept to the smooth
    fit   = ilink(eta),                    # response-scale fit
    lower = ilink(b0 + estimate - 1.96 * se),
    upper = ilink(b0 + estimate + 1.96 * se)
  )

## 4) Partial residuals for this smooth ----
pres <- partial_residuals(mod_gam1, select = "s(Altitude_scaled)") %>%
  mutate(
    eta_partial = b0 + .partial,     # add intercept
    y_partial   = ilink(eta_partial) # response-scale partial residuals
  )

## 5A) Plot on LINK scale (canonical for GAM smooths)
p_link <- ggplot() +
  geom_ribbon(data = sm,
              aes(x = Altitude_scaled, ymin = est - 1.96*se, ymax = est + 1.96*se),
              alpha = 0.35) +
  geom_line(data = sm, aes(x = Altitude_scaled, y = est), linewidth = 1.1) +
  geom_jitter(data = pres, aes(x = Altitude_scaled, y = .partial),
              width = 0.03, height = 0, alpha = 0.6, size = 1.6) +
  labs(x = "Altitude (scaled)", y = "Effect on link scale",
       subtitle = "gratia: s(Altitude_scaled) with partial residuals (link)") +
  theme_minimal(base_size = 12)

## 5B) Plot on RESPONSE scale (adds intercept, inverse link)
p_resp <- ggplot() +
  geom_ribbon(data = sm,
              aes(x = Altitude_scaled, ymin = pmax(0, lower), ymax = pmin(1, upper)),
              alpha = 0.35) +
  geom_line(data = sm, aes(x = Altitude_scaled, y = fit), linewidth = 1.1) +
  geom_jitter(data = pres, aes(x = Altitude_scaled, y = y_partial),
              width = 0.03, height = 0, alpha = 0.6, size = 1.6) +
  labs(x = "Altitude (scaled)", y = "Wings CWM (response scale)",
       subtitle = "gratia: s(Altitude_scaled) partial effect (response)") +
  theme_minimal(base_size = 12)

p_link
p_resp