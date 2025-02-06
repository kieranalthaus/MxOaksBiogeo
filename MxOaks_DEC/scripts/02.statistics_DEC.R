library(RevGadgets)

## Load trace files ------------------------------------------------------------
m0 = readTrace(paths = "MX_oaks_DEC/OUT/20240210_out_6states/m0/m0.params.log")
m1 = readTrace(paths = "MX_oaks_DEC/OUT/20240210_out_6states/m1/m1.params.log")

## Look at posterior
RevGadgets::summarizeTrace(trace = m0,
                           vars = "Posterior")
RevGadgets::summarizeTrace(trace = m1,
                           vars = "Posterior")