inhibitorsGCNN = [
    "BCRP_inhibitor",
    "BSEP_inhibitor",
    "MRP1_inhibitor",
    "OATP1B1_inhibitor",
    "OATP1B3_inhibitor",
    "PGP_inhibitor"
]
substratesGCNN = [
    "BCRP_substrate",
    "MRP1_substrate",
    "PGP_substrate"
]
proteinsGCNN = inhibitorsGCNN + substratesGCNN


inhibitorsSA = [
    "MATE1_inhibitor",
    "MATE2_inhibitor",
    "OAT1_inhibitor",
    "OAT3_inhibitor",
    "OCT2_inhibitor",
    "MRP2_inhibitor"

]
substratesSA = [
    "BSEP_substrate",
    "MATE1_substrate",
    "MATE2_substrate",
    "OAT1_substrate",
    "OAT3_substrate",
    "OATP1B1_substrate",
    "OATP1B3_substrate",
    "OCT2_substrate",
    "MRP2_substrate"
]
proteinsSA = inhibitorsSA + substratesSA


APPLICABILITY_DOMAIN_THRESHOLD = 0.2
#ROUND_DIGIT = 2
SA_THRESHOLD = 0.27