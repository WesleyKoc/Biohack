function toMmolPerL(value, unit, molarMass) {
    switch (unit) {
        case "g/L":
            return (value / molarMass) * 1000;
        case "% w/v":
            // % w/v = g per 100 mL = 10 g/L
            return ((value * 10) / molarMass) * 1000;
        case "mmol/L":
        default:
            return value;
    }
}

const MOLAR_MASSES = {
    // Carbon sources
    "corn steep liquor":        180, 
    "wheat bran hydrolysate":   162,  
    "cane molasses":            342,  
    "distillery stillage":      180,  
    "cheese whey (lactose)":    342,  
    "wheat straw hydrolysate":  150,   
    "potato waste":             162,   
    "glucose":                  180.16,
    "sucrose":                  342.30,
    "glycerol":                  92.09,
    "xylose":                   150.13,
    "fructose":                 180.16,
    // Nitrogen sources
    "soy protein hydrolysate":  119, 
    "distillers dried grains (DDGS)": 119,
    "rapeseed meal extract":    119,
    "ammonium sulphate":        132.14,
    "ammonium chloride":         53.49,
    "urea":                      60.06,
    "yeast extract":            119, 
    "peptone":                  119,
    // corn steep liquor can also serve as N source (same MW entry as C source)
};

// Carbon content (g C per g substrate) — used for C:N ratio calculation
const CARBON_FRACTION = {
    "corn steep liquor":         0.40,
    "wheat bran hydrolysate":    0.44,
    "cane molasses":             0.42,
    "distillery stillage":       0.40,
    "cheese whey (lactose)":     0.42,
    "wheat straw hydrolysate":   0.44,
    "potato waste":              0.44,
    "glucose":                   0.40,
    "sucrose":                   0.42,
    "glycerol":                  0.39,
    "xylose":                    0.40,
    "fructose":                  0.40,
};

// Nitrogen content (g N per g source)
const NITROGEN_FRACTION = {
    "corn steep liquor":         0.08,
    "soy protein hydrolysate":   0.16,
    "distillers dried grains (DDGS)": 0.13,
    "rapeseed meal extract":     0.12,
    "ammonium sulphate":         0.212,
    "ammonium chloride":         0.262,
    "urea":                      0.467,
    "yeast extract":             0.10,
    "peptone":                   0.14,
};

// Phosphate molar masses (g/mol)
const PHOSPHATE_MM = {
    "potassium dihydrogen phosphate":  136.09, 
    "dipotassium hydrogen phosphate":  174.18, 
    "na2hpo4":                         141.96,  
    "none":                            null,
};

/**
 * @param {number} cnRatio   Calculated C:N mass ratio
 * @returns {number}         Fraction between 0 and 1
 */
function cnRatioFactor(cnRatio) {
    const optimalLow  = 10;
    const optimalHigh = 25;
    if (cnRatio >= optimalLow && cnRatio <= optimalHigh) {
        return 1.0;
    }
    if (cnRatio < optimalLow) {
        // Excess nitrogen; mild penalty — microbe can always assimilate extra N
        return Math.max(0.5, cnRatio / optimalLow);
    }
    // Nitrogen limitation above the optimum window
    return Math.max(0, 1 - ((cnRatio - optimalHigh) / 100));
}

function nitrogenFluxBound(nConc_mmolL, cnFactor) {
    const BASE_N_UPTAKE = 5.0;
    return nConc_mmolL * BASE_N_UPTAKE * cnFactor / 1000;
}

function phosphateFluxBound(phosphateInput) {
    if (!phosphateInput || phosphateInput.source === "none" ||
        phosphateInput.concentration == null || isNaN(phosphateInput.concentration)) {
        return Infinity;   // leave exchange unconstrained
    }
    const mm = PHOSPHATE_MM[phosphateInput.source];
    if (!mm) return Infinity;
    const conc_mmolL = toMmolPerL(
        phosphateInput.concentration,
        phosphateInput.unit,
        mm
    );
    // Phosphate uptake kinetics: simple proportional bound (Vmax-like scaling)
    const BASE_P_UPTAKE = 3.0;
    return conc_mmolL * BASE_P_UPTAKE / 1000;
}

function organismsDefaults(organismName) {
    const name = (organismName || "").toLowerCase();

    if (/aspergill|penicill|trichoderm|fusari/.test(name)) {
        // Filamentous fungi: warm and mildly acidic
        return { criticalTemp: 30, phOptima: 6.0, maintenanceCoeff: 0.05 };
    }
    if (/yeast|saccharomyces|kluyveromyces|pichia|yarrowia/.test(name)) {
        // Yeasts: slightly higher optimum temperature
        return { criticalTemp: 32, phOptima: 5.5, maintenanceCoeff: 0.04 };
    }
    if (/bacill|corynebact|escherichia|e\.? ?coli/.test(name)) {
        // Common bacteria
        return { criticalTemp: 37, phOptima: 7.0, maintenanceCoeff: 0.03 };
    }
    if (/streptomyces|actinomycet/.test(name)) {
        return { criticalTemp: 28, phOptima: 7.2, maintenanceCoeff: 0.05 };
    }
    // Default (generic fungal, matches the original hard-coded values)
    return { criticalTemp: 30, phOptima: 6.2, maintenanceCoeff: 0.05 };
}

function range(start, stop, step) {
    const arr = [];
    for (let v = start; v <= stop; v += step) arr.push(parseFloat(v.toFixed(6)));
    return arr;
}

async function runScientificallyHardenedOptimization(modelXml, userInputs, bioMetadata = {}) {
    const solver = await glpk();
    const model  = escherFba.loadModel(modelXml);
    const MAX_KLA = 350; 

    const resolvedOrganismName = userInputs.organism
        || (bioMetadata.organism_type
                ? bioMetadata.organism_type.toLowerCase() 
                : "");
    const orgDefaults = organismsDefaults(resolvedOrganismName);

    const CRITICAL_TEMP = orgDefaults.criticalTemp;
    const PH_OPTIMA     = orgDefaults.phOptima;

    const MAINTENANCE_COEFFICIENT =
        (bioMetadata.maintenance_coefficient != null)
            ? bioMetadata.maintenance_coefficient / 30
            : orgDefaults.maintenanceCoeff;

    const maintenanceSource = (bioMetadata.maintenance_coefficient != null)
        ? "pipeline.py (BV-BRC)"
        : "heuristic fallback";
    console.log(
        `[FBA] organism="${resolvedOrganismName}" | ` +
        `maintenance_coefficient=${MAINTENANCE_COEFFICIENT} (source: ${maintenanceSource})`
    );

    if (Array.isArray(bioMetadata.required_supplements) && bioMetadata.required_supplements.length > 0) {
        console.warn(
            `[FBA] BV-BRC pipeline flagged required supplements for this organism: ` +
            `${bioMetadata.required_supplements.join(", ")}. ` +
            `Ensure these are present in your media formulation.`
        );
    }
    
    const cSource  = userInputs.carbonSource;
    const cMM      = MOLAR_MASSES[cSource.name] || 180;
    const cConc_mmolL = toMmolPerL(cSource.concentration, cSource.unit, cMM);
    const nSource  = userInputs.nitrogenSource;
    const nMM      = MOLAR_MASSES[nSource.name] || 119;
    const nConc_mmolL = toMmolPerL(nSource.concentration, nSource.unit, nMM);
    const cFraction = CARBON_FRACTION[cSource.name] || 0.40;
    const nFraction = NITROGEN_FRACTION[nSource.name] || 0.14;
    const cMass_gL  = (cConc_mmolL * cMM) / 1000;
    const nMass_gL  = (nConc_mmolL * nMM) / 1000;
    const carbonMass   = cMass_gL * cFraction;
    const nitrogenMass = nMass_gL * nFraction;
    const cnRatio = nitrogenMass > 0 ? carbonMass / nitrogenMass : Infinity;

    const cnFactor         = cnRatioFactor(cnRatio);
    const nFluxBound       = nitrogenFluxBound(nConc_mmolL, cnFactor);
    const phosphateBound   = phosphateFluxBound(userInputs.phosphate);

    const KS_CARBON          = 0.5;
    const VMAX_CARBON_UPTAKE = 10.0;  // mmol/g/h
    const cFluxBound = (VMAX_CARBON_UPTAKE * cConc_mmolL) / (KS_CARBON + cConc_mmolL);
    
    model.setUpperBound('EX_glc__D_e',  cFluxBound);     
    model.setUpperBound('EX_nh4_e',     nFluxBound);    
    if (isFinite(phosphateBound)) {
        model.setUpperBound('EX_pi_e',  phosphateBound); 
    }

    const inoculumFraction = (userInputs.inoculumPct > 0 && !isNaN(userInputs.inoculumPct))
        ? userInputs.inoculumPct / 100
        : 0.05;   // default 5 % v/v
    const fermentationHours = (userInputs.hours > 0 && !isNaN(userInputs.hours))
        ? userInputs.hours
        : 72;   // default 72 h

    let results = [];

    for (let temp of range(userInputs.temp.min, userInputs.temp.max, 2)) {
        for (let ph of range(userInputs.ph.min, userInputs.ph.max, 0.2)) {
            for (let rpm of range(userInputs.rpm.min, userInputs.rpm.max, 50)) {

                const tempFactor = Math.exp(-0.5 * Math.pow((temp - CRITICAL_TEMP) / 5, 2));
                const phFactor   = Math.exp(-0.5 * Math.pow((ph   - PH_OPTIMA)     / 1, 2));

                const o2Availability = (MAX_KLA * rpm) / (200 + rpm);
                model.setUpperBound('EX_o2_e', o2Availability);

                const fbaResult = escherFba.solve(model, solver);
                let rawGrowth   = fbaResult.objectiveValue;

                let netGrowth = (rawGrowth * tempFactor * phFactor * cnFactor)
                                - MAINTENANCE_COEFFICIENT;
                netGrowth = Math.max(0, netGrowth);

                // E. Industrial yield (uses inoculum-adjusted X0 and actual duration)
                const finalBiomass = calculateBiomassYield(
                    netGrowth,
                    userInputs.volume,
                    fermentationHours,
                    inoculumFraction
                );

                results.push({
                    temp,
                    ph,
                    rpm,
                    cnRatio:              cnRatio.toFixed(2),
                    cnFactor:             cnFactor.toFixed(3),
                    growthRate:           netGrowth.toFixed(4),
                    biomass_gL:           finalBiomass.toFixed(2),
                    o2_limited:           rawGrowth >= o2Availability * 0.95,
                    n_limited:            cnFactor < 0.8,
                    organism:             resolvedOrganismName || "unspecified",
                    tax_id:               bioMetadata.tax_id    ?? null,
                    pipeline_maintenance: bioMetadata.maintenance_coefficient ?? null,
                });
            }
        }
    }

    const optimum = results.reduce((prev, curr) =>
        parseFloat(curr.growthRate) > parseFloat(prev.growthRate) ? curr : prev
    );

    renderFullReport(results, optimum);
}

function calculateBiomassYield(mu, volume, time, inoculumFraction = 0.05) {
    // Starting concentration scales with inoculum fraction.
    // Assume the seed culture is at ≈ 10 g/L.
    const SEED_CONCENTRATION = 10; // g/L
    const X0 = inoculumFraction * SEED_CONCENTRATION;  // g/L
    return (X0 * Math.exp(mu * time)) * volume;         // g total
}
