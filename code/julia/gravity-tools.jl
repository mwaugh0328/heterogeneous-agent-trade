using FixedEffectModels

function gravity(tradedata)
    #function to perform basic gravity regression
    #assumes tradedata is a DataFrame, takes on the structure of EK dataset

    outreg = reg(tradedata, @formula(trade ~ fe(importer) + fe(exporter) +
         bin375 + bin750 + bin1500 + bin3000 + bin6000 + binmax  + border + 
                sharedlanguage + europeancom + efta + 0), save = true)

    norm_importer_fe = sum(unique(outreg.fe.fe_importer)) / 19

    dist_bins = outreg.coef[1:6] .+ norm_importer_fe

    return dist_bins

end