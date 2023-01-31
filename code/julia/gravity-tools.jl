using FixedEffectModels


function gravity!(tradedata, d, T, W, θ)
    #mulitple dispatch version, runs gravity regression
    # then fills matrix d with trade costs

    dist_coef, lang_coef, S, θm = gravity(tradedata)

    inv_θ = (1.0 / θ)

    for importer = 1:19

        foo = tradedata[tradedata.importer .== importer, :]

        T[importer] = exp( S[importer] + θ*log(W[importer]) )
        #equation (27) from EK, but with β = 1 (no round about)

        for exporter = 1:19

            if exporter != importer

                get_exporter = foo.exporter .== exporter

                distance_effect = exp(-inv_θ * dist_coef[Int(foo[get_exporter, :].distbin[1])]) 

                border_effect = exp(-inv_θ * lang_coef[1] * foo[get_exporter, :].border[1]) 

                language_effect = exp(-inv_θ * lang_coef[2] * foo[get_exporter, :].sharedlanguage[1]) 

                europeancom_effect = exp(-inv_θ * lang_coef[3] * foo[get_exporter, :].europeancom[1]) 

                efta_effect = exp(-inv_θ * lang_coef[3] * foo[get_exporter, :].efta[1])

                asym_effect = exp( -inv_θ * θm[importer] )           

                d[importer, exporter] =(distance_effect * border_effect * language_effect 
                                         * europeancom_effect * efta_effect * asym_effect)
                                         # equation (29) exponentiated

            elseif exporter == importer

                d[importer, exporter] = 1.0

            end

        end

    end

end

##########################################################################
##########################################################################

function gravity(tradedata; display = false)
    #function to perform basic gravity regression
    #assumes tradedata is a DataFrame, takes on the structure of EK dataset

    outreg = reg(tradedata, @formula(trade ~ fe(importer) + fe(exporter) +
         bin375 + bin750 + bin1500 + bin3000 + bin6000 + binmax  + border + 
                sharedlanguage + europeancom + efta + 0), save = true)

    lang_coef = outreg.coef[7:end]

    grp = groupby(outreg.fe, "exporter");

    S = get_first(grp, "fe_exporter")

    grp = groupby(outreg.fe, "importer");

    θm = get_first(grp, "fe_importer")
    
    θm = θm .+ S 

    norm_fe = sum(θm) / 19
    
    θm = θm .- norm_fe

    dist_bins = outreg.coef[1:6] .+ norm_fe

    if display == true

        println("Compare to Table III (1762)")
        println(" ")
        println("Distance Effects")
        dffoo = DataFrame(distance_effects = dist_bins);
        println(dffoo)
        println(" ")
        println("Border, language, Eupope, etc. Effects")
        dffoo = DataFrame(boder_lang_effects = lang_coef);
        println(dffoo)
        println(" ")
        println("Source and Destination Effects (The S's and θm's)")
        dffoo = DataFrame(source_effects = S, destination_effects = θm);
        println(dffoo)
    end

    return dist_bins, lang_coef, S, θm

end

##########################################################################
##########################################################################

function normalize_by_home_trade(df)
    
    home_trade = df[df.importer .== df.exporter,:].trade
    # grabs the home trade share by finding where importer = exporter
    
    df.trade .= df.trade / home_trade
    # normalize things
    
    return df
    
end 

##########################################################################
##########################################################################

function make_trade_shares(tradedata)
    # function to go from log, normalized trade data
    # back into the tradeshare matrix

    πdata = Array{Float64}(undef, 19, 19)
    fill!(πdata, 0.0)

    for importer = 1:19

        foo = tradedata[tradedata.importer .== importer, :]

        for exporter = 1:19

            if exporter != importer

            get_exporter = foo.exporter .== exporter

            πdata[importer, exporter] = exp(foo[get_exporter, : ].trade[1])

            end

        end

        hometrade = (1.0 + sum(πdata[importer, :]))^(-1.0)

        πdata[importer, :] = πdata[importer, :]*hometrade

        πdata[importer, importer] = hometrade

    end

    return πdata

end

##########################################################################
##########################################################################

function get_first(grp, variable)
       
    var = Array{Float64}(undef, length(grp))

    for i ∈ eachindex(grp)

        var[i[1]] = grp[i][1, variable]

    end
    
    return var

end 