#!/usr/local/bin/Rscript

regioni = c("Piemonte", "Valle d'Aosta", "Liguria", 
            "Lombardia", "P.A. Bolzano", "P.A. Trento", 
            "Veneto", "Friuli Venezia Giulia", "Emilia-Romagna", 
            "Toscana", "Umbria", "Marche", "Lazio", "Abruzzo", 
            "Molise", "Campania", "Puglia", "Basilicata", 
            "Calabria", "Sicilia", "Sardegna")

for(regione in regioni) 
{
  rmarkdown::render("Gompertz_model_fit.Rmd", 
                    params = list(regione = regione), 
                    output_file = paste0("Gompertz_model_fit", "_", regione, ".html"))
}
