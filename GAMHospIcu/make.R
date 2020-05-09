#!/usr/local/bin/Rscript

# date = Sys.Date()
date = "2020-05-08"

regioni = c("Piemonte", "Valle d'Aosta", "Liguria", 
            "Lombardia", "P.A. Bolzano", "P.A. Trento", 
            "Veneto", "Friuli Venezia Giulia", "Emilia-Romagna", 
            "Toscana", "Umbria", "Marche", "Lazio", "Abruzzo", 
            "Molise", "Campania", "Puglia", "Basilicata", 
            "Calabria", "Sicilia", "Sardegna")

# Fit by looping over regions ----
for(regione in regioni) 
{
  rmarkdown::render("HospitalisedICU.Rmd", 
                    params = list(regione = regione, date = date), 
                    output_file = paste0(regione, "_", date, ".html"))
}
