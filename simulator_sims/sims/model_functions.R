## @knitr models

R.utils::sourceDirectory("./data_generators")

make_my_model <- function(Nsamp, nullpi) {
  new_model(name = "pois_thin",
            label = sprintf("Poisson Thinned Sims with (Nsamp = %s, nullpi = %s)", Nsamp, nullpi),
            params = list(Nsamp = Nsamp, nullpi = nullpi,
                          path = "../../../data/gtex_tissue_gene_reads/"),
            simulate = rep_pois_thin(nsim, Nsamp, nullpi, path)
}
