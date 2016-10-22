sd_section("Package overview",
  "",
  c(
    "photobiologyInOut-package"
    )
)

sd_section("Import data files from spectrometers",
  "We provide functions for importing data from some common brands and types of instruments.",
  c(
    "read_avaspec_csv",
    "read_licor_prn",
    "read_macam_dta",
    "read_oo_jazirrad",
    "read_oo_pidata",
    "read_oo_ssirrad"
  )
)

sd_section("Import data files from dataloggers",
           "We provide functions for importing data from some common brands and types of instruments.",
           c(
             "read_csi_dat"
           )
)

sd_section("Import data from simulation models",
           "We provide functions for importing spectral data from models we use. As the file format may vary depending on settings, these are mainly provided as examples.",
  c(
    "read_uvspec_disort",
    "read_uvspec_disort_vesa",
    "read_fmi_cum",
    "read_tuv_usrout"
    )
)

sd_section("Exchange data with base R and other R packages.",
           "",
           c(
             "mat2mspct",
             "mspct2mat",
             "colorSpec2mspct",
             "hyperSpec2mspct",
             "rspec2mspct"
           )
)
