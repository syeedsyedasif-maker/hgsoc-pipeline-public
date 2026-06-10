# ============================================================================
# Dockerfile — fully reproducible environment for the simulated HGSOC pipeline.
# Based on the official Bioconductor image matching R 4.3.x / Bioconductor 3.18.
#
# Unlike the bare Windows host this was developed on, this Linux image HAS a
# C/C++ toolchain, so the optional "primary" tools (inferCNV, BayesPrism,
# spacexr) can actually be compiled here. We install JAGS (for inferCNV) and
# opt into the optional tools so the container runs the PRIMARY paths, not just
# the fallbacks. Drop the INSTALL_OPTIONAL line for a faster, core-only build.
#
# Build:  docker build -t hgsoc-sim .
# Run:    docker run --rm -v "$PWD/out:/sim/results" hgsoc-sim
#         (writes results/ back to ./out; figs/ and truth/ likewise if mounted)
# ============================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# JAGS is a system dependency of inferCNV (the GitHub tools compile fine here).
RUN apt-get update && apt-get install -y --no-install-recommends jags \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /sim
COPY . /sim

# Install packages. install.R pulls the core stack from CRAN/Bioconductor 3.18;
# INSTALL_OPTIONAL=1 also builds the primary tools (slower). For an EXACT
# version pin instead, swap the line below for: renv::restore(prompt = FALSE)
RUN R -e "Sys.setenv(INSTALL_OPTIONAL='1'); source('install.R')"

# Run the whole pipeline and print the 7-check PASS/FAIL summary.
CMD ["Rscript", "run_all.R"]
